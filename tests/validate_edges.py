"""Additional domain, sampling and plotting checks using saved original runs."""

import argparse
import json
import os
from pathlib import Path
import sys
import warnings

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def run(model):
    sys.path.insert(0, str(ROOT / model))
    from main_sim import Simulation

    out = ROOT / "test_results"
    os.environ["MPLCONFIGDIR"] = str(out / "mplconfig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from plotting_paper_comparison import PaperPlotter

    results = {}
    data = {mode: np.load(out / f"{model}_{mode}.npz") for mode in ["bCBF", "blending", "OI"]}
    envs = [Simulation(verbose=False) for _ in range(3)]
    envs[1].mu_array = data["blending"]["mu"]
    envs[2].mu_array = data["OI"]["mu"]
    # Files are written into a dedicated result folder, never over source assets.
    target = out / f"{model}_plots"
    (target / "plots").mkdir(parents=True, exist_ok=True)
    os.chdir(target)
    for latex in [True, False]:
        plt.rcdefaults()
        try:
            args = [data[m]["x"] for m in ["bCBF", "blending", "OI"]]
            args += [data[m]["u"] for m in ["bCBF", "blending", "OI"]]
            desired = np.zeros_like(data["bCBF"]["u"])
            for i in range(1, desired.shape[1]):
                x_prev = data["bCBF"]["x"][:, i - 1]
                desired[:, i] = envs[0].primaryControl(x_prev, (i - 1) * envs[0].del_t) if model == "double_integrator" else envs[0].primaryControl(x_prev)
            args += [desired]
            if model == "double_integrator":
                args += [envs[1].mu_array, envs[2].mu_array, envs[0]]
            else:
                args += envs
            with warnings.catch_warnings(record=True) as ws:
                warnings.simplefilter("always")
                PaperPlotter.plotter(*args, latex_plots=latex, save_plots=True, show_plots=False)
            results[f"plotting_latex_{latex}"] = {"completed": True, "warnings": sorted(set(str(w.message) for w in ws))}
        except Exception as exc:
            results[f"plotting_latex_{latex}"] = {"completed": False, "exception": type(exc).__name__, "message": str(exc)}
        finally:
            plt.close("all")

    if model == "double_integrator":
        env = Simulation(verbose=False, bCBF_QP=False, OI_QP_flag=True)
        # Exact continuous backup trajectory: p(t)=p0+v0*t-t*t/2.
        cases = []
        for p in [-.2, -.15, -.14, -.13, -.128, -.127]:
            x = np.array([p, .5])
            u, _, _ = env.OI_QP(x, 1.)
            exact_hI = min(-p - .5**2 / 2, 1 - .5)
            cases.append({"x": x.tolist(), "exact_hI": exact_hI, "tightening": env.del_t / 2 * env.Lh_const * env.sup_fcl, "mu": float(env.mu_array[-1]), "u": float(u), "within_bounds": bool(abs(u) <= 1 + 1e-8)})
        results["recoverable_boundary_states"] = cases
        dt_runs = []
        for dt in [.01, .04]:
            env = Simulation(verbose=False, bCBF_QP=False, OI_QP_flag=True)
            env.del_t = dt
            env.total_steps = round(7 / dt) + 1
            env.setupASIF()
            x, _, _, u, *_ = env.sim()
            # Exact maximum position within every zero-order-held control interval.
            p, v, a = x[0, :-1], x[1, :-1], u[0, 1:]
            t = np.divide(-v, a, out=np.zeros_like(v), where=a != 0)
            valid = (a < 0) & (t > 0) & (t < dt)
            peak = max(float(np.max(x[0])), float(np.max((p + v * t + .5 * a * t**2)[valid])) if valid.any() else -np.inf)
            dt_runs.append({"dt": dt, "finite": bool(np.isfinite(x).all()), "min_continuous_safety_margin": -peak, "max_abs_u": float(np.abs(u).max())})
        results["sampling_variation"] = dt_runs
    else:
        env = Simulation(verbose=False, bCBF_QP=False, OI_QP_flag=True)
        mismatches = []
        for idx in [0, 300, 600]:
            x = data["OI"]["x"][:, idx]
            captured = {}

            def capture(frame, event, arg):
                if event == "return" and frame.f_code.co_name == "OI_QP":
                    captured.update(frame.f_locals)

            sys.setprofile(capture)
            try:
                env.OI_QP(x, env.primaryControl(x))
            finally:
                sys.setprofile(None)
            used_rhs = captured["h_temp_i"]
            expected_rhs = -(captured["hdot_up"] + env.alpha_b(captured["hb_phi"]))
            g = np.array(captured["G"])[:, 0]
            h = np.array(captured["h"])
            h[-1] = expected_rhs
            corrected_mu = float(np.max(h[g > 0] / g[g > 0]))
            mismatches.append({"state_index": idx, "terminal_hb": float(captured["hb_phi"]), "used_rhs": float(used_rhs), "alpha_b_rhs": float(expected_rhs), "actual_mu": float(captured["mu_starKKT"]), "mu_with_alpha_b": corrected_mu})
        results["OI_terminal_constraint_alpha_mismatch"] = mismatches
        numeric_cases = []
        for label, state_index, value in [("P_minus_30", 6, -30.), ("H_12000", 5, 12000.)]:
            x = env.x0.astype(float).copy()
            x[state_index] = value
            with warnings.catch_warnings(record=True) as ws:
                warnings.simplefilter("always")
                u = env.backupControl(x)
                hb = env.hb_softmin(x, env.hb_funs, env.kappa)
                grad = env.grad_hb_softmin(x, env.hb_funs, env.gradhb_funs, env.kappa)
            numeric_cases.append({"case": label, "u": [str(v) for v in u], "hb": str(hb), "gradient_finite": bool(np.isfinite(grad).all()), "warnings": sorted(set(str(w.message) for w in ws))})
        results["off_nominal_numerics"] = numeric_cases

    (out / f"{model}_edges.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=["double_integrator", "aircraft"])
    run(parser.parse_args().model)
