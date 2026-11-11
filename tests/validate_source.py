"""Run original simulations and independent numerical checks without editing them.

Usage: .test-environment/Scripts/python.exe tests/validate_source.py double_integrator
       .test-environment/Scripts/python.exe tests/validate_source.py aircraft
"""

import argparse
import ast
import json
from pathlib import Path
import sys
import time
import traceback
import warnings

import numpy as np
import quadprog
from scipy.optimize._numdiff import approx_derivative

ROOT = Path(__file__).resolve().parents[1]


def run(model):
    sys.path.insert(0, str(ROOT / model))
    from main_sim import Simulation

    report = {"model": model, "python": sys.version, "simulations": [], "checks": {}}
    output = ROOT / "test_results"
    output.mkdir(exist_ok=True)

    def save():
        (output / f"{model}.json").write_text(
            json.dumps(report, indent=2, allow_nan=False), encoding="utf-8"
        )

    # Syntax validation without writing pyc files.
    files = list((ROOT / model).glob("*.py"))
    for file in files:
        ast.parse(file.read_text(encoding="utf-8"), filename=str(file))
    report["checks"]["syntax"] = {"files": len(files), "passed": True}

    trajectories = {}
    for mode, method in [("bCBF", "vanilla_bCBF"), ("blending", "vanilla_blending"), ("OI", "OI_QP")]:
        env = Simulation(verbose=False, bCBF_QP=mode == "bCBF", regFun=mode == "blending", OI_QP_flag=mode == "OI")
        original = getattr(env, method)
        counters = {"calls": 0, "QP_fallbacks": 0}

        def monitored(x, u_des):
            result = original(x, u_des)
            counters["calls"] += 1
            counters["QP_fallbacks"] += result[2] is None
            if counters["calls"] % 100 == 0:
                print(f"{model} {mode}: {counters['calls']}/{env.total_steps - 1}", flush=True)
            return result

        setattr(env, method, monitored)
        start = time.perf_counter()
        try:
            with warnings.catch_warnings(record=True) as ws:
                warnings.simplefilter("always")
                result = env.sim()
            x, _, _, u, interventions, avg, maximum = result
            margins = -x[0] if model == "double_integrator" else env.h_x_geo_array(x)
            bounds = np.asarray(env.u_bounds)
            ua = u[:, 1:]  # Column zero is an unused placeholder in the source.
            violation = max(0., float(np.max(bounds[:, 0, None] - ua)), float(np.max(ua - bounds[:, 1, None])))
            mu = np.asarray(env.mu_array)
            item = {
                "mode": mode, "completed": True, "seconds": time.perf_counter() - start,
                "steps": x.shape[1] - 1, "simulated_seconds": (x.shape[1] - 1) * env.del_t,
                "finite": bool(np.isfinite(x).all() and np.isfinite(u).all()),
                "min_safety_margin": float(np.min(margins)),
                "safety_violations_below_minus_1e_8": int(np.count_nonzero(margins < -1e-8)),
                "max_input_bound_violation": violation,
                "input_min": ua.min(axis=1).tolist(), "input_max": ua.max(axis=1).tolist(),
                "mu_range": [float(mu.min()), float(mu.max())] if mu.size else None,
                "interventions": len(interventions), "returned_avg_time": avg, "returned_max_time": maximum,
                "warnings": sorted(set(str(w.message) for w in ws)), **counters,
            }
            trajectories[mode] = (env, result)
            np.savez_compressed(output / f"{model}_{mode}.npz", x=x, u=u, mu=mu)
        except Exception:
            item = {"mode": mode, "completed": False, "error": traceback.format_exc(), **counters}
        report["simulations"].append(item)
        save()
        print(json.dumps(item), flush=True)

    # Rejecting an empty controller selection currently has no explicit guard.
    env = Simulation(verbose=False, bCBF_QP=False, regFun=False, OI_QP_flag=False)
    env.total_steps = 2
    try:
        env.sim()
        report["checks"]["no_controller_selected"] = {"exception": None}
    except Exception as exc:
        report["checks"]["no_controller_selected"] = {"exception": type(exc).__name__, "message": str(exc)}

    # Independent Jacobian check against central differences of actual backup dynamics.
    env = Simulation(verbose=False)
    rng = np.random.default_rng(42)
    jac_errors, grad_errors = [], []
    for _ in range(25):
        x = env.x0.astype(float).copy()
        if model == "aircraft":
            x += rng.uniform(-1, 1, 8) * np.array([.8, .25, 3., 2000., 2000., 300., .5, .4])
            analytic = env.computeJacobianCL(x)
            numeric = approx_derivative(lambda y: env.f_x(y) + env.g_x(y) @ env.backupControl(y), x)
            pairs = [(env.h_x_geo, env.grad_h_geo)] + list(zip(env.hb_funs, env.gradhb_funs))
        else:
            x += rng.uniform(-1, 1, 2)
            analytic = env.computeJacobianSTM(x)
            numeric = approx_derivative(lambda y: env.f_x(y) + env.g_x(y) * env.backupControl(y), x)
            pairs = [(env.h1_x, env.grad_h1), (env.hb_x, env.grad_hb)]
        jac_errors.append(float(np.max(np.abs(analytic - numeric))))
        for fun, grad in pairs:
            grad_errors.append(float(np.max(np.abs(grad(x) - approx_derivative(fun, x).ravel()))))
    report["checks"]["finite_difference"] = {"states": 25, "max_jacobian_absolute_error": max(jac_errors), "max_barrier_gradient_absolute_error": max(grad_errors)}

    if model == "double_integrator":
        report["checks"]["OI_vs_bCBF_max_state_difference"] = float(np.max(np.abs(trajectories["OI"][1][0] - trajectories["bCBF"][1][0])))
        # These states satisfy h>=0 but lie outside the recoverable set for T=1.
        cases = []
        env = Simulation(verbose=False, bCBF_QP=False, OI_QP_flag=True)
        for x in [np.array([-.01, .5]), np.array([-2., 1.1]), np.array([-.127, .5])]:
            u, _, _ = env.OI_QP(x, 1.)
            cases.append({"x": x.tolist(), "u": float(u), "mu": float(env.mu_array[-1]), "within_bounds": bool(-1 - 1e-8 <= u <= 1 + 1e-8)})
        report["checks"]["OI_initial_state_domain"] = cases
    else:
        # Current soft-min should agree with stable logsumexp on finite inputs.
        from scipy.special import logsumexp
        x = env.x0.astype(float).copy()
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            soft = env.hb_softmin(x, env.hb_funs, env.kappa)
            grad = env.grad_hb_softmin(x, env.hb_funs, env.gradhb_funs, env.kappa)
        stable = -logsumexp(-env.kappa * np.array([fun(x) for fun in env.hb_funs])) / env.kappa
        report["checks"]["softmin_at_default_initial_state"] = {"source_value": str(soft), "stable_value": float(stable), "gradient_finite": bool(np.isfinite(grad).all()), "warnings": sorted(set(str(w.message) for w in ws))}
        x[6] = 30.
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            u = env.backupControl(x)
        report["checks"]["backup_saturation_extreme_P_30"] = {"finite": bool(np.isfinite(u).all()), "u": [str(v) for v in u], "warnings": sorted(set(str(w.message) for w in ws))}

    # Reference QP: compare the actual closed form to a solver on exactly the same constraints.
    if "OI" in trajectories:
        oi_env, oi_result = trajectories["OI"]
        discrepancies, residuals, infeasible = [], [], 0
        # Capture locals at the return of the original method, without changing its arithmetic.
        for idx in np.linspace(0, oi_result[0].shape[1] - 2, 12, dtype=int):
            x = oi_result[0][:, idx]
            desired = oi_env.primaryControl(x, 0) if model == "double_integrator" else oi_env.primaryControl(x)
            captured = {}

            def trace(frame, event, arg):
                if event == "return" and frame.f_code.co_name == "OI_QP":
                    captured.update(frame.f_locals)
                return trace

            sys.setprofile(trace)
            try:
                # Class method bypasses the monitoring wrapper.
                Simulation.OI_QP(oi_env, x, desired)
            finally:
                sys.setprofile(None)
            g = np.array(captured["G"])[:, 0]
            h = np.array(captured["h"])
            mu = captured["mu_starKKT"]
            residuals.append(float(np.min(g * mu - h)))
            try:
                reference = quadprog.solve_qp(np.eye(2), np.zeros(2), np.column_stack([g, np.zeros_like(g)]).T, h.astype(float))[0][0]
                discrepancies.append(float(abs(mu - reference)))
            except ValueError:
                infeasible += 1
        report["checks"]["OI_closed_form_vs_reference_QP"] = {"states": 12, "infeasible_QPs": infeasible, "max_mu_difference": max(discrepancies) if discrepancies else None, "min_constraint_residual": min(residuals)}

    save()
    print(json.dumps(report["checks"], indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("model", choices=["double_integrator", "aircraft"])
    run(parser.parse_args().model)
