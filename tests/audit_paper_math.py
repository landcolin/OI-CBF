"""Audit original Python against independent paper expressions; export MATLAB fixtures.

Run after validate_source.py for both models. Original implementation is read-only.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.io import savemat

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "test_results"


def paper_aircraft_f(x, e):
    ph, th, ps, _, _, _, P, Nz = x
    return np.array([P + Nz * e.g / e.VT * np.sin(ph) * np.tan(th),
                     e.g / e.VT * (Nz * np.cos(ph) - np.cos(th)),
                     Nz * e.g * np.sin(ph) / (e.VT * np.cos(th)),
                     e.VT * np.cos(th) * np.cos(ps), e.VT * np.cos(th) * np.sin(ps),
                     e.VT * np.sin(th), -P / e.tp, -Nz / e.tz])


def capture(method, *args):
    local = {}
    code = method.__func__.__code__

    def profiler(frame, event, result):
        if event == "return" and frame.f_code is code:
            local.update(frame.f_locals)

    sys.setprofile(profiler)
    try:
        result = method(*args)
    finally:
        sys.setprofile(None)
    return result, local


def audit(model):
    sys.path.insert(0, str(ROOT / model))
    Simulation = importlib.import_module("main_sim").Simulation
    e = Simulation(verbose=False, bCBF_QP=False, OI_QP_flag=True)
    state = np.load(OUT / f"{model}_OI.npz")["x"]
    indices = np.linspace(0, state.shape[1] - 2, 12, dtype=int)
    X = state[:, indices]
    n = X.shape[0]
    f, g, J, primary, backup, hs, hb, dh, dhb, hi, dhi = [], [], [], [], [], [], [], [], [], [], []
    for x in X.T:
        f.append(e.f_x(x)); g.append(np.atleast_2d(e.g_x(x)).reshape(n, -1))
        J.append(e.computeJacobianSTM(x) if n == 2 else e.computeJacobianCL(x))
        primary.append(np.atleast_1d(e.primaryControl(x, 0) if n == 2 else e.primaryControl(x)))
        backup.append(np.atleast_1d(e.backupControl(x)))
        if n == 2:
            hs.append(e.h1_x(x)); hb.append(e.hb_x(x)); dh.append(e.grad_h1(x)); dhb.append(e.grad_hb(x))
            hi.append([e.hb_x(x)]); dhi.append(e.grad_hb(x)[None, :])
        else:
            hs.append(e.h_x_geo(x)); hb.append(e.hb_softmin(x, e.hb_funs, e.kappa))
            dh.append(e.grad_h_geo(x)); dhb.append(e.grad_hb_softmin(x, e.hb_funs, e.gradhb_funs, e.kappa))
            hi.append([fun(x) for fun in e.hb_funs]); dhi.append([fun(x) for fun in e.gradhb_funs])

    fixture = dict(x=X, f=np.array(f).T, g=np.stack(g, axis=2), J=np.stack(J, axis=2),
                   primary=np.array(primary).T, backup=np.array(backup).T, h=np.array(hs), hb=np.array(hb),
                   dh=np.array(dh), dhb=np.array(dhb), hi=np.array(hi).T, dhi=np.stack(dhi, axis=2))
    flows, Cs, ds, controls, qp_C, qp_d, qp_u = [], [], [], [], [], [], []
    flow_errors = []
    for idx in [0, 5, 10]:
        x = X[:, idx]
        desired = fixture["primary"][:, idx]
        if n == 2:
            desired = float(desired[0])
        S0 = np.eye(n)
        flow = e.integrateStateBackup(np.r_[x, S0.ravel()], e.tspan_b, e.int_options)
        S = flow[n:].reshape(n, n, -1)
        fx, gx = e.f_x(x), e.g_x(x)
        fp = fx + (gx * desired if n == 2 else gx @ desired)
        fb = fx + (gx * e.backupControl(x) if n == 2 else gx @ e.backupControl(x))
        directional = e.integrateStateBackupDirectional(np.r_[x, fp, fb], e.tspan_b, e.int_options)
        qp_stm = np.einsum("ijt,j->it", S, fp)
        qb_stm = np.einsum("ijt,j->it", S, fb)
        flow_errors.append(float(max(np.abs(qp_stm - directional[n:2*n]).max(), np.abs(qb_stm - directional[2*n:]).max())))
        flows.append(dict(x=x, phi=flow[:n], S=S, qp=directional[n:2*n], qb=directional[2*n:]))
        result, local = capture(e.OI_QP, x, desired)
        Cs.append(np.array(local["G"])[:, 0]); ds.append(np.array(local["h"])); controls.append(np.atleast_1d(result[0]))
        result, local = capture(e.vanilla_bCBF, x, desired)
        qp_C.append(np.array(local["G"])[:, :e.u_bounds.__len__()]); qp_d.append(np.array(local["h"]).ravel()); qp_u.append(np.atleast_1d(result[0]))
    fixture["flow"] = np.array(flows, dtype=object)
    fixture["OI_C"] = np.array(Cs).T; fixture["OI_d"] = np.array(ds).T; fixture["OI_u"] = np.array(controls).T
    fixture["QP_C"] = np.stack(qp_C, axis=2); fixture["QP_d"] = np.array(qp_d).T; fixture["QP_u"] = np.array(qp_u).T
    for mode in ["bCBF", "blending", "OI"]:
        data = np.load(OUT / f"{model}_{mode}.npz")
        fixture["trajectory_"+mode] = {key: data[key] for key in data.files}
    savemat(OUT / f"{model}_reference.mat", {"ref": fixture}, do_compression=True)

    rng = np.random.default_rng(20261108)
    dyn_errors = []
    for _ in range(40):
        x = e.x0.astype(float).copy()
        if n == 2:
            x += rng.uniform(-1, 1, 2)
            reference = np.array([x[1], 0])
        else:
            x += rng.uniform(-1, 1, 8) * [0.7, 0.3, 3, 4000, 4000, 500, 0.5, 0.5]
            reference = paper_aircraft_f(x, e)
        dyn_errors.append(float(np.max(np.abs(e.f_x(x) - reference))))
    checks = {"dynamics_against_paper_max_error": max(dyn_errors), "directional_vs_STM_max_absolute_errors": flow_errors}
    x_initial = e.x0.astype(float)
    ub_initial = e.backupControl(x_initial)
    fcl_initial = e.f_x(x_initial) + (e.g_x(x_initial) * ub_initial if n == 2 else e.g_x(x_initial) @ ub_initial)
    grad_initial = e.grad_h1(x_initial) if n == 2 else e.grad_h_geo(x_initial)
    checks["tightening_bound_at_default_state"] = {"configured_sup_fcl": float(e.sup_fcl),
        "actual_backup_dynamics_norm": float(np.linalg.norm(fcl_initial)),
        "actual_abs_safety_derivative": float(abs(grad_initial @ fcl_initial)),
        "is_upper_bound_at_this_state": bool(e.sup_fcl >= np.linalg.norm(fcl_initial))}
    if n == 8:
        x = e.x0.astype(float).copy(); x[5] += 150
        source_h3 = e.hb3_x(x); printed_h3 = e.c_H**2 - (e.H_star - x[5])**2
        checks["altitude_barrier_scaling"] = {"source": float(source_h3), "printed": float(printed_h3), "same_individual_zero_set": True, "same_softmin": False}
        x = e.x0.astype(float).copy(); x[2] = .5
        los = np.arctan2(e.pE_primary-x[4], e.pN_primary-x[3]); err = np.arctan2(np.sin(los-x[2]), np.cos(los-x[2]))
        printed_phi = e.phimax * np.tanh(e.VT/e.g * e.kpsi*err/e.phimax)
        source_phi = e.phimax * np.tanh(e.VT/e.g * e.kpsi*abs(err)/e.phimax)
        checks["nominal_heading_absolute_value"] = {"heading_error": float(err), "source_phi_desired": float(source_phi), "signed_phi_desired": float(printed_phi)}
        x = np.array([e.phi_star, 0, 0, 26000, 0, e.H_star, 0, e.Nz_star])
        checks["printed_backup_set_not_subset_of_safe_set"] = {"x": x.tolist(), "rho": float(e.rho), "h": float(e.h_x_geo(x)), "h6": float(e.hb6_x(x)), "exact_equalities_of_eq33": True}
        x[3] = 24990
        backup_path = e.integrateStateBackup(x, e.tspan_b, e.int_options)
        checks["safe_backup_set_state_leaves_geofence"] = {"x": x.tolist(), "initial_h": float(e.h_x_geo(x)), "initial_h6": float(e.hb6_x(x)), "initial_softmin_hb": float(e.hb_softmin(x,e.hb_funs,e.kappa)), "minimum_h_over_backup": float(e.h_x_geo_array(backup_path).min())}
        u = e.backupControl(x)
        checks["soft_saturation_equilibrium_residual"] = {"command_minus_turn_equilibrium": (u-np.array([0,e.Nz_star])).tolist()}
    else:
        x = np.array([-.128, .5]); result, local = capture(e.OI_QP,x,1.)
        C = np.array(local["G"])[:,0]; d = np.array(local["h"])
        checks["recoverable_but_QP_infeasible"] = {"x": x.tolist(), "continuous_backup_hI": .003, "tightening": .002,
            "mu": float(local["mu_starKKT"]), "u": float(result[0]), "minimum_constraint_residual": float(np.min(C*local["mu_starKKT"]-d)),
            "backup_mu1_minimum_constraint_residual": float(np.min(C-d))}
    checks["source_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((ROOT / model).glob("*.py"))}
    (OUT / f"{model}_paper_audit.json").write_text(json.dumps(checks,indent=2),encoding="utf-8")
    print(json.dumps(checks,indent=2),flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("model", choices=["aircraft","double_integrator"])
    audit(parser.parse_args().model)
