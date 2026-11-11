"""Independent quadprog oracle for the MATLAB 1D/2D projection solver."""
from pathlib import Path
import numpy as np
import quadprog
from scipy.io import savemat

root = Path(__file__).resolve().parents[1]
rng = np.random.default_rng(20261108)
cases = []
for n in [1, 2]:
    for idx in range(40):
        center = rng.uniform(-.5, .5, n)
        normals = rng.normal(size=(15, n))
        C = np.vstack([np.eye(n), -np.eye(n), normals])
        d = np.r_[-np.ones(2*n), normals @ center - rng.uniform(.05, 1, 15)]
        if idx % 10 == 0:
            # Deliberately impossible relative to the upper input bound.
            C = np.vstack([C, np.eye(n)[0]])
            d = np.r_[d, 2.]
        desired = rng.uniform(-4, 4, n)
        # quadprog requires at least two coordinates, as in the source DI case.
        Cs = C if n == 2 else np.column_stack([C, np.zeros(C.shape[0])])
        qs = desired if n == 2 else np.r_[desired, 0.]
        try:
            u = quadprog.solve_qp(np.eye(2), qs, Cs.T, d)[0][:n]
            feasible = True
        except ValueError:
            u = np.full(n, np.nan)
            feasible = False
        cases.append(dict(C=C, d=d, up=desired, u=u, feasible=feasible, n=n))
savemat(root / "test_results" / "projection_reference.mat", {"cases": np.array(cases, dtype=object)})
print(f"Exported {len(cases)} seeded projection cases; {sum(c['feasible'] for c in cases)} feasible.")
