"""
capture.py — independent mpmath oracle for the Noumi-Yamada A_4^(1) GENERIC
(non-rational) transcendent.  test/corpus_noumi_yamada_test.jl (epic
padetaylor-p1v0, bead padetaylor-lgpc).  Run: python3 capture.py

Catalogue ny-a4-generic-mpmath (01_corpus_catalogue.md:752-767):
  A_4^(1): f_j' = f_j*(f_{j+1}-f_{j+2}+f_{j+3}-f_{j+4}) + alpha_j  (mod 5)
  alpha = (0.1, 0.15, 0.2, 0.25, 0.3)   (sum = 1)
  IC at t0=1: f(1) = (0.3, 0.25, 0.15, 0.2, 0.1)   (sum = 1 = t0)

The RHS bracket matches src/NoumiYamada.jl noumi_yamada_rhs (n=2):
  bracket_j = sum_{r=1}^{n} f_{j+2r-1} - sum_{r=1}^{n} f_{j+2r}
            = f_{j+1} - f_{j+2} + f_{j+3} - f_{j+4}.
This is the FIRST external (non-rational) validation of the NY vector walk.
We integrate from t0=1 to a DOWNSTREAM t with mpmath.odefun at high dps and
pin f(t) for the Julia cross-check.  Constraint sum(f)=t is self-checked.
"""
from mpmath import mp, mpf, odefun

mp.dps = 40
n = 2
d = 2*n + 1
alpha = [mpf('0.1'), mpf('0.15'), mpf('0.2'), mpf('0.25'), mpf('0.3')]
f0 = [mpf('0.3'), mpf('0.25'), mpf('0.15'), mpf('0.2'), mpf('0.1')]

assert sum(alpha) == 1, "alpha must sum to 1"
assert sum(f0) == 1, "f0 must sum to t0=1 (constraint)"

def rhs(t, f):
    out = []
    for j in range(d):
        bracket = sum(f[(j + 2*r - 1) % d] for r in range(1, n+1)) \
                  - sum(f[(j + 2*r) % d] for r in range(1, n+1))
        out.append(f[j]*bracket + alpha[j])
    return out

sol = odefun(rhs, mpf('1'), f0, tol=mpf('1e-30'))

print("# === ny-a4-generic-mpmath: A_4^(1) transcendent, mpmath odefun dps=40 ===")
print("# alpha =", [float(a) for a in alpha], " sum =", float(sum(alpha)))
print("# f0(t0=1) =", [float(x) for x in f0])
# Pin at a few downstream points; keep |t-t0| modest so the walk stays
# pole-free (the generic transcendent has movable poles for large |t-t0|).
for tq in ('1.2', '1.4', '1.5'):
    t = mpf(tq)
    f = sol(t)
    constraint = sum(f) - t
    print(f"# --- t = {tq} ---")
    for j in range(d):
        print(f"#   f[{j}] = {f[j]}")
    print(f"#   sum(f)-t = {constraint}   (constraint, must be ~0)")
    print(f"#   sum(f)-t (float) = {float(constraint):.3e}")
print("# === done ===")
