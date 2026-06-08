#!/usr/bin/env python3
"""capture.py — mpmath + sympy oracle for corpus_pathnet_walls_rows_test.jl.

Family I (docs/test_corpus/02_corpus_extension_plan.md:345-385): PATH-NETWORK
routing targets that stress the FW 2011 §3.1 spanning-tree walk around
structured pole sets.  Four backgrounds, each with an exact closed form so both
the EVALUATED value and the CHOSEN ROUTE (no visited node lands on a known pole;
no spurious Froissart doublet) are independently pinnable.

This script is the residual==0 GATE plus the numerical-pin source.  No value is
emitted before sympy proves the recast residual ≡ 0.  mpmath at dps=50.

  CPN.1  logistic vertical wall.   u' = u(1-u),  recast  u'' = u(1-u)(1-2u),
         u = 1/(1 + c·e^{-z}),  c = 2.  Pole WALL at Re z = log 2, Im z = π(2k+1).
         IC z0 = 0:  u(0) = 1/3,  u'(0) = 2/9 (exact rationals, inline in test).
         Targets straddling the wall {2±2i, -1±2i}.

  CPN.4  near-coalescent pair, δ-sweep.  Exact u = 1/(z-a) + 1/(z-(a+δ)), a = 1.
         The 2nd-order forcing is the explicit z-dependent
         f(z,u,up) = 2/(z-a)^3 + 2/(z-(a+δ))^3  (= u'' for this u; residual ≡ 0).
         IC z0 = 0.2.  δ ∈ {0.4, 0.1, 0.02, 0.005}.  Past-pole EXACT rationals
         (e.g. δ=0.1: u(0.5) = -11/3, u(1.5) = 9/2) — written inline in the test.

  CPN.8  off-node dense-eval (reuse logistic c=2).  Spot values strictly between
         visited tree nodes, pinned to the closed form.

  CPN.2  Airy-Riccati semi-infinite real row.  f = 2·u·up - 1,  u = -Ai'(z)/Ai(z)
         (recast of the Riccati u' = u^2 - z).  IC z0 = 0.5.  Poles at the Airy
         zeros a_k on the NEGATIVE real axis; between zeros u has its own zero
         (small |u|) next to a pole spike — the adversarial case for :min_u.
         Targets {-1, -1.5, -3, -5, -2+0.5i}.

Run:  python3 external/probes/corpus-oracles/pathnet-routing/capture.py
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 50


def emit(label, val):
    if isinstance(val, mp.mpc):
        print(f"const {label} = complex({mp.nstr(val.real, 45)}, "
              f"{mp.nstr(val.imag, 45)})")
    else:
        print(f"const {label} = {mp.nstr(val, 45)}")


# ===========================================================================
# (1) RESIDUAL==0 GATES — sympy, symbolic.  Run BEFORE any value is emitted.
# ===========================================================================
z, c, a, d = sp.symbols('z c a d')

# CPN.1 logistic recast.  u = 1/(1 + c e^{-z}) must satisfy u'' = u(1-u)(1-2u)
# AND the first-order u' = u(1-u).
u_log = 1 / (1 + c * sp.exp(-z))
res1 = sp.simplify(sp.diff(u_log, z) - u_log * (1 - u_log))
res2 = sp.simplify(sp.diff(u_log, z, 2) - u_log * (1 - u_log) * (1 - 2 * u_log))
assert res1 == 0, f"CPN.1 1st-order residual != 0: {res1}"
assert res2 == 0, f"CPN.1 2nd-order recast residual != 0: {res2}"

# CPN.4 rational two-pole.  u = 1/(z-a) + 1/(z-(a+d)) must satisfy the explicit
# forcing  u'' = 2/(z-a)^3 + 2/(z-(a+d))^3.  (Equivalently f(z,u,up) reduces to
# the z-only forcing since u'' is the exact 2nd derivative.)
b = a + d
u_rat = 1 / (z - a) + 1 / (z - b)
f_rat = 2 / (z - a) ** 3 + 2 / (z - b) ** 3
res_rat = sp.simplify(sp.diff(u_rat, z, 2) - f_rat)
assert res_rat == 0, f"CPN.4 forcing residual != 0: {res_rat}"

# CPN.2 Airy-Riccati.  u = -Ai'(z)/Ai(z) solves u' = u^2 - z (Ai'' = z·Ai), so
# u'' = 2 u u' - 1 = f(z,u,up).  Verify symbolically via the Airy DE: with
# w(z) = Ai(z), u = -w'/w, u' = -(w''/w) + (w'/w)^2 = -z + u^2 (w'' = z w).
w = sp.Function('w')
u_sym = -w(z).diff(z) / w(z)
up_sym = sp.diff(u_sym, z).subs(w(z).diff(z, 2), z * w(z))
res_airy = sp.simplify(up_sym - (u_sym ** 2 - z))
assert res_airy == 0, f"CPN.2 Riccati residual != 0: {res_airy}"
# and the 2nd-order recast u'' = 2 u u' - 1 follows by differentiating u'=u^2-z:
#   u'' = 2 u u' - 1.  (trivially exact)

print("# residual==0 GATES PASSED (CPN.1 logistic, CPN.4 rational, CPN.2 Airy)")
print()


# ===========================================================================
# (2) CPN.1 — logistic vertical wall, c = 2.  u(z) = 1/(1 + 2 e^{-z}).
#     Wall at Re z = log 2, Im z = π(2k+1).  Targets straddling the wall.
# ===========================================================================
C1 = mp.mpf(2)


def u_log_cf(zz):
    return 1 / (1 + C1 * mp.e ** (-zz))


print("# CPN.1  logistic c=2:  u(z) = 1/(1 + 2 e^{-z}).  Wall at Re z = log 2.")
emit("CPN1_LOG2", mp.log(2))                        # wall real-part
emit("CPN1_U_2p2i",   u_log_cf(mp.mpc(2,  2)))
emit("CPN1_U_2m2i",   u_log_cf(mp.mpc(2, -2)))
emit("CPN1_U_m1p2i",  u_log_cf(mp.mpc(-1, 2)))
emit("CPN1_U_m1m2i",  u_log_cf(mp.mpc(-1, -2)))
print()


# ===========================================================================
# (8) CPN.8 — off-node dense-eval (reuse logistic c=2).  Spot values strictly
#     between visited tree nodes, pinned to the closed form.  Two interior
#     reals on the real axis (pole-free) used as off-node probes.
# ===========================================================================
print("# CPN.8  logistic c=2 off-node dense-eval (interior, real-axis).")
emit("CPN8_U_1p17",  u_log_cf(mp.mpf('1.17')))      # ~ near a visited node + h/3
emit("CPN8_U_0p83",  u_log_cf(mp.mpf('0.83')))
print()


# ===========================================================================
# (2) CPN.2 — Airy-Riccati.  u(z) = -Ai'(z)/Ai(z).  IC z0 = 0.5.
#     Poles at the Airy zeros a_k (negative axis); targets march past the row.
# ===========================================================================
def u_airy_cf(zz):
    return -mp.airyai(zz, 1) / mp.airyai(zz)        # airyai(z,1) = Ai'(z)


z0_airy = mp.mpf('0.5')
u0_airy = u_airy_cf(z0_airy)
up0_airy = u0_airy ** 2 - z0_airy                   # u' = u^2 - z

print("# CPN.2  Airy-Riccati:  u(z) = -Ai'(z)/Ai(z),  IC z0 = 0.5.")
emit("CPN2_U0",  u0_airy)                            # IC u(0.5)
emit("CPN2_UP0", up0_airy)                           # IC u'(0.5) = u0^2 - 0.5
emit("CPN2_U_m1",   u_airy_cf(mp.mpf('-1')))
emit("CPN2_U_m1p5", u_airy_cf(mp.mpf('-1.5')))
emit("CPN2_U_m3",   u_airy_cf(mp.mpf('-3')))
emit("CPN2_U_m5",   u_airy_cf(mp.mpf('-5')))
emit("CPN2_U_m2p0p5i", u_airy_cf(mp.mpc(-2, '0.5')))
print()

# Airy zeros a_k (the negative-axis pole row).  The walk must NOT land a visited
# node within eps of any of these.
print("# CPN.2  Airy zeros a_k (negative-axis pole row); k = 1..5.")
for k in range(1, 6):
    emit(f"CPN2_AIRY_ZERO_{k}", mp.airyaizero(k))
print()
print("# u'(0.5) sanity:  u0^2 - 0.5 =", mp.nstr(up0_airy, 30))
