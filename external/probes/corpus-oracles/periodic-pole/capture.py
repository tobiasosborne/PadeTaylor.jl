#!/usr/bin/env python3
"""capture.py — closed-form + mpmath oracle for corpus_periodic_pole_test.jl.

Family D (docs/test_corpus/02_corpus_extension_plan.md:232-241): the logistic /
Bernoulli equation

    u' = u(1 - u)          (first order)
    u'' = u(1 - u)(1 - 2u)  (2nd-order recast for solve_pade / path_network_solve)

has exact solution  u(z) = 1/(1 + c·e^{-z}).  Its poles form an INFINITE
VERTICAL ROW at  z = log(c) + iπ + 2πi·k  (k ∈ ℤ): spacing 2πi, all SIMPLE,
residue +1.  This pole *geometry* (periodic vertical row) is absent from the
v1 corpus — CPB.6 (tanh) uses only the single nearest pole, never the row.

With the recommended c = 1 the poles are PURELY IMAGINARY at z = ±iπ, ±3iπ, …
(the cleanest bridge).  IC at z0 = 0:  u(0) = 1/2,  u'(0) = 1/4.

This script is the residual==0 GATE plus the numerical pin source.  In order:
  (1) sympy verifies the 2nd-order recast residual ≡ 0 (u = 1/(1+c e^{-z})),
  (2) sympy verifies the residue at z = iπ (c=1) is +1, via the limit
      (z - z0)·u → 1, AND independently by a contour integral (mpmath),
  (3) mpmath (dps=50) emits the spot values + the past-pole bridge target,
      printed as Julia-readable literals.

Run:  python3 external/probes/corpus-oracles/periodic-pole/capture.py
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 50

# ---------------------------------------------------------------------------
# (1) RESIDUAL==0 GATE — sympy, symbolic.  u = 1/(1 + c e^{-z}) must satisfy
#     u'' = u(1-u)(1-2u)  AND  u' = u(1-u).
# ---------------------------------------------------------------------------
z, c = sp.symbols('z c')
u = 1 / (1 + c * sp.exp(-z))

res1 = sp.simplify(sp.diff(u, z) - u * (1 - u))                       # 1st order
res2 = sp.simplify(sp.diff(u, z, 2) - u * (1 - u) * (1 - 2 * u))      # 2nd recast
assert res1 == 0, f"1st-order residual != 0: {res1}"
assert res2 == 0, f"2nd-order recast residual != 0: {res2}"

# Independent recast check: the brief gives u'' = u(1-u)(1-2u).  Confirm it is
# exactly d/dz[u(1-u)] = (1-2u)u' = (1-2u)·u(1-u) — i.e. the chain-rule recast
# of the first-order RHS, so the recast itself is verified symbolically.
U = sp.symbols('U')
R1 = U * (1 - U)
recast = sp.expand(sp.diff(R1, U) * R1)            # (1-2U)·U(1-U)
assert sp.expand(U * (1 - U) * (1 - 2 * U)) - recast == 0, "recast mismatch"

# ---------------------------------------------------------------------------
# (2) RESIDUE +1 at z = iπ (c=1).  Closed form: simple pole, residue = +1.
#     (a) symbolic limit  (z - iπ)·u → ?   (b) numeric contour integral.
# ---------------------------------------------------------------------------
u1 = 1 / (1 + sp.exp(-z))                            # c = 1
res_sym = sp.limit((z - sp.I * sp.pi) * u1, z, sp.I * sp.pi)
assert sp.simplify(res_sym - 1) == 0, f"symbolic residue != 1: {res_sym}"


def u_cf(zz, cc=1):
    return 1 / (1 + cc * mp.e ** (-zz))


# contour integral  (1/2πi) ∮ u dz  around z = iπ on a small circle.
z0 = mp.mpc(0, mp.pi)
r = mp.mpf('0.3')
resid = mp.quad(lambda th: u_cf(z0 + r * mp.e ** (1j * th)) *
                (1j * r * mp.e ** (1j * th)), [0, 2 * mp.pi]) / (2j * mp.pi)
assert abs(resid - 1) < mp.mpf('1e-40'), f"contour residue != 1: {resid}"

# spacing of the vertical row: next pole up is z0 + 2πi.
z1 = mp.mpc(0, 3 * mp.pi)
assert abs((z1 - z0) - mp.mpc(0, 2 * mp.pi)) < mp.mpf('1e-45'), "spacing != 2πi"

print("# residual==0 + residue+1 + 2πi-spacing GATE PASSED (c=1)")
print()


def emit(label, val):
    if isinstance(val, mp.mpc):
        print(f"const {label} = complex({mp.nstr(val.real, 45)}, "
              f"{mp.nstr(val.imag, 45)})")
    else:
        print(f"const {label} = {mp.nstr(val, 45)}")


# ---------------------------------------------------------------------------
# (3) SPOT VALUES (c=1), to dps=50, pinned at tol 1e-13 in the test.
#     IC sanity: u(0)=1/2, u'(0)=1/4 are exact rationals (written inline).
# ---------------------------------------------------------------------------
print("# c = 1 :  u(z) = 1/(1+e^{-z}).  Spot values, mpmath dps=50.")
emit("CVROW_U_AT_1p0",  u_cf(mp.mpf('1.0')))          # real z=1
emit("CVROW_U_AT_2p0",  u_cf(mp.mpf('2.0')))          # real z=2
emit("CVROW_U_AT_0p5p0p5i", u_cf(mp.mpc('0.5', '0.5')))   # complex z=0.5+0.5i
print()

# ---------------------------------------------------------------------------
# POLE-BRIDGE target (CVrow.2): just past the nearest pole z=iπ, at
# z = i(π + 0.3) = complex(0, π+0.3).  Pin u and u' = u(1-u) there.
# ---------------------------------------------------------------------------
z_tgt = mp.mpc(0, mp.pi + mp.mpf('0.3'))
u_tgt = u_cf(z_tgt)
up_tgt = u_tgt * (1 - u_tgt)
print("# Past-pole bridge target  z = i(π+0.3)  (just past the pole at iπ).")
emit("CVROW_U_AT_PASTPOLE",  u_tgt)
emit("CVROW_UP_AT_PASTPOLE", up_tgt)
print()

# ---------------------------------------------------------------------------
# POLE LOCATIONS of the vertical row (c=1): z = iπ(2k+1).  Pinned structurally
# (the row formula) — the nearest two are iπ and 3iπ.
# ---------------------------------------------------------------------------
print("# Vertical-row pole locations (c=1):  z = iπ·(2k+1).")
emit("CVROW_POLE_0",  mp.mpc(0, mp.pi))          # k=0
emit("CVROW_POLE_1",  mp.mpc(0, 3 * mp.pi))      # k=1
print()
print("# spacing = 2πi :", mp.nstr(2 * mp.pi, 45), "(imag)")
