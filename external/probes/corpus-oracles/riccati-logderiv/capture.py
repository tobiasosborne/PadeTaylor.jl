#!/usr/bin/env python3
"""capture.py — exact-rational oracle for corpus_riccati_rational_test.jl.

Family A (docs/test_corpus/02_corpus_extension_plan.md:134-162): a Riccati
u = w'/w is the log-derivative of a linear ODE w''+p w'+q w = 0, with a SIMPLE
pole of residue +1 at every zero of w.  Two RATIONAL members (the solution is a
ratio of polynomials, so every pinned value is an EXACT rational at machine
precision — no transcendental floor):

  CRic.3  Hermite  u = H_4'/H_4 = 2n H_{n-1}/H_n   (physicists' H, n=4)
  CRic.4  Chebyshev u = T_4'/T_4                    (first kind, n=4)

Both are recast to 2nd order for solve_pade via  u'' = ∂R/∂z + (∂R/∂u)·R  where
R(z,u) is the first-order Riccati RHS.  This script is the residual==0 GATE: it
verifies (a) the first-order RHS, (b) the 2nd-order recast against the exact
u'' (both symbolic residuals ≡ 0) BEFORE emitting any number, then prints the
Julia-readable exact rationals and pole locations the test pins inline.

Run:  python3 external/probes/corpus-oracles/riccati-logderiv/capture.py
"""
import sympy as sp

z, U = sp.symbols('z U')


def riccati_recast(R1_of_zU):
    """u'' = ∂R/∂z + (∂R/∂u)·R  (chain rule through u' = R)."""
    return sp.diff(R1_of_zU, z) + sp.diff(R1_of_zU, U) * R1_of_zU


# ----------------------------------------------------------------------------
# CRic.3  Hermite  (n = 4)
# ----------------------------------------------------------------------------
n = 4
H4 = sp.expand(sp.hermite(4, z))                 # 16 z^4 - 48 z^2 + 12
u3 = sp.diff(H4, z) / H4
up3 = sp.diff(u3, z)

R1_h = -U**2 + 2*z*U - 2*n                        # u' = -u^2 + 2 z u - 2n
assert sp.simplify(up3 - R1_h.subs(U, u3)) == 0, "CRic.3 1st-order RHS"

f2_h = sp.expand(riccati_recast(R1_h))            # n=4: 2U^3-6U^2 z+4U z^2+18U-16z
brief_h = 2*U**3 - 6*U**2*z + 4*U*z**2 + 18*U - 16*z
assert sp.expand(f2_h - brief_h) == 0, "CRic.3 recast vs brief"
assert sp.simplify(sp.diff(u3, z, 2) - f2_h.subs(U, u3)) == 0, "CRic.3 2nd-order residual"
assert sp.simplify(u3.subs(z, 0)) == 0 and sp.simplify(up3.subs(z, 0)) == -8

# ----------------------------------------------------------------------------
# CRic.4  Chebyshev T  (n = 4)
# ----------------------------------------------------------------------------
T4 = sp.expand(sp.chebyshevt(4, z))               # 8 z^4 - 8 z^2 + 1
u4 = sp.diff(T4, z) / T4
up4 = sp.diff(u4, z)

R1_c = -U**2 + (z*U - n**2)/(1 - z**2)            # u' = -u^2 + (z u - n^2)/(1-z^2)
assert sp.simplify(up4 - R1_c.subs(U, u4)) == 0, "CRic.4 1st-order RHS"

f2_c = riccati_recast(R1_c)
# clean form: numerator over (1-z^2)^2 (DERIVED, not hand-copied):
num_c = sp.expand(sp.simplify(f2_c * (1 - z**2)**2))
assert sp.expand(num_c -
                 (2*U**3*(1 - z**2)**2 + 3*U**2*z*(z**2 - 1)
                  + U*(33 - 30*z**2) - 48*z)) == 0, "CRic.4 numerator form"
assert sp.simplify(sp.diff(u4, z, 2) - f2_c.subs(U, u4)) == 0, "CRic.4 2nd-order residual"
assert sp.simplify(u4.subs(z, 0)) == 0 and sp.simplify(up4.subs(z, 0)) == -16

print("# residual==0 GATE PASSED — emitting oracle values")
print("# CRic.4 recast numerator over (1-z^2)^2:")
print("#   2U^3(1-z^2)^2 + 3U^2 z (z^2-1) + U(33-30 z^2) - 48 z")


def rat(expr):
    return sp.nsimplify(sp.simplify(expr))


# Hermite zeros (poles of u3), residue +1 at each:
z3_in = sp.sqrt(6 - 2*sp.sqrt(6))/2               # inner pair ±
z3_out = sp.sqrt(2*sp.sqrt(6) + 6)/2              # outer pair ±
print("\n# CRic.3 Hermite poles (= H_4 zeros), residue +1:")
print(f"#   inner ±{float(z3_in):.16g}   outer ±{float(z3_out):.16g}")
print(f"#   residue inner = {sp.simplify(sp.limit((z - z3_in)*u3, z, z3_in))}")
print("# CRic.3 past/around-pole EXACT rationals (z : u, u'):")
for pt in [sp.Rational(3, 10), sp.Rational(4, 5), sp.Integer(1), sp.Rational(6, 5)]:
    print(f"#   z={pt}: u={rat(u3.subs(z, pt))}  up={rat(up3.subs(z, pt))}")

# Chebyshev zeros (poles of u4), residue +1; ±1 are RHS coeff-singularities:
z4_in = sp.cos(3*sp.pi/8)                         # inner pair ±0.38268...
z4_out = sp.cos(sp.pi/8)                          # outer pair ±0.92388...
print("\n# CRic.4 Chebyshev poles (= T_4 zeros = cos((2k-1)π/8)), residue +1:")
print(f"#   inner ±{float(z4_in):.16g}   outer ±{float(z4_out):.16g}")
print(f"#   residue inner = {sp.simplify(sp.limit((z - z4_in)*u4, z, z4_in))}")
print("# CRic.4 past/around-pole EXACT rationals inside (-1,1) (z : u, u'):")
for pt in [sp.Rational(1, 5), sp.Rational(1, 2), sp.Rational(3, 5), sp.Rational(7, 10)]:
    print(f"#   z={pt}: u={rat(u4.subs(z, pt))}  up={rat(up4.subs(z, pt))}")
