#!/usr/bin/env python3
# ============================================================================
# Oracle capture for the orthogonal-polynomial BVP corpus (CBvx.1).
#
# GROUND TRUTH (Rule 5): for each named polynomial P(z) we (a) verify the
# ODE residual  u'' - f(z, u, u')  is IDENTICALLY ZERO as a symbolic rational
# function (the gate below refuses to emit unless residual.simplify() == 0),
# and (b) print the EXACT RATIONAL boundary + interior oracle values.  No
# transcendental constants appear anywhere: every pinned number is a Python
# Fraction / sympy Rational.  These are independent of the Julia routine.
#
# References:  docs/test_corpus/02_corpus_extension_plan.md Family H (CBvx.1);
#              DLMF 18.8 (Hermite/Legendre/Chebyshev ODEs, table of σ, τ).
#
# The three ODEs (n = 5, physicists'/standard normalisations):
#   Hermite    H_5 :  u'' = 2 z u'           - 2 n u           (regular ∀ z)
#   Legendre   P_5 :  u'' = (2 z u' - n(n+1) u) / (1 - z^2)    (sing. at ±1)
#   Chebyshev  T_5 :  u'' = (  z u' -   n^2 u) / (1 - z^2)     (sing. at ±1)
# ============================================================================
import sympy as sp

z = sp.symbols('z')
n = 5

# Named degree-5 polynomials (exact integer/rational coefficients).
H5 = 32*z**5 - 160*z**3 + 120*z                       # physicists' Hermite
P5 = (63*z**5 - 70*z**3 + 15*z) / 8                    # Legendre
T5 = 16*z**5 - 20*z**3 + 5*z                           # Chebyshev (first kind)


def rhs_hermite(u):
    return 2*z*sp.diff(u, z) - 2*n*u


def rhs_legendre(u):
    return (2*z*sp.diff(u, z) - n*(n + 1)*u) / (1 - z**2)


def rhs_chebyshev(u):
    return (z*sp.diff(u, z) - n**2*u) / (1 - z**2)


def emit(label, poly, rhs):
    """Gate on residual == 0, then print exact-rational BC + interior pins."""
    upp = sp.diff(poly, z, 2)
    residual = sp.simplify(upp - rhs(poly))
    assert residual == 0, f"{label}: ODE residual NOT identically 0: {residual}"
    print(f"# {label}: ODE residual u''-f == 0  VERIFIED symbolically")
    print(f"#   poly = {sp.expand(poly)}")
    return poly


def val(poly, q):
    """Exact value of poly at rational q, as a sympy Rational."""
    return sp.nsimplify(poly.subs(z, sp.Rational(q)))


# --- Hermite: segment [-1, 1], BC u(-1)=8, u(1)=-8 --------------------------
H = emit("Hermite H_5", H5, rhs_hermite)
print(f"#   BC  u(-1) = {val(H,-1)}   u(1) = {val(H,1)}")
print(f"#   int u(-1/2) = {val(H,sp.Rational(-1,2))}"
      f"   u(1/4) = {val(H,sp.Rational(1,4))}"
      f"   u(1/2) = {val(H,sp.Rational(1,2))}")
print()

# --- Legendre: segment [-4/5, 4/5] (strictly inside ±1) ---------------------
P = emit("Legendre P_5", P5, rhs_legendre)
print(f"#   BC  u(-4/5) = {val(P,sp.Rational(-4,5))}   u(4/5) = {val(P,sp.Rational(4,5))}")
print(f"#   int u(-2/5) = {val(P,sp.Rational(-2,5))}"
      f"   u(0) = {val(P,0)}"
      f"   u(2/5) = {val(P,sp.Rational(2,5))}")
print()

# --- Chebyshev: segment [-9/10, 9/10] (strictly inside ±1) ------------------
C = emit("Chebyshev T_5", T5, rhs_chebyshev)
print(f"#   BC  u(-9/10) = {val(C,sp.Rational(-9,10))}   u(9/10) = {val(C,sp.Rational(9,10))}")
print(f"#   int u(-1/2) = {val(C,sp.Rational(-1,2))}"
      f"   u(1/4) = {val(C,sp.Rational(1,4))}"
      f"   u(1/2) = {val(C,sp.Rational(1,2))}")
print()
print("# ALL THREE RESIDUALS == 0; oracle rationals above are independent ground truth.")
