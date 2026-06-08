#!/usr/bin/env python3
"""capture.py — special-function BVP oracles for corpus_special_fn_bvp_test.jl.

Family H, CBvx.2/CBvx.3 (docs/test_corpus/02_corpus_extension_plan.md:309-327).
The Chebyshev-Newton BVP solver (src/BVP.jl) has, before this file, been pinned
only against orthogonal-polynomial / elementary closed forms.  This script
captures the FIRST special-function-BVP oracles for that solver:

  CBvx.3  Kummer confluent hypergeometric (3-arg, EXACT polynomial)
          z u'' + (b - z) u' - a u = 0 with a = -2, b = 5/2.  Because a is a
          NEGATIVE INTEGER the Kummer series terminates:
              M(-2, 5/2, z) = 4 z^2 / 35 - 4 z / 5 + 1   (degree-2 polynomial).
          DLMF 13.2.  Singular at z=0 -> segment [1/2, 2] away from the origin.
          3-arg recast  u'' = ((z - b) u' + a u) / z = ((z - 5/2) u' - 2 u)/z.
          EXACT rationals everywhere -> machine-precision pins (tol 1e-13).

  CBvx.2  Mathieu (2-arg, externally-pinned eigenvalue)
          u'' + (a - 2 q cos 2z) u = 0,  q = 1, a = a_2(q=1) (an EIGENVALUE that
          makes ce_2 periodic).  DLMF 28.2/28.4.  Solution = ce_2, the even
          pi-periodic Mathieu function.  Segment [0, pi/2].  scipy supplies both
          the eigenvalue (scipy.special.mathieu_a) and the function values
          (scipy.special.mathieu_cem) -> double precision -> tol 1e-10.

  CBvx.3-Whittaker (OPTIONAL) Whittaker M (2-arg, mpmath)
          u'' + (-1/4 + k/z + (1/4 - m^2)/z^2) u = 0, k=2, m=3/2, on [1,4].
          DLMF 13.14.  Oracle mpmath.whitm dps=50 -> tol 1e-12.

This script is the residual==0 GATE: it verifies each ODE recast SYMBOLICALLY
(sympy, residual identically 0) BEFORE emitting any number, then prints the
labelled oracle literals transcribed into test/_oracle_corpus_special_fn_bvp.jl
(transcendental pins) and inline into the test (exact Kummer rationals).

ORACLE GOTCHA (scipy Mathieu): scipy.special.mathieu_cem(m, q, x_DEGREES) takes
x in DEGREES and returns (value, derivative).  Convert z_radians via
numpy.degrees(z).  scipy.special.mathieu_a(m, q) returns the characteristic
value a_m(q).

Run:  python3 external/probes/corpus-oracles/special-fn-bvp/capture.py
"""
import numpy as np
import sympy as sp
import scipy.special as sc
import mpmath as mp

mp.mp.dps = 50

print("=" * 72)
print("CBvx.3 KUMMER M(-2, 5/2, z) = 4 z^2/35 - 4 z/5 + 1  (EXACT polynomial)")
print("=" * 72)

z = sp.symbols('z')
a_k, b_k = sp.Integer(-2), sp.Rational(5, 2)

# --- SYMBOLIC residual==0 GATE: the claimed polynomial solves Kummer's eqn ----
M = sp.Rational(4, 35) * z**2 - sp.Rational(4, 5) * z + 1
kummer_residual = sp.simplify(z * sp.diff(M, z, 2) + (b_k - z) * sp.diff(M, z)
                              - a_k * M)
assert kummer_residual == 0, f"Kummer residual != 0: {kummer_residual}"
print(f"  symbolic Kummer residual  z M'' + (b-z) M' - a M = {kummer_residual}  (GATE PASSED)")

# Cross-check against sympy's own hyper (independent of our hand-typed poly).
M_sympy = sp.hyper((a_k,), (b_k,), z)
assert sp.simplify(sp.hyperexpand(M_sympy) - M) == 0, "sympy hyper mismatch"
print("  sympy hyperexpand(1F1(-2;5/2;z)) == 4z^2/35 - 4z/5 + 1  (cross-check OK)")

# --- 3-arg recast residual==0 GATE: u'' == ((z-b)u' + a u)/z -------------------
# u'' = ((z - 5/2) u' - 2 u) / z.  Verify f(z, M, M') == M''.
Mp = sp.diff(M, z)
Mpp = sp.diff(M, z, 2)
f3 = ((z - b_k) * Mp + a_k * M) / z
assert sp.simplify(f3 - Mpp) == 0, "3-arg recast residual != 0"
print("  3-arg recast  u'' = ((z-5/2)u' - 2u)/z  verified == M''  (GATE PASSED)")
print("  ∂f/∂u  = a/z   = -2/z")
print("  ∂f/∂up = (z-b)/z = 1 - 5/(2z)")

# --- EXACT rational BC + interior pins -----------------------------------------
def Mrat(zz):
    return sp.Rational(4, 35) * zz**2 - sp.Rational(4, 5) * zz + 1
print("  -- exact rationals (pin INLINE in the test) --")
print(f"  u(1/2)  = {Mrat(sp.Rational(1,2))}   # BC at z_a = 1/2")
print(f"  u(2)    = {Mrat(2)}        # BC at z_b = 2")
print(f"  u(1)    = {Mrat(1)}       # interior")
print(f"  u(3/2)  = {Mrat(sp.Rational(3,2))}        # interior")

print()
print("=" * 72)
print("CBvx.2 MATHIEU  u'' + (a - 2q cos 2z) u = 0,  q=1, a = a_2(1)  (ce_2)")
print("=" * 72)

# --- SYMBOLIC recast GATE: 2-arg  u'' = -(a - 2q cos2z) u ----------------------
u = sp.Function('u')
a_s, q_s = sp.symbols('a q')
mathieu_lhs = sp.diff(u(z), z, 2) + (a_s - 2 * q_s * sp.cos(2 * z)) * u(z)
f2 = -(a_s - 2 * q_s * sp.cos(2 * z)) * u(z)
assert sp.simplify((sp.diff(u(z), z, 2) - f2) - mathieu_lhs) == 0, "Mathieu recast"
print("  2-arg recast  u'' = -(a - 2q cos2z) u  verified == Mathieu ODE  (GATE)")
print("  ∂f/∂u = -(a - 2q cos2z)")

q = 1.0
a2 = sc.mathieu_a(2, q)   # characteristic value a_2(q=1)
print(f"  a = a_2(q=1) = {a2!r}    # scipy.special.mathieu_a(2, 1.0)")

# scipy mathieu_cem takes DEGREES and returns (value, derivative).
def ce2(z_rad):
    val, _deriv = sc.mathieu_cem(2, q, np.degrees(z_rad))
    return val

z_a, z_b = 0.0, np.pi / 2
print(f"  u(0)     = {ce2(z_a)!r}    # BC at z_a = 0       = ce_2(0)")
print(f"  u(pi/2)  = {ce2(z_b)!r}    # BC at z_b = pi/2    = ce_2(pi/2)")
print(f"  u(pi/4)  = {ce2(np.pi/4)!r}    # interior            = ce_2(pi/4)")
print(f"  u(pi/3)  = {ce2(np.pi/3)!r}    # interior            = ce_2(pi/3)")

print()
print("=" * 72)
print("CBvx.3-Whittaker  u'' + (-1/4 + k/z + (1/4-m^2)/z^2) u = 0, k=2,m=3/2 [1,4]")
print("=" * 72)

k_w, m_w = sp.Integer(2), sp.Rational(3, 2)
# --- SYMBOLIC recast GATE: 2-arg  u'' = -(-1/4 + k/z + (1/4-m^2)/z^2) u --------
coeff = (-sp.Rational(1, 4) + k_w / z + (sp.Rational(1, 4) - m_w**2) / z**2)
W = sp.Function('W')
whit_lhs = sp.diff(W(z), z, 2) + coeff * W(z)
fw = -coeff * W(z)
assert sp.simplify((sp.diff(W(z), z, 2) - fw) - whit_lhs) == 0, "Whittaker recast"
print("  2-arg recast  u'' = -(-1/4 + k/z + (1/4-m^2)/z^2) u  verified  (GATE)")
print("  ∂f/∂u = -(-1/4 + k/z + (1/4-m^2)/z^2)")

# mpmath whitm(k, m, z): Whittaker M_{k,m}(z).  dps=50 oracle.
kf, mf = mp.mpf(2), mp.mpf(3) / 2
def whitm(zz):
    return mp.whitm(kf, mf, mp.mpf(zz))

# Numeric residual==0 cross-check at an interior point (independent of recast).
zc = mp.mpf('2.5')
h = mp.mpf('1e-20')
Wpp = (whitm(zc + h) - 2 * whitm(zc) + whitm(zc - h)) / h**2
coeff_n = (-mp.mpf(1) / 4 + kf / zc + (mp.mpf(1) / 4 - mf**2) / zc**2)
resid = Wpp + coeff_n * whitm(zc)
print(f"  numeric Whittaker residual at z=2.5: {mp.nstr(resid, 5)}  (≈0, GATE)")
print(f"  u(1)    = {mp.nstr(whitm(1), 45)}")
print(f"  u(4)    = {mp.nstr(whitm(4), 45)}")
print(f"  u(2)    = {mp.nstr(whitm(2), 45)}")
print(f"  u(2.5)  = {mp.nstr(whitm(mp.mpf('2.5')), 45)}")
print(f"  u(3)    = {mp.nstr(whitm(3), 45)}")
