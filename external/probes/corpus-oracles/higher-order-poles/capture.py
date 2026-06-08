#!/usr/bin/env python3
# capture.py -- sympy/mpmath oracle for the HIGHER-ORDER-POLE corpus
# (bead padetaylor-gf8i, epic padetaylor-25og; plan Family B,
#  docs/test_corpus/02_corpus_extension_plan.md:177-200).
#
# Emits, for the three Family-B cases, the INDEPENDENT ground truth the
# test pins against.  Every value is gated on an exact symbolic identity
# (sympy residual == 0) or a high-precision Laurent/invariant check
# (mpmath) BEFORE it is printed -- never copied from the package.
#
#   CHop.1  order-k fixed pole   u'' = k(k+1) u / (1-z)^2,   u = (1-z)^-k.
#           Residual must be EXACTLY 0 (sympy).  Past-pole values are
#           EXACT RATIONALS: u(1.05) = (-0.05)^-k, u(0.5) = 2^k, and
#           u'(z) = k (1-z)^-(k+1).
#   CHop.3  mixed-order shared pole vector  f = [2 y2, 3 y1^2], IC (1,1),
#           y = ((1-z)^-2, (1-z)^-3).  Residual must be EXACTLY 0 (sympy,
#           both components).  Past-pole values are EXACT RATIONALS.
#   CHop.2  P-prime triple pole   f = [y2, 6 y1^2 - g2/2], (y1,y2)=(P,P').
#           LEMNISCATIC lattice g2=1, g3=0 (real e_i => dps=50-clean).
#           P via the validated Jacobi-sn recipe; P' = mp.diff(P).  Gated
#           on (P')^2 == 4 P^3 - g2 P - g3 and on the triple-pole Laurent
#           leading coefficient (-2/w^3) before any value is emitted.
#
# Run:  python3 capture.py    (sympy 1.14, mpmath 1.3)

from fractions import Fraction
import sympy as sp
import mpmath as mp


def jl_rational(fr):
    """Render a Python Fraction as an exact Julia rational literal."""
    return f"{fr.numerator}//{fr.denominator}"


print("# ====================================================================")
print("# CHop.1  order-k fixed pole:  u'' = k(k+1) u/(1-z)^2,  u=(1-z)^-k")
print("# ====================================================================")
z = sp.symbols('z')
for k in (2, 3, 4, 5):
    u = (1 - z) ** (-k)
    upp = sp.diff(u, z, 2)
    rhs = sp.Integer(k * (k + 1)) * u / (1 - z) ** 2
    resid = sp.simplify(upp - rhs)
    assert resid == 0, f"CHop.1 k={k} recast residual != 0: {resid}"
    # IC at z0=0: u(0)=1, u'(0)=k.
    assert sp.nsimplify(u.subs(z, 0)) == 1
    assert sp.nsimplify(sp.diff(u, z).subs(z, 0)) == k
    # Exact past-pole rationals.
    u_1p05 = Fraction(-5, 100) ** (-k)          # (-0.05)^-k
    u_0p5 = Fraction(1, 2) ** (-k)              # 2^k
    up_1p05 = k * Fraction(-5, 100) ** (-(k + 1))  # u'(z)=k(1-z)^-(k+1) at 1.05
    # Cross-check against sympy substitution (exact).
    assert u.subs(z, sp.Rational(105, 100)) == sp.Rational(
        u_1p05.numerator, u_1p05.denominator)
    assert u.subs(z, sp.Rational(1, 2)) == sp.Rational(
        u_0p5.numerator, u_0p5.denominator)
    assert sp.diff(u, z).subs(z, sp.Rational(105, 100)) == sp.Rational(
        up_1p05.numerator, up_1p05.denominator)
    print(f"#  k={k}: residual==0 OK; residue 0 (pure order-{k})")
    print(f"const CHOP1_U_1p05_K{k}  = {jl_rational(u_1p05)}")
    print(f"const CHOP1_U_0p5_K{k}   = {jl_rational(u_0p5)}")
    print(f"const CHOP1_UP_1p05_K{k} = {jl_rational(up_1p05)}")
print()

print("# ====================================================================")
print("# CHop.3  mixed-order shared pole vector:  f=[2 y2, 3 y1^2], IC (1,1)")
print("#         y = ((1-z)^-2, (1-z)^-3),  shared pole z=1, orders 2 & 3")
print("# ====================================================================")
y1 = (1 - z) ** (-2)
y2 = (1 - z) ** (-3)
# Verify the first-order system: y1' = 2 y2, y2' = 3 y1^2.
r1 = sp.simplify(sp.diff(y1, z) - 2 * y2)
r2 = sp.simplify(sp.diff(y2, z) - 3 * y1 ** 2)
assert r1 == 0, f"CHop.3 y1' = 2 y2 residual != 0: {r1}"
assert r2 == 0, f"CHop.3 y2' = 3 y1^2 residual != 0: {r2}"
# IC at z0=0: (1,1).
assert y1.subs(z, 0) == 1 and y2.subs(z, 0) == 1
print("#  residual==0 OK for BOTH components (y1'=2y2, y2'=3y1^2); IC (1,1)")
for zq in (2, sp.Rational(5, 2)):
    v1 = y1.subs(z, zq)
    v2 = y2.subs(z, zq)
    f1 = Fraction(int(sp.numer(v1)), int(sp.denom(v1)))
    f2 = Fraction(int(sp.numer(v2)), int(sp.denom(v2)))
    tag = "2p0" if zq == 2 else "2p5"
    print(f"const CHOP3_Y1_{tag} = {jl_rational(f1)}")
    print(f"const CHOP3_Y2_{tag} = {jl_rational(f2)}")
print()

print("# ====================================================================")
print("# CHop.2  P' triple pole (LEMNISCATIC g2=1, g3=0):  f=[y2, 6 y1^2-g2/2]")
print("#         (y1,y2)=(P, P').  Oracle: P via Jacobi-sn, P' via mp.diff.")
print("# ====================================================================")
mp.mp.dps = 50
g2 = mp.mpf(1)
g3 = mp.mpf(0)
# Roots e_i of 4 t^3 - g2 t - g3 = 4 t^3 - t = t(4t^2-1):  e = {1/2, 0, -1/2}.
# Order e1 > e2 > e3 (DLMF 23.x convention) => e1=1/2, e2=0, e3=-1/2.
e1 = mp.mpf(1) / 2
e2 = mp.mpf(0)
e3 = -mp.mpf(1) / 2
for e in (e1, e2, e3):
    assert abs(4 * e ** 3 - g2 * e - g3) < mp.mpf(10) ** (-45), "e_i not a root"
m = (e2 - e3) / (e1 - e3)        # Jacobi modulus^2 parameter


def wp(zz):
    # P(z; g2,g3) = e3 + (e1-e3)/sn^2( z*sqrt(e1-e3), m ).  DLMF 23.6.
    arg = zz * mp.sqrt(e1 - e3)
    return e3 + (e1 - e3) / mp.ellipfun('sn', arg, m) ** 2


def wpp(zz):
    return mp.diff(wp, zz)


# GATE 1: the elliptic invariant (P')^2 = 4 P^3 - g2 P - g3 at several z.
for zz in (mp.mpf('0.4'), mp.mpf('0.9'), mp.mpf('1.3'),
           mp.mpc('0.7', '0.2')):
    P = wp(zz)
    Pp = wpp(zz)
    inv = Pp ** 2 - (4 * P ** 3 - g2 * P - g3)
    assert abs(inv) < mp.mpf(10) ** (-40), f"invariant fails at z={zt}: {inv}"
# GATE 2: the triple pole.  Real lattice half-period (period 2w with
# w = real quarter related to K).  Nearest real pole is at z = 0 (lattice
# point); P ~ 1/z^2, P' ~ -2/z^3.  Check the Laurent leading coeff of P'
# near a small z: z^3 * P'(z) -> -2.
for eps in ('1e-3', '1e-4'):
    e = mp.mpf(eps)
    lead = e ** 3 * wpp(e)
    assert abs(lead + 2) < mp.mpf('1e-5'), f"P' Laurent lead != -2: {lead}"
print("#  invariant (P')^2=4P^3-g2 P-g3 OK (dps=50); P' Laurent lead -2/w^3 OK")

# Real half-period 2*omega of the lemniscatic lattice: poles at z = 2*omega*k.
# omega = K(m)/sqrt(e1-e3) where K is the complete elliptic integral.  The
# REAL pole nearest a chosen IC point z0 is the bridging target.
Kc = mp.ellipk(m)
omega = Kc / mp.sqrt(e1 - e3)
pole = 2 * omega           # first real pole at z = 2*omega (z=0 is also one)
print(f"const CHOP2_POLE_REAL = {mp.nstr(2 * omega, 40)}")
# IC point z0 and a past-pole evaluation target straddling z=2*omega.
z0 = mp.mpf('0.9')
ztgt = 2 * omega + mp.mpf('0.5')   # past the pole at 2*omega
for label, zz in (('Z0', z0), ('ZTGT', ztgt)):
    print(f"const CHOP2_{label}      = {mp.nstr(zz, 40)}")
    print(f"const CHOP2_P_{label}    = {mp.nstr(wp(zz).real, 40)}")
    print(f"const CHOP2_PP_{label}   = {mp.nstr(wpp(zz).real, 40)}")
