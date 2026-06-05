"""
capture.py — independent mpmath/sympy oracle values for test/corpus_vector_*_test.jl
(epic padetaylor-p1v0, bead padetaylor-jmcf).  Run: python3 capture.py

Each block re-derives the oracle from canonical ground truth (closed form /
mpmath), NOT from the Julia routine under test (CLAUDE.md Rule 5).  The
emitted constants are pinned verbatim in the test files.
"""
import sympy as sp
from mpmath import mp, mpf, mpc, tan as mtan, cos as mcos, sec as msec, pi as mpi

# =====================================================================
# (1) coupled-riccati-shared-denom-canonical  y=(1/(z-z0),-1/(z-z0)^2), z0=1.5
#     Verify the ODE y1'=y2, y2'=2 y2^2/y1 symbolically (sympy), and emit
#     the EXACT Taylor jets at z=0 and the closed-form values past the pole.
# =====================================================================
print("# === (1) coupled-riccati-shared-denom-canonical (z0=3/2) ===")
z = sp.symbols('z')
z0 = sp.Rational(3,2)
y1 = 1/(z-z0); y2 = -1/(z-z0)**2
# ODE checks
r1 = sp.simplify(sp.diff(y1,z) - y2)                 # y1' - y2 == 0
r2 = sp.simplify(sp.diff(y2,z) - 2*y2**2/y1)         # y2' - 2 y2^2/y1 == 0
print("#   ODE residuals (must be 0):", r1, r2)
# Exact closed-form values at evaluation points (past the pole z0=1.5)
for zq in (sp.Rational(1,2), sp.Rational(2), sp.Rational(5,2)):
    v1 = sp.nsimplify(y1.subs(z,zq)); v2 = sp.nsimplify(y2.subs(z,zq))
    print(f"#   z={float(zq)}: y1={v1}={float(v1):.17g}  y2={v2}={float(v2):.17g}")
# Jets at z=0 (catalogue: c1[k]=(-1)^k/z0^(k+1); c2[k]=(-1)^(k+1)(k+1)/z0^(k+2))
print("#   Taylor jet at z=0 (first 6 of each, exact rational):")
s1 = sp.series(y1, z, 0, 6).removeO()
s2 = sp.series(y2, z, 0, 6).removeO()
print("#   y1:", [sp.nsimplify(s1.coeff(z,k)) for k in range(6)])
print("#   y2:", [sp.nsimplify(s2.coeff(z,k)) for k in range(6)])
print()

# =====================================================================
# (2) coupled-riccati-tan-sec2  y'=(y2, 2 y1 y2), y=(tan z, sec^2 z); IC (0,1)
#     Bridge pole at pi/2.  mpmath values past the pole.
# =====================================================================
print("# === (2) coupled-riccati-tan-sec2 (shared pole pi/2) ===")
# sympy ODE check
zt = sp.symbols('zt')
t1 = sp.tan(zt); t2 = sp.sec(zt)**2
print("#   y1'-y2  =", sp.simplify(sp.diff(t1,zt)-t2))
print("#   y2'-2y1y2=", sp.simplify(sp.diff(t2,zt)-2*t1*t2))
mp.dps = 40
def tanval(x): return (mtan(x), msec(x)**2)
for xq in (mpf('2.0'), mpf('2.5')):
    a,b = tanval(xq)
    print(f"#   z={float(xq)}: tan={mp.nstr(a,30)}  sec^2={mp.nstr(b,30)}")
# Float64-ready strings
print("CRTS_TAN_2p0 =", repr(float(mtan(mpf('2.0')))))
print("CRTS_SEC2_2p0 =", repr(float(msec(mpf('2.0'))**2)))
print("CRTS_TAN_2p5 =", repr(float(mtan(mpf('2.5')))))
print("CRTS_SEC2_2p5 =", repr(float(msec(mpf('2.5'))**2)))
print()

# =====================================================================
# (3) tan-sec2-vector-jet-b13  y1'=y2, y2'=2 y1 y2; y=(tan(z+pi/4), sec^2(z+pi/4))
#     IC y1(0)=1, y2(0)=2.  Independent jet via mpmath.diff AND the recurrence.
# =====================================================================
print("# === (3) tan-sec2-vector-jet-b13 (jet at z=0, a=pi/4) ===")
mp.dps = 50
a = mpi/4
order = 30
# (a) recurrence c1[n+1]=c2[n]/(n+1); c2[n+1]=2 sum_k c1[k]c2[n-k]/(n+1)
c1 = [mpf(0)]*(order+1); c2 = [mpf(0)]*(order+1)
c1[0] = mtan(a)            # = 1
c2[0] = 1/mcos(a)**2       # = 2
for n in range(order):
    c1[n+1] = c2[n]/(n+1)
    conv = sum(c1[k]*c2[n-k] for k in range(n+1))
    c2[n+1] = 2*conv/(n+1)
# (b) independent mpmath.diff cross-check of a few coefficients
from mpmath import diff, factorial
f1 = lambda x: mtan(x+a)
f2 = lambda x: 1/mcos(x+a)**2
ok = True
for k in (0,1,2,3,5,10,29,30):
    d1 = diff(f1, 0, k)/factorial(k)
    d2 = diff(f2, 0, k)/factorial(k)
    e1 = abs(d1-c1[k]); e2 = abs(d2-c2[k])
    if e1>mpf('1e-40') or e2>mpf('1e-40'): ok=False
    print(f"#   k={k:2d}: c1={mp.nstr(c1[k],20)} (diff err {float(e1):.1e})  "
          f"c2={mp.nstr(c2[k],20)} (diff err {float(e2):.1e})")
print("#   recurrence-vs-diff agree to 1e-40:", ok)
# Emit the full 50-dps coefficient strings for the test
print("CTSJ_C1 = [")
for k in range(order+1): print("  ", repr(mp.nstr(c1[k], 48)), ",", sep="")
print("]")
print("CTSJ_C2 = [")
for k in range(order+1): print("  ", repr(mp.nstr(c2[k], 48)), ",", sep="")
print("]")
print()

# =====================================================================
# (4) pii-u1-companion-2vector  PII(alpha=1): u=-1/t, companion y=(-1/t, 1/t^2)
#     IC y(1)=(-1,1). ODE: y1'=y2, y2'=2 y1^3 + t y1 + 1.  Exact rational.
# =====================================================================
print("# === (4) pii-u1-companion-2vector (u=-1/t, alpha=1) ===")
tt = sp.symbols('tt')
u1 = -1/tt; u2 = 1/tt**2
print("#   y1'-y2     =", sp.simplify(sp.diff(u1,tt)-u2))
print("#   y2'-(2y1^3+t y1+1) =", sp.simplify(sp.diff(u2,tt) - (2*u1**3+tt*u1+1)))
for tq in (sp.Rational(1,2), sp.Rational(3,2), sp.Rational(2)):
    print(f"#   t={float(tq)}: y1={float((-1/tq)):.17g}  y2={float((1/tq**2)):.17g}")
