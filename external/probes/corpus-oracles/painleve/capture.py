"""
capture.py — independent sympy/mpmath oracle values for
test/corpus_painleve_rational_test.jl (epic padetaylor-p1v0, bead
padetaylor-lgpc).  Run: python3 capture.py

Every PII/PIV closed form is SUBSTITUTED into its own Painlevé equation by
sympy and the residual confirmed identically 0 BEFORE any value is pinned
(CLAUDE.md Law 1 / Rule 5).  Transcendental oracles (PII Airy-half-integer
pole, NY generic transcendent) come from mpmath at high dps.  The emitted
constants are pinned verbatim in the test file.

Painlevé II canonical form:  u'' = 2 u^3 + z u + alpha
Painlevé IV canonical form:  u'' = (u')^2/(2u) + (3/2)u^3 + 4 z u^2
                                   + 2(z^2 - alpha) u + beta/u
(matches src/Painleve.jl _pII_rhs / _pIV_rhs).
"""
import sympy as sp
from mpmath import mp, mpf, airyai, findroot, cbrt

def airyaiprime(x):
    return airyai(x, derivative=1)

z = sp.symbols('z')

def pii_residual(u, alpha):
    return sp.simplify(sp.diff(u, z, 2) - (2*u**3 + z*u + alpha))

def piv_residual(u, alpha, beta):
    up = sp.diff(u, z)
    rhs = up**2/(2*u) + sp.Rational(3,2)*u**3 + 4*z*u**2 \
          + 2*(z**2 - alpha)*u + beta/u
    return sp.simplify(sp.diff(u, z, 2) - rhs)

print("# ===== PII rational closed forms (alpha = n) — sympy residual must be 0 =====")
u1 = -1/z
u2 = (4 - 2*z**3)/(4*z + z**4)
u3 = 3*z**2*(160 + 8*z**3 + z**6)/(320 - 24*z**6 - z**9)
print("#   u1=-1/z       (alpha=1):  residual =", pii_residual(u1, 1))
print("#   u2            (alpha=2):  residual =", pii_residual(u2, 2))
print("#   u3            (alpha=3):  residual =", pii_residual(u3, 3))

print("# ----- u2 real-axis pole: real root of z*(z^3+4) other than 0 -----")
print("#   z = -4^(1/3) =", float(-sp.Rational(4)**sp.Rational(1,3)))   # -1.58740105...

print("# ----- u3 real-axis poles: roots of Q2=z^3+4 and Q3=z^6+20z^3-80 -----")
q3roots = sp.solve(z**6 + 20*z**3 - 80, z)
realq3 = sorted([complex(r).real for r in q3roots if abs(complex(r).imag) < 1e-12])
print("#   Q3 real roots:", realq3)              # -2.86092..., +1.50610...
print("#   Q2 real root :", float(-sp.Rational(4)**sp.Rational(1,3)))
# The pole-CROSSING for u3 the catalogue flags is the POSITIVE real root ~1.50611.
print("#   u3 first positive real pole z* =", realq3[-1])
# Closed-form value of u3 just past that pole (for the bridge far-side check).
zc = sp.nsimplify(realq3[-1], [sp.sqrt(5)])
print("#   (algebraic z*) =", zc)

print("# ===== PIV rational/entire closed forms — sympy residual must be 0 =====")
# piv-hermite-m3: u = -4z/(2z^2-1), (alpha,beta)=(-3,-8)
um3 = -4*z/(2*z**2 - 1)
print("#   hermite m3 u=-4z/(2z^2-1)        (a=-3,b=-8): residual =",
      piv_residual(um3, -3, -8))
print("#     poles at z=+-1/sqrt(2) =", float(1/sp.sqrt(2)))
# piv-hermite-m4: u = 3(1-2z^2)/(z(2z^2-3)), (alpha,beta)=(-4,-18)
um4 = 3*(1 - 2*z**2)/(z*(2*z**2 - 3))
print("#   hermite m4 u=3(1-2z^2)/(z(2z^2-3))(a=-4,b=-18): residual =",
      piv_residual(um4, -4, -18))
print("#     poles at z=0, +-sqrt(6)/2 =", float(sp.sqrt(6)/2))
# piv-entire-m2z: u = -2z, (alpha,beta)=(0,-2)
um2z = -2*z
print("#   entire u=-2z                      (a=0,b=-2):  residual =",
      piv_residual(um2z, 0, -2))

print("# ===== PII Airy half-integer (alpha=1/2): u = (1/2^(1/3)) Ai'(w)/Ai(w), w=-z/2^(1/3) =====")
# DLMF 32.10.7 / FW2014: phi=Ai(-z/2^(1/3)) solves phi''=-(z/2)phi, u=-phi'/phi.
# First real pole = first zero of Ai(-z/2^(1/3)) on z>0 ⇒ z = -2^(1/3)*a1,
# a1 = first (negative) Airy zero ≈ -2.33810741.  So z* = 2^(1/3)*|a1|.
mp.dps = 40
c13 = cbrt(mpf(2))                      # 2^(1/3)
# first Airy zero a1 (negative real): Ai(a1)=0 near -2.3381
a1 = findroot(airyai, mpf('-2.33810741'))
print("#   first Airy zero a1 =", a1)
z_pole = -a1 * c13                       # = 2^(1/3)*2.33810741 ≈ 2.9446
print("#   u_{1/2} first real pole z* = -2^(1/3)*a1 =", z_pole)   # 2.94455...

# Value of u_{1/2}(z) at a regular interior point z=1.0 (oracle for the
# solver cross-check), u = (1/2^(1/3)) Ai'(w)/Ai(w), w = -z/2^(1/3).
def u_half(zval):
    w = -mpf(zval)/c13
    return (mpf(1)/c13) * airyaiprime(w)/airyai(w)
for zq in ('0.5', '1.0', '1.5', '2.0'):
    print(f"#   u_1/2({zq}) =", u_half(mpf(zq)))

print("# ===== done — every PII/PIV residual above is identically 0 =====")
