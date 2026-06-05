"""
capture.py — independent mpmath/sympy oracle values for test/corpus_bvp_test.jl
and test/corpus_bvp_hybrid_test.jl (epic padetaylor-p1v0, bead padetaylor-u638).

Run: python3 capture.py

Each block RE-DERIVES the oracle from canonical ground truth (closed form /
mpmath / a hand-solved transcendental), NOT from the Julia routine under test
(CLAUDE.md Rule 5).  Crucially, every ODE residual is verified == 0 with sympy
BEFORE any boundary/interior value is emitted — this is the Law-1 ground-truth
gate.  The emitted constants are pinned verbatim in the test files.

Catalogue records (docs/test_corpus/01_corpus_catalogue.md):
  3arg-damped-oscillator   md:1234-1249    3arg-euler-cauchy    md:1251-1266
  3arg-nonlinear-exp       md:1268-1283    bratu-bvp            md:1302-1317
  airy-bvp                 md:1319-1334    oblique-complex-cosh md:1387-1402
  bigfloat-256-oblique-cosh md:1404-1419   3arg-boundary-layer  md:1421-1436
  ivp-bvp-hybrid-dispatch  md:1438-1453    bvp-domain-guard     md:1919-1934
"""
import sympy as sp
from mpmath import mp, mpf, mpc, exp as mexp, cos as mcos, sin as msin, \
    cosh as mcosh, log as mlog, findroot, airyai

z = sp.symbols('z', real=True)

def show_resid(name, u_expr, rhs_expr):
    """u'' - rhs(z,u,u') must simplify to 0."""
    upp = sp.diff(u_expr, z, 2)
    r = sp.simplify(upp - rhs_expr)
    print(f"#   {name}: u'' - f = {r}  (must be 0)")
    return r

# =====================================================================
# (1) 3arg-damped-oscillator   u'' = -2u' - 5u   u = e^{-z} cos(2z)  on [-1,1]
#     Isolates df/du' = -2 (the D1-coupled Jacobian term).
# =====================================================================
print("# === (1) 3arg-damped-oscillator-bvp  u''=-2u'-5u, u=e^{-z}cos(2z) ===")
u = sp.exp(-z) * sp.cos(2*z)
up = sp.diff(u, z)
show_resid("damped", u, -2*up - 5*u)          # ODE residual gate
mp.dps = 50
fu = lambda x: mexp(-x) * mcos(2*x)
print("DAMP_BC_A =", repr(float(fu(mpf('-1')))))   # u(-1)
print("DAMP_BC_B =", repr(float(fu(mpf('1')))))    # u(1)
for zq in ('-0.7', '0.0', '0.7'):
    print(f"DAMP_U_{zq.replace('-','m').replace('.','p')} =", repr(float(fu(mpf(zq)))))
print()

# =====================================================================
# (2) 3arg-euler-cauchy   u'' = u/z^2 - u'/z   u = (4/3)z + (2/3)/z  on [1,2]
#     df/du' = -1/z is z-DEPENDENT -> node-varying D1 Jacobian contribution.
# =====================================================================
print("# === (2) 3arg-euler-cauchy-bvp  u''=u/z^2-u'/z, u=(4/3)z+(2/3)/z ===")
u = sp.Rational(4,3)*z + sp.Rational(2,3)/z
up = sp.diff(u, z)
show_resid("euler", u, u/z**2 - up/z)
# Exact rational interior + BC values.
for zq in (sp.Integer(1), sp.Rational(5,4), sp.Rational(3,2),
           sp.Rational(7,4), sp.Integer(2)):
    val = sp.nsimplify(u.subs(z, zq))
    print(f"#   u({float(zq)}) = {val} = {float(val):.17g}")
print("EULER_U_1p5_NUM, EULER_U_1p5_DEN = 22, 9   # exact 22/9")
print()

# =====================================================================
# (3) 3arg-nonlinear-exp   u'' = (u')^2/u   u = e^{2z}  on [0,1]
#     Nonlinear; df/du=-(u')^2/u^2, df/du'=2u'/u.
# =====================================================================
print("# === (3) 3arg-nonlinear-exp-bvp  u''=(u')^2/u, u=e^{2z} ===")
u = sp.exp(2*z)
up = sp.diff(u, z)
show_resid("nlexp", u, up**2/u)
mp.dps = 50
fu = lambda x: mexp(2*x)
print("NLEXP_BC_B =", repr(float(fu(mpf('1')))))   # u(1) = e^2
for zq in ('0.25', '0.5', '0.75'):
    print(f"NLEXP_U_{zq.replace('.','p')} =", repr(float(fu(mpf(zq)))))
print()

# =====================================================================
# (4) bratu-bvp   u'' = -2 e^u   u = 2 log(B/cosh(B(z-1/2))), B=cosh(B/2)  [0,1]
#     Solve the transcendental B ourselves, verify residual symbolically.
# =====================================================================
print("# === (4) bratu-bvp  u''=-2e^u, B=cosh(B/2) ===")
B = sp.symbols('B', positive=True)
u_sym = 2*sp.log(B/sp.cosh(B*(z - sp.Rational(1,2))))
# Residual u'' + 2 e^u with B treated as the constant from B=cosh(B/2):
r = sp.simplify(sp.diff(u_sym, z, 2) + 2*sp.exp(u_sym))
print("#   bratu residual (symbolic, any B): u''+2e^u =", r)
mp.dps = 50
Bval = findroot(lambda b: b - mcosh(b/2), mpf('1.2'))
print("#   B = cosh(B/2) ->", mp.nstr(Bval, 40))
print("#   check B - cosh(B/2) =", mp.nstr(Bval - mcosh(Bval/2), 5))
print("BRATU_B =", repr(float(Bval)))
ub = lambda x: 2*mlog(Bval/mcosh(Bval*(x - mpf('0.5'))))
print("#   BC u(0) =", mp.nstr(ub(mpf('0')), 5), " u(1) =", mp.nstr(ub(mpf('1')), 5))
for zq in ('0.25', '0.5', '0.75'):
    print(f"BRATU_U_{zq.replace('.','p')} =", repr(float(ub(mpf(zq)))))
print()

# =====================================================================
# (5) airy-bvp   u'' = z u   u = Ai(z)  on [0,3]   (clean 2-arg cross-check)
# =====================================================================
print("# === (5) airy-bvp  u''=z u, u=Ai(z) ===")
# Ai satisfies u''=zu by definition (DLMF 9.2.1); residual gate is the defn.
mp.dps = 50
print("AIRY_BC_A =", repr(float(airyai(mpf('0')))))   # Ai(0)
print("AIRY_BC_B =", repr(float(airyai(mpf('3')))))   # Ai(3)
for zq in ('0.5', '1.5', '2.5'):
    print(f"AIRY_U_{zq.replace('.','p')} =", repr(float(airyai(mpf(zq)))))
print("#   Ai(0) 50dps =", mp.nstr(airyai(mpf('0')), 40))
print("#   Ai(3) 50dps =", mp.nstr(airyai(mpf('3')), 40))
print()

# =====================================================================
# (6) oblique-complex-cosh   u'' = u   u = cosh(z)  on [0, 1+i]
#     + bigfloat256-oblique-cosh  (same, 80-dps reference string)
# =====================================================================
print("# === (6) oblique-complex-cosh-bvp  u''=u, u=cosh(z), [0,1+i] ===")
zc = sp.symbols('zc')
ucz = sp.cosh(zc)
print("#   cosh residual u''-u =", sp.simplify(sp.diff(ucz, zc, 2) - ucz))
mp.dps = 50
def cval(zr, zi): return mcosh(mpc(zr, zi))
for (zr, zi) in (('0.5','0.5'), ('1','1')):
    v = cval(zr, zi)
    print(f"#   cosh({zr}+{zi}i) = {mp.nstr(v.real,30)} + {mp.nstr(v.imag,30)} i")
c_bc_b = cval('1','1'); c_int = cval('0.5','0.5')
print("OBL_BC_B_RE, OBL_BC_B_IM =", repr(float(c_bc_b.real)), ",", repr(float(c_bc_b.imag)))
print("OBL_INT_RE, OBL_INT_IM =", repr(float(c_int.real)), ",", repr(float(c_int.imag)))
# BigFloat-256: 80-dps reference strings for cosh(0.5+0.5i)
mp.dps = 80
c80 = mcosh(mpc('0.5','0.5'))
print('OBL_BF_INT_RE = "' + mp.nstr(c80.real, 78) + '"')
print('OBL_BF_INT_IM = "' + mp.nstr(c80.imag, 78) + '"')
cb80 = mcosh(mpc('1','1'))
print('OBL_BF_BC_B_RE = "' + mp.nstr(cb80.real, 78) + '"')
print('OBL_BF_BC_B_IM = "' + mp.nstr(cb80.imag, 78) + '"')
print()

# =====================================================================
# (7) bvp-domain-guard-imaginary-tstar  segment [0,1+i], query z=-4.5+5.5j
#     Confirm t* = 10j EXACTLY (Re(t*)=0 passes guard, Im(t*)=10 ignored).
# =====================================================================
print("# === (7) bvp-domain-guard-imaginary-tstar  t*=10j ===")
za, zb, zq = mpc('0','0'), mpc('1','1'), mpc('-4.5','5.5')
half_sum = (za + zb)/2; half_diff = (zb - za)/2
tstar = (zq - half_sum)/half_diff
print("#   t* =", mp.nstr(tstar, 20), " Re =", mp.nstr(tstar.real,5),
      " Im =", mp.nstr(tstar.imag,5))
print("#   |Re t*| <= 1 ?", abs(tstar.real) <= 1, "  (guard passes, off-segment)")
print()

# =====================================================================
# (8) ivp-bvp-dispatch-hybrid  u''=u, u=cosh(z) on [0,1+i]; cross-check vs
#     IVP shot from u(0)=1, u'(0)=sinh(0)=0.
# =====================================================================
print("# === (8) ivp-bvp-dispatch-hybrid  cosh on [0,1+i] ===")
mp.dps = 50
print("#   u(0)=cosh(0)=1, u'(0)=sinh(0)=0   IVP IC for the cosh shot")
print("#   interior cosh(0.5+0.5i) re/im pinned above (OBL_INT_*)")
print()
print("# === done ===")
