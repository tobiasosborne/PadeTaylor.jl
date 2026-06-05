"""
Independent verification of the cm-n2-imaginary-collision closed form
(catalogue confidence = MEDIUM, "not yet independently verified").

Catalogue claims (01_corpus_catalogue.md:574):
  ODE: x1'' = 4/(x1-x2)^3, x2'' = 4/(x2-x1)^3  (repulsive rational CM, EOM const 4)
  IC : x1(0)=+i*r0, x2(0)=-i*r0, v1(0)=v2(0)=0, r0=1
  closed form: x1(t) = (i/2)*sqrt(4 r0^2 - 2 t^2 / r0^2),  x2 = -x1
               collision (shared pole) at t* = r0^2*sqrt(2) ~ 1.414

Two independent checks (Law 1):
  (a) SYMBOLIC substitution in sympy: does x1(t)=(I/2)*sqrt(4-2t^2) satisfy
      x1'' = 4/(x1-x2)^3 with x2=-x1 ?
  (b) HIGH-DPS RK4 in mpmath of the first-order CM system from the imaginary IC,
      compared to the closed form at several real t in [0, t*).
"""
import sympy as sp

print("="*72)
print("(a) SYMBOLIC sympy check of cm-n2-imaginary-collision")
print("="*72)
t, r0 = sp.symbols('t r0', positive=True)
# Catalogue closed form (r0 generic):
x1 = (sp.I/2) * sp.sqrt(4*r0**2 - 2*t**2/r0**2)
x2 = -x1
# x1 - x2 = 2 x1
diff = sp.simplify(x1 - x2)
# LHS: x1''
x1pp = sp.simplify(sp.diff(x1, t, 2))
# RHS: 4/(x1-x2)^3
rhs1 = sp.simplify(4/diff**3)
resid1 = sp.simplify(x1pp - rhs1)
print("x1            =", x1)
print("x1''          =", x1pp)
print("4/(x1-x2)^3   =", rhs1)
print("x1'' - RHS    =", resid1, " (must be 0)")
# x2'' = 4/(x2-x1)^3
x2pp = sp.simplify(sp.diff(x2, t, 2))
rhs2 = sp.simplify(4/(x2-x1)**3)
resid2 = sp.simplify(x2pp - rhs2)
print("x2'' - RHS2   =", resid2, " (must be 0)")
# Collision: x1 = x2  <=>  x1 = 0  <=>  4r0^2 - 2 t^2/r0^2 = 0  =>  t = r0^2 sqrt(2)
tstar = sp.solve(sp.Eq(4*r0**2 - 2*t**2/r0**2, 0), t)
print("collision t*  =", tstar, " (catalogue: r0^2*sqrt(2))")
# IC check
print("x1(0)         =", sp.simplify(x1.subs(t,0)), " (catalogue: I*r0)")
v1 = sp.diff(x1,t)
print("v1(0)         =", sp.simplify(v1.subs(t,0)), " (catalogue: 0)")
SYM_OK = (resid1 == 0 and resid2 == 0)
print(">>> SYMBOLIC VERIFIED:", SYM_OK)

print()
print("="*72)
print("(b) HIGH-DPS mpmath RK4 of the first-order CM system from imaginary IC")
print("="*72)
from mpmath import mp, mpf, mpc, sqrt as msqrt, j as mj
mp.dps = 50

def cm_rhs(y):
    # y = (x1, x2, v1, v2); first-order; EOM const 4
    x1_, x2_, v1_, v2_ = y
    d = x1_ - x2_
    a1 = 4/d**3
    return [v1_, v2_, a1, -a1]

def rk4(y0, t1, n):
    h = t1/n
    y = list(y0)
    for _ in range(n):
        k1 = cm_rhs(y)
        k2 = cm_rhs([y[i]+(h/2)*k1[i] for i in range(4)])
        k3 = cm_rhs([y[i]+(h/2)*k2[i] for i in range(4)])
        k4 = cm_rhs([y[i]+h*k3[i] for i in range(4)])
        y = [y[i]+(h/6)*(k1[i]+2*k2[i]+2*k3[i]+k4[i]) for i in range(4)]
    return y

def closed_x1(tt):       # r0 = 1
    return (mj/2)*msqrt(4 - 2*tt**2)
def closed_v1(tt):
    # d/dt (I/2 sqrt(4-2t^2)) = (I/2)*(-2t)/sqrt(4-2t^2) = -I t / sqrt(4-2t^2)
    return -mj*tt/msqrt(4-2*tt**2)

y0 = [mj, -mj, mpf(0), mpf(0)]   # x1=+i, x2=-i, v=0
RK_OK = True
for tt in (mpf('0.3'), mpf('0.7'), mpf('1.0'), mpf('1.2')):
    yr = rk4(y0, tt, 200000)
    x1e = closed_x1(tt); v1e = closed_v1(tt)
    ex1 = abs(yr[0]-x1e); ev1 = abs(yr[2]-v1e)
    print(f"t={float(tt):.2f}  RK4 x1={complex(yr[0])}  closed={complex(x1e)}  |err|={float(ex1):.2e}")
    print(f"          RK4 v1={complex(yr[2])}  closed={complex(v1e)}  |err|={float(ev1):.2e}")
    if ex1 > mpf('1e-15') or ev1 > mpf('1e-15'):
        RK_OK = False
print(">>> RK4 VERIFIED (all |err|<1e-15):", RK_OK)

print()
print("="*72)
print("OVERALL cm-n2-imaginary-collision VERIFIED:", SYM_OK and RK_OK)
print("="*72)
# Emit pinned reference values at the test points (r0=1, dps=40)
mp.dps = 40
for tt in (mpf('0.5'), mpf('0.9'), mpf('1.2')):
    print(f"x1({float(tt)}) = {mp.nstr(closed_x1(tt).imag, 36)} * I    "
          f"v1 = {mp.nstr(closed_v1(tt).imag, 36)} * I")
