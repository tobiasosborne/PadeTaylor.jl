#!/usr/bin/env python3
# =============================================================================
# Oracle capture — CBvx.4 oblique-complex BVP/VectorBVP domain-guard (53tu).
#
# GROUND TRUTH (independent of the Julia routine under test — CLAUDE.md Rule 5):
#   ODE  u'' = u,  exact  u = cosh(z),  on the OBLIQUE segment z_a=0, z_b=1+i
#   (both Re and Im of z_b - z_a nonzero).  Vector mirror y'=(y2,y1) with the
#   closed form y = (cosh z, sinh z): y1'=y2=sinh=cosh' ✓, y2'=y1=cosh=sinh' ✓.
#
# THE t* MAPPING (read verbatim from src/BVP.jl:489 / src/VectorBVP.jl:316):
#   half_sum  = (z_a + z_b)/2 = 0.5 + 0.5i
#   half_diff = (z_b - z_a)/2 = 0.5 + 0.5i
#   t*(z*) = (z* - half_sum) / half_diff
#   z*(t*) = half_sum + half_diff * t*
#
# POINTS PINNED:
#   - BC endpoints:  u(z_a)=cosh(0)=1,  u(z_b)=cosh(1+i).
#   - on-segment midpoint z_mid = 0.5+0.5i  (t* = 0, real ⇒ guard correct).
#   - off-segment z_off with t* = 0.5 + 3i: INTERIOR Re(t*)=0.5 in (-1,1) but
#     large Im(t*)=3.  The real(t*)-only guard PASSES it, the barycentric
#     interpolant silently EXTRAPOLATES.  This is the INTERIOR-Re complement
#     to CBV.7's PERPENDICULAR t*=10i (Re=0) case — a genuinely distinct,
#     more insidious trigger for bug padetaylor-53tu.
#       z_off = half_sum + half_diff*(0.5+3i)
#             = (0.5+0.5i) + (0.5+0.5i)(0.5+3i) = -0.75 + 2.25i  (exact).
#
# RESIDUAL GATE: sympy verifies u''-u == 0 for u=cosh(z) symbolically, and the
# vector system residual y'-(y2,y1) == 0 for (cosh,sinh) — so the oracle ODE
# is certified before any numerical value is emitted.
# =============================================================================
import mpmath as mp
import sympy as sp

mp.mp.dps = 30

# ---- sympy residual gate (== 0 exactly) -------------------------------------
z = sp.symbols('z')
u = sp.cosh(z)
assert sp.simplify(sp.diff(u, z, 2) - u) == 0, "scalar ODE residual u''-u != 0"
y1, y2 = sp.cosh(z), sp.sinh(z)
assert sp.simplify(sp.diff(y1, z) - y2) == 0, "vector residual y1'-y2 != 0"
assert sp.simplify(sp.diff(y2, z) - y1) == 0, "vector residual y2'-y1 != 0"

# ---- the segment + the t* mapping -------------------------------------------
z_a = mp.mpc(0, 0)
z_b = mp.mpc(1, 1)
half_sum  = (z_a + z_b) / 2          # 0.5 + 0.5i
half_diff = (z_b - z_a) / 2          # 0.5 + 0.5i

def zof(tstar):
    return half_sum + half_diff * tstar

def tof(zstar):
    return (zstar - half_sum) / half_diff

z_mid = zof(mp.mpc(0, 0))            # t* = 0      (on-segment NODE for even N)
z_qrt = zof(mp.mpc("0.3", 0))        # t* = 0.3    (on-segment OFF-NODE: exercises
                                     #              the barycentric interpolant)
t_off = mp.mpc("0.5", "3")          # interior Re, large Im
z_off = zof(t_off)                  # = -0.75 + 2.25i

# round-trip check: t*(z_off) reproduces 0.5+3i
assert abs(tof(z_off) - t_off) < mp.mpf("1e-25"), "t* round-trip failed"
assert abs(z_off - mp.mpc("-0.75", "2.25")) < mp.mpf("1e-25"), "z_off mismatch"

def cx(name, val):
    print(f"const {name} = {mp.nstr(val.real, 25)} + {mp.nstr(val.imag, 25)}im")

print("# --- CBvx.4 oblique-guard oracle (mpmath dps=30; sympy residual==0) ---")
print(f"# z_a=0, z_b=1+i; t*(z*)=(z*-({mp.nstr(half_sum.real,3)}+{mp.nstr(half_sum.imag,3)}i))/({mp.nstr(half_diff.real,3)}+{mp.nstr(half_diff.imag,3)}i)")
print(f"# z_mid={mp.nstr(z_mid,6)} (t*=0); z_off={mp.nstr(z_off,6)} (t*=0.5+3i)")
cx("OG_BC_A",   mp.cosh(z_a))        # u(z_a)=cosh(0)=1
cx("OG_BC_B",   mp.cosh(z_b))        # u(z_b)=cosh(1+i)
cx("OG_COSH_MID", mp.cosh(z_mid))    # on-segment NODE value (t*=0)
cx("OG_COSH_QRT", mp.cosh(z_qrt))    # on-segment OFF-NODE value (t*=0.3)
cx("OG_COSH_OFF", mp.cosh(z_off))    # true cosh at z_off (vs the extrapolant)
cx("OG_SINH_A",   mp.sinh(z_a))      # y2(z_a)=sinh(0)=0 (vector BC)
cx("OG_SINH_B",   mp.sinh(z_b))      # sinh(1+i)
cx("OG_SINH_MID", mp.sinh(z_mid))    # vector on-segment 2nd comp (node)
cx("OG_SINH_QRT", mp.sinh(z_qrt))    # vector on-segment OFF-NODE 2nd comp
cx("OG_SINH_OFF", mp.sinh(z_off))    # vector true value at z_off
print(f"# z_qrt literal: {mp.nstr(z_qrt.real,6)} + {mp.nstr(z_qrt.imag,6)}im  (t* = 0.3, off-node)")
print("# z_off literal:  -0.75 + 2.25im   (t* = 0.5 + 3.0im)")
