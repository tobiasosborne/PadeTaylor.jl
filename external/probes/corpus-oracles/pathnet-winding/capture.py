#!/usr/bin/env python3
"""capture.py -- sympy/mpmath oracle for corpus_pathnet_winding_test.jl (CPN.7).

Family I (docs/test_corpus/02_corpus_extension_plan.md:374-381): the HEADLINE
padetaylor-61um catch under a REAL `path_network_solve(cross_branch=true)`
walk -- the FIRST time the walker-side sheet bookkeeping (BranchTracker
step_sheet_update via SheetTracker winding_delta) is exercised by an actual
wedge step that subtends |dtheta| >= pi about a tracked branch point.  The
existing CWD.5 (corpus_winding_test.jl) and CBr.3 (corpus_elementary_branch
_test.jl) are only BARE winding_delta unit fixtures; CPN.7 drives the bug
through the full Stage-1 walk.

Background ODE (residual==0 gated below before any value is emitted):

  CPN.7  y = sqrt(1 - z^2), C=1.  recast for path_network_solve as the
         2nd-order RHS  f(z,y,yp) = (-1 - yp^2)/y .
    y'  = -z/y                       (first order)
    y'' = (-1 - (y')^2)/y            (2nd-order recast: y*y'' + (y')^2 = -1,
                                      i.e. d/dz(y*y') = -1, the recast of
                                      y^2 = 1 - z^2 differentiated twice)
  Branch points at z = +-1 (zeros of 1 - z^2); a FINITE cut between them.
  IC z0 = 0:  y0 = sqrt(1) = 1,  y'0 = -0/1 = 0.

THE +-sqrt SHEET ORACLE (the SEMANTIC meaning of a sheet flip).  At the
plan's verified test point z = 1.5 + 0.5i:
    principal  y = +sqrt(1 - z^2) = 0.6335517... - 1.1838022...i
    after one CCW loop about z=1 the sqrt sheet FLIPS SIGN:
               y = -sqrt(1 - z^2) = -0.6335517... + 1.1838022...i
This closed-form pair is the INDEPENDENT oracle for what `sheet_index = +1`
*means*; it is NOT the per-node Padre output (the single-step Padre across a
near-branch leap is numerically poor -- value-level pins there would be
dishonest; we pin the BOOKKEEPING, which is exactly what 61um corrupts).

THE 61um REAL-WALK TRIGGER GEOMETRY (verified in Julia, re-confirmed here).
Two single-wedge-step `cross_branch=true` walks about branch=1 on a ring of
radius r=0.45, differing ONLY in step coarseness:

  COARSE leg -- IC at arg +0.82*pi, target at arg +1.18*pi relative to the
    branch, h=0.48.  The single wedge step straddles the arg=pi cut with a
    TRUE subtended angle of +0.555*pi (< pi, CCW).  winding_delta is EXACT
    here; step_sheet_update bumps +1; visited_sheet=[1] is CORRECT.

  FINE leg -- IC at arg -0.47*pi, target at arg +0.59*pi relative to the
    branch, h=0.70.  The single wedge step straddles the cut with a TRUE
    subtended angle of +1.083*pi (> pi, CCW).  winding_delta CLAMPS the raw
    +1.083*pi-ish arg-difference into (-pi,pi] and returns a NEGATIVE value
    (~ -2.88 rad), so step_sheet_update bumps -1 instead of +1:
    visited_sheet=[-1] is the 61um-CORRUPTED sheet.  The TRUE crossing is
    CCW (+1).

The two legs prove the corruption is STEP-SIZE-TRIGGERED: same ODE, same
branch, same cut, same CCW crossing direction; only the step's subtended
angle differs across the pi threshold.

Run:  python3 external/probes/corpus-oracles/pathnet-winding/capture.py
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 50


def emit(label, val):
    if isinstance(val, mp.mpc):
        print(f"const {label} = complex({mp.nstr(val.real, 45)}, "
              f"{mp.nstr(val.imag, 45)})")
    else:
        print(f"const {label} = {mp.nstr(val, 45)}")


# ===========================================================================
# CPN.7  background residual GATE.  y = sqrt(1 - z^2) satisfies
#   y'  = -z / y
#   y'' = (-1 - (y')^2)/y           (the path_network_solve recast RHS)
# ===========================================================================
z = sp.symbols('z')
y = sp.sqrt(1 - z**2)
yp = sp.diff(y, z)
ypp = sp.diff(y, z, 2)
# 1st-order: y' = -z/y.
assert sp.simplify(yp - (-z / y)) == 0, "CPN.7 1st-order residual y'=-z/y"
# 2nd-order recast RHS: f(z,y,yp) = (-1 - yp^2)/y must equal y''.
rhs_recast = (-1 - yp**2) / y
assert sp.simplify(ypp - rhs_recast) == 0, \
    "CPN.7 2nd-order recast f=(-1-yp^2)/y != y''"
# Branch points are the zeros of 1 - z^2: z = +-1.
bp = sp.solve(1 - z**2, z)
assert set(bp) == {-1, 1}, f"CPN.7 branch points != +-1: {bp}"

print("# CPN.7 sqrt: residual==0 GATE PASSED  (y'=-z/y, y''=(-1-yp^2)/y).")
print("# branch points = +-1 (zeros of 1-z^2); FINITE cut between them.")
print()

# IC at z0=0: y0 = sqrt(1) = 1, y'0 = -0/1 = 0.
y0 = mp.sqrt(mp.mpf('1.0'))
assert abs(y0 - 1) < mp.mpf('1e-45'), f"CPN.7 y0 != 1: {y0}"
print("# CPN.7 IC at z0=0: y0=1, y'0=0.")
emit("CPN7_Y0", mp.mpc(y0, 0))
emit("CPN7_YP0", mp.mpc(0, 0))
print()


# ===========================================================================
# THE +-sqrt SHEET ORACLE at z = 1.5 + 0.5i.  mpmath sqrt is the PRINCIPAL
# branch; the flipped sheet is exactly the negation.
# ===========================================================================
zt = mp.mpc('1.5', '0.5')
y_principal = mp.sqrt(1 - zt**2)
y_flipped = -y_principal
print("# +-sqrt sheet oracle at z = 1.5 + 0.5i (principal vs sheet-flipped).")
emit("CPN7_YPRIN_AT_1p5_0p5i", y_principal)   # 0.6335... - 1.1838...i
emit("CPN7_YFLIP_AT_1p5_0p5i", y_flipped)     # -0.6335... + 1.1838...i
# Cross-check: the flip is exactly a sign change.
assert abs((y_flipped + y_principal)) < mp.mpf('1e-45'), "flip != -principal"
print()


# ===========================================================================
# THE 61um REAL-WALK TRIGGER GEOMETRY.  branch = 1; ring radius r = 0.45.
# Endpoints at arg*pi relative to the branch.  Report, for each leg:
#   * the IC and target ring points (full precision),
#   * the TRUE (continuous) subtended angle of the IDEAL chord (arg_new -
#     arg_old, no wrap),
#   * the winding_delta the package returns (the same difference normalised
#     into (-pi,pi]).
# NOTE: the package's ACTUAL wedge step lands slightly OFF the target (the
# +-22.5/+-45 deg wedge offsets perturb it), so the realised subtended angle
# is measured in Julia at test time; these literals pin the RING GEOMETRY and
# the IDEAL-chord winding that motivate the legs.  The load-bearing oracle is
# the TRUE CROSSING DIRECTION (CCW => +1) which is robust to the wedge offset.
# ===========================================================================
branch = mp.mpc(1, 0)
r = mp.mpf('0.45')


def ring_pt(arg_frac):
    return branch + r * mp.e ** (1j * mp.mpf(arg_frac) * mp.pi)


# COARSE leg: arg +0.82*pi -> +1.18*pi  (ideal gap +0.36*pi < pi).
coarse_P = ring_pt('0.82')
coarse_Q = ring_pt('1.18')
coarse_ideal_sub = mp.mpf('1.18') * mp.pi - mp.mpf('0.82') * mp.pi  # +0.36*pi
assert coarse_ideal_sub < mp.pi, "coarse ideal subtended must be < pi"
print("# COARSE leg ring points (arg +0.82*pi -> +1.18*pi rel branch=1, r=0.45).")
emit("CPN7_COARSE_P", coarse_P)
emit("CPN7_COARSE_Q", coarse_Q)
emit("CPN7_COARSE_IDEAL_SUBTENDED", coarse_ideal_sub)   # +0.36*pi
print()

# FINE leg: arg -0.47*pi -> +0.59*pi  (ideal gap +1.06*pi > pi -- the trigger).
fine_P = ring_pt('-0.47')
fine_Q = ring_pt('0.59')
fine_ideal_sub = mp.mpf('0.59') * mp.pi - mp.mpf('-0.47') * mp.pi  # +1.06*pi
assert fine_ideal_sub > mp.pi, "fine ideal subtended must be > pi (the trigger)"
# The wrapped winding_delta of the ideal chord: normalise +1.06*pi into (-pi,pi].
fine_wrapped = fine_ideal_sub
while fine_wrapped > mp.pi:
    fine_wrapped -= 2 * mp.pi
while fine_wrapped <= -mp.pi:
    fine_wrapped += 2 * mp.pi
assert fine_wrapped < 0, "fine wrapped winding_delta must be NEGATIVE (the flip)"
# The error is exactly one missed CCW revolution (-2*pi below the truth).
assert abs((fine_wrapped - fine_ideal_sub) - (-2 * mp.pi)) < mp.mpf('1e-40'), \
    "fine wrap error != -2*pi"
print("# FINE leg ring points (arg -0.47*pi -> +0.59*pi rel branch=1, r=0.45).")
print("#   ideal subtended = +1.06*pi (> pi, CCW); winding_delta wraps to a")
print("#   NEGATIVE value (loses a +2*pi revolution) => sheet bumps -1 not +1.")
emit("CPN7_FINE_P", fine_P)
emit("CPN7_FINE_Q", fine_Q)
emit("CPN7_FINE_IDEAL_SUBTENDED", fine_ideal_sub)        # +1.06*pi
emit("CPN7_FINE_IDEAL_WRAPPED", fine_wrapped)            # -0.94*pi
print()
print("# All residual==0 / geometry GATES PASSED; values above are dps=50 pins.")
