#!/usr/bin/env python3
"""capture.py — sympy/mpmath oracle for corpus_elementary_branch_test.jl.

Family E (docs/test_corpus/02_corpus_extension_plan.md:245-270): elementary
closed-form branch-point ODEs that bind the SheetTracker / BranchTracker
winding machinery to genuinely multivalued ground truth (v1 only ever
exercised the sheet primitives through the PVI/PIII transforms or the entire
ODE u''=u', where branches are pure routing constraints).

Three cases (each gated by a sympy residual==0 check BEFORE any value is
emitted, then pinned numerically by mpmath at dps=50):

  CBr.2  log single branch.   u' = 1/z, recast u'' = -1/z^2, u = log z.
         Branch at z=0; cut along arg=pi (BranchTracker default).  Real probe
         from z0=1 (u(1)=0, u'(1)=1).  Winding: CCW unit circle, total 2*pi,
         sheet +1, log jumps exactly +2*pi*i per loop.

  CBr.3  atanh two-branch (the padetaylor-61um catcher).
         u' = 2/(1-z^2), recast u'' = 4z/(1-z^2)^2, u = 2*atanh(z) =
         log((1+z)/(1-z)).  Branch points +-1.  Real probe inside (-1,1) from
         z0=0 (u(0)=0, u'(0)=2): u(1/2) = 2*atanh(1/2) = log 3.  THE 61um BUG
         STEP: a single step about branch z=1 subtending |dtheta| ~ 1.1*pi CCW
         (z_old = 1 + 0.1*exp(-0.45*pi*i), z_new = 1 + 0.1*exp(+0.65*pi*i)).
         TRUE subtended angle = +1.1*pi; winding_delta clamps to (-pi,pi] and
         returns the WRAPPED -0.9*pi (loses a full revolution).

  CBr.1  sqrt finite cut.   y' = -z/y, recast y'' = -C/y^3, y = sqrt(C-z^2),
         C=2.  Branch points +-sqrt2, a FINITE cut between them.  Principal-
         sheet value via path_network_solve + eval_at_sheet([0,0]) vs
         +sqrt(2-z^2); cross-sheet -sqrt via eval_at_sheet([0,1]) if clean.

Run:  python3 external/probes/corpus-oracles/elementary-branch/capture.py
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
# CBr.2  log single branch — RESIDUAL GATE.  u = log z  satisfies
#   u'  = 1/z       (first order)
#   u'' = -1/z^2    (2nd-order recast for solve_pade)
# ===========================================================================
z = sp.symbols('z')
u_log = sp.log(z)
assert sp.simplify(sp.diff(u_log, z) - 1 / z) == 0, "CBr.2 1st-order residual"
assert sp.simplify(sp.diff(u_log, z, 2) + 1 / z**2) == 0, "CBr.2 2nd-order recast"

# Monodromy: log jumps by +2*pi*i per CCW loop about z=0 (the branch).
# log(e^{i*2pi}) - log(1) = 2*pi*i analytically continued.  Pin spot values.
print("# CBr.2 log: residual==0 GATE PASSED (u'=1/z, u''=-1/z^2).")
print()
print("# CBr.2 real probe from z0=1 (u(1)=0, u'(1)=1).  u=log z, u'=1/z.")
emit("CBR2_U_AT_2", mp.log(mp.mpf('2.0')))            # log 2
emit("CBR2_U_AT_E", mp.log(mp.e))                     # = 1 exactly
emit("CBR2_UP_AT_2", 1 / mp.mpf('2.0'))               # u'(2) = 1/2
# The per-loop log jump is exactly 2*pi*i (purely imaginary).
emit("CBR2_LOG_JUMP_PER_LOOP_IMAG", 2 * mp.pi)        # imag part of 2*pi*i
print()


# ===========================================================================
# CBr.3  atanh two-branch — RESIDUAL GATE.  u = 2*atanh(z) = log((1+z)/(1-z))
#   u'  = 2/(1-z^2)        (first order)
#   u'' = 4z/(1-z^2)^2     (2nd-order recast)
# ===========================================================================
u_at = 2 * sp.atanh(z)
assert sp.simplify(sp.diff(u_at, z) - 2 / (1 - z**2)) == 0, "CBr.3 1st-order"
assert sp.simplify(sp.diff(u_at, z, 2) - 4 * z / (1 - z**2)**2) == 0, \
    "CBr.3 2nd-order recast"
# The two closed forms agree on (-1,1) (cosmetic principal-branch relabel
# elsewhere — per the plan we pin the DERIVATIVE identity + the numeric value,
# NOT the raw log-vs-atanh STRING form, which sympy refuses to collapse without
# force=True precisely because of that branch relabel).  We GATE on:
#   (i) the canonical log expansion of 2*atanh agreeing with log(1+z)-log(1-z),
#   (ii) the numeric value 2*atanh(1/2) == log 3 to >60 digits.
atanh_as_log = sp.simplify(
    sp.expand_log(2 * sp.atanh(z).rewrite(sp.log), force=True)
    - (sp.log(1 + z) - sp.log(1 - z)))
assert atanh_as_log == 0, f"CBr.3 atanh==log identity (expand_log): {atanh_as_log}"
num_gap = abs(sp.N(2 * sp.atanh(sp.Rational(1, 2)) - sp.log(3), 80))
assert num_gap < sp.Float('1e-60'), f"CBr.3 atanh==log numeric gap: {num_gap}"
# Branch points are the zeros of (1-z^2): z = +-1.
bp = sp.solve(1 - z**2, z)
assert set(bp) == {-1, 1}, f"CBr.3 branch points != +-1: {bp}"

print("# CBr.3 atanh: residual==0 GATE PASSED (u'=2/(1-z^2), u''=4z/(1-z^2)^2).")
print("# branch points = +-1 (zeros of 1-z^2).")
print()
print("# CBr.3 real probe inside (-1,1) from z0=0 (u(0)=0, u'(0)=2).")
emit("CBR3_U_AT_HALF", 2 * mp.atanh(mp.mpf('0.5')))   # = log 3
emit("CBR3_LOG3", mp.log(mp.mpf('3.0')))              # closed-form cross-check
emit("CBR3_UP_AT_HALF", 2 / (1 - mp.mpf('0.5')**2))   # u'(1/2) = 8/3
print()

# ---------------------------------------------------------------------------
# THE 61um BUG STEP.  Branch at z=1.  Both endpoints OFF the real axis at
# NONZERO arg relative to the branch (distinct from CWD.5, where the branch
# is at 0 and both endpoints lie ON the real axis).  Subtended 1.1*pi CCW.
#   z_old - 1 = 0.1 * exp(-0.45*pi*i)  ->  arg = -0.45*pi
#   z_new - 1 = 0.1 * exp(+0.65*pi*i)  ->  arg = +0.65*pi
#   TRUE subtended dtheta = 0.65*pi - (-0.45*pi) = +1.1*pi  (CCW, > pi)
# winding_delta forms angle(z_new-1) - angle(z_old-1) = +1.1*pi, then the
# (-pi,pi] normalisation subtracts 2*pi -> WRAPPED -0.9*pi (a lost revolution).
# ---------------------------------------------------------------------------
branch = mp.mpc(1, 0)
# Full-precision exponents (avoid Float64 contamination from 0.45/0.65 literals).
ang_old = -mp.mpf('0.45') * mp.pi
ang_new = +mp.mpf('0.65') * mp.pi
z_old = 1 + mp.mpf('0.1') * mp.e ** (1j * ang_old)
z_new = 1 + mp.mpf('0.1') * mp.e ** (1j * ang_new)
arg_old = mp.arg(z_old - branch)
arg_new = mp.arg(z_new - branch)
true_subtended = arg_new - arg_old                     # +1.1*pi (no wrap)
assert abs(true_subtended - mp.mpf('1.1') * mp.pi) < mp.mpf('1e-40'), \
    f"CBr.3 true subtended != 1.1*pi: {true_subtended}"
# What winding_delta returns: the same difference normalised into (-pi,pi].
wrapped = true_subtended
while wrapped > mp.pi:
    wrapped -= 2 * mp.pi
while wrapped <= -mp.pi:
    wrapped += 2 * mp.pi
assert abs(wrapped - (-mp.mpf('0.9') * mp.pi)) < mp.mpf('1e-40'), \
    f"CBr.3 wrapped != -0.9*pi: {wrapped}"
# The error is EXACTLY one missed revolution (-2*pi).
assert abs((wrapped - true_subtended) - (-2 * mp.pi)) < mp.mpf('1e-40'), \
    "CBr.3 wrap error != -2*pi"
print("# CBr.3 61um step about branch z=1, endpoints at arg -0.45*pi / +0.65*pi.")
print("#   true subtended = +1.1*pi ; winding_delta wraps to -0.9*pi "
      "(loses 2*pi).")
emit("CBR3_61UM_TRUE_SUBTENDED", true_subtended)       # +1.1*pi
emit("CBR3_61UM_WRAPPED", wrapped)                     # -0.9*pi
emit("CBR3_61UM_Z_OLD", z_old)
emit("CBR3_61UM_Z_NEW", z_new)
print()


# ===========================================================================
# CBr.1  sqrt finite cut — RESIDUAL GATE.  y = sqrt(C - z^2), C=2 satisfies
#   y'  = -z / y          (first order)
#   y'' = -C / y^3        (2nd-order recast)
# Branch points +-sqrt(C) = +-sqrt(2); a FINITE cut between them.
# ===========================================================================
C = sp.Integer(2)
y_sqrt = sp.sqrt(C - z**2)
assert sp.simplify(sp.diff(y_sqrt, z) - (-z / y_sqrt)) == 0, "CBr.1 1st-order"
assert sp.simplify(sp.diff(y_sqrt, z, 2) - (-C / y_sqrt**3)) == 0, \
    "CBr.1 2nd-order recast"
bp1 = sp.solve(C - z**2, z)
assert set(bp1) == {-sp.sqrt(2), sp.sqrt(2)}, f"CBr.1 branch pts: {bp1}"

print("# CBr.1 sqrt: residual==0 GATE PASSED (y'=-z/y, y''=-2/y^3, C=2).")
print("# branch points = +-sqrt(2); FINITE cut between them.")
print()

# IC off the cut on the imaginary axis: z0 = 0 + 0.5i.
#   y0 = sqrt(2 - (0.5i)^2) = sqrt(2 + 0.25) = sqrt(2.25) = 1.5  (principal +).
z0 = mp.mpc(0, '0.5')
C_m = mp.mpf('2.0')
y0 = mp.sqrt(C_m - z0**2)
assert abs(y0 - mp.mpf('1.5')) < mp.mpf('1e-45'), f"CBr.1 y0 != 1.5: {y0}"
yp0 = -z0 / y0
print("# CBr.1 IC off the cut: z0 = 0.5i, y0 = +sqrt(2-(0.5i)^2) = 1.5.")
emit("CBR1_Y0", y0)                                    # 1.5 + 0i (principal +)
emit("CBR1_YP0", yp0)                                  # -z0/y0

# Principal-sheet (+sqrt) probe values along the imaginary axis (off the cut,
# safe for the path-network walk).  These pin eval_at_sheet([0,0]) vs +sqrt.
print("# CBr.1 principal-sheet (+sqrt) probe values y = +sqrt(2-z^2):")
for zz, lab in [(mp.mpc(0, '0.5'), "AT_0p5i"),
                (mp.mpc(0, '1.0'), "AT_1i"),
                (mp.mpc('0.3', '0.5'), "AT_0p3_0p5i")]:
    yv = mp.sqrt(C_m - zz**2)              # mpmath sqrt = principal (+) branch
    emit(f"CBR1_YPRIN_{lab}", yv)
# The opposite (-sqrt) sheet is exactly the negation.
emit("CBR1_YMINUS_AT_0p5i", -mp.sqrt(C_m - mp.mpc(0, '0.5')**2))
print()
emit("CBR1_BRANCH_POS", mp.sqrt(2))                    # +sqrt(2)
emit("CBR1_BRANCH_NEG", -mp.sqrt(2))                   # -sqrt(2)
print()
print("# All residual==0 GATES PASSED; values above are mpmath dps=50 pins.")
