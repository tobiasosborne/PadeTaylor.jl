# test/_oracle_corpus_elementary_branch.jl
#
# Pinned sympy/mpmath oracle constants for the elementary branch-point corpus
# (bead padetaylor-g9lg, epic padetaylor-25og).  Family E of
# docs/test_corpus/02_corpus_extension_plan.md:245-270.
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/elementary-branch/capture.py  (mpmath 1.3,
#   mp.dps = 50; sympy residual==0 GATE on every recast + the atanh==log
#   derivative identity + the ±√2 / ±1 branch-point sets, BEFORE any value is
#   emitted).  Each transcendental constant below is a verbatim copy of that
#   script's labelled output.  To update: re-run `python3 capture.py`, re-pin.
#
# Exact rationals (the ICs u(1)=0/u'(1)=1, u(0)=0/u'(0)=2, y0=1.5; the value
# u'(1/2)=8/3, u'(2)=1/2) are written inline in the test where they are exact.
#
# CLOSED FORMS:
#   CBr.2  u = log z              ; u' = 1/z,         u'' = -1/z^2.
#   CBr.3  u = 2*atanh z          ; u' = 2/(1-z^2),   u'' = 4z/(1-z^2)^2.
#   CBr.1  y = sqrt(C-z^2), C=2   ; y' = -z/y,        y'' = -2/y^3.

# ---------------------------------------------------------------------------
# CBr.2 log: real probe from z0=1 (u(1)=0, u'(1)=1).  u=log z, u'=1/z.
# The per-loop monodromy jump of log is exactly 2*pi*i (purely imaginary);
# the imag part = 2*pi is pinned here, the realness inline.
# ---------------------------------------------------------------------------
const CBR2_U_AT_2 = 0.693147180559945309417232121458176568075500134
const CBR2_U_AT_E = 1.0
const CBR2_UP_AT_2 = 0.5
const CBR2_LOG_JUMP_PER_LOOP_IMAG = 6.2831853071795864769252867665590057683943388

# ---------------------------------------------------------------------------
# CBr.3 atanh: real probe inside (-1,1) from z0=0 (u(0)=0, u'(0)=2).
# u(1/2) = 2*atanh(1/2) = log 3 (both pins agree to dps=50); u'(1/2) = 8/3.
# ---------------------------------------------------------------------------
const CBR3_U_AT_HALF = 1.09861228866810969139524523692252570464749056
const CBR3_LOG3 = 1.09861228866810969139524523692252570464749056
const CBR3_UP_AT_HALF = 2.66666666666666666666666666666666666666666667

# ---------------------------------------------------------------------------
# CBr.3 — THE padetaylor-61um STEP.  Branch at z=1; both endpoints OFF the real
# axis at NONZERO arg relative to the branch (z_old-1 at arg -0.45*pi, z_new-1
# at arg +0.65*pi).  TRUE subtended angle = +1.1*pi (CCW, > pi); winding_delta
# normalises into (-pi,pi] and returns the WRAPPED -0.9*pi (loses one full
# revolution = exactly -2*pi).  See capture.py for the mpmath verification.
# ---------------------------------------------------------------------------
const CBR3_61UM_TRUE_SUBTENDED = 3.45575191894877256230890772160745317261688634   # +1.1*pi
const CBR3_61UM_WRAPPED = -2.82743338823081391461637904495155259577745246          # -0.9*pi

# ---------------------------------------------------------------------------
# CBr.1 sqrt finite cut: IC off the cut at z0=0.5i (principal +sqrt sheet),
# branch points ±√2.  Principal-sheet (+sqrt) probe values y=+sqrt(2-z^2);
# the opposite (-sqrt) sheet is exactly the negation.
# ---------------------------------------------------------------------------
const CBR1_Y0 = complex(1.5, 0.0)
const CBR1_YP0 = complex(0.0, -0.333333333333333333333333333333333333333333333)
const CBR1_YPRIN_AT_0p5i = complex(1.5, 0.0)
const CBR1_YPRIN_AT_1i = complex(1.73205080756887729352744634150587236694280525, 0.0)
const CBR1_YPRIN_AT_0p3_0p5i = complex(1.47321651861604337335288261555953634983658773, -0.101818027495993415730376477039859439389017948)
const CBR1_YMINUS_AT_0p5i = complex(-1.5, 0.0)
const CBR1_BRANCH_POS = 1.41421356237309504880168872420969807856967188
const CBR1_BRANCH_NEG = -1.41421356237309504880168872420969807856967188
