# test/_oracle_corpus_special_fn_bvp.jl
#
# Pinned TRANSCENDENTAL oracle constants for the special-function BVP corpus
# (bead padetaylor-evlo, epic padetaylor-25og).  Tagged CBvx.2 (Mathieu) and
# CBvx.3-Whittaker.  The CBvx.3 Kummer case is an EXACT degree-2 polynomial and
# is pinned INLINE as rationals in the test (no constant lives here).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/special-fn-bvp/capture.py  (python3 +
#   sympy + scipy.special + mpmath 1.3 dps=50).  Every recast's ODE residual is
#   verified IDENTICALLY ZERO symbolically (sympy) — a GATE that refuses to emit
#   otherwise — before any number below is printed.  Each constant is a verbatim
#   copy of that script's labelled output.  To update: re-run `python3
#   capture.py` and re-pin.  Ground-truth recipe is
#   docs/test_corpus/02_corpus_extension_plan.md:309-327 (Family H).
#
# These are package-INDEPENDENT (Rule 5): scipy.special supplies BOTH the
# Mathieu characteristic value a_2(1) AND the ce_2 function values; mpmath
# supplies the Whittaker M values.  No PadeTaylor output is pinned.

# ---------------------------------------------------------------------------
# CBvx.2  Mathieu  u'' + (a - 2q cos2z) u = 0,  q = 1.  The eigenvalue
# a = a_2(q=1) is what makes ce_2 (the even pi-periodic Mathieu function) the
# solution.  DLMF 28.2/28.4.  scipy double precision -> the values carry only a
# ~1e-16 floor, hence the test tol 1e-10.
# scipy GOTCHA: mathieu_cem(m, q, x_DEGREES) -> (value, deriv); convert via
# numpy.degrees(z_radians).
# ---------------------------------------------------------------------------
const CBvx2_A      = 4.371300982735086     # scipy.special.mathieu_a(2, 1.0)
const CBvx2_Q      = 1.0
const CBvx2_BC_A   = 1.0859619368890114    # u(0)    = ce_2(0)
const CBvx2_BC_B   = -0.8157268390500979   # u(pi/2) = ce_2(pi/2)
const CBvx2_U_pi4  = 0.2986515729070586    # interior ce_2(pi/4)
const CBvx2_U_pi3  = -0.21370917207980372  # interior ce_2(pi/3)

# ---------------------------------------------------------------------------
# CBvx.3-Whittaker  u'' + (-1/4 + k/z + (1/4 - m^2)/z^2) u = 0, k=2, m=3/2, on
# [1,4].  Solution = Whittaker M_{k,m}(z).  DLMF 13.14.  mpmath.whitm dps=50 ->
# 45 significant digits captured; the BVP solver's Float64 spectral floor sets
# the test tol 1e-12.
# ---------------------------------------------------------------------------
const CBvxW_BC_A  = 0.606530659712633423603799534991180453441918135  # u(1)
const CBvxW_BC_B  = 2.16536453178580307030399191955975045452210473   # u(4)
const CBvxW_U_2   = 1.47151776468576928638209508064584346978324452   # interior
const CBvxW_U_2p5 = 1.79065498037618812703053391654898501745719245   # interior
const CBvxW_U_3   = 2.00817144133586846039952423687611269207954466   # interior
