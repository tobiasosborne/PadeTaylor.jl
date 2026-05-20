# =============================================================================
# _oracle_shared_pade.jl — pinned Calgo-766 type-II shared-Q oracle literals.
#
# These are the double-precision outputs of ACM TOMS Algorithm 766
# (Cabay–Jones–Labahn 1997, "Algorithm 766: Experiments with a Weakly Stable
# Algorithm for Computing Padé–Hermite and Simultaneous Padé Approximants",
# ACM TOMS 23(1), 91–110, DOI 10.1145/244768.244790) — an INDEPENDENT FORTRAN
# reference implementation of the type-II Hermite–Padé (simultaneous-Padé /
# shared-denominator) approximant, using the Beckermann–Labahn striped-Sylvester
# recurrence — a different algorithm from `SharedPade`'s block-Toeplitz SVD.
#
# Probe (capture script, FORTRAN sources, build/sanity record):
#   external/probes/calgo766-oracle/   (oracles.txt, capture.jl, README.md)
# Captured 2026-05-20 via libcalgo766_dp.so (gfortran 13.3.0, -fdefault-real-8 →
# IEEE double precision; accuracy floor ~1e-15).
#
# CONVENTION CAVEAT — read external/probes/calgo766-oracle/README.md §"Convention
# mapping" before using these.  Calgo's type-II degrees are `||n|| - n(j)`, set
# by the per-series degree vector `n`; they CANNOT in general be made equal to
# `SharedPade`'s single shared `m`.  For the d=2 case `SP_1_2` the oracle was
# run at `n = [2,3,2]`, so its `Q` is degree 5 — the true degree-2
# `Q = [1,-0.4,0.2]` times a degree-3 SPURIOUS common factor that also divides
# `P_1, P_2`.  The pinned coefficients therefore look nothing like the true
# low-degree coefficients BY DESIGN.  Compare oracle and `SharedPade` ONLY by
#   (a) rational-function VALUES at sample points (the spurious factor cancels
#       in the reduced ratio P_i/Q), and
#   (b) POLE LOCATIONS — the true poles are the SUBSET of the roots of the
#       degree-5 oracle `Q` that coincide with the roots of the reduced `Q`.
# NEVER compare raw coefficient vectors against these literals.
#
# The d=1 case `SP_1_1` was run at `n = [1,2]`, giving an honest degree-2 `Q`;
# there the coefficients DO match (no spurious factor) — that case is also a
# direct coefficient cross-check.
#
# Q, P_i coefficients are low-to-high (z⁰ first), renormalised so Q(0)=1, matching
# `SharedPade.shared_denominator_pade`'s convention.  Transcribed verbatim from
# external/probes/calgo766-oracle/oracles.txt.  The test suite asserts against
# THESE LITERALS — it does NOT `ccall` FORTRAN, so the suite has no gfortran
# dependency (the v0.1 `_oracles.jl` pinning convention).
# =============================================================================

# --- Case SP_1_1_d1_scalar — f(z) = (1+2z)/(1-0.5z+0.3z²) ---------------------
# Calgo n = [1, 2], k = 1, flag = 0.  Honest degree-2 Q (no spurious factor).
const CALGO_SP_1_1_Q  = [1.000000000000000000e+00,
                         -5.000000000000002220e-01,
                          2.999999999999999334e-01]
const CALGO_SP_1_1_P1 = [1.000000000000000000e+00,
                          2.000000000000000000e+00]
const CALGO_SP_1_1_flag = 0

# --- Case SP_1_2_d2_shared_Q — f₁=(1+0.3z)/Q, f₂=(2-0.7z+0.5z²)/Q, ------------
# --- Q = 1-0.4z+0.2z²  -------------------------------------------------------
# Calgo n = [2, 3, 2], k = 2, flag = 2.  Degree-5 Q: true degree-2 Q times a
# degree-3 SPURIOUS factor (see the convention caveat above).  flag = 2 is
# benign for the type-II row-0 result — see calgo766-oracle/README.md
# §"Sanity-check result".  Compare by VALUES and POLE-SUBSET only.
const CALGO_SP_1_2_Q  = [ 1.000000000000000000e+00,
                         -1.931237386841950787e+00,
                          9.350227757596708544e-01,
                         -3.982334371216044167e-01,
                          4.169549674220130592e-02,
                         -8.594966268811616086e-03]
const CALGO_SP_1_2_P1 = [ 1.000000000000000000e+00,
                         -1.231237386841950832e+00,
                         -3.368433950296947366e-01,
                         -6.216485037190882211e-03,
                         -1.289244940321745275e-02]
const CALGO_SP_1_2_P2 = [ 2.000000000000000000e+00,
                         -3.762474773683901486e+00,
                          1.816921812835146754e+00,
                         -9.373378308251147928e-01,
                          9.134629245228582939e-02,
                         -2.148741567202904976e-02]
const CALGO_SP_1_2_flag = 2
