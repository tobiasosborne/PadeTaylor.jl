# test/_oracle_corpus_robust_pade.jl
#
# Pinned oracle constants for corpus_robust_pade_test.jl (bead
# padetaylor-7pav, epic padetaylor-p1v0).  Split out of the test body to
# keep the test file under the 200-code-LOC cap (CLAUDE.md Rule 6).
#
# PROVENANCE of every constant in this file:
#
#   * CRP_EXP77_A / _B — exact rational Padé(7,7) of exp, b[1]=1 convention,
#     a_k = C(14-k,7)·7!·7! / (14!·k!), b_k = (-1)^k a_k.  Closed form;
#     regenerate with the python3 one-liner in the catalogue record
#     `exp-77-exact-rational` (docs/test_corpus/01_corpus_catalogue.md:1468).
#     Cross-checked here by `fractions.Fraction`: a[1]=1/2, b[1]=-1/2 exactly.
#
#   * CRP_TAN33_A / _B — closed-form (3,3) Padé of tan: p = z - z^3/15,
#     q = 1 - 2z^2/5.  Catalogue `tan-33-exact-pade` (md:1482).  The
#     denominator root is ±sqrt(5/2) = ±1.5811388300841898.
#
#   * CRP_EXPI_R1 — exp(iz) Padé(7,7) evaluated at z=1; the exact target is
#     exp(i) itself (the rational interpolates exp(iz) through O(z^15)).
#     Catalogue `exp-iz-complex-pade` (md:132, exp-iz-77-complex-bigfloat
#     md:1544): |r(1) - exp(i)| ~ 2.5e-16, near machine ε.
#
#   * CRP_SVD_SIGMA2 — singular value σ_2 of the 10×11 lower-triangular
#     Toeplitz of c_k = 1 + 1e-13·(1/2)^k.  mpmath dps=50, captured by
#     external/probes/corpus-oracles/robust-pade/capture.py:
#         σ_1 = 10.48808848170153451124194
#         σ_2 = 4.759563572881758889598953e-14
#         σ_3 = 4.85e-51  (exact-arithmetic zero)
#     Matches catalogue `adr0002-rank-boundary-svd` (md:1533) to 20 sig figs.
#
#   * CRP_NOISY_C — deterministic tiny-noise coefficient vector
#     c_k = 1 + 1e-9·sin(k+1), captured by the same capture.py at dps=50
#     and rounded to Float64.  Used to pin the geom-noisy-1on1z structural
#     reduction (tol > noise ⇒ type (0,1)) WITHOUT an RNG dependence — the
#     catalogue's frozen Octave randn(state,42) oracle (test/_oracles.jl
#     2.1.6) is unreproducible locally, so we regenerate the *structural*
#     claim deterministically instead (catalogue note md:1508).

# --- exp-77-exact-rational (closed form, b[1]=1 normalisation) --------------
const CRP_EXP77_A = [1.0, 0.5, 3/26, 5/312, 5/3432, 1/11440, 1/308880, 1/17297280]
const CRP_EXP77_B = [1.0, -0.5, 3/26, -5/312, 5/3432, -1/11440, 1/308880, -1/17297280]

# --- tan-33-exact-pade (closed form) ---------------------------------------
# p = z - z^3/15  ->  a = [0, 1, 0, -1/15];  q = 1 - 2z^2/5 -> b = [1, 0, -2/5].
const CRP_TAN33_A = [0.0, 1.0, 0.0, -1/15]
const CRP_TAN33_B = [1.0, 0.0, -2/5]
const CRP_TAN33_ROOT = sqrt(5 / 2)            # 1.5811388300841898

# --- exp-iz-complex-pade ----------------------------------------------------
const CRP_EXPI_R1 = exp(im)                   # exact target for r(1)

# --- adr0002-rank-boundary-svd (mpmath dps=50, capture.py) ------------------
const CRP_SVD_SIGMA1 = 10.48808848170153451124194
const CRP_SVD_SIGMA2 = 4.759563572881758889598953e-14
# Same σ_2 as a string, for full-precision BigFloat comparison (the Float64
# literal above loses the digits past the 16th; the BF256 Jacobi assertion
# needs all 25 to bite at rtol 1e-18).
const CRP_SVD_SIGMA2_STR = "4.759563572881758889598953e-14"
# Float64 GGT threshold 1e-14·‖c‖_2 = 4.582575694955883e-14 (σ_2/thr = 1.039).
const CRP_SVD_THRESHOLD_F64 = 4.582575694955883e-14

# --- geom-noisy-1on1z (deterministic mpmath noise, capture.py) -------------
const CRP_NOISY_C = [1.00000000084147098, 1.00000000090929743, 1.00000000014112001,
    0.999999999243197505, 0.999999999041075725, 0.999999999720584502,
    1.0000000006569866, 1.00000000098935825, 1.00000000041211849,
    0.999999999455978889, 0.999999999000009793, 0.999999999463427082,
    1.00000000042016704, 1.00000000099060736, 1.00000000065028784,
    0.999999999712096683, 0.999999999038602508, 0.999999999249012753,
    1.00000000014987721, 1.00000000091294525, 1.00000000083665564]
