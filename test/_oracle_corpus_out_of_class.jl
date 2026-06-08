# test/_oracle_corpus_out_of_class.jl
#
# Pinned closed-form / mpmath oracle constants for the out-of-class / fail-loud
# corpus (bead padetaylor-nzcj, epic padetaylor-25og; bug padetaylor-v1ub).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/out-of-class/capture.py  (mpmath 1.3,
#   mp.dps = 50; sympy residual==0 GATE + ESSENTIAL-singularity confirmation
#   before any value is emitted).  Every transcendental constant below is a
#   verbatim copy of that script's labelled output.  To update: re-run
#   `python3 capture.py` and re-pin.  Ground-truth recipe is from
#   docs/test_corpus/02_corpus_extension_plan.md:292-306 (Family G).
#
# CLOSED FORM:  u(z) = e^{1/z}  solves  u'' = u·(1+2z)/z^4  (sympy residual ≡ 0).
# z = 0 is an ESSENTIAL singularity (Laurent series 1 + 1/z + 1/(2z²) + … has
# infinitely many negative-power terms — NOT a pole).  The package is
# meromorphic-only (Rule 1); this is the out-of-class input that documents the
# candidate fail-loud gap padetaylor-v1ub.  On z<0, u = e^{1/z} ∈ (0,1) and
# → 0 as z→0⁻; on z>0 it blows up as z→0⁺.

# ---------------------------------------------------------------------------
# Negative-axis reference values u(z) = e^{1/z} (approach z=0⁻ from z0 = -1).
# These are exactly e^{-1/|z|}; the relerr-vs-distance curve is measured here.
# ---------------------------------------------------------------------------
const CFAIL_U_AT_NEG0p5  = 0.135335283236612691893999494972484403407631546   # e^{-2}
const CFAIL_U_AT_NEG0p2  = 0.00673794699908546709663604842314842424884958503 # e^{-5}
const CFAIL_U_AT_NEG0p1  = 0.0000453999297624848515355915155605506102379180889 # e^{-10}
const CFAIL_U_AT_NEG0p05 = 0.00000000206115362243855782796594038015582097637580728 # e^{-20}

# ---------------------------------------------------------------------------
# Across-zero target: integrate z0 = -0.2 (u = e^{-5}) → z = +0.2 (u = e^{+5}).
# Crossing z=0 crosses the essential singularity; e^{1/z} → +∞ as z→0⁺.
# ---------------------------------------------------------------------------
const CFAIL_U_AT_POS0p2 = 148.413159102576603421115580040552279623487668      # e^{+5}
