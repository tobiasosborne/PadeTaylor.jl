# test/_oracle_corpus_periodic_pole.jl
#
# Pinned closed-form / mpmath oracle constants for the periodic-vertical-pole
# corpus (bead padetaylor-5rks, epic padetaylor-25og).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/periodic-pole/capture.py  (mpmath 1.3,
#   mp.dps = 50; sympy residual==0 + residue+1 + 2πi-spacing GATE before any
#   value is emitted).  Every transcendental constant below is a verbatim copy
#   of that script's labelled output.  To update: re-run `python3 capture.py`
#   and re-pin.  Ground-truth recipe is from
#   docs/test_corpus/02_corpus_extension_plan.md:232-241 (Family D).
#
# Exact rationals (the IC u(0)=1/2, u'(0)=1/4; the past-pole real part Re u=1/2)
# are written inline in the test, not here — only transcendental pins live here.
#
# CLOSED FORM:  c = 1  ⇒  u(z) = 1/(1 + e^{-z}).  Poles z = iπ·(2k+1), an
# infinite vertical row spaced 2πi, all simple, residue +1 (capture.py GATE).

# ---------------------------------------------------------------------------
# Spot values, c = 1, u(z) = 1/(1+e^{-z}).  Real z and one complex z.
# ---------------------------------------------------------------------------
const CVROW_U_AT_1p0 = 0.73105857863000487925115924182183627436514464
const CVROW_U_AT_2p0 = 0.880797077977882444059729141302396795206384299
const CVROW_U_AT_0p5p0p5i = complex(0.629935440242770308157204804956469613032098034, 0.119545057814595607504560270188972458607672303)

# ---------------------------------------------------------------------------
# Past-pole bridge target  z = i(π+0.3)  (just past the nearest pole at iπ).
# Re u = 1/2 EXACTLY (logistic reflection u(z)+u(-z)=1 about the pole); the
# imaginary part is the transcendental pin.  u' = u(1-u) is REAL there
# (imag part ~1e-50, i.e. numerically zero) — pinned inline.
# ---------------------------------------------------------------------------
const CVROW_U_AT_PASTPOLE  = complex(0.5, -3.308295752794975049337198542995000070038868)
const CVROW_UP_AT_PASTPOLE = 11.1948207879612706619684740398105774354926394

# ---------------------------------------------------------------------------
# Vertical-row pole locations (c=1):  z = iπ·(2k+1).  Nearest two.
# spacing 2πi = 6.283185307179586... (purely imaginary).
# ---------------------------------------------------------------------------
const CVROW_POLE_0 = complex(0.0, 3.1415926535897932384626433832795028841971694)
const CVROW_POLE_1 = complex(0.0, 9.4247779607693797153879301498385086525915082)
const CVROW_SPACING_2PI = 6.2831853071795864769252867665590057683943388
