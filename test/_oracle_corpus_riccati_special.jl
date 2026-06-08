# test/_oracle_corpus_riccati_special.jl
#
# Pinned mpmath-dps=50 oracle constants for the special-function log-derivative
# Riccati corpus, TRANSCENDENTAL subfamily (bead padetaylor-tyef, epic
# padetaylor-25og).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/riccati-special/capture.py  (mpmath 1.3,
#   mp.dps = 50; sympy + numeric residual==0 GATE — CRic.1/2/5 recasts verified
#   identically 0 symbolically and <1e-30 against the closed-form special
#   function — before any value is emitted).  Every constant below is a verbatim
#   copy of that script's labelled output.  To update: re-run `python3
#   capture.py` and re-pin.  Ground-truth recipe is
#   docs/test_corpus/02_corpus_extension_plan.md:134-166 (Family A).
#
# Residues are NOT pinned as numbers (they are exact integers ±1 by the
# log-derivative theorem — capture verifies them); the tests assert the
# pole-bridging SIGN-FLIP that a residue-±1 simple pole forces.

# ---------------------------------------------------------------------------
# CRic.1  Airy  u(z) = -Ai'(z)/Ai(z).  IC z0=0; u0 = -Ai'(0)/Ai(0), up0 = u0^2.
# Poles at the Airy zeros a_k on the NEGATIVE real axis, residue -1.
# ---------------------------------------------------------------------------
const CRic1_U0  = 0.729011132947226981418636264703935975972769969
const CRic1_UP0 = 0.531457231960999452867124932783597512879376261
const CRic1_A1 = -2.33810741045976703848919725244673544063854015
const CRic1_A2 = -4.0879494441309706166369887014573910602247647
const CRic1_A3 = -5.52055982809555105912985551293129357379721428
const CRic1_U_m0p5 = 0.428988058385715955920865714014064631235275612
const CRic1_U_m1 = 0.0189718245555636599559750467452062630588004464
const CRic1_U_m1p5 = -0.665982954300023111374311204206832333010023949
const CRic1_U_m2 = -2.71872834423580371936441409033060593400755086
const CRic1_U_m2p5 = 6.043644114423309776075097315779501533372407
const CRic1_U_p0p5 = 0.97072394910167399736949306240730576569846182
const CRic1_U_p1 = 1.17632196714370102308934626278734984847391164

# ---------------------------------------------------------------------------
# CRic.2  Bessel  u(z) = -J1(z)/J0(z) = J0'/J0.  Seed z0=1 (1/z, 1/z^2 RHS
# coeff-singularity at the origin).  IC u0 = -J1(1)/J0(1), up0 = -u0^2-u0-1.
# Poles at j_{0,k} on the POSITIVE real axis, residue +1.
# ---------------------------------------------------------------------------
const CRic2_U0  = -0.575080915004305960499443395318508359290982728
const CRic2_UP0 = -0.755637143797883815908522675518263414915886543
const CRic2_J1 = 2.40482555769577276862163187932645464312424491
const CRic2_J2 = 5.52007811028631064959660411281302742522186548
const CRic2_U_p1p5 = -1.09008664189217921419197948654486324774393116
const CRic2_U_p2 = -2.57592032136822195685749678231504449061298195
const CRic2_U_p3 = 1.30381238108283730118996432576765514140432764
const CRic2_U_p3p5 = 0.361398321961270650121208883640227644952131419

# ---------------------------------------------------------------------------
# CRic.5  parabolic cylinder  u(z) = -U'(0,sqrt2 z)/U(0,sqrt2 z), w''=z^2 w.
# BRIEF CORRECTION (plan 02:63-68): NO real poles (recessive PCF U(0,.) is
# non-oscillatory).  Poles are a COMPLEX-conjugate pair near arg~=3pi/4,
# residue -1.  Pole oracle = root-find of w(z) = pcfu(0, sqrt2 z) = 0.
# ---------------------------------------------------------------------------
const CRic5_U0  = 0.675978240067284728995447685126052948591996091
const CRic5_UP0 = 0.456946581044463625374966623163157739126844321
const CRic5_P1 = complex(-1.4925974108469686254059802458433553078872661, 1.60304589241592527451540078382591864804891987)
const CRic5_U_a = complex(0.409473021516156916871843310856592541402087361, 0.0411884824551727349837268255373468549235369407)
const CRic5_U_b = complex(-0.172657164730854156308925912490337319165901944, -0.541491846505941151863289130791967693784851081)
const CRic5_U_c = complex(2.11942973654678410160049258174519915112698838, 2.01972229727285599626889794762420288104912393)
