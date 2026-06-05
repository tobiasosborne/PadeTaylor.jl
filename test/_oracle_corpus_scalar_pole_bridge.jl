# test/_oracle_corpus_scalar_pole_bridge.jl
#
# Pinned closed-form / mpmath oracle constants for the scalar pole-bridge
# corpus (bead padetaylor-2tj5, epic padetaylor-p1v0).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/scalar-pole-bridge/capture.py  (mpmath
#   1.3, mp.dps = 50).  Every transcendental / Bessel / elliptic constant
#   below is a verbatim copy of that script's labelled output.  To update:
#   re-run `python3 capture.py` and re-pin.  Ground-truth recipes are from
#   docs/test_corpus/01_corpus_catalogue.md "Full candidate records".
#
# Closed-form constants (rational / exact algebraic) are written inline in
# the test, not here — only transcendental pins live in this file.
#
# CROSS-VALIDATION of the WP recipe against FW2011 Table 5.1 (eq 5.3,
# references/markdown/FW2011_painleve_methodology_JCP230/
# FW2011_painleve_methodology_JCP230.md:281-302):
#   capture.py WP_CHECK_U30  = 1.09509825595974419974...  ==  FW 1.095098255959744
#   capture.py WP_CHECK_U1e4 = 21.0253033947105499877...  ==  FW 21.02530339471055
# The recipe is therefore the published keystone oracle, not a self-pin.

# ---------------------------------------------------------------------------
# wp-equianharmonic-g0g2 : u'' = 6u^2, u(z) = WP(z-1; 0, 2)   (catalogue:281)
# IC at z=0 (catalogue ic_bc) and pole-bridge eval points around pole z=1.
# ---------------------------------------------------------------------------
const WP_U0       = 1.07182251641691742068161429481250292922453268
const WP_UP0      = 1.71033735317678620846433189149624950839231872
const WP_U_AT_0p5 = 4.00446466900308745041894518926888663277715402
const WP_U_AT_1p05  = 400.000000446428571466898057302094343353962712   # PAST pole
const WP_UP_AT_1p05 = -15999.9999642857142780489599678824003935939971
const WP_U_AT_1p4   = 6.25182861258210936858821123830667476852571992
# 80-dps strings for the BigFloat-256 path (catalogue precision: arbitrary).
const WP_U_AT_1p05_STR =
    "400.00000044642857146689805730209434335396271221358299108294803421544061930158100"

# ---------------------------------------------------------------------------
# tan : u' = 1+u^2 -> tan(z); recast 2nd-order u'' = 2u + 2u^3 (catalogue:417)
# IC u(0)=0, u'(0)=1; simple pole at z=pi/2; bridge crosses it.
# Exact at multiples of pi (tan(pi)=0); transcendental elsewhere.
# ---------------------------------------------------------------------------
const TAN_U_AT_2p0  = -2.18503986326151899164330610231368254343201775   # PAST pole
const TAN_UP_AT_2p0 = 5.77439920404191761241276760432323738242865200
const TAN_U_AT_2p5  = -0.747022297238660279355352687825274557904116957

# ---------------------------------------------------------------------------
# coth : u' = 1-u^2 -> coth(z); recast 2nd-order u'' = -2u + 2u^3 (cat:449)
# IC at z=-1; simple pole at z=0; bridge to z=+1 (odd symmetry u(1)=-u(-1)).
# ---------------------------------------------------------------------------
const COTH_U0_M1  = -1.31303528549933130363616124693084783291201394
const COTH_UP0_M1 = -0.724061660966310466407988789130590577259388341
const COTH_U_AT_1p0  = 1.31303528549933130363616124693084783291201394   # PAST pole
const COTH_UP_AT_1p0 = -0.724061660966310466407988789130590577259388341

# ---------------------------------------------------------------------------
# fw-riccati-bessel : u' = t^2+u^2 -> t*J_{3/4}(t^2/2)/J_{-1/4}(t^2/2);
# recast 2nd-order u'' = 2t + 2t^2 u + 2u^3 (catalogue:499).  IC u(0)=0,
# u'(0)=0; non-algebraic simple pole at t1 = 2.0031...; FW2011 §2.2.1 demo.
# ---------------------------------------------------------------------------
const FW_FIRST_POLE_T = 2.00314735942688470800461097905429922381014482
const FW_U_AT_0p5 = 0.0417911461546818632207688068491774981077532562
const FW_U_AT_1p5 = 1.51744754388000185172957981149904220354506183   # near pole
const FW_U_AT_2p1 = -10.1854755476848350240823121476236161797863928  # PAST pole

# ---------------------------------------------------------------------------
# tanh-complex-pole : u' = 1-u^2 -> tanh(z); recast u'' = -2u + 2u^3
# (catalogue:516).  IC u(0)=0, u'(0)=1; poles at i*(pi/2+k pi) OFF the real
# axis; imaginary-direction bridge crosses pole at i*pi/2.  tanh(i*pi/4)=i
# and tanh(i*3pi/4)=-i are EXACT (written inline as ±im in the test).
# ---------------------------------------------------------------------------
# (no transcendental pins; the bridge targets are exact ±im.)

# ---------------------------------------------------------------------------
# jacobi-sn-2ndorder : u'' = -(5/4)u + (1/2)u^3, u(z) = sn(z, m=1/4)
# (catalogue:465).  Exact algebraic ICs at z=i*Kp/2: u=i*sqrt(2),
# u'=3/sqrt(2); complex bridge through pole at i*Kp to z=i*3Kp/2 where
# u=-i*sqrt(2) (EXACT).  Kp = K'(1/4) = K(3/4).  Real-axis sn values smooth.
# ---------------------------------------------------------------------------
const SN_KP        = 2.15651564749964323543867499880032202886411022
const SN_U_AT_1p0  = 0.822635578129862359676230338653976488440647117  # real, no pole
const SN_UP_AT_1p0 = 0.518246096435056696352842932212043893602178190
