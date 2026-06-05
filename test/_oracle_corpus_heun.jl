# test/_oracle_corpus_heun.jl
#
# Pinned oracle constants for the Heun + special-function corpus
# (bead padetaylor-h5ho, epic padetaylor-p1v0).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/heun/capture.py  (mpmath 1.3, mp.dps = 50)
#   — an INDEPENDENT from-scratch re-implementation of the HeunG / HeunC
#   Frobenius three-term recurrences (DLMF 31.3.2-31.3.4 / Motygin 2015
#   eqs 3-4).  It does NOT call the Julia heun_g / heun_c (that would be
#   self-referential).  Every Heun NUMERIC value below was additionally
#   cross-checked against Mathematica's native HeunG / HeunC via
#   wolframscript (commands printed at the bottom of capture.py and run
#   on 2026-06-05; see the per-block notes).
#
# CLOSED-FORM REDUCTIONS verified BEFORE use (Law 1 — verify the reduction
# identity first, never trust the catalogue's numeric):
#   * HeunG(2,2,1,2,3/2,1;z) == 2F1(1/2,1;3/2;2z-z^2)   (|diff|<1e-44 @ z<=0.7)
#   * HeunG(2,4,1,4,2,2;z)   == 1/(1-z), all c_n == 1   (|diff|=0 @ z<=0.5)
#   * HeunG(2,q*,-1,2,3,4;z) == 1 + (q*/6) z,  q* = -6+-2sqrt(6),  c_2 == 0
#   * HeunG(a,0,0,b,g,d;z)   == 1   (c_1 = q/(a g) = 0; all higher c_n = 0)
#   * F(1/2,3;3;z)           == (1-z)^(-1/2)            (DLMF 15.4.6; |diff|=0)
# All five identities reproduced by the recurrence to >=44 digits — they are
# genuine reductions, not catalogue lore.
#
# ERRATA corrected here (docs/test_corpus/ERRATA.md), independently re-derived:
#   * heunc-epsilon-neg3 @ z=0.9: catalogue -10.40476 is a TRUNCATION artefact;
#     converged value (N~300, Wolfram-confirmed) = -13.5052005053212.
#   * heung-complex-a @ z=0.3: catalogue real part wrong in the 7th figure;
#     exact-rational a=1/2+i/2 gives 1.10598130983985... (Wolfram-confirmed).
#
# Constants are 40-digit strings (Float64 parses exact); the Float64 literal
# is what the test compares against (heun_g/heun_c are Float64-internally).

# --- CHN.1  heung-2f1-reduction : HeunG(2,2,1,2,3/2,1;z) = 2F1(1/2,1;3/2;2z-z^2)
const HEUNG_2F1_Z0p3 = 1.254074179264320151832734154415091038561
const HEUNG_2F1_Z0p5 = 1.520691992601892695062188509760811648403
const HEUNG_2F1_Z0p7 = 1.964297346984075972059722556526677935491

# --- CHN.3a heung-polynomial-degree1 : Hp_{1,m}(z) = 1 + (q*/6) z
const HEUNG_POLY_QP     = -1.101020514433643803605431850588217216068   # -6+2sqrt6
const HEUNG_POLY_QM     = -10.89897948556635619639456814941178278393    # -6-2sqrt6
const HEUNG_POLY_P_Z0p5 = 0.908248290463863016366214012450981898661
const HEUNG_POLY_P_Z0p9 = 0.8348469228349534294591852224117674175898
const HEUNG_POLY_M_Z0p5 = 0.09175170953613698363378598754901810133901

# --- CHN.4a heung-complex-a : HeunG(1/2+i/2,1,1,2,3,4;0.3) [ERRATA-corrected]
const HEUNG_COMPLEX_A_Z0p3 =
    1.105981309839851012572693554038350726328 -
    0.1059215738336126137955457572831439266118im

# --- CHN.4b heung-small-a : HeunG(1/2,1,1,2,3,4;z)  (|a|<1, R = 1/2)
const HEUNG_SMALL_A_Z0p1  = 1.068937255689352422622185520186435937921
const HEUNG_SMALL_A_Z0p3  = 1.210685750907479745041848259922301720356
const HEUNG_SMALL_A_Z0p45 = 1.28554731854205577005880510232445487146

# --- CHN.4c heung-large-q-a3 : HeunG(3,10,1,2,2,3;z)  (|a|>1, large q)
const HEUNG_LARGE_Q_Z0p25 = 1.625058391255755947175928954899294589489
const HEUNG_LARGE_Q_Z0p5  = 3.28667048162199769691378416516780620348

# --- CHN.6  heunc-epsilon-neg3 / -pos2 / -fractional-gamma
const HEUNC_NEG3_Z0p1  = 0.945254571984878077364473848951811165683
const HEUNC_NEG3_Z0p5  = 0.4856076424570184331703870753157375054546
const HEUNC_NEG3_Z0p9  = -13.50520050532119415919845708316748397631   # ERRATA
const HEUNC_POS2_Z0p1  = 0.9498651796735502930663803794116258596603
const HEUNC_POS2_Z0p5  = 0.7227641496996699014792923292954837401356
const HEUNC_FRACG_Z0p1 = 0.776948521236988019842405997138506560981
const HEUNC_FRACG_Z0p5 = -1.34971152247169554130517865974550709304

# --- CSF special-function IVP oracle constants (mpmath natives, dps=50) ---
# Airy: u'' = z u, u = Ai(z)
const AIRY_AI0   = 0.355028053887817239260063186004183176398
const AIRY_AIP0  = -0.2588194037928067984051835601892039634791
const AIRY_Z1    = 0.1352924163128814155241474235154663061749
const AIRY_Z2    = 0.03492413042327437913532208079180760976106
const AIRY_Z3    = 0.006591139357460719144257448407961351071733
const AIRY_ZM1   = 0.5355608832923521187995165656388747074669
const AIRY_ZM2   = 0.2274074282016855759919244360378737994608
# Bessel J0: u'' + u'/z + u = 0, u = J_0(z); IC at z=1
const BESSEL_J0_1 = 0.7651976865579665514497175261026632209093
const BESSEL_J1_1 = 0.4400505857449335159596822037189149131274   # -J0'(1)
const BESSEL_J0_2 = 0.2238907791412356680518274546499486258252
const BESSEL_J0_3 = -0.2600519549019334376241546959773314368196
const BESSEL_J0_5 = -0.1775967713143383043473970130747587110711
# Gauss 2F1 IVP: z(1-z)u'' + [5/4-(11/6)z]u' - (1/6)u = 0, u = 2F1(1/2,1/3;5/4;z)
const GAUSS2F1_Z0p4 = 1.065978090961767260771336016200299294642
const GAUSS2F1_Z0p7 = 1.147425312518162185641645278974814224702
const GAUSS2F1_Z0p9 = 1.254059730476018748482012530394667385029
# IC at z0=0.05 (z=0 is a regular singular point of the divided-through ODE):
const GAUSS2F1_U0_005  = 1.006819404546595452616772785168014896833
const GAUSS2F1_UP0_005 = 0.139537924530179305922656121376982143256
# 2F1 elementary: F(1/2,3;3;z) = (1-z)^(-1/2); IC at z0=0.05
const F2F1ELEM_U0_005  = 1.025978352085154095456675073673642136763
const F2F1ELEM_UP0_005 = 0.5399886063606074186614079335124432298755
