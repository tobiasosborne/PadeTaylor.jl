# test/_oracle_corpus_pathnet_walls_rows.jl
#
# Pinned closed-form / mpmath oracle constants for the path-network
# walls/rows corpus (bead padetaylor-vt02, epic padetaylor-25og).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/pathnet-routing/capture.py  (mpmath 1.3,
#   mp.dps = 50; sympy residual==0 GATE — CPN.1 logistic recast, CPN.4 rational
#   forcing, CPN.2 Airy-Riccati — before any value is emitted).  Every
#   transcendental constant below is a verbatim copy of that script's labelled
#   output.  To update: re-run `python3 capture.py` and re-pin.  Ground-truth
#   recipe is docs/test_corpus/02_corpus_extension_plan.md:345-385 (Family I).
#
# Exact rationals (CPN.1 IC u(0)=1/3, u'(0)=2/9; CPN.4 ICs and past-pole
# rationals like u(0.5)=-11/3, u(1.5)=9/2; CPN.4 forcing constants) are written
# INLINE in the test, not here — only transcendental pins live here.

# ---------------------------------------------------------------------------
# CPN.1  logistic c=2:  u(z) = 1/(1 + 2 e^{-z}).  Vertical pole WALL at
# Re z = log 2, Im z = π(2k+1).  Targets {2±2i, -1±2i} straddle the wall.
# ---------------------------------------------------------------------------
const CPN1_LOG2 = 0.693147180559945309417232121458176568075500134
const CPN1_U_2p2i = complex(1.04643494543496863314996197483848007447504426, 0.290240988288702816778851119133311471159050082)
const CPN1_U_2m2i = complex(1.04643494543496863314996197483848007447504426, -0.290240988288702816778851119133311471159050082)
const CPN1_U_m1p2i = complex(-0.0484956028027716837817707134828209769879100136, 0.189903425906408701884139198011199811547985954)
const CPN1_U_m1m2i = complex(-0.0484956028027716837817707134828209769879100136, -0.189903425906408701884139198011199811547985954)

# ---------------------------------------------------------------------------
# CPN.8  logistic c=2 off-node dense-eval (interior, real-axis, pole-free).
# Spot values strictly between visited tree nodes.
# ---------------------------------------------------------------------------
const CPN8_U_1p17 = 0.617004439025104027009434010993258075720073155
const CPN8_U_0p83 = 0.534159907289591776473498279978214325253986481

# ---------------------------------------------------------------------------
# CPN.2  Airy-Riccati:  u(z) = -Ai'(z)/Ai(z),  IC z0 = 0.5.  Poles at the Airy
# zeros a_k on the NEGATIVE real axis.  IC pins + target values + zero row.
# ---------------------------------------------------------------------------
const CPN2_U0 = 0.97072394910167399736949306240730576569846182
const CPN2_UP0 = 0.442304985359549369484597437077338280832526764
const CPN2_U_m1 = 0.0189718245555636599559750467452062630588004464
const CPN2_U_m1p5 = -0.665982954300023111374311204206832333010023949
const CPN2_U_m3 = 0.830443239515891886122050086436360142483619195
const CPN2_U_m5 = -0.932808408393958909309433016790640458742659615
const CPN2_U_m2p0p5i = complex(-0.650657279527084317158529945153808032135614015, 1.68684317245392454504042461788955540507544529)

# Airy zeros a_k (the negative-axis pole row); k = 1..5.  z = 0.5 IC walking
# toward the negative axis must NOT land a visited node within eps of any.
const CPN2_AIRY_ZERO_1 = -2.33810741045976703848919725244673544063854015
const CPN2_AIRY_ZERO_2 = -4.0879494441309706166369887014573910602247647
const CPN2_AIRY_ZERO_3 = -5.52055982809555105912985551293129357379721428
const CPN2_AIRY_ZERO_4 = -6.78670809007175899878024638449617696605388248
const CPN2_AIRY_ZERO_5 = -7.9441335871208531231382805557982685321406744
