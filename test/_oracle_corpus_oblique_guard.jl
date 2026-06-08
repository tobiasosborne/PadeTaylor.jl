# test/_oracle_corpus_oblique_guard.jl
#
# Pinned transcendental oracle for CBvx.4 (oblique-complex BVP/VectorBVP
# domain guard, bug padetaylor-53tu).  Every constant below is emitted by
# external/probes/corpus-oracles/oblique-guard/capture.py (mpmath dps=30,
# sympy residual==0 gate) — INDEPENDENT of the Julia routine under test
# (CLAUDE.md Rule 5).  Re-run: `python3 capture.py` in that probe dir.
#
# Segment z_a=0, z_b=1+i.  t*(z*) = (z* - (0.5+0.5i)) / (0.5+0.5i).
#   on-segment NODE      z_mid = 0.5 + 0.5i      ⇒ t* = 0       (real; for even
#                                                   N, t*=0 IS a node, so this
#                                                   point hits the coincident-
#                                                   node guard, not barycentry)
#   on-segment OFF-NODE  z_qrt = 0.65 + 0.65i    ⇒ t* = 0.3     (real; exercises
#                                                   the barycentric interpolant)
#   off-segment query    z_off = -0.75 + 2.25i   ⇒ t* = 0.5 + 3i
#                                                   (Re interior, Im large)
#
# ERRATUM (reported to caller): the bead text quoted
# cosh(1+i)=0.833…−0.988…i; the correct value has a POSITIVE imaginary part
# (cosh(1+i)=cosh1·cos1 + i·sinh1·sin1, both factors positive).  mpmath,
# Python cmath, and Julia all agree on +0.988….  OG_BC_B below is correct.

# scalar u = cosh(z)
const OG_BC_A     = 1.0 + 0.0im                                                   # u(z_a)=cosh(0)
const OG_BC_B     = 0.8337300251311490488838854 + 0.9888977057628650963821295im  # u(z_b)=cosh(1+i)
const OG_COSH_MID = 0.9895848833999199364440705 + 0.2498263975004615314895597im  # cosh(z_mid) node
const OG_COSH_QRT = 0.970261602208795953219371 + 0.4216621310419375486022782im   # cosh(z_qrt) off-node
const OG_COSH_OFF = -0.8132858892140288446917218 - 0.6398226084717389934521876im # cosh(z_off) (true)

# vector y = (cosh z, sinh z)  for the y'=(y2,y1) mirror
const OG_SINH_A   = 0.0 + 0.0im                                                   # y2(z_a)=sinh(0)
const OG_SINH_B   = 0.6349639147847361082550822 + 1.298457581415977294826042im   # sinh(1+i)
const OG_SINH_MID = 0.4573041531842492216075127 + 0.5406126857131533803537029im  # sinh(z_mid) node
const OG_SINH_QRT = 0.554669417228393882780027 + 0.7375971383096157387964337im   # sinh(z_qrt) off-node
const OG_SINH_OFF = 0.5165576805256538361923466 + 1.007358362265867162828912im   # sinh(z_off) (true)

# geometry literals
const OG_Z_A   = 0.0 + 0.0im
const OG_Z_B   = 1.0 + 1.0im
const OG_Z_MID = 0.5 + 0.5im    # t* = 0     (node)
const OG_Z_QRT = 0.65 + 0.65im  # t* = 0.3   (off-node, on-segment)
const OG_Z_OFF = -0.75 + 2.25im # t* = 0.5 + 3.0im
const OG_T_OFF = 0.5 + 3.0im
