# test/_oracle_corpus_pathnet_winding.jl
#
# Pinned sympy/mpmath oracle constants for the path-network winding corpus
# (bead padetaylor-nkw8, epic padetaylor-25og).  Family I of
# docs/test_corpus/02_corpus_extension_plan.md:374-381 -- the HEADLINE
# padetaylor-61um catch under a REAL path_network_solve(cross_branch=true) walk.
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/pathnet-winding/capture.py  (mpmath 1.3,
#   mp.dps = 50; sympy residual==0 GATE on the y=sqrt(1-z^2) recast
#   f=(-1-yp^2)/y and the +-1 branch-point set, BEFORE any value is emitted;
#   a >pi / <pi geometry GATE on the two trigger legs).  Each transcendental
#   constant below is a verbatim copy of that script's labelled output.  To
#   update: re-run `python3 capture.py`, re-pin.
#
# Exact quantities (the IC y0=1 / y'0=0; the +-1 branch points) are written
# inline in the test where they are exact rationals/integers.
#
# CLOSED FORM:
#   CPN.7  y = sqrt(1 - z^2), C=1 ; y' = -z/y ; y'' = (-1 - (y')^2)/y.

# ---------------------------------------------------------------------------
# CPN.7 IC at z0=0: y0 = sqrt(1) = 1 (principal +), y'0 = -0/1 = 0.
# ---------------------------------------------------------------------------
const CPN7_Y0 = complex(1.0, 0.0)
const CPN7_YP0 = complex(0.0, 0.0)

# ---------------------------------------------------------------------------
# THE +-sqrt SHEET ORACLE at z = 1.5 + 0.5i (the SEMANTIC meaning of a sheet
# flip).  Principal y = +sqrt(1-z^2); after one CCW loop about z=1 the sqrt
# sheet FLIPS SIGN to -sqrt(1-z^2) (exactly the negation).  This is the
# INDEPENDENT closed-form context for what `sheet_index = +1` MEANS; it is
# NOT the per-node Pade output (single-step Pade across a near-branch leap is
# numerically poor -- 61um corrupts the BOOKKEEPING, which is what we pin).
# ---------------------------------------------------------------------------
const CPN7_YPRIN_AT_1p5_0p5i = complex(0.633551749161816554509756960322285244696908042, -1.1838022718621540655654375257135366817724956)
const CPN7_YFLIP_AT_1p5_0p5i = complex(-0.633551749161816554509756960322285244696908042, 1.1838022718621540655654375257135366817724956)

# ---------------------------------------------------------------------------
# THE 61um REAL-WALK TRIGGER GEOMETRY.  branch = 1; ring radius r = 0.45.
# Two single-wedge-step cross_branch=true legs differing only in step
# coarseness (the ideal-chord subtended angle straddles the pi threshold):
#
#   COARSE leg -- IC at arg +0.82*pi, target at arg +1.18*pi (ideal gap
#     +0.36*pi < pi).  winding_delta EXACT; visited_sheet=[1] CORRECT.
#   FINE leg -- IC at arg -0.47*pi, target at arg +0.59*pi (ideal gap
#     +1.06*pi > pi).  winding_delta WRAPS the +1.06*pi into (-pi,pi] and
#     returns a NEGATIVE value (loses a +2*pi revolution) => visited_sheet
#     bumps -1 not +1 -- the 61um corruption.  TRUE crossing is CCW (+1).
#
# The IDEAL-chord literals pin the ring geometry + the wrap arithmetic; the
# package's realised wedge step lands slightly off-target (wedge offsets), so
# the realised subtended angle is re-measured in Julia at test time.  The
# load-bearing oracle is the TRUE CROSSING DIRECTION (CCW => +1), robust to
# the wedge offset.
# ---------------------------------------------------------------------------
const CPN7_COARSE_P = complex(0.620052433524093214653148871214993322578253341, 0.241122057740548478222088945540437990128608906)
const CPN7_COARSE_Q = complex(0.620052433524093214653148871214993322578253341, -0.241122057740548478222088945540437990128608906)
const CPN7_COARSE_IDEAL_SUBTENDED = 1.13097335529232556584655161798062103831098098      # +0.36*pi

const CPN7_FINE_P = complex(1.04234874099333144331297008199846415038548268, -0.448002884071386005803955119896578738425669941)
const CPN7_FINE_Q = complex(0.874454002282346836666037097223697775627102488, 0.432132158554624382288430960220047881987620207)
const CPN7_FINE_IDEAL_SUBTENDED = 3.33008821280518083277040198627627305724899956        # +1.06*pi
const CPN7_FINE_IDEAL_WRAPPED = -2.95309709437440564415488478028273271114533924         # -0.94*pi
