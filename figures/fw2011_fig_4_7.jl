# figures/fw2011_fig_4_7.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 4.7 — six Painlevé-I
# pole-location plots, one per choice of initial condition at z = 0,
# each displayed over the square z = x + iy, -50 <= x,y <= 50.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:265-308, and the panel
#         IC labels read from the figure image _page_11_Figure_2.jpeg.
#
# FW md:265-271 (§4.5, "Transitions in pole field patterns"): far from
# the origin every PI solution settles into the five-fold sectorial
# structure of Fig. 1.1, but the local pattern near z = 0 can differ
# (e.g. a six-fold pattern for large ICs), so transitions must occur.
# The six panels illustrate the range:
#
#   (a)  u(0) =  0,                  u'(0) = 0
#   (b)  u(0) =  0,                  u'(0) = 1.8518          (= Fig. 4.3)
#   (c)  u(0) =  2,                  u'(0) = 2
#   (d)  u(0) =  4,                  u'(0) = 0
#   (e)  u(0) = -0.1875,             u'(0) = 0.3049          (= Fig. 3.1)
#   (f)  tritronquée: u(0) = -0.1875543083404949,
#                     u'(0) =  0.3049055602612289           (FW eq. 4.1)
#
# Panel (f) is the Boutroux tritronquée — pole-free in four of its five
# sectors (FW §2 lines 52-54), populated only in the wedge around the
# positive real axis.
#
# ----------------------------------------------------------------------
# Method — and why panel (f) is computed differently
#
# A pole-location plot is exactly what `PoleField.extract_poles` reads
# back out of a solved path-network: the roots of every visited node's
# stored Padé denominator, mapped to the z-plane (bead `padetaylor-xvf`,
# worklog 027).
#
# But the underlying solve cannot be a plain `path_network_solve` for
# every panel.  FW md:401 is explicit: "smooth regions are unstable
# regions for any IVP solver, and the associated loss of accuracy ought
# not be carried back into the pole field ... force the path selection
# algorithm to complete pole fields before stepping into smooth
# regions".  The tritronquée (f) has four *unbounded* pole-free
# sectors; a plain IVP integrated across them accumulates error that
# perturbs the solution off the tritronquée manifold, blooming spurious
# poles where there should be none (confirmed empirically — over
# [-50,50]² the plain solve's pole field is nearly angle-uniform).
#
# So every panel here uses `edge_gated_pole_field_solve` (bead
# `padetaylor-dmb`, worklog 028): the region-growing solve that starts
# from a seed at the IC and expands the IVP's target set only into
# cells the FW §3.2.2 edge detector confirms as pole-field.  For the
# generic panels (a, c, d) the pole field covers essentially the whole
# window, so the gated solve grows to fill it — equivalent to the plain
# solve, just reached safely.  For the panels with smooth regions (b,
# e, f) the gate keeps the IVP out of them, and the pole-free regions
# stay correctly empty.  This is FW's own "fully automated" workflow
# (md:401), uniform across all six panels.
#
# Step length h = 0.5 and Taylor order = 30 are FW 2011's working
# defaults (md:279).

using PadeTaylor
using CairoMakie
using Printf
include(joinpath(@__DIR__, "figutil.jl"))

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
const HALF   = 50.0              # window is [-HALF, HALF]^2 (FW md:308)
const N      = 101               # 101 x 101 lattice — spacing 1.0, fine
                                 # enough for the §3.2.2 edge detector to
                                 # separate pole-field from smooth (see
                                 # EdgeGatedSolve "Grid resolution matters")
const H_STEP = 0.5               # FW 2011 working step length (md:279)
const ORDER  = 30                # FW 2011 standard Taylor order (md:279)
const GROW   = 4                 # edge-gated region-growing ring width
const OUTPNG = joinpath(@__DIR__, "output", "fw2011_fig_4_7.png")

# PI right-hand side: u'' = 6 u^2 + z.
pI(z, u, up) = 6 * u^2 + z

# The six panels: (label, u(0), u'(0)).  Order matches FW's 3x2 layout.
const PANELS = [
    ("(a)  u(0) = 0,  u'(0) = 0",                 0.0,                 0.0),
    ("(b)  u(0) = 0,  u'(0) = 1.8518",            0.0,                 1.8518),
    ("(c)  u(0) = 2,  u'(0) = 2",                 2.0,                 2.0),
    ("(d)  u(0) = 4,  u'(0) = 0",                 4.0,                 0.0),
    ("(e)  u(0) = -0.1875,  u'(0) = 0.3049",      -0.1875,             0.3049),
    ("(f)  tritronquée (FW eq. 4.1)",             -0.1875543083404949, 0.3049055602612289),
]

# ----------------------------------------------------------------------
# Solve + extract poles, one panel at a time
# ----------------------------------------------------------------------
xs = range(-HALF, HALF; length = N)
ys = range(-HALF, HALF; length = N)

fig = Figure(; size = (980, 1380))

for (k, (label, u0, up0)) in enumerate(PANELS)
    @printf("FW 2011 Fig 4.7 panel %d/%d — %s\n", k, length(PANELS), label)
    prob = PadeTaylorProblem(pI, (u0, up0), (0.0, HALF); order = ORDER)

    t0  = time()
    egs = edge_gated_pole_field_solve(prob, xs, ys;
                                      h = H_STEP, order = ORDER, grow_rings = GROW)
    poles = extract_poles(egs.pn_solution)
    @printf("  %d region-growing passes, %d pole-field cells, %d visited nodes; %d poles extracted; %.1f s\n",
            egs.iterations, count(egs.field_mask),
            length(egs.pn_solution.visited_z), length(poles), time() - t0)

    row = (k + 1) ÷ 2          # 1,1,2,2,3,3
    col = isodd(k) ? 1 : 2     # 1,2,1,2,1,2
    pole_scatter_axis(fig[row, col], poles;
                      title = label, xlim = (-HALF, HALF), ylim = (-HALF, HALF),
                      markersize = 1.7)
end

Label(fig[0, :], "FW 2011 Fig. 4.7 — Painlevé-I pole fields over [-50, 50]²";
      fontsize = 17, font = :bold)

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
