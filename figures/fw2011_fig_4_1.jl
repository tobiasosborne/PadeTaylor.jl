# figures/fw2011_fig_4_1.jl
#
# Reproduces Fornberg & Weideman 2011, Fig. 4.1 — a near-tritronquée
# Painlevé-I solution `|u(z)|` over z = x + iy, x ∈ [-10,10],
# y ∈ [-10,30], assembled by FW's three-step BVP+IVP composition.
#
# Source: references/markdown/FW2011_painleve_methodology_JCP230/
#         FW2011_painleve_methodology_JCP230.md:214-229, and the figure
#         image _page_7_Figure_2.jpeg.
#
# FW md:220 — "Fig. 4.1.  A near-tritronquée case calculated by the
# steps (i) use a BVP solver over [-20i,20i] to obtain values for u
# and u' at z = 0 and z = 20i  (ii) use this data to separately run
# out the two pole fields  (iii) fill in the area between pole field
# edges by further BVP solutions."
#
# The figure is a `|u(z)|` surface: a vast smooth plain with two
# localised pole fields — one run out from z = 0 (the tritronquée
# wedge around the positive real axis), one from z = 20i (the "upper
# pole field", which for the exact tritronquée would sit at infinity,
# FW md:216) — and the solid imaginary-axis line marking the step-(i)
# BVP segment, with the points z = 0 and z = 20i where u, u' are read
# off.
#
# ----------------------------------------------------------------------
# The three steps, mapped to this package
#
#   (i)  `bvp_solve` on the segment [-20i, 20i] with the leading-term
#        asymptote u(z) = -√(-z/6) (FW eq. 1.2, negative branch — the
#        tritronquée selector, see test/fw_fig_41_test.jl) as Dirichlet
#        BCs.  This is the recipe pinned in `test/fw_fig_41_test.jl`
#        (bead `padetaylor-0c3`).  It yields u(0), u'(0) — essentially
#        the exact tritronquée ICs (FW eq. 4.1) — and u(20i), u'(20i),
#        plus the smooth solution along the whole imaginary axis.
#
#   (ii) Two `edge_gated_pole_field_solve` runs (bead `padetaylor-dmb`,
#        worklog 028): the IVP path-network confined to the pole field
#        so it never wanders into — and is corrupted by — the smooth
#        plain (FW md:401).  One run from the IC (0, u(0), u'(0)), one
#        from (20i, u(20i), u'(20i)).  Each returns the pole field's
#        `|u|` on its cells plus a `field_mask`.
#
#   (iii) The smooth region — every cell in neither pole field — is
#        filled per grid row by `bvp_solve`, FW's "one BVP per grid
#        line" (md:190).  Each smooth run is bridged between its
#        flanking anchors: a pole-field-edge cell's IVP value, the
#        step-(i) BVP value where the row crosses the imaginary axis,
#        or the leading-term asymptote at a grid edge.  A run whose
#        BVP fails to converge falls back to the leading term (counted
#        and reported — there should be none on this geometry).
#
# Step length h = 0.5, Taylor order = 30 — FW 2011 working defaults
# (md:279); BVP node count N tuned per segment length.

using PadeTaylor
using CairoMakie
using Printf
include(joinpath(@__DIR__, "figutil.jl"))

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
const XHALF  = 10.0              # x ∈ [-XHALF, XHALF]            (FW md:220)
const YLO    = -10.0             # y ∈ [YLO, YHI]                 (FW md:220)
const YHI    =  30.0
const STEP   = 0.5               # lattice step (fine enough for the
                                 # §3.2.2 edge detector — EdgeGatedSolve
                                 # "Grid resolution matters")
const H_STEP = 0.5               # FW 2011 working step length (md:279)
const ORDER  = 30                # FW 2011 standard Taylor order (md:279)
const Z_BVP  = 20.0              # step-(i) BVP segment is [-Z_BVP·i, Z_BVP·i]
const N_AXIS = 240               # step-(i) BVP node count (worklog 016)
const U_CAP  = 40.0              # |u| z-axis cap, matching FW's Fig 4.1 framing
const OUTPNG = joinpath(@__DIR__, "output", "fw2011_fig_4_1.png")

# Painlevé I.  IVP form u'' = f(z, u, u'); BVP form u'' = f(z, u) with
# analytic ∂f/∂u (the Chebyshev–Newton solver takes the Jacobian
# directly — no autodiff).
pI_ivp(z, u, up) = 6 * u^2 + z
bvp_f(z, u)      = 6 * u^2 + z
bvp_∂f_∂u(z, u)  = 12 * u

# Leading-term asymptote, negative branch (FW eq. 1.2; the tritronquée
# selector — see test/fw_fig_41_test.jl FF.1.6).
leading(z) = -sqrt(-z / 6)

# Index of the lattice axis entry nearest a coordinate.
nearest_idx(v, x) = argmin(abs.(collect(v) .- x))

# ----------------------------------------------------------------------
# Step (i) — the imaginary-axis BVP
# ----------------------------------------------------------------------
@printf("FW 2011 Fig 4.1 — step (i): BVP on [-%gi, %gi] (N = %d)\n",
        Z_BVP, Z_BVP, N_AXIS)
z_a = -Z_BVP * im
z_b =  Z_BVP * im
spine = bvp_solve(bvp_f, bvp_∂f_∂u, z_a, z_b, leading(z_a), leading(z_b);
                  N = N_AXIS, initial_guess = leading, tol = 1e-13, maxiter = 20)

u0,  up0  = spine(0.0 + 0.0im)
u20, up20 = spine(Z_BVP * im)
@printf("  u(0)  = %.12f %+.12fim   u'(0)  = %.12f %+.12fim\n",
        real(u0), imag(u0), real(up0), imag(up0))
@printf("  u(20i)= %.6f %+.6fim    u'(20i)= %.6f %+.6fim\n",
        real(u20), imag(u20), real(up20), imag(up20))

# ----------------------------------------------------------------------
# The shared lattice
# ----------------------------------------------------------------------
xs = range(-XHALF, XHALF; step = STEP)
ys = range(YLO,    YHI;   step = STEP)
nx, ny = length(xs), length(ys)
i_axis = findfirst(x -> abs(x) < STEP / 4, xs)        # column with x ≈ 0
@assert i_axis !== nothing "lattice must contain the x = 0 column"
@printf("  lattice %d × %d, step %g\n", nx, ny, STEP)

# ----------------------------------------------------------------------
# Step (ii) — run out the two pole fields with the edge-gated IVP
# ----------------------------------------------------------------------
@printf("FW 2011 Fig 4.1 — step (ii): edge-gated pole-field run-outs\n")
prob0 = PadeTaylorProblem(pI_ivp, (ComplexF64(u0),  ComplexF64(up0)),
                          (0.0 + 0.0im, ComplexF64(XHALF)); order = ORDER)
prob20 = PadeTaylorProblem(pI_ivp, (ComplexF64(u20), ComplexF64(up20)),
                           (Z_BVP * im, ComplexF64(XHALF)); order = ORDER)

t0 = time()
egs0 = edge_gated_pole_field_solve(prob0, xs, ys; h = H_STEP, order = ORDER,
                                   grow_rings = 4)
@printf("  pole field from z=0  : %d passes, %d field cells, %.1f s\n",
        egs0.iterations, count(egs0.field_mask), time() - t0)
t0 = time()
egs20 = edge_gated_pole_field_solve(prob20, xs, ys; h = H_STEP, order = ORDER,
                                    grow_rings = 4)
@printf("  pole field from z=20i: %d passes, %d field cells, %.1f s\n",
        egs20.iterations, count(egs20.field_mask), time() - t0)

# ----------------------------------------------------------------------
# Step (iii) — per-row BVP fill of the smooth region
# ----------------------------------------------------------------------
# `uval[i,j]` accumulates the stitched solution.  Seed it with the two
# IVP run-outs (whichever is finite at a cell); everything still NaN is
# a smooth cell to be BVP-bridged.
uval = fill(ComplexF64(NaN, NaN), nx, ny)
for j in 1:ny, i in 1:nx
    if isfinite(egs0.u_grid[i, j])
        uval[i, j] = egs0.u_grid[i, j]
    elseif isfinite(egs20.u_grid[i, j])
        uval[i, j] = egs20.u_grid[i, j]
    end
end

@printf("FW 2011 Fig 4.1 — step (iii): per-row smooth-band BVP fill\n")
t0 = time()
# Wrapped in `let` so the run counters are loop-local — a top-level
# `for` would hit Julia's soft-scope ambiguity.
n_runs, n_fallback = let nr = 0, nf = 0
    for j in 1:ny
        y = ys[j]
        # Row anchors: the IVP run-outs above, plus the step-(i) BVP
        # spine value where this row crosses the imaginary axis
        # (|y| ≤ Z_BVP).
        if isnan(real(uval[i_axis, j])) && abs(y) ≤ Z_BVP
            uval[i_axis, j] = spine(im * y)[1]
        end
        known(i) = isfinite(uval[i, j])

        i = 1
        while i ≤ nx
            if known(i)
                i += 1
                continue
            end
            # Maximal run of unknown (smooth) cells [run_lo, run_hi].
            run_lo = i
            while i ≤ nx && !known(i)
                i += 1
            end
            run_hi = i - 1

            # Bracket the run with endpoints: a flanking known cell, or
            # the grid edge with the leading-term BC.
            iL = run_lo > 1  ? run_lo - 1 : run_lo
            iR = run_hi < nx ? run_hi + 1 : run_hi
            z_L = xs[iL] + im * y
            z_R = xs[iR] + im * y
            u_L = known(iL) ? uval[iL, j] : leading(z_L)
            u_R = known(iR) ? uval[iR, j] : leading(z_R)

            nr += 1
            seglen = abs(z_R - z_L)
            Nrun = clamp(round(Int, 8 * seglen), 16, 200)
            filled = false
            try
                rsol = bvp_solve(bvp_f, bvp_∂f_∂u, z_L, z_R, u_L, u_R;
                                 N = Nrun, initial_guess = leading, maxiter = 25)
                for k in run_lo:run_hi
                    uval[k, j] = rsol(xs[k] + im * y)[1]
                end
                filled = true
            catch err
                err isa ErrorException || rethrow(err)
            end
            if !filled
                # BVP did not converge on this run — fall back to the
                # leading-term asymptote (FW md:229: "crude" but smooth).
                nf += 1
                for k in run_lo:run_hi
                    uval[k, j] = leading(xs[k] + im * y)
                end
            end
        end
    end
    (nr, nf)
end
@printf("  %d smooth runs BVP-filled (%d leading-term fallbacks); %.1f s\n",
        n_runs, n_fallback, time() - t0)

# ----------------------------------------------------------------------
# Render — |u(z)| surface + the step-(i) BVP segment line
# ----------------------------------------------------------------------
Z = clamp.(abs.(uval), 0.0, U_CAP)
# Domain is twice as long in y as in x; stretch the y footprint to match.
fig, ax = pole_field_figure(xs, ys, Z;
    title = "FW 2011 Fig. 4.1 — near-tritronquée PI via BVP+IVP+BVP composition",
    ucap = U_CAP, aspect = (1.0, 2.0, 0.5), size = (900, 1100))

# The solid imaginary-axis line marks the step-(i) BVP segment; the two
# crosses sit at z = 0 and z = 20i where u, u' are read off (FW md:218).
yline = collect(range(max(YLO, -Z_BVP), min(YHI, Z_BVP); length = 200))
zline = [Z[i_axis, nearest_idx(ys, y)] for y in yline]
lines!(ax, zeros(length(yline)), yline, zline; color = :black, linewidth = 3)
for ymark in (0.0, Z_BVP)
    j = nearest_idx(ys, ymark)
    scatter!(ax, [0.0], [ymark], [Z[i_axis, j]];
             color = :black, marker = :cross, markersize = 14)
end

mkpath(dirname(OUTPNG))
save(OUTPNG, fig)
@printf("  wrote %s\n", OUTPNG)
