"""
F2 tests for `PadeTaylor.VectorPathNetwork`'s **Stage-2 fine-grid fill**
(bead `padetaylor-0ln.20`; v0.2 plan row F2) and the **B1 true-radius
validity gate** (bead `padetaylor-0ln.37.5`; ADR-0025 Amendment 1).

V7 (`vector_path_network_test.jl`) covers the Stage-1 visited-tree
walk.  This file covers the Stage-2 fine-grid fill: when
`vector_path_network_solve` is called with a `fine_grid` kwarg, the
solve evaluates each node's stored shared-`Q` Padé densely over that
grid, producing the `grid_z`/`grid_y` solution fields — a filled
heatmap of `y(z)`, the substrate for the `P_I⁽²⁾` tritronquée surface
figure.  Stage 2 is the direct vector lift of v0.1's
`PathNetwork.path_network_solve` Stage-2 loop (`src/PathNetwork.jl:669-693`).

## The B1 true-radius gate (ADR-0025 Amendment 1)

Through V0.1 the Stage-2 fill gated each visited node's validity at the
*step size* `h_v` (`abs(z_f−z_v) > h_v ⇒ NaN`).  The Phase-A1 spike
(`external/probes/pade-validity-radius/REPORT.md`) measured that `h_v`
is **not** an honesty guarantee — it over-estimates the true validity
radius for 36 % of nodes.  ADR-0025 Lever 1 replaces it with the
per-node **true validity radius**:

    R_gate = min( s(tol)·h_JZ·h_v ,  0.5·h_v·min|t*| )

where `h_JZ = vector_step_jorba_zou(rescaled order-`order` jet, tol)` is
the project's own per-node Jorba–Zou truncation-radius estimator and
`min|t*|` is the smallest-magnitude root of the node's shared `Q`.  The
gate semantics in `VPN.2.*` below have been **updated** from the old
`h_v` gate to this true-radius gate — legitimate semantic evolution per
ADR-0025, not tolerance-relaxation; the new assertions test against
*computed* reference radii.

These tests (`VPN.2.*`) assert an invariant against a known-correct
value (Rule 5):

    VPN.2.1  closed-form oracle  — the d=3 Riccati system y_i' = -y_i²/c_i
                                   has the EXACT solution y_i(z) = c_i/(z−p).
                                   Stage-2-fill a fine grid; assert every
                                   filled cell inside a node's disc matches
                                   the closed form to local Padé accuracy.
    VPN.2.2  harmonic oscillator — the d=2 pole-free system y' = (y₂, −y₁)
                                   with IC (0,1) has the exact solution
                                   (sin z, cos z).  Stage-2-fill a real
                                   sub-interval; assert the fill matches
                                   (sin, cos).
    VPN.2.3  fail-soft / extrapolate — a grid point far from every visited
                                   node returns a d-vector of NaN (no
                                   throw); with extrapolate=true the same
                                   point evaluates to a finite vector.
    VPN.2.4  node-consistency    — Stage-2 at a grid point equal to a
                                   visited node reproduces that node's
                                   `visited_y` (t = 0 ⇒ Pᵢ(0)/Q(0)).
    VPN.2.5  fail-fast / empties — empty `fine_grid` throws; no `fine_grid`
                                   leaves the Stage-2 fields empty.
    VPN.2.6  B1 true-radius gate — the per-node `R_gate` is honoured: a
                                   query just inside vs just outside
                                   `R_gate` is finite vs NaN; the
                                   pole-adjacency clamp binds when a
                                   Q-root is near; the cached `visited_jets`
                                   match the no-cache fallback; degenerate
                                   inputs throw (Rule 1).
    VPN.2.7  B4 Voronoi honesty   — the retained nearest-node Voronoi fill
                                   (ADR-0025 Amendment 2, bead B4) is
                                   honest by construction: every covered
                                   (non-NaN) fill cell is evaluated by a
                                   node whose verified disc `R_gate`
                                   STRICTLY contains it — the Padé is
                                   never extrapolated.

## The B4 audition — Voronoi retained over a blend (ADR-0025 Amendment 2)

Bead `padetaylor-0ln.37.8` auditioned the nearest-node Voronoi fill
against a distance-weighted blend of all overlapping in-disc Padé
evaluations.  The audition
(`external/probes/stage2-blend-audition/audition.jl`) measured, on the
*actual* `P_I⁽²⁾` tritronquée wedge underlay, that the B1 true-radius
discs barely overlap (98.2 % of covered pixels lie in exactly one
node's disc; only 14.9 % of Voronoi seams lie inside a real overlap;
the inter-disc disagreement at the overlaps is `~6·10⁻¹⁰` — machine
precision).  A blend would change `~1.8 %` of a `~5 %` underlay by
`~10⁻⁹` — a no-op.  The audition therefore **retains** the Voronoi
fill (a well-evidenced retain is a legitimate audition outcome,
Rule 9).  `VPN.2.7` pins the honesty invariant that motivated the
audition's "the blend cannot be more honest than what we have"
conclusion: the retained fill never evaluates a Padé outside its disc.

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/vector_path_network_stage2_test.jl`) and under
`runtests.jl`.
"""

using Test
using PadeTaylor
using PadeTaylor.VectorPathNetwork: vector_path_network_solve,
                                    VectorPathNetworkSolution
using PadeTaylor.VectorPathNetworkStage2: _validity_radius, _safety_factor
using PadeTaylor.VectorCoefficients: vector_taylor_coefficients

# -----------------------------------------------------------------------------
# Oracle A — the d-component Riccati pole system (the V7 oracle, reused).
#
# y_i' = -y_i²/c_i  has the exact solution  y_i(z) = c_i/(z-p).
# Given IC y_i(z0) = c_i/(z0-p), the pole is p = z0 - c_i/y_i(z0).
# This is the figure-relevant case: a pole-rich vector meromorphic
# solution whose components all blow up at the SAME z = p.
# -----------------------------------------------------------------------------
function riccati_pole_problem_s2(p, z0, cs; order = 24)
    f  = (z, y) -> [-(y[i]^2) / cs[i] for i in eachindex(y)]
    y0 = [cs[i] / (z0 - p) for i in eachindex(cs)]
    return VectorPadeTaylorProblem(f, y0, (z0, z0 + 1.0 + 0im); order = order)
end
riccati_exact(z, p, cs) = [cs[i] / (z - p) for i in eachindex(cs)]

# -----------------------------------------------------------------------------
# Oracle B — the d=2 harmonic-oscillator system y' = (y₂, -y₁).
#
# With IC y(0) = (0, 1) the exact solution is y(z) = (sin z, cos z):
#   y₁' = y₂ = cos z = (sin z)',   y₂' = -y₁ = -sin z = (cos z)'.
# Pole-free — the whole plane is in some node's disc of validity, so
# this oracle has no fail-soft gaps to reason around.
# -----------------------------------------------------------------------------
function harmonic_problem(z0, z1; order = 24)
    f  = (z, y) -> [y[2], -y[1]]
    y0 = ComplexF64[sin(z0), cos(z0)]
    return VectorPadeTaylorProblem(f, y0, (z0, z1); order = order)
end
harmonic_exact(z) = ComplexF64[sin(z), cos(z)]

@testset "VectorPathNetwork Stage-2 fill (VPN.2.*)" begin

    # -------------------------------------------------------------------------
    # VPN.2.1 — closed-form oracle.  Riccati y_i(z) = c_i/(z-p); the
    # Stage-2 fill must reproduce it to local Padé accuracy on every
    # grid cell that lands inside some visited node's disc.
    # -------------------------------------------------------------------------
    @testset "VPN.2.1 closed-form Riccati oracle" begin
        p   = 1.0 + 0.6im
        z0  = -0.4 + 0.0im
        cs  = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem_s2(p, z0, cs; order = 24)

        # Targets BRACKET the pole `p` without sitting on it: the B2
        # adaptive walk caps `h` toward zero near a pole and cannot
        # honestly walk *onto* one (a target equal to `p` collapses `h`
        # below `h_min`, Rule 1).  The offset grid still surrounds `p`.
        targets = ComplexF64[x + y * im
                             for x in 0.0:0.4:1.6 for y in -0.4:0.4:1.2]
        # The fine grid: a denser lattice over the same region.
        fine = ComplexF64[x + y * im
                          for x in 0.0:0.2:1.6 for y in -0.2:0.2:1.2]

        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.25,
                                        fine_grid = fine)

        # The Stage-2 fields are populated and shaped correctly.
        @test length(sol.grid_z) == length(fine)
        @test length(sol.grid_y) == length(fine)
        @test sol.grid_z == fine

        # Every covered cell (a finite, non-NaN fill) must match the
        # closed form.  The min-‖y‖ walk steers around the pole p, so
        # cells very close to p may be NaN (uncovered) or large-error;
        # we assert accuracy only on covered cells comfortably away
        # from the singularity.
        n_checked = 0
        for (zf, yf) in zip(sol.grid_z, sol.grid_y)
            all(isfinite, yf) || continue          # skip fail-soft NaN cells
            abs(zf - p) > 0.35 || continue         # skip near-pole cells
            exact = riccati_exact(zf, p, cs)
            err = maximum(abs.(yf .- exact))
            @test err < 1.0e-6
            n_checked += 1
        end
        # The fill must actually have covered a substantial number of
        # cells — otherwise "every covered cell matches" is vacuous.
        @test n_checked ≥ 20
    end

    # -------------------------------------------------------------------------
    # VPN.2.2 — harmonic-oscillator oracle.  Pole-free: (sin, cos) is the
    # exact solution; every Stage-2 cell is covered and must match it.
    # -------------------------------------------------------------------------
    @testset "VPN.2.2 harmonic-oscillator oracle" begin
        z0  = 0.0 + 0.0im
        z1  = 2.0 + 0.0im
        prob = harmonic_problem(z0, z1; order = 24)

        # Targets + fine grid along the real axis (the solution is
        # entire, so a real-axis walk is well-posed).
        targets = ComplexF64[x + 0.0im for x in 0.3:0.3:1.8]
        fine    = ComplexF64[x + 0.0im for x in 0.1:0.1:1.7]

        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.3,
                                        fine_grid = fine)

        @test length(sol.grid_y) == length(fine)
        # Pole-free ⇒ every cell within h of a node ⇒ no NaN gaps for
        # this h-and-target choice; assert the fill matches (sin, cos).
        n_checked = 0
        for (zf, yf) in zip(sol.grid_z, sol.grid_y)
            all(isfinite, yf) || continue
            exact = harmonic_exact(zf)
            @test abs(yf[1] - exact[1]) < 1.0e-8       # sin component
            @test abs(yf[2] - exact[2]) < 1.0e-8       # cos component
            n_checked += 1
        end
        @test n_checked ≥ 14
    end

    # -------------------------------------------------------------------------
    # VPN.2.3 — fail-soft NaN vs extrapolate=true.
    # -------------------------------------------------------------------------
    @testset "VPN.2.3 fail-soft / extrapolate" begin
        z0  = 0.0 + 0.0im
        z1  = 2.0 + 0.0im
        prob = harmonic_problem(z0, z1; order = 24)
        targets = ComplexF64[x + 0.0im for x in 0.3:0.3:1.5]

        # A grid point planted FAR from any plausible visited node — the
        # walk never goes near 12+12im, so this cell is uncovered.
        far  = 12.0 + 12.0im
        near = 0.6 + 0.0im                  # comfortably inside a node disc
        fine = ComplexF64[near, far]

        # Default (extrapolate=false): the far cell is a d-vector of NaN,
        # the near cell is finite — fail-soft, NOT a throw (ADR-0015 / Rule 1).
        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.3,
                                        fine_grid = fine)
        @test all(isfinite, sol.grid_y[1])             # near: covered
        @test all(z -> isnan(real(z)) && isnan(imag(z)), sol.grid_y[2])
        @test length(sol.grid_y[2]) == 2               # still a d-vector

        # extrapolate=true: the B1 true-radius gate is skipped entirely,
        # so the far cell now evaluates to a finite (if inaccurate)
        # vector — the ADR-0015 escape hatch, retained by ADR-0025 B1.
        sol_x = vector_path_network_solve(prob, targets; order = 24, h = 0.3,
                                          fine_grid = fine, extrapolate = true)
        @test all(isfinite, sol_x.grid_y[2])
        # Inside-disc cells are identical between the two modes — the
        # gate only changes behaviour for out-of-disc cells.
        @test sol.grid_y[1] == sol_x.grid_y[1]
    end

    # -------------------------------------------------------------------------
    # VPN.2.4 — node-consistency.  A grid point placed exactly on a
    # visited node must reproduce that node's `visited_y`: t = 0, and
    # the canonical shared-Q approximant satisfies Pᵢ(0)/Q(0) = yᵢ.
    # -------------------------------------------------------------------------
    @testset "VPN.2.4 node-consistency" begin
        z0  = 0.0 + 0.0im
        z1  = 2.0 + 0.0im
        prob = harmonic_problem(z0, z1; order = 24)
        targets = ComplexF64[x + 0.0im for x in 0.3:0.3:1.5]

        # First walk with no fine grid to learn the visited nodes.
        sol0 = vector_path_network_solve(prob, targets; order = 24, h = 0.3)
        @test isempty(sol0.grid_z)                     # Stage-2 fields empty
        @test isempty(sol0.grid_y)

        # Re-solve with the visited-node positions themselves as the
        # fine grid.  Each grid point z = visited_z[k] is its own
        # nearest node ⇒ t = 0 ⇒ the fill must equal visited_y[k].
        node_grid = copy(sol0.visited_z)
        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.3,
                                        fine_grid = node_grid)
        @test length(sol.grid_y) == length(sol.visited_z)
        for k in eachindex(sol.visited_z)
            yk    = sol.visited_y[k]
            fill_k = sol.grid_y[k]
            @test maximum(abs.(fill_k .- yk)) < 1.0e-10
        end
    end

    # -------------------------------------------------------------------------
    # VPN.2.5 — fail-fast on an empty fine_grid; no-fine_grid invariant.
    # -------------------------------------------------------------------------
    @testset "VPN.2.5 fail-fast / no-grid invariant" begin
        z0  = 0.0 + 0.0im
        z1  = 2.0 + 0.0im
        prob = harmonic_problem(z0, z1; order = 24)
        targets = ComplexF64[x + 0.0im for x in 0.3:0.3:1.5]

        # An empty fine_grid is a caller mistake — pass nothing to skip.
        @test_throws ArgumentError vector_path_network_solve(
            prob, targets; order = 24, h = 0.3, fine_grid = ComplexF64[])

        # No fine_grid: the Stage-2 fields stay empty, the V7 scatter
        # behaviour is byte-unaffected.
        sol = vector_path_network_solve(prob, targets; order = 24, h = 0.3)
        @test sol.grid_z == ComplexF64[]
        @test sol.grid_y == Vector{ComplexF64}[]

        # The 6-arg backward-compat constructor also yields empty Stage-2
        # fields (the hand-built-fixture path) AND an empty `visited_jets`
        # (the pre-B1 jet cache — ADR-0025 Amendment 1).
        sol6 = VectorPathNetworkSolution{Float64}(
            ComplexF64[1.0], [ComplexF64[1.0]], ComplexF64[0.5],
            [[ComplexF64[1.0]]], [ComplexF64[1.0]], [0])
        @test isempty(sol6.grid_z)
        @test isempty(sol6.grid_y)
        @test isempty(sol6.visited_jets)
    end

    # -------------------------------------------------------------------------
    # VPN.2.6 — the B1 true-radius validity gate (ADR-0025 Amendment 1).
    #
    # The gate `R_gate = min(s(tol)·h_JZ·h_v, 0.5·h_v·min|t*|)` replaces
    # the old `h_v` step gate.  These assertions test the gate against
    # computed reference radii: an in-disc vs out-of-disc query (finite
    # vs NaN), the pole-adjacency clamp, the jet-cache / fallback match,
    # and the Rule-1 degenerate-input throws.
    # -------------------------------------------------------------------------
    @testset "VPN.2.6 B1 true-radius gate" begin

        # --- the safety factor s(tol) — calibrated, ADR-0025 Amendment 1 ---
        @test _safety_factor(1.0e-6)  == 0.34
        @test _safety_factor(1.0e-8)  == 0.36          # production default
        @test _safety_factor(1.0e-10) == 0.30

        # --- the jet cache is populated and matches the no-cache fallback ---
        # The Stage-1 walk caches each node's order-`order` Taylor jet in
        # `visited_jets`; recomputing it from the node `(z,y)` state via
        # `vector_taylor_coefficients` must give the bit-identical jet —
        # the cache is a pure speed-up, not a different computation.
        prob_h  = harmonic_problem(0.0 + 0.0im, 2.0 + 0.0im; order = 24)
        tgts_h  = ComplexF64[x + 0.0im for x in 0.3:0.3:1.8]
        sol_h   = vector_path_network_solve(prob_h, tgts_h; order = 24, h = 0.3)
        @test !isempty(sol_h.visited_jets)
        @test length(sol_h.visited_jets) == length(sol_h.visited_z)
        for k in eachindex(sol_h.visited_z)
            recomputed = vector_taylor_coefficients(
                prob_h.f, sol_h.visited_z[k], sol_h.visited_y[k], 24)
            @test sol_h.visited_jets[k] == recomputed
        end

        # --- in-disc vs out-of-disc: a query just inside R_gate is finite,
        #     just outside is NaN.  R_gate is computed directly from the
        #     gate helper, so the boundary is a known reference value.
        #     The probe points are planted along the IMAGINARY axis: the
        #     harmonic walk's nodes lie on the real axis, so a point
        #     i·r is closest to node 1 (z = 0) for r up to ≈ R_g1 — node
        #     1 is unambiguously the gating node for both probes. ---
        h_v1   = real(sol_h.visited_h[1])
        R_g1   = _validity_radius(sol_h.visited_jets[1],
                                  sol_h.visited_denominator[1], h_v1, 1.0e-8)
        @test R_g1 > 0
        z_in   = (0.97 * R_g1) * im                 # just inside the disc
        z_out  = (1.03 * R_g1) * im                 # just outside the disc
        sol_io = vector_path_network_solve(prob_h, tgts_h; order = 24,
                                           h = 0.3, tol = 1.0e-8,
                                           fine_grid = ComplexF64[z_in, z_out])
        @test all(isfinite, sol_io.grid_y[1])        # inside  ⇒ finite
        @test all(z -> isnan(real(z)) && isnan(imag(z)), sol_io.grid_y[2])
        # extrapolate=true bypasses the gate even for the out-of-disc cell.
        sol_iox = vector_path_network_solve(prob_h, tgts_h; order = 24,
                                            h = 0.3, extrapolate = true,
                                            fine_grid = ComplexF64[z_out])
        @test all(isfinite, sol_iox.grid_y[1])

        # --- the gate is tol-dependent: a stricter tol shrinks R_gate ---
        # h_JZ = (ε/‖c_k‖)^(1/k) is increasing in ε, and s(1e-10) < s(1e-8),
        # so the 1e-10 truncation radius is strictly the smaller of the two.
        R_strict = _validity_radius(sol_h.visited_jets[1],
                                    sol_h.visited_denominator[1], h_v1, 1.0e-10)
        @test R_strict < R_g1

        # --- the pole-adjacency clamp.  Build a synthetic node by hand: a
        #     fast-decaying jet (the exp jet — huge truncation radius) over
        #     a shared Q with a root planted very close in the t-variable.
        #     Q = [1, -1/t0] has its single root at t = t0; with t0 small
        #     the clamp 0.5·h_v·t0 must bind below the (huge) truncation
        #     radius.  This is the deep-wedge pole-adjacent regime. ---
        exp_jet = [ComplexF64[1 / factorial(big(k)) for k in 0:24]]
        h_v     = 0.3
        t0      = 0.05                              # Q-root very close
        Q_near  = ComplexF64[1.0, -1.0 / t0]        # root at t = t0
        R_clamped = _validity_radius(exp_jet, Q_near, h_v, 1.0e-8)
        @test R_clamped ≈ 0.5 * h_v * t0  rtol=1.0e-10   # clamp binds exactly
        # Same jet, a Q with its root far away ⇒ the clamp is inert and the
        # truncation term (much larger than the clamp above) governs.
        Q_far   = ComplexF64[1.0, -1.0 / 50.0]      # root at t = 50
        R_trunc = _validity_radius(exp_jet, Q_far, h_v, 1.0e-8)
        @test R_trunc > R_clamped                   # truncation term binds
        @test R_trunc < 0.5 * h_v * 50.0            # ... below the far clamp

        # --- a constant Q (no roots) ⇒ the clamp is simply absent; the
        #     radius is the pure truncation term, not an error. ---
        R_noclamp = _validity_radius(exp_jet, ComplexF64[1.0], h_v, 1.0e-8)
        @test R_noclamp == R_trunc                  # Q_far root was inert too
        @test isfinite(R_noclamp) && R_noclamp > 0

        # --- Rule 1: degenerate inputs throw with a suggestion/detail ---
        # Empty shared-Q denominator — a malformed node.
        @test_throws ArgumentError _validity_radius(
            exp_jet, ComplexF64[], h_v, 1.0e-8)
        # Empty jet — no component series to estimate a radius from.
        @test_throws ArgumentError _validity_radius(
            Vector{ComplexF64}[], ComplexF64[1.0], h_v, 1.0e-8)
        # A component jet too short for the two-coefficient JZ formula.
        @test_throws ArgumentError _validity_radius(
            [ComplexF64[1.0]], ComplexF64[1.0], h_v, 1.0e-8)
    end

    # -------------------------------------------------------------------------
    # VPN.2.7 — B4 Voronoi honesty (ADR-0025 Amendment 2, bead 0ln.37.8).
    #
    # The B4 audition retained the nearest-node Voronoi fill.  Its
    # honesty contract — the load-bearing property the audition relied
    # on to conclude "a blend cannot be more honest" — is: every covered
    # (finite) Stage-2 cell is evaluated by SOME visited node whose B1
    # true-radius disc `R_gate` strictly contains that cell; the Padé is
    # never evaluated outside its verified disc.  We assert this
    # directly: for each finite fill cell, at least one node satisfies
    # `abs(z_f − z_v) ≤ R_gate(node)`.  A Riccati pole problem is used
    # so the grid genuinely straddles covered and NaN regions.
    # -------------------------------------------------------------------------
    @testset "VPN.2.7 B4 Voronoi honesty" begin
        p    = 1.0 + 0.6im
        z0   = -0.4 + 0.0im
        cs   = ComplexF64[1.0, 0.7, -1.3]
        prob = riccati_pole_problem_s2(p, z0, cs; order = 24)
        targets = ComplexF64[x + y * im
                             for x in 0.0:0.4:1.6 for y in -0.4:0.4:1.2]
        fine = ComplexF64[x + y * im
                          for x in 0.0:0.2:1.6 for y in -0.2:0.2:1.2]
        tol  = 1.0e-8
        sol  = vector_path_network_solve(prob, targets; order = 24, h = 0.25,
                                         tol = tol, fine_grid = fine)

        # Recompute every node's B1 R_gate — the exact gate `_stage2_fill`
        # applies (same `_validity_radius`, same cached `visited_jets`).
        n_nodes = length(sol.visited_z)
        R_gate  = [_validity_radius(sol.visited_jets[k],
                                    sol.visited_denominator[k],
                                    real(sol.visited_h[k]), tol)
                   for k in 1:n_nodes]

        # Every finite fill cell must lie strictly inside SOME node's
        # verified disc — the honesty contract of the retained Voronoi
        # fill (no Padé evaluated outside its disc).
        n_covered = 0
        for (zf, yf) in zip(sol.grid_z, sol.grid_y)
            all(isfinite, yf) || continue
            n_covered += 1
            in_some_disc = any(abs(zf - sol.visited_z[k]) <= R_gate[k]
                               for k in 1:n_nodes)
            @test in_some_disc
        end
        # The fill must actually cover cells — otherwise the test is vacuous.
        @test n_covered ≥ 20

        # The complementary half: a cell OUTSIDE every node's disc must be
        # NaN (fail-soft), never silently evaluated.  Plant one far away.
        far  = 40.0 + 40.0im
        sol_f = vector_path_network_solve(prob, targets; order = 24, h = 0.25,
                                          tol = tol,
                                          fine_grid = ComplexF64[far])
        @test all(z -> isnan(real(z)) && isnan(imag(z)), sol_f.grid_y[1])
        @test all(abs(far - sol_f.visited_z[k]) > R_gate[k]
                  for k in 1:length(sol_f.visited_z))
    end

end

# =============================================================================
# Mutation-proof record (CLAUDE.md Rule 4).
#
# Procedure: perturb `src/VectorPathNetworkStage2.jl`, rerun this file,
# confirm RED, restore.  Four mutations bit; 1 false-trail mutation is
# recorded because the "why it doesn't bite" is itself instructive.
#
#   M1 — drop the 1/h_v rescale.  In `_stage2_fill`, change
#          t = (z_f - z_v) / h_v
#        to
#          t = (z_f - z_v)
#        (evaluate the canonical Padé at the raw z-displacement instead
#        of the rescaled variable t = (z - z_node)/h the approximant
#        actually lives in).
#        Expected: every covered off-node Stage-2 cell evaluates the
#        Padé at the wrong argument, so the closed-form-oracle
#        assertions fail.
#        Result: RED — VPN.2.1 (Riccati err < 1e-6) and VPN.2.2
#        (harmonic sin/cos err < 1e-8) bit on essentially every covered
#        cell.  VPN.2.4 node-consistency correctly STILL passed: at a
#        node z_f = z_v the displacement is 0 either way (t = 0
#        unchanged), so M1 specifically corrupts only off-node cells —
#        exactly the load-bearing 1/h_v factor's job.  Restored to GREEN.
#
#   M2-false — double every Horner coefficient.  In `_eval_poly`, change
#          s = s * t + c[k]
#        to
#          s = s * t + 2 * c[k]
#        Result: GREEN — NO test bit.  Doubling every coefficient
#        computes 2·Pᵢ(t) for the numerator AND 2·Q(t) for the
#        denominator; the factor of 2 cancels exactly in the ratio
#        2·Pᵢ(t)/(2·Q(t)) = Pᵢ(t)/Q(t).  A scalar multiple of `_eval_poly`
#        is invisible to a shared-Q quotient.  Recorded as a false trail.
#
#   M2 — reverse the Horner loop direction.  In `_eval_poly`, change
#          for k in length(c):-1:1
#        to
#          for k in 1:length(c)
#        (Horner requires high-to-low coefficient order; iterating
#        low-to-high evaluates the *reversed* polynomial tⁿ·P(1/t)).
#        Result: RED — VPN.2.1, VPN.2.2 and VPN.2.4 all failed (the
#        reversed evaluation breaks the oracle match AND Pᵢ(0)/Q(0)=yᵢ
#        at the node).  Restored to GREEN.
#
# -- B1 true-radius gate mutations (bead padetaylor-0ln.37.5) --
#
#   M3 — neuter the gate to the old `h_v` step gate.  In
#        `_validity_radius`, replace the whole body's return with
#          return float(real(C))(h_v)
#        (the pre-B1 behaviour: gate at the bare step size).
#        Expected: the in-disc/out-of-disc boundary moves from the
#        computed true radius `R_gate` to `h_v`.  For the harmonic walk
#        `R_g1 ≈ 1.52 ≫ h_v = 0.3`, so the VPN.2.6 `z_in`/`z_out` cells
#        (both at ~`R_g1`, hence both > h_v) flip: `z_in` becomes NaN.
#        Result: RED — VPN.2.6 `@test all(isfinite, sol_io.grid_y[1])`
#        bit (z_in, an honest in-disc point, wrongly NaN-masked).  The
#        tol-dependence `R_strict < R_g1` also bit (the mutant returns a
#        tol-independent `h_v`).  Restored to GREEN.
#
#   M4 — drop the pole-adjacency clamp.  In `_validity_radius`, delete
#        the `if length(denominator) ≥ 2 … end` clamp block, returning
#        `R_trunc` unconditionally.
#        Expected: the clamp no longer binds for a near Q-root, so the
#        synthetic clamp node's radius jumps from `0.5·h_v·t0 = 0.0075`
#        to the huge exp-jet truncation radius.
#        Result: RED — VPN.2.6 `R_clamped ≈ 0.5·h_v·t0` bit (the mutant
#        returns the truncation radius instead of the clamp).  Restored
#        to GREEN.
#
#   M5 — corrupt the safety factor.  In `_safety_factor`, change the
#        `tol = 1e-8` branch `return 0.36` to `return 1.20` (s > 1, a
#        dishonest over-cover that evaluates Padé well outside its
#        verified disc).
#        Result: RED — VPN.2.6 `_safety_factor(1e-8) == 0.36` bit
#        directly (1 assertion).  Note the in/out-of-disc cells did NOT
#        flip: VPN.2.6 computes its `R_g1` reference through the *same*
#        mutated `_validity_radius`, so `z_in`/`z_out` track the
#        inflated radius and the in/out test stays self-consistent.
#        The direct `_safety_factor` assertion is the load-bearing
#        catch — it pins the calibrated constant against drift.
#        Restored to GREEN.
#
#   M6 — drop the h_v rescale of the jet inside the gate.  In
#        `_validity_radius`, change
#          rescaled = [_rescale_by_powers(j, C(h_v)) for j in jet]
#        to
#          rescaled = jet
#        (feed `vector_step_jorba_zou` the raw, un-rescaled jet — h_JZ
#        is then a step in the wrong variable, not the t-variable the
#        canonical Padé lives in).
#        Result: RED — VPN.2.1 closed-form Riccati oracle bit (the
#        un-rescaled jet of the slowly-decaying 1/(z−p) tail yields a
#        too-large h_JZ, so the gate OVER-covers; a near-pole cell that
#        should be NaN-masked is instead evaluated and its `err < 1e-6`
#        assertion fails).  This is the honesty failure the rescale
#        prevents — caught by the oracle, not just a gate unit test.
#        Restored to GREEN.
#
# -- B4 Voronoi-honesty mutation (bead padetaylor-0ln.37.8) --
#
#   M7 — relax the Stage-2 gate comparison to over-cover.  In
#        `_stage2_fill`, change the fail-soft gate test
#          if abs(z_f - z_v) > R_gate
#        to
#          if abs(z_f - z_v) > 100 * R_gate
#        (the fill now evaluates the canonical Padé at grid cells far
#        OUTSIDE the node's verified disc — the dishonest extrapolation
#        the Voronoi honesty contract forbids).  VPN.2.7 recomputes its
#        OWN `R_gate` reference through the *untouched* `_validity_radius`,
#        so its in-some-disc check is the honest one and correctly fails
#        for the now-extrapolated cells.
#        Expected: VPN.2.7's per-cell `@test in_some_disc` bites on every
#        cell the relaxed gate wrongly evaluates between `R_gate` and
#        `100·R_gate`, and the `far = 40+40im` fail-soft assertion bites
#        (the far cell is now finite, not NaN).
#        Result: RED — VPN.2.7 went 41 failures (the in-some-disc honesty
#        assertions for the extrapolated cells + the fail-soft `far`
#        cell).  VPN.2.3 and VPN.2.6 also bit, as expected — they share
#        the same `_stage2_fill` gate comparison.  VPN.2.7 is the
#        load-bearing B4 catch: it pins the honesty contract the
#        audition relied on (no Padé evaluated outside its disc).
#        Restored to GREEN.
#
# Certified bites: M1 → VPN.2.1+2.2 RED; M2 → VPN.2.1+2.2+2.4 RED;
# M3 → VPN.2.6 (3 assertions: tol-dependence + clamp regime) RED;
# M4 → VPN.2.6 (2 assertions: clamp + no-clamp consistency) RED;
# M5 → VPN.2.6 s-factor RED; M6 → VPN.2.1 oracle RED;
# M7 → VPN.2.7 Voronoi-honesty RED (41 assertions).  M2-false is a
# recorded non-bite.  Restored to GREEN after each mutation.
# =============================================================================
