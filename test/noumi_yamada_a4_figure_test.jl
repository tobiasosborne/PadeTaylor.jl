# test/noumi_yamada_a4_figure_test.jl
#
# Quantitative acceptance test for `figures/noumi_yamada_a4_pole_
# field.jl` — the A_4^(1) Noumi–Yamada pole-field figure (bead
# `padetaylor-0ln.17`, v0.2 plan row V8a; the PRD `n ≥ 4` acceptance
# item).
#
# ## Why this test exists / what it pins
#
# There is no published A_4^(1) pole-field reference, so this test
# does NOT pin eyeballed pole coordinates (CLAUDE.md Rule 5: a test
# asserts an invariant, never "the number I happened to see").  It
# reproduces the figure's *exact* computation — the
# `vector_path_network_solve` + `extract_poles_shared_q` pipeline —
# by `include`-ing the same helper kernel the figure script uses
# (`figures/_noumi_yamada_a4_helpers.jl`), then asserts genuine
# structural invariants:
#
#   NYF.1.1 — IC / α validity: Σα = 1, Σf0 = t0 exactly, every f0
#             component non-zero (the V5b all-nonzero requirement).
#   NYF.1.2 — the flow invariant `Σf_j(t) = t` holds at every visited
#             node of the walk (pillar B §1.2).
#   NYF.1.3 — the extracted pole field is finite and non-empty over
#             the chosen region.
#   NYF.1.4 — pole genuineness: at each extracted pole the shared-Q
#             denominator of the supporting node genuinely
#             (near-)vanishes there.
#   NYF.1.5 — equation-level conjugate symmetry: the A_4^(1) RHS is
#             *exactly* conjugate-symmetric (real α), and the
#             canonical shared-Q built at a conjugate node `(t̄, f̄)`
#             is the *exact* conjugate of the one built at `(t, f)` —
#             so the *system's* pole field is conjugate-symmetric.
#             (The wedge *walk* tree is not conjugate-symmetric — it
#             is a deterministic but asymmetric traversal — so we pin
#             the symmetry at the equation level where it is exact,
#             not on the extracted set; see the note at NYF.1.5.)
#   NYF.1.6 — reproducibility: a fixed walk produces a deterministic
#             visited tree and a deterministic pole set (the PRD's
#             "reproducible from figures/").
#
# Self-contained: `using Test, PadeTaylor` only; runnable standalone
# (`julia --project=. test/noumi_yamada_a4_figure_test.jl`) and under
# `runtests.jl`.
#
# ## Helper-include path note
#
# The figure helper lives under `figures/` and `using`s `PadeTaylor`
# (no CairoMakie — that is kept in the figure script).  We resolve
# the path relative to this test file's `@__DIR__` so the include
# works from either the repo root or the `test/` directory.

using Test
using PadeTaylor
using PadeTaylor.VectorStepper: VectorPadeStepperState,
                                vector_pade_step_with_pade!
using PadeTaylor.NoumiYamada:   noumi_yamada_rhs
using Polynomials: Polynomial, roots

# The shared figure kernel — the EXACT computation the figure renders.
include(joinpath(@__DIR__, "..", "figures", "_noumi_yamada_a4_helpers.jl"))

@testset "A_4^(1) Noumi–Yamada pole-field figure acceptance (NYF.*)" begin

    # ---- NYF.1.1 — IC / α validity ----------------------------------
    @testset "NYF.1.1: generic α + all-nonzero IC" begin
        @test length(NY_ALPHA) == NY_D
        @test length(NY_F0)    == NY_D
        # Σα = 1, the k = 1 normalisation.
        @test sum(NY_ALPHA) ≈ 1.0 + 0.0im atol = 1.0e-14
        # Σf0 = t0 exactly — the constraint Σf_j = t at the start point.
        @test sum(NY_F0) ≈ NY_T0 atol = 1.0e-14
        # Every IC component non-zero — the V5b shared-Q degeneracy
        # is triggered by an identically-zero component.
        @test all(c -> abs(c) > 1.0e-3, NY_F0)
        # The α are genuinely *generic*: five distinct values (not the
        # Type-C 1/5,…,1/5 rational-solution parameters, which would
        # give a pole-free figure).
        @test length(unique(NY_ALPHA)) == NY_D
        # Building the problem must not throw — the constructor
        # re-validates Σα and the constraint.
        @test a4_problem() isa NoumiYamadaProblem
    end

    # ---- Run the figure's walk ONCE; reuse below --------------------
    sol = a4_walk()
    poles = a4_pole_field(sol)

    # ---- NYF.1.2 — the Σf_j = t constraint along the walk -----------
    # The constraint is an invariant of the A_4^(1) flow (pillar B
    # §1.2): summing the system, every f_a f_b product cancels, so
    # d/dt(Σf − t) = 0.  A faithful walk keeps |Σf − t| near zero at
    # every visited node.  Tolerance: the walk takes order-24 Padé
    # steps of magnitude h = 0.3 over a region of radius ~3.5, so a
    # few e-9 of accumulated truncation drift is expected; 1e-6 is a
    # generous-but-genuine ceiling that the actual ~2e-9 clears by
    # three orders.
    @testset "NYF.1.2: Σf_j(t) = t holds along the walk" begin
        resid = a4_constraint_residual(sol)
        @info "NYF.1.2 max |Σf_j − t| over $(length(sol.visited_z)) nodes = $resid"
        @test resid ≤ 1.0e-6
        # And it really is *small* — pin the order of magnitude so a
        # regression that loosens it to e.g. 1e-3 is caught.
        @test resid ≤ 1.0e-7
    end

    # ---- NYF.1.3 — finite, non-empty pole field ---------------------
    @testset "NYF.1.3: pole field finite + non-empty" begin
        @test !isempty(poles)
        @test all(isfinite, poles)
        # The generic problem has a genuinely non-trivial pole cloud
        # — a rational (Type A/B/C) solution would give zero poles.
        # ≥ 10 is a structural floor well below the observed ~26.
        @info "NYF.1.3 extracted $(length(poles)) poles"
        @test length(poles) ≥ 10
        # Poles lie in a sensible neighbourhood of the sampled region
        # (radius_t = 5 · h = 0.3 caps how far a local Padé reaches;
        # combined with the [-2.5,2.5]² region, all poles fall well
        # inside |t| ≤ 8).
        @test all(p -> abs(p) ≤ 8.0, poles)
    end

    # ---- NYF.1.4 — pole genuineness ---------------------------------
    # Each extracted pole is the cluster representative of a shared-Q
    # root.  Pole genuineness: the *supporting* node's shared
    # denominator Q must genuinely (near-)vanish at the rescaled root
    # t* = (p − z_node)/h_node.  We scan every node, root its Q, map
    # roots to the t-plane, and confirm that for each extracted pole
    # there is a node whose Q has a root mapping to it AND whose Q
    # genuinely vanishes there (|Q(t*)| tiny).
    @testset "NYF.1.4: extracted poles are genuine Q-roots" begin
        # Gather (z-plane pole, |Q(root)|) pairs from every node.
        genuine = Tuple{ComplexF64, Float64}[]
        for k in eachindex(sol.visited_denominator)
            Q = sol.visited_denominator[k]
            length(Q) ≥ 2 || continue
            qp = Polynomial(Q)
            z_node = sol.visited_z[k]
            h_node = sol.visited_h[k]
            for t in roots(qp)
                abs(t) ≤ NY_RADIUS_T || continue
                push!(genuine,
                      (ComplexF64(z_node + h_node * t), Float64(abs(qp(t)))))
            end
        end
        @test !isempty(genuine)
        # Every Q-root is a genuine vanishing point of its Q
        # (root-finder residual): |Q(t*)| ≈ 0.
        worst_resid = maximum(g[2] for g in genuine)
        @info "NYF.1.4 worst |Q(root)| over all nodes = $worst_resid"
        @test worst_resid ≤ 1.0e-6
        # Every *extracted* pole coincides with at least one node's
        # Q-root (it is the cluster representative of such roots).
        for p in poles
            d = minimum(abs(p - g[1]) for g in genuine)
            @test d ≤ NY_CLUSTER_ATOL
        end
    end

    # ---- NYF.1.5 — equation-level conjugate symmetry ----------------
    # The chosen α are real, so the A_4^(1) RHS is conjugate-symmetric:
    # `rhs(t̄, f̄) = conj(rhs(t, f))`.  Consequently, if `f(t)` solves
    # the system so does `conj(f(conj t))`, and the *system's* pole
    # field is conjugate-symmetric.  We pin this where it is *exact*:
    #
    #   (a) the RHS satisfies rhs(t̄, f̄) = conj(rhs(t,f)) to 0;
    #   (b) the canonical shared-Q built at a node `(t̄, f̄)` is the
    #       *exact* conjugate of the canonical-Q built at `(t, f)` —
    #       hence its roots are the conjugate roots.
    #
    # NOTE — why we do NOT assert conjugate symmetry of the *extracted
    # pole set*: the wedge walk is a deterministic but conjugate-
    # *asymmetric* tree traversal (it grows greedily from t0 toward
    # each target, so the visited-node set is not closed under
    # conjugation).  The extracted set therefore samples the
    # conjugate-symmetric pole field asymmetrically.  Pinning
    # set-level symmetry would be pinning a walk artefact, not a real
    # invariant; the equation-level invariant (b) is the genuine,
    # exact structural property — and it is exactly what guarantees
    # the underlying pole field is symmetric.
    @testset "NYF.1.5: A_4^(1) RHS + canonical-Q conjugate symmetry" begin
        rhs = noumi_yamada_rhs(NY_N; α = NY_ALPHA)
        # (a) RHS conjugate symmetry at a generic complex sample.
        t_s = 0.4 + 0.3im
        y_s = ComplexF64[0.7 + 0.1im, -0.3 - 0.2im, 0.5 + 0.05im,
                         -0.55 + 0.1im, -0.35 - 0.05im]
        a = rhs(t_s, y_s)
        b = conj.(rhs(conj(t_s), conj.(y_s)))
        @test maximum(abs.(a .- b)) ≤ 1.0e-14

        # (b) canonical shared-Q conjugate symmetry.  Build the Padé
        # step's shared denominator at (z, y) and at (z̄, ȳ); the
        # latter must be the exact conjugate of the former.
        z = 0.4 + 0.5im
        y = ComplexF64[0.7 + 0.1im, -0.3 - 0.15im, 0.5 + 0.2im,
                       -0.55 - 0.05im, -0.35 - 0.1im]
        st1 = VectorPadeStepperState{ComplexF64}(z, y)
        _, _, Q1 = vector_pade_step_with_pade!(st1, rhs, NY_ORDER,
                                               ComplexF64(NY_H))
        st2 = VectorPadeStepperState{ComplexF64}(conj(z), conj.(y))
        _, _, Q2 = vector_pade_step_with_pade!(st2, rhs, NY_ORDER,
                                               ComplexF64(NY_H))
        @test length(Q1) == length(Q2)
        @test maximum(abs.(Q1 .- conj.(Q2))) ≤ 1.0e-12
        # Hence the root sets are exact conjugates: the system's pole
        # field is conjugate-symmetric.
        r1 = sort(roots(Polynomial(Q1)); by = x -> (real(x), imag(x)))
        r2 = sort(conj.(roots(Polynomial(Q2)));
                  by = x -> (real(x), imag(x)))
        @test maximum(abs.(r1 .- r2)) ≤ 1.0e-10
    end

    # ---- NYF.1.6 — reproducibility ----------------------------------
    # The PRD acceptance: "reproducible from figures/".  A fixed walk
    # must produce the identical visited tree and the identical pole
    # set on every call — the figure is a function of the (fixed)
    # helper constants, nothing else.
    @testset "NYF.1.6: deterministic walk + pole set" begin
        sol2 = a4_walk()
        @test length(sol2.visited_z) == length(sol.visited_z)
        @test sol2.visited_z      == sol.visited_z
        @test sol2.visited_parent == sol.visited_parent
        poles2 = a4_pole_field(sol2)
        @test length(poles2) == length(poles)
        @test poles2 == poles
    end

end

# =============================================================================
# Mutation-proof record (CLAUDE.md Rule 4 — "mutation-proving replaces
# the literal RED-first step").
#
# Procedure: perturb a load-bearing input/source, rerun this file,
# confirm RED, restore.  2 meaningful mutations applied + verified
# 2026-05-20.
#
#   M1 — corrupt α so Σα ≠ 1.  In `figures/_noumi_yamada_a4_helpers.jl`
#        change `NY_ALPHA`'s last entry `0.20` → `0.25` (now Σα = 1.05).
#        Expected: `NoumiYamadaProblem` / `noumi_yamada_rhs` throw the
#        `_check_alpha_sum` ArgumentError; the figure's `a4_problem()`
#        / `a4_walk()` calls fail before any pole is extracted.
#        Result: RED — `a4_problem()` in NYF.1.1 throws ArgumentError
#        ("the parameters must satisfy Σα_j = 1"); the shared `sol`
#        construction at the top of the testset also throws, so every
#        subsequent NYF.1.* testset errors.  Strongest possible bite
#        (the figure cannot even be built).  Restored to GREEN.
#
#   M2 — flip a sign in the A_4^(1) RHS.  In `src/NoumiYamada.jl`
#        change the bracket term
#          `fvec[slot(j + 2r - 1)] - fvec[slot(j + 2r)]`
#        →  `fvec[slot(j + 2r - 1)] + fvec[slot(j + 2r)]`
#        (the even-offset partial sum sign).  This is no longer the
#        Noumi–Yamada system: the `Σf_j = t` constraint is no longer
#        an invariant of the flow (the f_a f_b products no longer
#        cancel in the total sum).
#        Result: RED — NYF.1.2 bites hard: `max |Σf_j − t|` along the
#        walk jumps from ~2e-9 to O(1) (the perturbed flow does not
#        conserve Σf − t), so both the 1e-6 and the 1e-7 ceilings
#        fail.  NYF.1.5(a) RHS-conjugate-symmetry still holds (a sign
#        flip preserves realness) but NYF.1.5(b)'s exact-conjugate
#        canonical-Q comparison still passes too — the genuine
#        load-bearing detector for M2 is the constraint test NYF.1.2,
#        which is precisely the invariant the Noumi–Yamada cancellation
#        structure produces.  Restored to GREEN.
#
# Certified bites: M1 → all NYF.1.* error (problem cannot be built),
# M2 → NYF.1.2 (2 assertions RED).  Restored to GREEN after each.
# =============================================================================
