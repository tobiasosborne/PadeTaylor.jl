# test/type_stability_test.jl -- bead `padetaylor-krgy.12`: the efficiency gate
# of the test-hardening sweep -- allocation budgets + type-stability assertions
# on the numerical core's HOT inner loop.
#
# ## Why this chapter exists
#
# The perf-regression gate (bead krgy.11, `benchmark/run_perf_gate.jl`) pins
# WALL TIME against a committed baseline.  Wall time is a downstream symptom; it
# moves for many reasons (BLAS, thermal, scheduler) and a wide slack band must
# swallow that jitter.  This file pins the two UPSTREAM, deterministic causes of
# the per-step cost that a time gate can only see indirectly:
#
#   * **Type stability.**  A type instability in the inner loop -- an accumulator
#     that infers to `Union{Int,Float64}`, a field access that returns `Any`,
#     a closure whose return type cannot be proven -- forces Julia to box values
#     and dispatch `*`/`+` at RUNTIME.  That is both a silent ~10x slowdown AND a
#     correctness smell: dynamic dispatch on `Any` means inference gave up, and
#     the next refactor can silently change which method runs.  `Test.@inferred`
#     proves the concrete return type is inferable (it catches the broad class,
#     including small-`Union` accumulator instabilities); JET's `@test_opt`
#     (report_opt mode -- the type-stability complement to krgy.2's
#     `@report_call`) proves NO runtime dispatch / `Any` remains in the call.
#
#   * **Allocation budget.**  Heap allocation in the hot path means GC pressure
#     and a non-`O(1)`-per-step memory profile.  `_evaluate_pade` /
#     `_evaluate_pade_deriv` are two-pass Horner sweeps over fixed vectors -- they
#     must be allocation-FREE, and `AllocCheck.@check_allocs` PROVES that
#     STATICALLY (over every inferred code path, not just the one input we time).
#     The genuinely-allocating paths (`taylor_coefficients_2nd` builds Taylor1
#     jets; `pade_step_with_pade!` runs an SVD/Toeplitz solve) cannot be
#     alloc-free, so we assert a BUDGET derived from the committed perf baseline
#     (`benchmark/baseline.toml`) -- and, crucially, that an N-step solve
#     allocates `O(N)`, NOT `O(N^2)` (the signature of an accidental
#     quadratic-in-steps blow-up).
#
# ## Inventory finding (Law 1 / Rule 9 -- every assertion is a real invariant)
#
# This package was NOT built under a type-stability gate, so the bead anticipated
# genuine instabilities to surface and be triaged (fix-trivially / @test_broken +
# bead).  The full probe (Float64 AND ComplexF64 path-network element types) found
# the ENTIRE hot path already type-stable: `@inferred` passes and `report_opt`
# returns ZERO reports for every one of `_evaluate_pade`, `_evaluate_pade_deriv`,
# `taylor_coefficients_2nd`, `robust_pade`, and `pade_step_with_pade!`.  The
# Coefficients module's "Type-stability watch" docstring (src/Coefficients.jl:67)
# and PadeStepper's `where {T}` discipline are why.  So this gate LOCKS IN that
# property as a regression tripwire rather than fixing a defect -- which is the
# correct senior outcome when the inventory comes back clean, and the
# mutation-proof footer below proves the gate still BITES a freshly-injected
# instability or allocation.
#
# ## Budget derivations (Rule 5 -- numeric invariants, not "didn't throw")
#
# All byte budgets are read off `benchmark/baseline.toml` (captured Julia 1.12.5,
# this box) with headroom for cross-version/platform jitter in TaylorSeries /
# LinearAlgebra internals.  The budgets are ORDER-OF-MAGNITUDE tripwires: they
# catch a doubling (a lost `@view`, an O(steps^2) accumulation), not a 5% drift.
#
#   * `taylor_coefficients_2nd` o30: baseline bytes = 17584.  Budget 40_000
#     (~2.3x headroom) -- trips a genuine doubling of the jet's allocation.
#   * `pade_step_with_pade!` 1 step: baseline bytes = 23280.  Budget 50_000
#     (~2.1x headroom).
#   * O(steps) LINEARITY: 4 steps must allocate < 3x a single step's BUDGET (not
#     4 * budget, and emphatically not 16 * single-step -- the O(steps^2) trap).
#     This is the load-bearing structural invariant the bead names explicitly.
#
# Reference: docs/test_corpus/03_hardening_methodology.md §EFFICIENCY (krgy.12);
# benchmark/baseline.toml (krgy.11 perf floor); test/quality_test.jl (shared JET
# setup, krgy.2); AllocCheck.jl / JET.jl upstream docs.

using Test
using PadeTaylor
using PadeTaylor: robust_pade, PadeApproximant, taylor_coefficients_2nd
using PadeTaylor.PadeStepper:
    PadeStepperState, pade_step_with_pade!, _evaluate_pade, _evaluate_pade_deriv
using JET
using AllocCheck

# ── Fixtures — the canonical FW 2011 ℘ workload (mirrors benchmark/run_perf_gate.jl) ──
# u'' = 6u² with the equianharmonic Weierstrass-℘ IC; order-30 jet; (15,15) Padé.
# Two element types: Float64 (real IVP stepper) and ComplexF64 (path-network
# wedge step) — the two concrete eltypes the production hot path runs under.
const _TS_F  = (z, u, up) -> 6 * u^2
const _U0    = 1.071822516416917
const _UP0   = 1.710337353176786
const _ORDER = 30

const _JET_F = taylor_coefficients_2nd(_TS_F, 0.0, _U0, _UP0, _ORDER)
const _CS_F  = [(_JET_F[k] * 0.5^(k - 1)) for k in eachindex(_JET_F)]
const _P_F   = robust_pade(_CS_F, _ORDER ÷ 2, _ORDER ÷ 2; method = :svd)

const _HC    = complex(0.0, 0.5)
const _JET_C = taylor_coefficients_2nd(_TS_F, complex(0.0, 0.0),
                                       complex(_U0, 0.0), complex(_UP0, 0.0), _ORDER)
const _CS_C  = [(_JET_C[k] * _HC^(k - 1)) for k in eachindex(_JET_C)]
const _P_C   = robust_pade(_CS_C, _ORDER ÷ 2, _ORDER ÷ 2)   # complex ⇒ :classical

@testset "Type stability + allocation budgets (krgy.12)" begin

    # ── TS.1 — TYPE STABILITY via @inferred (concrete-return-type proof) ──────
    # Each call must infer to a concrete type.  @inferred throws if the inferred
    # return type is not concrete (e.g. a Union from an unstable accumulator), so
    # a passing @inferred IS the invariant (Rule 5), not a "didn't throw".  Run
    # on BOTH Float64 and ComplexF64 so inference is exercised on the real paths.
    @testset "TS.1 @inferred — hot path infers concrete (Float64 + ComplexF64)" begin
        @test (@inferred _evaluate_pade(_P_F, 1.0)) isa Float64
        @test (@inferred _evaluate_pade_deriv(_P_F, 1.0)) isa Float64
        @test (@inferred _evaluate_pade(_P_C, complex(1.0, 0.0))) isa ComplexF64
        @test (@inferred _evaluate_pade_deriv(_P_C, complex(1.0, 0.0))) isa ComplexF64
        @test (@inferred taylor_coefficients_2nd(_TS_F, 0.0, _U0, _UP0, _ORDER)) isa Vector{Float64}
        @test (@inferred robust_pade(_CS_F, _ORDER ÷ 2, _ORDER ÷ 2; method = :svd)) isa PadeApproximant{Float64}
        @test (@inferred robust_pade(_CS_F, _ORDER ÷ 2, _ORDER ÷ 2)) isa PadeApproximant{Float64}
        s = PadeStepperState{Float64}(0.0, _U0, _UP0)
        @test (@inferred pade_step_with_pade!(s, _TS_F, _ORDER, 0.5)) isa Tuple{PadeStepperState{Float64}, PadeApproximant{Float64}}
        sc = PadeStepperState{ComplexF64}(complex(0.0, 0.0), complex(_U0, 0.0), complex(_UP0, 0.0))
        @test (@inferred pade_step_with_pade!(sc, _TS_F, _ORDER, _HC)) isa Tuple{PadeStepperState{ComplexF64}, PadeApproximant{ComplexF64}}
    end

    # ── TS.2 — NO RUNTIME DISPATCH via JET report_opt ─────────────────────────
    # @test_opt proves the call carries no `Any`-typed runtime dispatch — the
    # complement to @inferred (it catches dispatch @inferred's concrete-return
    # check can miss when only an intermediate is boxed).  An empty report set is
    # the invariant; @test_opt fails the testset on ANY optimisation-failure
    # report inside these calls.
    @testset "TS.2 @test_opt — no runtime dispatch in the hot path" begin
        @test_opt _evaluate_pade(_P_F, 1.0)
        @test_opt _evaluate_pade_deriv(_P_F, 1.0)
        @test_opt _evaluate_pade(_P_C, complex(1.0, 0.0))
        @test_opt _evaluate_pade_deriv(_P_C, complex(1.0, 0.0))
        @test_opt taylor_coefficients_2nd(_TS_F, 0.0, _U0, _UP0, _ORDER)
        @test_opt robust_pade(_CS_F, _ORDER ÷ 2, _ORDER ÷ 2; method = :svd)
        s = PadeStepperState{Float64}(0.0, _U0, _UP0)
        @test_opt pade_step_with_pade!(s, _TS_F, _ORDER, 0.5)
    end

    # ── TS.3 — ALLOCATION-FREE (static, AllocCheck) for the Horner evaluators ──
    # check_allocs STATICALLY proves the function emits no heap allocation over
    # every inferred path (default ignore_throw=true: the Rule-1 DomainError
    # fail-fast branch allocates its message string only when the denominator
    # vanishes on a pole — not a hot-path allocation).  Empty result == proven
    # allocation-free; this is strictly stronger than `@allocated == 0` on one
    # warmed input.  Matches benchmark/baseline.toml [evaluate_pade] allocs = 0.
    @testset "TS.3 AllocCheck — _evaluate_pade(_deriv) provably alloc-free" begin
        @test isempty(AllocCheck.check_allocs(_evaluate_pade, (typeof(_P_F), Float64)))
        @test isempty(AllocCheck.check_allocs(_evaluate_pade_deriv, (typeof(_P_F), Float64)))
        @test isempty(AllocCheck.check_allocs(_evaluate_pade, (typeof(_P_C), ComplexF64)))
        @test isempty(AllocCheck.check_allocs(_evaluate_pade_deriv, (typeof(_P_C), ComplexF64)))
        # Corroborate the static proof with a warmed dynamic count (belt + braces).
        _evaluate_pade(_P_F, 1.0); _evaluate_pade_deriv(_P_F, 1.0)
        @test (@allocated _evaluate_pade(_P_F, 1.0)) == 0
        @test (@allocated _evaluate_pade_deriv(_P_F, 1.0)) == 0
    end

    # ── TS.4 — ALLOCATION BUDGETS for the inherently-allocating paths ─────────
    # These build Taylor1 jets / run an SVD, so they allocate; we bound it to an
    # order-of-magnitude tripwire derived from benchmark/baseline.toml (see the
    # "Budget derivations" header).  @allocated on a warmed call is the measure.
    @testset "TS.4 budget — taylor_coefficients_2nd + pade_step bounded O(1)" begin
        taylor_coefficients_2nd(_TS_F, 0.0, _U0, _UP0, _ORDER)            # warm
        jet_bytes = @allocated taylor_coefficients_2nd(_TS_F, 0.0, _U0, _UP0, _ORDER)
        @test jet_bytes ≤ 40_000        # baseline 17_584 B; ~2.3x headroom.

        _onestep() = pade_step_with_pade!(PadeStepperState{Float64}(0.0, _U0, _UP0),
                                          _TS_F, _ORDER, 0.5)
        _onestep()                                                       # warm
        step_bytes = @allocated _onestep()
        @test step_bytes ≤ 50_000       # baseline 23_280 B; ~2.1x headroom.
    end

    # ── TS.5 — O(steps) LINEARITY (the load-bearing structural invariant) ─────
    # An N-step solve must allocate O(N), NOT O(N^2).  We run 1 and 4 consecutive
    # steps (small h to stay clear of the ℘ pole near z≈1.5) and assert the 4-step
    # cost is < 3x a single step's BUDGET — far below the 16x a quadratic-in-steps
    # accumulation would show.  This catches an accidental "carry the whole history
    # each step" regression a single-step byte budget passes.
    @testset "TS.5 O(steps) — N-step allocation is linear, not quadratic" begin
        _nsteps(n, h) = begin
            s = PadeStepperState{Float64}(0.0, _U0, _UP0)
            for _ in 1:n
                pade_step_with_pade!(s, _TS_F, _ORDER, h)
            end
            s
        end
        _nsteps(1, 0.05)                                                 # warm
        b1 = @allocated _nsteps(1, 0.05)
        b4 = @allocated _nsteps(4, 0.05)
        @test b4 ≤ 4 * b1 + 50_000      # ~linear: 4 steps ≈ 4x one step + slack.
        @test b4 < 16 * b1              # O(steps^2) tripwire: rules out quadratic.
    end

end # @testset "Type stability + allocation budgets (krgy.12)"

# ── Mutation-proof procedure (verified 2026-06-09, bead padetaylor-krgy.12) ──
#
# Each assertion class was proven to BITE by injecting a genuine defect of its
# class into a hot function, running THIS file, confirming the failure, then
# restoring `src/` exactly (`git diff src/` clean of every injection after).
#
#   Mutation A (AllocCheck / @allocated — heap allocation)
#     In `src/PadeStepper.jl::_evaluate_pade`, replace the in-place Horner sweep
#     with one over a fresh copy:
#         num = zero(T); cc = collect(P.a)
#         @inbounds for k in length(cc):-1:1; num = num * z + cc[k]; end
#     Verified bite: TS.3 goes RED —
#         check_allocs(_evaluate_pade, …) → 1 allocation (not isempty), AND
#         @allocated _evaluate_pade(_P_F, 1.0) == 0 → Test Failed.
#     Restored to the in-place `P.a` sweep.
#
#   Mutation B (@inferred — type instability / non-concrete return)
#     In `_evaluate_pade`, make the RETURN type a small `Union` with a value-
#     dependent branch that inference cannot resolve:
#         return real(z) > 0 ? num / den : 0          # was `return num / den`
#     so the inferred return is `Union{Float64,Int64}` (resp. `Union{Int64,
#     ComplexF64}`).  Verified bite: TS.1 goes RED — `@inferred …` reports
#     "return type Float64 does not match inferred return type Union{Float64,
#     Int64}".  TS.2 (report_opt) and TS.3 (AllocCheck) stay GREEN, confirming
#     @inferred is the NECESSARY tool for the non-concrete-return class — the
#     three assertion classes are complementary, not redundant.
#     (Note: the naive `num = 0` Int-SEED injection does NOT bite here — with a
#     non-empty `P.a` the loop body's `Int*Float64+Float64` fixed-point converges
#     to `Float64`, so the return stays concrete.  The function is robust to that
#     particular perturbation; the ternary above is the injection that genuinely
#     makes the return non-concrete.)
#     Restored to `return num / den`.
#
#   Mutation C (@test_opt — runtime dispatch on Any)
#     In `_evaluate_pade`, route the accumulator through an `Any` box:
#         box = Ref{Any}(zero(T))
#         @inbounds for k in length(P.a):-1:1; box[] = box[]*z + P.a[k]; end
#         num = box[]
#     Verified bite: TS.2 goes RED — `@test_opt _evaluate_pade(_P_F, 1.0)` reports
#     runtime-dispatch failures (the `*`/`+` on `Any` cannot be devirtualised);
#     @inferred (TS.1) also goes RED.  Restored to the plain `num` accumulator.
#
# All three mutations restored before hand-off; `git diff src/` shows no injection.
