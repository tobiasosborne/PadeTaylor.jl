"""
V2 spec-from-scratch tests for `PadeTaylor.VectorCoefficients`
(bead `padetaylor-0ln.7`; v0.2 plan row V2).

`vector_taylor_coefficients(f, z0, y0, order)` generates the Taylor jet
of a first-order **vector** ODE `y' = f(z, y)`, `y ∈ ℂ^d`.  These tests
(`VC.*`) are written RED-first per CLAUDE.md Rule 4; each asserts an
invariant against a known-correct value (Rule 5):

    VC.1.1  d=1 reduction         — bit-identical to `taylor_coefficients_1st`
    VC.1.2  decoupled known sys   — y₁'=y₁, y₂'=2y₂ ⇒ c⁽¹⁾_k=1/k!, c⁽²⁾_k=2ᵏ/k!
    VC.1.3  genuinely coupled sys — y₁'=y₂, y₂'=−y₁ ⇒ cos h, −sin h
    VC.1.4  BigFloat type stability — eltype exactly BigFloat, no widening
    VC.1.5  failure modes         — order<1, empty y0 throw informatively
    VC.1.6  mutation-proof        — recorded at end of file

Self-contained: `using Test, PadeTaylor` only — runnable standalone
(`julia --project=. test/vector_coefficients_test.jl`) and under
`runtests.jl`.

Oracle for VC.1.1 is the v0.1 scalar `taylor_coefficients_1st`; oracles
for VC.1.2/1.3/1.4 are closed-form Taylor coefficients of exp / cos / sin.
"""

using Test
using PadeTaylor
using PadeTaylor.VectorCoefficients: vector_taylor_coefficients
using PadeTaylor.Coefficients: taylor_coefficients_1st

@testset "VectorCoefficients (V2)" begin

    @testset "VC.1.1 d=1 reduces to scalar taylor_coefficients_1st" begin
        # A 1-component vector RHS must give a jet bit-identical to the
        # scalar v0.1 routine on the same scalar ODE.  Scalar oracle:
        # dy/dz = z² + y², y(0) = 0, order 14 (the FW 2011 §2.2.1 demo,
        # also covered by scalar test 3.1.2).
        scalar_jet = taylor_coefficients_1st((z, y) -> z^2 + y^2, 0.0, 0.0, 14)
        vec_jet = vector_taylor_coefficients(
            (z, y) -> [z^2 + y[1]^2], 0.0, [0.0], 14)
        @test length(vec_jet) == 1
        @test vec_jet[1] == scalar_jet          # bit-identical, not isapprox
        @test eltype(vec_jet[1]) == Float64

        # A second, unrelated scalar ODE for good measure: exp(z).
        scalar_exp = taylor_coefficients_1st((z, y) -> y, 0.0, 1.0, 10)
        vec_exp = vector_taylor_coefficients((z, y) -> [y[1]], 0.0, [1.0], 10)
        @test vec_exp[1] == scalar_exp
    end

    @testset "VC.1.2 decoupled system y₁'=y₁, y₂'=2y₂" begin
        # Decoupled: y₁(h) = exp(h), y₂(h) = exp(2h).
        # ⇒ c⁽¹⁾_k = 1/k!, c⁽²⁾_k = 2ᵏ/k!.
        order = 12
        jet = vector_taylor_coefficients(
            (z, y) -> [y[1], 2 * y[2]], 0.0, [1.0, 1.0], order)
        @test length(jet) == 2
        @test length(jet[1]) == order + 1
        @test length(jet[2]) == order + 1
        fact = 1.0
        for k in 0:order
            k > 0 && (fact *= k)
            @test isapprox(jet[1][k+1], 1 / fact;        rtol = 1e-13, atol = 1e-300)
            @test isapprox(jet[2][k+1], 2.0^k / fact;    rtol = 1e-13, atol = 1e-300)
        end
    end

    @testset "VC.1.3 coupled system y₁'=y₂, y₂'=−y₁" begin
        # Genuinely coupled.  With y0 = [1, 0]: y₁(h) = cos h, y₂(h) = −sin h.
        # cos:  c_k = 0 (odd k), (−1)^(k/2)/k! (even k).
        # −sin: c_k = 0 (even k), −(−1)^((k−1)/2)/k! (odd k).
        order = 16
        jet = vector_taylor_coefficients(
            (z, y) -> [y[2], -y[1]], 0.0, [1.0, 0.0], order)
        @test length(jet) == 2
        fact = 1.0
        for k in 0:order
            k > 0 && (fact *= k)
            cos_k = iseven(k) ? (-1.0)^(k ÷ 2) / fact : 0.0
            msin_k = isodd(k) ? -(-1.0)^((k - 1) ÷ 2) / fact : 0.0
            @test isapprox(jet[1][k+1], cos_k;  rtol = 1e-12, atol = 1e-14)
            @test isapprox(jet[2][k+1], msin_k; rtol = 1e-12, atol = 1e-14)
        end
    end

    @testset "VC.1.4 BigFloat type stability, no widening" begin
        # Repeat the coupled cos/−sin system at 256-bit BigFloat; assert the
        # element type of every returned vector is exactly BigFloat and the
        # high-order coefficients are accurate to 2⁻²⁰⁰.
        setprecision(BigFloat, 256) do
            order = 40
            y0 = BigFloat[1, 0]
            jet = vector_taylor_coefficients(
                (z, y) -> [y[2], -y[1]], BigFloat(0), y0, order)
            @test length(jet) == 2
            @test eltype(jet[1]) == BigFloat
            @test eltype(jet[2]) == BigFloat
            tol = BigFloat(2)^(-200)
            fact = BigFloat(1)
            for k in 0:order
                k > 0 && (fact *= k)
                cos_k = iseven(k) ? BigFloat(-1)^(k ÷ 2) / fact : BigFloat(0)
                msin_k = isodd(k) ? -BigFloat(-1)^((k - 1) ÷ 2) / fact : BigFloat(0)
                scale = BigFloat(2)^(-256)
                @test abs(jet[1][k+1] - cos_k) / max(abs(cos_k), scale) < tol
                @test abs(jet[2][k+1] - msin_k) / max(abs(msin_k), scale) < tol
            end
        end
    end

    @testset "VC.1.5 failure modes throw informatively" begin
        rhs = (z, y) -> [y[1], 2 * y[2]]
        # order < 1 — a zero-order expansion carries no info about f.
        @test_throws ArgumentError vector_taylor_coefficients(rhs, 0.0, [1.0, 1.0], 0)
        @test_throws ArgumentError vector_taylor_coefficients(rhs, 0.0, [1.0, 1.0], -3)
        # empty y0 — no components, nothing to integrate.
        @test_throws ArgumentError vector_taylor_coefficients(
            (z, y) -> Float64[], 0.0, Float64[], 5)
        # messages are informative (mention the offending input).
        err_ord = try
            vector_taylor_coefficients(rhs, 0.0, [1.0, 1.0], 0)
        catch e
            e
        end
        @test occursin("order", err_ord.msg)
        err_emp = try
            vector_taylor_coefficients((z, y) -> Float64[], 0.0, Float64[], 5)
        catch e
            e
        end
        @test occursin("y0", err_emp.msg) || occursin("empty", err_emp.msg)
    end

end

#=
VC.1.6 mutation-proof procedure — VERIFIED 2026-05-20
─────────────────────────────────────────────────────
Two independent mutations to `src/VectorCoefficients.jl`, each confirmed
to RED the suite, then restored.  Baseline: 158 pass / 158.

  Mutation A — off-by-one in the integration step:
      y[i][k] = f_t[i][k-1] / k   →   y[i][k] = f_t[i][k-1] / (k + 1)
    VERIFIED: 82 of 158 assertions RED (76 pass).  REDs VC.1.1 (vs the
    scalar oracle), VC.1.2, VC.1.3, VC.1.4 — every coefficient at k ≥ 1
    is wrong by the factor k/(k+1).  Confirms the integration recurrence
    is genuinely tested.

  Mutation B — update only the first component (loop `for i in 1:1`
  instead of `for i in 1:d`):
    VERIFIED: 68 of 158 assertions RED (90 pass).  VC.1.1 stays fully
    GREEN (d = 1, single component — Mutation A alone could be satisfied
    by a one-component impl); VC.1.2 / VC.1.3 / VC.1.4 go RED on the
    second component (jet[2] never advances past its constant term).
    Proves the component-wise bootstrap is genuinely exercised.
=#
