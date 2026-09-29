"""
Component order must preserve shared-denominator Padé values (padetaylor-brbv).

The six known-correct counterexamples from worklog 083:88-108 use the exact
geometric function 1/(1-t/2) and the rotated Noumi–Yamada Type A seed (0,t,0,0,0).
They exercise the primitive, one step, both step policies and dense output.
Every cyclic rotation of several rational, polynomial and full-column-rank
jets must preserve the normalised denominator and component values to roundoff.
Literal rational formulas supply an oracle independent of SVD/QR (Rules 4–5).
See ADR-0036 and worklog 088 for the measured QR residual and mutation proof.
"""

using Test, PadeTaylor
using PadeTaylor.SharedPade: shared_denominator_pade
using PadeTaylor.NoumiYamada: NoumiYamadaProblem
using PadeTaylor.VectorStepper: VectorPadeStepperState, vector_pade_step_with_pade!
using PadeTaylor.VectorProblems: vector_solve_pade

# Formal division of literal polynomials, independent of the Padé construction.
function order_invariance_jet(p, q, order)
    c = zeros(eltype(q), order + 1)
    for k in 0:order
        c[k+1] = k < length(p) ? p[k+1] : zero(eltype(q))
        for j in 1:min(k, length(q)-1)
            c[k+1] -= q[j+1] * c[k-j+1]
        end
        c[k+1] /= q[1]
    end
    c
end

@testset "SharedPade component-order invariance (brbv)" begin
    @testset "Worklog 083 six known-correct counterexamples" begin
        live = [0.5^k for k in 0:6]
        nums, den = shared_denominator_pade([zeros(7), live], 3)
        @test sum(nums[2]) / sum(den) ≈ 2.0 atol=1e-14 rtol=0
        prob = NoumiYamadaProblem(2; α=[0.,1.,0.,0.,0.], f0=zeros(5),
                                 tspan=(0.,1.), order=6)
        st = VectorPadeStepperState{Float64}(0.,zeros(5))
        vector_pade_step_with_pade!(st,prob.problem.f,6,.125)
        @test st.y ≈ [0.,.125,0.,0.,0.] atol=1e-14 rtol=0
        for policy in (:fixed,:jorba_zou)
            sol = vector_solve_pade(prob;h=.125,step_policy=policy)
            @test sol.y[end] ≈ [0.,1.,0.,0.,0.] atol=1e-14 rtol=0
            @test sol(.0625) ≈ [0.,.0625,0.,0.,0.] atol=1e-14 rtol=0
        end
    end

    @testset "Every rotation preserves Q and component values" begin
        for T in (Float64, ComplexF64)
            max_den_error = max_value_error = max_oracle_error = 0.0
            # One/two poles, zero blocks, duplicate constraints, polynomials,
            # and a well-conditioned least-squares stack with no exact null.
            q1 = T[1, -0.5]
            q2 = T[1, -0.75, 0.125]
            T === ComplexF64 && (q2 = T[1, -0.5-0.25im, 0.125im])
            cases = ((q1, [T[1], T[0]], 3),
                     (q1, [T[1], T[0], T[2], T[0], T[-1, 0.25]], 3),
                     (q2, [T[1, 0.25], T[0], T[2, -0.5]], 4),
                     (T[1], [T[0, 1], T[0], T[2, -0.5]], 3))
            for (q, ps, m) in cases
                jets = [order_invariance_jet(p, q, 2m) for p in ps]
                base_nums, base_den = shared_denominator_pade(jets, m)
                @test base_den ≈ q atol=64eps() rtol=64eps()
                for rotation in 0:length(jets)-1
                    perm = circshift(collect(eachindex(jets)), rotation)
                    nums, den = shared_denominator_pade(jets[perm], m)
                    @test den ≈ base_den atol=64eps() rtol=64eps()
                    length(den) == length(base_den) &&
                        (max_den_error = max(max_den_error, maximum(abs, den .- base_den)))
                    for (i, original) in enumerate(perm), t in (-0.5, 0.25, 1.0)
                        value = evalpoly(t, nums[i]) / evalpoly(t, den)
                        baseline = evalpoly(t, base_nums[original]) / evalpoly(t, base_den)
                        oracle = evalpoly(t, ps[original]) / evalpoly(t, q)
                        @test value ≈ baseline atol=64eps() rtol=64eps()
                        @test value ≈ oracle atol=64eps() rtol=64eps()
                        max_value_error = max(max_value_error, abs(value-baseline))
                        max_oracle_error = max(max_oracle_error, abs(value-oracle))
                    end
                end
            end
            jets = [T[1, 2, 3, 4, 5], T[2, -1, 1, -2, 3], T[0, 1, -1, 2, -3]]
            base_nums, base_den = shared_denominator_pade(jets, 2)
            for rotation in 0:length(jets)-1
                perm = circshift(collect(eachindex(jets)), rotation)
                nums, den = shared_denominator_pade(jets[perm], 2)
                @test den ≈ base_den atol=64eps() rtol=64eps()
                length(den) == length(base_den) &&
                    (max_den_error = max(max_den_error, maximum(abs, den .- base_den)))
                for (i, original) in enumerate(perm), t in (-0.5, 0.25, 1.0)
                    value = evalpoly(t, nums[i]) / evalpoly(t, den)
                    baseline = evalpoly(t, base_nums[original]) / evalpoly(t, base_den)
                    @test value ≈ baseline atol=64eps() rtol=64eps()
                    max_value_error = max(max_value_error, abs(value-baseline))
                end
            end
            @info "brbv rotation errors" T rotations=16 max_den_error max_value_error max_oracle_error
        end
    end
end
