"""
Zero components through the vector Padé–Taylor pipeline (padetaylor-0o9).

Matsuda's A_4 rational seeds are Type A `(t,0,0,0,0)` and Type B
`(t/3,t/3,t/3,0,0)` with the same constant vectors as their parameters:
`references/tex/noumi_yamada/Matsuda2012_rational_A4_NoumiYamada_JMP53/main.tex:324-332`.
Type A also satisfies every even-parity system by direct substitution into
NY1998's equation (`references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex:85-88`).
These literal formulas, independent of the rational-oracle implementation,
pin canonical Taylor jets, direct SharedPade construction, the dispatched
VectorStepper, and vector_solve_pade's nodes and dense output. At order 30,
cyclic rotations also put zero components before and after live siblings.
Starting at t=0 distinguishes an all-zero initial state from an all-zero
Taylor jet: the seed has a live slope.

Canonical seeds at orders 6 and 30 cover cell-B's d>m fallback and the
production (15,15) dispatch. Fixed and Jorba–Zou policies preserve these seeds.
A separate live-first meromorphic mixed jet pins a genuine shared pole.
This is a measured scope: zero-first blocks can defeat unpivoted QR in direct
SharedPade, and order-6 rotated seeds can fail end to end. Those counterexamples
and the unrelated ComplexF64 driver failure are recorded with reproduction
code in the worklog. Mutation evidence is in
`docs/worklog/083-zero-component-end-to-end.md` (Rules 4–5).
"""

using Test
using PadeTaylor
using PadeTaylor.VectorCoefficients: vector_taylor_coefficients
using PadeTaylor.SharedPade: shared_denominator_pade
using PadeTaylor.VectorStepper: VectorPadeStepperState, vector_pade_step_with_pade!
using PadeTaylor.VectorProblems: vector_solve_pade
using PadeTaylor.NoumiYamada: NoumiYamadaProblem
using PadeTaylor.NoumiYamadaSymmetry: noumi_yamada_rational

@testset "Zero components end to end (0o9)" begin
    @testset "Mixed meromorphic jets retain the known shared pole" begin
        live = [0.5^k for k in 0:6]  # 1/(1-t/2), exact binary coefficients
        for d in (2, 5)
            slot = 1
            jets = [i == slot ? copy(live) : zeros(7) for i in 1:d]
            nums, den = shared_denominator_pade(jets, 3)
            @test den ≈ [1.0, -0.5] atol = 1e-14 rtol = 0
            @test nums[slot] ≈ [1.0] atol = 1e-14 rtol = 0
            for i in setdiff(1:d, [slot])
                @test nums[i] == [0.0]
            end
            @info "0o9 mixed jet" d slot denominator=den
        end
    end

    seeds = ((1, :A, [1//1, 0//1, 0//1]),
             (2, :A, [1//1, 0//1, 0//1, 0//1, 0//1]),
             (2, :B, [1//3, 1//3, 1//3, 0//1, 0//1]))
    @testset "Type A/B Taylor, SharedPade, stepper, solver and dense values" begin
        for (n, kind, weights) in seeds
            α_oracle, seed_oracle = noumi_yamada_rational(n, kind)
            @test α_oracle == weights
            @test seed_oracle(7//3) == (7//3) .* weights
            max_step_error = max_node_error = max_dense_error = 0.0
            trajectories = segments = 0
            for T in (Float64, ComplexF64), order in (6, 30), t0 in (0.0, 1.0)
                slopes = T.(weights)
                zero_slots = findall(iszero, slopes)
                z0, h = T(t0), T(0.125)
                prob = NoumiYamadaProblem(n; α = slopes, f0 = z0 .* slopes,
                                         tspan = (z0, z0 + one(T)), order = order)
                jets = vector_taylor_coefficients(prob.problem.f, z0,
                                                  prob.problem.y0, order)
                # The full known Taylor jet includes the identically-zero legs.
                @test jets == [[z0 * c, c, zeros(T, order-1)...] for c in slopes]
                scaled = [[jet[k+1] * h^k for k in 0:order] for jet in jets]
                nums, den = shared_denominator_pade(scaled, order ÷ 2)
                @test den == T[1]
                @test nums == [iszero(c) ? T[0] : T[z0*c, h*c] for c in slopes]

                st = VectorPadeStepperState{T}(z0, z0 .* slopes)
                _, step_nums, step_den = vector_pade_step_with_pade!(
                    st, prob.problem.f, order, h)
                @test st.z == z0 + h
                @test st.y ≈ (z0 + h) .* slopes atol = 1e-14 rtol = 0
                @test step_den == T[1]
                @test step_nums[zero_slots] == [T[0] for _ in zero_slots]
                max_step_error = max(max_step_error,
                                     maximum(abs, st.y .- (z0+h) .* slopes))

                T === ComplexF64 && continue  # Driver bug is outside 0o9; see worklog.
                for policy in (:fixed, :jorba_zou)
                    sol = vector_solve_pade(prob; h = 0.125, step_policy = policy)
                    @test sol.z[end] == z0 + one(T)
                    @test length(sol.h) == 8
                    for (z, y) in zip(sol.z, sol.y)
                        @test y ≈ z .* slopes atol = 1e-14 rtol = 0
                        @test y[zero_slots] == zeros(T, length(zero_slots))
                        max_node_error = max(max_node_error,
                                             maximum(abs, y .- z .* slopes))
                    end
                    for k in eachindex(sol.h)
                        z = sol.z[k] + sol.h[k] / 2
                        y = sol(z)
                        @test y ≈ z .* slopes atol = 1e-14 rtol = 0
                        @test y[zero_slots] == zeros(T, length(zero_slots))
                        max_dense_error = max(max_dense_error,
                                              maximum(abs, y .- z .* slopes))
                    end
                    trajectories += 1
                    segments += length(sol.h)
                end
            end
            # Default-order dispatch recovers the rotated cases, although direct
            # cell A can lose their slope. Compare values, not a nonunique Q.
            for t0 in (0.0, 1.0), rotation in 1:(length(weights)-1)
                slopes = Float64.(circshift(weights, rotation))
                prob = NoumiYamadaProblem(n; α = slopes, f0 = t0 .* slopes,
                                         tspan = (t0, t0+1.0), order = 30)
                st = VectorPadeStepperState{Float64}(t0, t0 .* slopes)
                vector_pade_step_with_pade!(st, prob.problem.f, 30, 0.125)
                @test st.y ≈ (t0+0.125) .* slopes atol = 1e-14 rtol = 0
                max_step_error = max(max_step_error,
                                     maximum(abs, st.y .- (t0+0.125) .* slopes))
                for policy in (:fixed, :jorba_zou)
                    sol = vector_solve_pade(prob; h = 0.125, step_policy = policy)
                    @test sol.z[end] == t0+1.0
                    @test length(sol.h) == 8
                    for (z, y) in zip(sol.z, sol.y)
                        @test y ≈ z .* slopes atol = 1e-14 rtol = 0
                        @test y[iszero.(slopes)] == zeros(count(iszero, slopes))
                        max_node_error = max(max_node_error,
                                             maximum(abs, y .- z .* slopes))
                    end
                    for k in eachindex(sol.h)
                        z = sol.z[k] + sol.h[k]/2
                        @test sol(z) ≈ z .* slopes atol = 1e-14 rtol = 0
                        @test sol(z)[iszero.(slopes)] == zeros(count(iszero, slopes))
                        max_dense_error = max(max_dense_error,
                                              maximum(abs, sol(z) .- z .* slopes))
                    end
                    trajectories += 1
                    segments += length(sol.h)
                end
            end
            @info "0o9 rational seed" n kind trajectories segments max_step_error max_node_error max_dense_error
        end
    end
end
