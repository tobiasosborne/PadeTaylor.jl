# Investigate why VPN.3.1's fixed-h negative control no longer throws under the
# ADR-0028 dispatch.  The test @test_throws on the fixed-h :min_y walk (the A2
# pole-landing failure mode).  Did the dispatch genuinely make it thread (better
# cell-B values ⇒ different :min_y trajectory ⇒ avoids the degenerate pole), or
# does it silently return a pole-degenerate lie (a Rule-1 break)?  Reproduce both
# walks and inspect.
using PadeTaylor
using PadeTaylor: painleve_hierarchy, pI2_tritronquee_ic, vector_path_network_solve,
                  VectorPathNetworkSolution

f = painleve_hierarchy(:I, 2; t = 0.0)
y_seed = ComplexF64.(pI2_tritronquee_ic(-3.0 + 0.0im; t = 0.0, n_terms = 2))
prob = VectorPadeTaylorProblem(f, y_seed, (-3.0 + 0im, 20.0 + 0im); order = 24)
fan = ComplexF64[r * cis(a) for r in 2.0:2.0:18.0 for a in range(-0.4, 0.4; length = 5)]

function report(label, walkf)
    print(label, ": ")
    try
        sol = walkf()
        nz = length(sol.visited_z)
        maxx = maximum(real, sol.visited_z)
        finite = all(all(isfinite, y) for y in sol.visited_y)
        # how many fan targets reached (a node within 0.3 of the target)?
        reached = count(t -> any(abs(z - t) < 0.3 for z in sol.visited_z), fan)
        println("RETURNED  nodes=$nz  maxRe(x)=", round(maxx, digits = 2),
                "  reached=$reached/$(length(fan))  all_finite=$finite")
    catch e
        println("THREW  ", typeof(e), "  ", first(sprint(showerror, e), 90))
    end
end

report("fixed-h  :min_y  (the negative control)",
       () -> vector_path_network_solve(prob, fan; order = 24, h = 0.1,
                                       adaptive = false, step_policy = :min_y))
report("adaptive :max_q_root (the B2 walk)",
       () -> vector_path_network_solve(prob, fan; order = 24, h = 0.3,
                                       adaptive = true, step_policy = :max_q_root))
