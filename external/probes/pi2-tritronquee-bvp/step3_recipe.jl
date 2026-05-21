# Step 3 — find a BC / continuation recipe that selects the tritronquee.
#
# Step 2 finding: a 2+2 split pinning only (u,u') at each end has the
# tritronquee in its solution set but ALSO a wrong-branch solution; plain
# Newton from the asymptotic seed lands on the wrong branch (case B) or
# diverges (cases A,C,D). The BC residual of the wrong branch is 0 — it is
# a genuine P_I^(2) solution, just not the tritronquee.
#
# The cure is to give the BC MORE information so the tritronquee is the
# unique solution in a wide Newton basin: pin all 4 asymptotic components
# (u,u',u'',u''') at the OUTER (most negative) end where the asymptotic
# series is most accurate, plus enough at the inner end. We test:
#   R1: 4+0  — pin full asymptotic state at z_a only (B_a=I, B_b=0)
#   R2: 3+1  — pin u,u',u'' at z_a and u at z_b
#   R3: continuation in interval length, warm-starting Newton
# and cross-check against the KKG asymptotic.

using PadeTaylor.PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP: vector_bvp_solve, VectorBVPSolution
using LinearAlgebra

t   = 0.0
f   = painleve_hierarchy(:I, 2; t = t)
Jf  = painleve_hierarchy_jacobian(:I, 2; t = t)
u_kkg(x)   = -cbrt(6.0) * abs(x)^(1//3)
seed(z; n) = pI2_tritronquee_ic(z; t = t, n_terms = n)

trunc_msg(e) = first(sprint(showerror, e), 140)

function report(sol, za, zb; label = "")
    xs = range(za, zb; length = 13)
    unum = [real(sol(x)[1]) for x in xs]
    ukkg = [u_kkg(x) for x in xs]
    err  = maximum(abs.(unum .- ukkg))
    println("[$label] iters=$(sol.iterations) ODEres=$(sol.residual_inf)")
    println("  u num: ", round.(unum, digits=4))
    println("  u KKG: ", round.(ukkg, digits=4))
    println("  max|u_num - u_KKG| = $err  -> ",
            err < 5e-2 ? "TRITRONQUEE BRANCH" : "WRONG BRANCH")
    return err
end

println("="^70)
println("STEP 3 — recipe search for the P_I^(2) tritronquee")
println("="^70)

# --- R1: 4+0, pin full asymptotic state at z_a (outer end) ---------------
println("\n### R1: 4+0 split, B_a = I, B_b = 0, g = seed(z_a)")
for (za, zb, N, nt) in [(-20.0,-2.0,80,2), (-20.0,-2.0,128,2),
                        (-40.0,-2.0,128,2), (-20.0,-10.0,80,2)]
    try
        ya = seed(za; n = nt)
        Ba = Matrix{Float64}(I,4,4); Bb = zeros(4,4)
        sol = vector_bvp_solve(f, za, zb, Ba, Bb, ya;
                               N=N, maxiter=60, jacobian=Jf,
                               initial_guess = z->seed(z;n=nt))
        report(sol, za, zb; label="R1 [$za,$zb] N=$N")
    catch e
        println("[R1 [$za,$zb] N=$N] threw: ", trunc_msg(e))
    end
end

# --- R2: 3+1, pin u,u',u'' at z_a and u at z_b ---------------------------
println("\n### R2: 3+1 split, pin (u,u',u'') at z_a, u at z_b")
for (za, zb, N, nt) in [(-20.0,-2.0,80,2), (-20.0,-2.0,128,2),
                        (-20.0,-10.0,80,2)]
    try
        ya = seed(za; n = nt); yb = seed(zb; n = nt)
        Ba = zeros(4,4); Bb = zeros(4,4)
        Ba[1,1]=1; Ba[2,2]=1; Ba[3,3]=1; Bb[4,1]=1
        g = [ya[1], ya[2], ya[3], yb[1]]
        sol = vector_bvp_solve(f, za, zb, Ba, Bb, g;
                               N=N, maxiter=60, jacobian=Jf,
                               initial_guess = z->seed(z;n=nt))
        report(sol, za, zb; label="R2 [$za,$zb] N=$N")
    catch e
        println("[R2 [$za,$zb] N=$N] threw: ", trunc_msg(e))
    end
end

# --- R3: continuation in interval length, 2+2 split, warm-started --------
# Start on a SHORT segment near the outer end (where the tritronquee basin
# is wide because the seed is accurate), then lengthen the inner endpoint
# progressively, warm-starting Newton from the previous solution.
println("\n### R3: interval-lengthening continuation, 2+2 (u,u') split")
let
    za  = -20.0
    nt  = 2
    zbs = [-19.0, -18.0, -16.0, -14.0, -12.0, -10.0, -8.0, -6.0, -4.0, -2.0]
    N   = 80
    Ba = zeros(4,4); Bb = zeros(4,4)
    Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
    prevsol = nothing
    for zb in zbs
        ya = seed(za; n=nt); yb = seed(zb; n=nt)
        g  = [ya[1], ya[2], yb[1], yb[2]]
        # warm start: previous solution where defined, seed elsewhere
        guess = prevsol === nothing ? (z->seed(z;n=nt)) :
            (z -> (real(z) >= real(prevsol.z_a) && real(z) <= real(prevsol.z_b)) ?
                   prevsol(z) : seed(z;n=nt))
        try
            sol = vector_bvp_solve(f, za, zb, Ba, Bb, g;
                                   N=N, maxiter=60, jacobian=Jf,
                                   initial_guess = guess)
            e = report(sol, za, zb; label="R3 cont -> zb=$zb")
            prevsol = sol
        catch ex
            println("[R3 zb=$zb] threw: ", trunc_msg(ex))
            break
        end
    end
end
