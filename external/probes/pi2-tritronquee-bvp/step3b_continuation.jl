# Step 3b — continuation with relaxed tol + finer schedule.
#
# Step 3 finding: R3 continuation NEARLY converged on the shortest step
# (zb=-19: ‖ΔY‖ stalled at 6e-11, just above the default tol 1.8e-12).
# That is a near-converged tritronquee, not a divergence. The default
# step-norm tol = eps^(3/4) ~ 1.8e-12 is unreachable here because the
# collocation residual floors at cond(D1)*eps for this stiff operator.
#
# Recipe under test: 2+2 (u,u') split, relaxed tol ~ 1e-9, generous
# maxiter, fine interval-lengthening continuation warm-started from the
# previous VectorBVPSolution. Check each converged solution against KKG.

using PadeTaylor.PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP: vector_bvp_solve
using LinearAlgebra

t   = 0.0
f   = painleve_hierarchy(:I, 2; t = t)
Jf  = painleve_hierarchy_jacobian(:I, 2; t = t)
u_kkg(x)   = -cbrt(6.0) * abs(x)^(1//3)
seed(z; n) = pI2_tritronquee_ic(z; t = t, n_terms = n)

function kkg_err(sol, za, zb)
    xs = range(za, zb; length = 25)
    maximum(abs(real(sol(x)[1]) - u_kkg(x)) for x in xs)
end

println("="^70)
println("STEP 3b — relaxed-tol continuation for the P_I^(2) tritronquee")
println("="^70)

# --- first, a single short-interval solve with relaxed tol --------------
println("\n### Single short solve [-20,-19], 2+2 (u,u') split, tol=1e-9")
let
    za, zb, nt = -20.0, -19.0, 2
    Ba = zeros(4,4); Bb = zeros(4,4)
    Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
    ya = seed(za;n=nt); yb = seed(zb;n=nt)
    g  = [ya[1], ya[2], yb[1], yb[2]]
    for (N, tol) in [(60,1e-9),(60,1e-8),(100,1e-9),(40,1e-9)]
        try
            sol = vector_bvp_solve(f, za, zb, Ba, Bb, g;
                                   N=N, tol=tol, maxiter=100, jacobian=Jf,
                                   initial_guess = z->seed(z;n=nt))
            e = kkg_err(sol, za, zb)
            println("  N=$N tol=$tol -> iters=$(sol.iterations) ",
                    "ODEres=$(sol.residual_inf) KKGerr=$e ",
                    e < 5e-3 ? "TRITRONQUEE" : "WRONG")
        catch ex
            println("  N=$N tol=$tol threw: ", first(sprint(showerror,ex),120))
        end
    end
end

# --- continuation: lengthen inner endpoint, warm-started ----------------
function continuation(zbs; za=-20.0, N=80, tol=1e-9, nt=2, maxiter=120)
    Ba = zeros(4,4); Bb = zeros(4,4)
    Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
    prevsol = nothing
    finalsol = nothing
    for zb in zbs
        ya = seed(za;n=nt); yb = seed(zb;n=nt)
        g  = [ya[1], ya[2], yb[1], yb[2]]
        guess = prevsol === nothing ? (z->seed(z;n=nt)) :
            (z -> (real(z) >= -20.0 && real(z) <= real(prevsol.z_b)) ?
                   prevsol(z) : seed(z;n=nt))
        try
            sol = vector_bvp_solve(f, za, zb, Ba, Bb, g;
                                   N=N, tol=tol, maxiter=maxiter, jacobian=Jf,
                                   initial_guess = guess)
            e = kkg_err(sol, za, zb)
            println("  zb=$zb iters=$(sol.iterations) ODEres=",
                    round(sol.residual_inf,sigdigits=3),
                    " KKGerr=", round(e,sigdigits=4), "  ",
                    e < 5e-3 ? "TRITRONQUEE" : "WRONG-BRANCH")
            prevsol = sol; finalsol = sol
        catch ex
            println("  zb=$zb threw: ", first(sprint(showerror,ex),120))
            return finalsol
        end
    end
    return finalsol
end

println("\n### Continuation: fine schedule, N=80, tol=1e-9")
sched = collect(-20.0:0.5:-2.0)[2:end]
sol = continuation(sched; N=80, tol=1e-9)

if sol !== nothing
    println("\n### Final converged solution diagnostics")
    println("  segment [$(sol.z_a), $(sol.z_b)], N=$(sol.N)")
    xs = -20.0:2.0:-2.0
    for x in xs
        un = real(sol(x)[1]); uk = u_kkg(x)
        println("  x=$x  u_num=", round(un,digits=5),
                "  u_KKG=", round(uk,digits=5),
                "  |Δ|=", round(abs(un-uk),sigdigits=3))
    end
end
