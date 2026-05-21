# Step 5 — characterise the wrong branch + test the asymptotic-regime recipe.
#
# Every 2+2 split converges to the same wrong branch; 3+1 splits give a
# singular Jacobian (the separatrix growing mode leaves J rank-deficient).
# This step:
#   (a) examines the wrong branch in full -- is it a real P_I^(2) solution
#       or a pole-crossing artifact? print all 4 components along the seg.
#   (b) tests the IVP from the seed -- does forward integration diverge
#       (confirming the separatrix instability)?
#   (c) tests whether seeding the Newton guess with the EXACT tritronquee
#       (a high-N self-consistent solve) keeps it on-branch -- i.e. is the
#       tritronquee a STABLE fixed point of the collocation Newton map, or
#       does even an on-branch start slide off?

using PadeTaylor.PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic, PainleveHierarchyProblem
using PadeTaylor.VectorProblems: vector_solve_pade
using PadeTaylor.VectorBVP: vector_bvp_solve
using LinearAlgebra

t   = 0.0
f   = painleve_hierarchy(:I, 2; t = t)
Jf  = painleve_hierarchy_jacobian(:I, 2; t = t)
u_kkg(x)   = -cbrt(6.0) * abs(x)^(1//3)
seed(z; n) = pI2_tritronquee_ic(z; t = t, n_terms = n)

println("="^70)
println("STEP 5 — wrong-branch characterisation + recipe tests")
println("="^70)

# --- (a) the wrong branch, full state -----------------------------------
println("\n### (a) Wrong branch full state, [-20,-2] N=80, (u,u') 2+2")
let
    za, zb, nt = -20.0, -2.0, 2
    Ba=zeros(4,4); Bb=zeros(4,4)
    Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
    ya=seed(za;n=nt); yb=seed(zb;n=nt)
    g=[ya[1],ya[2],yb[1],yb[2]]
    sol = vector_bvp_solve(f,za,zb,Ba,Bb,g; N=80,tol=1e-9,maxiter=120,
                           jacobian=Jf, initial_guess=z->seed(z;n=nt))
    for x in -20.0:2.0:-2.0
        y = real.(sol(x))
        println("  x=$x  y=", round.(y,sigdigits=5))
    end
    # ODE residual check at midpoint: does f(x,y) match y' from D-matrix?
    println("  -> wrong branch ODEres=", sol.residual_inf,
            " (genuine solution if ~1e-9)")
end

# --- (b) forward IVP from the seed --------------------------------------
println("\n### (b) Forward IVP from the n_terms=2 seed at x0=-20")
let
    x0 = -20.0
    y0 = seed(x0; n=2)
    prob = PainleveHierarchyProblem(2; y0=y0, xspan=(x0,-2.0), t=t, order=30)
    try
        s = vector_solve_pade(prob; h=0.25, max_steps=2000)
        # report u at a few x
        println("  IVP solved; sampling u(x):")
        xs = s.x_values
        for k in 1:max(1,div(length(xs),8)):length(xs)
            println("    x=", round(xs[k],digits=3),
                    "  u=", round(real(s.y_values[k][1]),sigdigits=5))
        end
    catch ex
        println("  IVP threw (separatrix divergence?): ",
                first(sprint(showerror,ex),200))
    end
end

# --- (c) self-consistent restart: does an on-branch guess stay on? ------
# Build the KKG tritronquee numerically on [-20,-2] by GUESSING the exact
# asymptotic u_kkg at every node (not just endpoints) and see if Newton
# from that interior-correct guess converges to the tritronquee or slides.
println("\n### (c) Newton from an interior-correct (KKG asymptotic) guess")
let
    za, zb, nt = -20.0, -2.0, 2
    Ba=zeros(4,4); Bb=zeros(4,4)
    Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
    ya=seed(za;n=nt); yb=seed(zb;n=nt)
    g=[ya[1],ya[2],yb[1],yb[2]]
    # interior guess: the FULL asymptotic seed at every x (n_terms=2)
    # this is close to the tritronquee EVERYWHERE, not just endpoints
    guess = z -> seed(z; n=2)
    for N in (40, 80, 160, 256)
        try
            sol = vector_bvp_solve(f,za,zb,Ba,Bb,g; N=N,tol=1e-9,maxiter=200,
                                   jacobian=Jf, initial_guess=guess)
            e = maximum(abs(real(sol(x)[1])-u_kkg(x)) for x in range(za,zb,length=25))
            println("  N=$N iters=$(sol.iterations) ODEres=",
                    round(sol.residual_inf,sigdigits=3),
                    " KKGerr=", round(e,sigdigits=4), "  ",
                    e<5e-3 ? "*** TRITRONQUEE ***" : "wrong")
        catch ex
            println("  N=$N threw: ", first(sprint(showerror,ex),110))
        end
    end
end

# --- (d) single Newton step from the seed: which direction? -------------
println("\n### (d) Newton with maxiter=1 from the seed -- does step 1 go off?")
let
    za, zb, nt = -20.0, -19.0, 2
    Ba=zeros(4,4); Bb=zeros(4,4)
    Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
    ya=seed(za;n=nt); yb=seed(zb;n=nt)
    g=[ya[1],ya[2],yb[1],yb[2]]
    for mi in (1,2,3,5)
        try
            sol = vector_bvp_solve(f,za,zb,Ba,Bb,g; N=40,tol=1e-9,maxiter=mi,
                                   jacobian=Jf, initial_guess=z->seed(z;n=2))
            mid = real(sol((za+zb)/2)[1])
            println("  maxiter=$mi converged: u(mid)=", round(mid,sigdigits=5))
        catch ex
            # extract the final u(mid) is impossible after throw; just note
            println("  maxiter=$mi did not converge")
        end
    end
end
