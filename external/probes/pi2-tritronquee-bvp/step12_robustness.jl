# Step 12 — robustness of the correct-sign tritronquee recipe.
#
# Step 11 found the recipe: plain 2+2 (u,u') endpoint BVP with the
# CORRECT-SIGN asymptotic seed (u = +cbrt(6)|x|^{1/3}+r^{-2}/36 on x<0).
# Converges in 3 Newton iters, KKGerr ~ 2e-4, ODEres ~ 1e-11.
#
# This step stresses it:
#  (a) wide N sweep + segment sweep, autodiff Jacobian (no analytic Jf),
#  (b) ZERO interior guess (only endpoints seeded) -- does it still land
#      on the tritronquee, or does the basin need the interior guess?
#  (c) a linear-interpolation interior guess (KKG's alternative),
#  (d) cross-check the converged u against more asymptotic series terms,
#  (e) does it survive a 3+1 / 4+0 BC split now that the sign is right?

using PadeTaylor.PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian
using PadeTaylor.VectorBVP: vector_bvp_solve
using LinearAlgebra

t  = 0.0
f  = painleve_hierarchy(:I, 2; t = t)
Jf = painleve_hierarchy_jacobian(:I, 2; t = t)
const C6 = cbrt(6.0)

# correct-sign seed (x<0, r=-x): u = C6 r^{1/3} + (1/36) r^{-2}, n_terms=2
function tri_seed(x; n_terms=2)
    r = -x
    g1 = -(1/3)*r^(-2/3)
    g2 = -(-1)*((1/3)*(-2/3)*r^(-5/3))
    g3 = -(-1)*((1/3)*(-2/3)*(-5/3)*r^(-8/3))
    u   = C6*r^(1/3); up = C6*g1; upp = C6*g2; uppp = C6*g3
    if n_terms >= 2
        h0=(1/36)*r^(-2); h1=-(1/36)*(-2)*r^(-3)
        h2=-(-1)*((1/36)*(-2)*(-3)*r^(-4))
        h3=-(-1)*((1/36)*(-2)*(-3)*(-4)*r^(-5))
        u+=h0; up+=h1; upp+=h2; uppp+=h3
    end
    [u,up,upp,uppp]
end
u_tri(x) = C6*abs(x)^(1/3) + (1/36)*abs(x)^(-2)

println("="^70)
println("STEP 12 — robustness of the correct-sign tritronquee recipe")
println("="^70)

# --- (a) wide N + segment sweep, autodiff Jacobian ----------------------
println("\n### (a) N/segment sweep, AUTODIFF Jacobian (no analytic Jf)")
for (za,zb) in [(-20.0,-2.0),(-40.0,-3.0),(-60.0,-5.0),(-15.0,-2.0),
                (-25.0,-8.0)]
    for N in (40,80,160,256)
        Ba=zeros(4,4);Bb=zeros(4,4)
        Ba[1,1]=1;Ba[2,2]=1;Bb[3,1]=1;Bb[4,2]=1
        ya=tri_seed(za);yb=tri_seed(zb)
        g=[ya[1],ya[2],yb[1],yb[2]]
        try
            sol=vector_bvp_solve(f,za,zb,Ba,Bb,g; N=N,tol=1e-9,maxiter=40,
                                 initial_guess=z->tri_seed(z))  # autodiff
            e=maximum(abs(real(sol(x)[1])-u_tri(x)) for x in range(za,zb,length=25))
            print("  [$za,$zb] N=$N it=$(sol.iterations) KKGerr=",
                  round(e,sigdigits=3), e<5e-3 ? " OK" : " ** WRONG **","   ")
        catch ex
            print("  [$za,$zb] N=$N THREW   ")
        end
    end
    println()
end

# --- (b) ZERO interior guess --------------------------------------------
println("\n### (b) ZERO interior guess (endpoints seeded only)")
for (za,zb,N) in [(-20.0,-2.0,80),(-20.0,-2.0,128),(-40.0,-3.0,160)]
    Ba=zeros(4,4);Bb=zeros(4,4)
    Ba[1,1]=1;Ba[2,2]=1;Bb[3,1]=1;Bb[4,2]=1
    ya=tri_seed(za);yb=tri_seed(zb)
    g=[ya[1],ya[2],yb[1],yb[2]]
    try
        sol=vector_bvp_solve(f,za,zb,Ba,Bb,g; N=N,tol=1e-9,maxiter=80,jacobian=Jf)
        e=maximum(abs(real(sol(x)[1])-u_tri(x)) for x in range(za,zb,length=25))
        umin=minimum(real(sol(x)[1]) for x in range(za,zb,length=25))
        println("  [$za,$zb] N=$N it=$(sol.iterations) KKGerr=",round(e,sigdigits=4),
                " umin=",round(umin,digits=3),"  ",
                e<5e-3 ? "*** TRITRONQUEE ***" : "WRONG (basin needs interior guess)")
    catch ex
        println("  [$za,$zb] N=$N threw: ", first(sprint(showerror,ex),100))
    end
end

# --- (c) linear-interpolation interior guess (KKG's alternative) --------
println("\n### (c) linear-interpolation interior guess")
for (za,zb,N) in [(-20.0,-2.0,128),(-40.0,-3.0,160)]
    Ba=zeros(4,4);Bb=zeros(4,4)
    Ba[1,1]=1;Ba[2,2]=1;Bb[3,1]=1;Bb[4,2]=1
    ya=tri_seed(za);yb=tri_seed(zb)
    g=[ya[1],ya[2],yb[1],yb[2]]
    lin = z -> ya .+ (yb.-ya).*((real(z)-za)/(zb-za))
    try
        sol=vector_bvp_solve(f,za,zb,Ba,Bb,g; N=N,tol=1e-9,maxiter=80,
                             jacobian=Jf, initial_guess=lin)
        e=maximum(abs(real(sol(x)[1])-u_tri(x)) for x in range(za,zb,length=25))
        println("  [$za,$zb] N=$N it=$(sol.iterations) KKGerr=",round(e,sigdigits=4),
                "  ", e<5e-3 ? "*** TRITRONQUEE ***" : "WRONG")
    catch ex
        println("  [$za,$zb] N=$N threw: ", first(sprint(showerror,ex),100))
    end
end

# --- (d) BigFloat solve -------------------------------------------------
println("\n### (d) BigFloat solve, [-20,-2] N=80")
let
    setprecision(BigFloat, 113)
    za,zb = big"-20.0", big"-2.0"
    tb = big"0.0"
    fb  = painleve_hierarchy(:I,2;t=tb)
    Jfb = painleve_hierarchy_jacobian(:I,2;t=tb)
    C6b = cbrt(big"6.0")
    function tri_seed_b(x)
        r=-x
        g1=-(big(1)/3)*r^(big(-2)/3)
        g2= ((big(1)/3)*(big(2)/3)*r^(big(-5)/3))
        g3=-((big(1)/3)*(big(2)/3)*(big(5)/3)*r^(big(-8)/3))
        u=C6b*r^(big(1)/3)+(big(1)/36)*r^(big(-2))
        up=C6b*g1+(-(big(1)/36)*(-2)*r^(big(-3)))
        upp=C6b*g2+((big(1)/36)*(big(2))*(big(3))*r^(big(-4)))
        uppp=C6b*g3+(-(big(1)/36)*(big(2))*(big(3))*(big(4))*r^(big(-5)))
        [u,up,upp,uppp]
    end
    Ba=zeros(BigFloat,4,4);Bb=zeros(BigFloat,4,4)
    Ba[1,1]=1;Ba[2,2]=1;Bb[3,1]=1;Bb[4,2]=1
    ya=tri_seed_b(za);yb=tri_seed_b(zb)
    g=[ya[1],ya[2],yb[1],yb[2]]
    try
        sol=vector_bvp_solve(fb,za,zb,Ba,Bb,g; N=80,maxiter=40,jacobian=Jfb,
                             initial_guess=z->tri_seed_b(z))
        utri_b(x)=C6b*abs(x)^(big(1)/3)+(big(1)/36)*abs(x)^(big(-2))
        e=maximum(abs(real(sol(x)[1])-utri_b(BigFloat(x))) for x in range(-20.0,-2.0,length=25))
        println("  BigFloat: it=$(sol.iterations) ODEres=",Float64(sol.residual_inf),
                " KKGerr=",Float64(e),"  ", e<5e-3 ? "*** TRITRONQUEE ***" : "WRONG")
    catch ex
        println("  BigFloat threw: ", first(sprint(showerror,ex),120))
    end
end

# --- (e) 3+1 and 4+0 BC splits ------------------------------------------
println("\n### (e) 3+1 and 4+0 BC splits with correct-sign seed")
let
    za,zb,N = -20.0,-2.0,128
    ya=tri_seed(za);yb=tri_seed(zb)
    # 3+1: (u,u',u'') at z_a, u at z_b
    Ba=zeros(4,4);Bb=zeros(4,4)
    Ba[1,1]=1;Ba[2,2]=1;Ba[3,3]=1;Bb[4,1]=1
    g=[ya[1],ya[2],ya[3],yb[1]]
    try
        sol=vector_bvp_solve(f,za,zb,Ba,Bb,g;N=N,tol=1e-9,maxiter=80,jacobian=Jf,
                             initial_guess=z->tri_seed(z))
        e=maximum(abs(real(sol(x)[1])-u_tri(x)) for x in range(za,zb,length=25))
        println("  3+1: it=$(sol.iterations) KKGerr=",round(e,sigdigits=4),
                "  ", e<5e-3 ? "*** TRITRONQUEE ***" : "WRONG")
    catch ex
        println("  3+1 threw: ", first(sprint(showerror,ex),100))
    end
    # 4+0: full state at z_a
    Ba=Matrix{Float64}(I,4,4);Bb=zeros(4,4)
    try
        sol=vector_bvp_solve(f,za,zb,Ba,Bb,ya;N=N,tol=1e-9,maxiter=80,jacobian=Jf,
                             initial_guess=z->tri_seed(z))
        e=maximum(abs(real(sol(x)[1])-u_tri(x)) for x in range(za,zb,length=25))
        println("  4+0: it=$(sol.iterations) KKGerr=",round(e,sigdigits=4),
                "  ", e<5e-3 ? "*** TRITRONQUEE ***" : "WRONG")
    catch ex
        println("  4+0 threw: ", first(sprint(showerror,ex),100))
    end
end
