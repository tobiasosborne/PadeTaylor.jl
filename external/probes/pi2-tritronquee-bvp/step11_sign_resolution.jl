# Step 11 — RESOLVE THE SIGN. The prior probe's premise is inverted.
#
# The P_I^(2) leading balance is  u'''' + 10u'^2 + 20uu'' + 40(u^3+6x) = 0
# => to leading order 40u^3 + 240x = 0  => u^3 = -6x.
# For x<0:  u^3 = -6x = 6|x| > 0  =>  the REAL root is  u = +cbrt(6)|x|^{1/3}
# i.e. on the negative real axis the genuine real P_I^(2) solution is
# POSITIVE. The "-branch" u = -cbrt(6)|x|^{1/3} has ODE residual ~ -4800
# at x=-10 (step10 check) -- it is NOT a solution at all.
#
# `pI2_tritronquee_ic` returns -cbrt(6)*|x|^{1/3} (negative). KKG's
# asymptotic is u ~ -cbrt(6)*x^{1/3}; for x<0 the real cube root x^{1/3}
# is NEGATIVE, so -cbrt(6)*x^{1/3} = +cbrt(6)*|x|^{1/3} -- POSITIVE.
# The src builder used |x|^{1/3} in place of x^{1/3} => sign slip for x<0.
#
# This step: build the CORRECT-SIGN asymptotic seed (u>0 on x<0) and its
# companion derivatives, run the plain 2+2 (u,u') BVP, and check it lands
# on a CLEAN monotone positive tritronquee with small ODE residual.
#
# Correct asymptotic at t=0:  u_f = sum b_n x^{-(7n-1)/3},
#   b_0=-cbrt(6), b_1=1/36  (tritronquee_coeff.tex:2686-2688).
# On x<0 write x = -r, r=|x|>0. Using the real branch consistent with the
# ODE, the leading term magnitude is cbrt(6) r^{1/3} and SIGN +:
#   u_a(x) = +cbrt(6) * r^{1/3}  + correction.
# Correction: u = Y + Y^{-6} with Y = +cbrt(6) r^{1/3} (so that Y^3=6r,
# matching u^3=6r). Y^{-6} = 6^{-2} r^{-2} = r^{-2}/36 -- consistent with
# b_1 = 1/36. d/dx = -d/dr.

using PadeTaylor.PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian
using PadeTaylor.VectorBVP: vector_bvp_solve
using LinearAlgebra

t   = 0.0
f   = painleve_hierarchy(:I, 2; t = t)
Jf  = painleve_hierarchy_jacobian(:I, 2; t = t)
const C6 = cbrt(6.0)

# correct-sign tritronquee asymptotic on x<0 (r=-x):
#   u = C6 r^{1/3} + (1/36) r^{-2}
# derivatives wrt x (d/dx = -d/dr):
#   u'   = -(d/dr)u
#   etc.
function tri_seed(x; n_terms=2)
    r = -x                       # > 0
    # leading Y = C6 r^{1/3}
    u   =  C6*r^(1//3)
    up  = -C6*(1//3)*r^(-2//3)            # dY/dx = -dY/dr
    upp = -C6*(1//3)*(-2//3)*r^(-5//3) * (-1)   # see below; compute cleanly
    # compute cleanly via repeated d/dx = -d/dr on  g(r)=r^{1/3}
    # g    = r^{1/3}
    # dg/dr= (1/3) r^{-2/3} ; d/dx g = -(1/3) r^{-2/3}
    # d2/dx2 g = d/dx(-(1/3)r^{-2/3}) = -(-1)*(1/3)(-2/3) r^{-5/3} = (1/3)(-2/3)(-1)... be explicit:
    g1 = -(1/3)*r^(-2/3)                              # d/dx r^{1/3}
    g2 = -(-1)*( (1/3)*(-2/3)*r^(-5/3) )              # d/dx g1 = -d/dr g1
    g3 = -(-1)*( (1/3)*(-2/3)*(-5/3)*r^(-8/3) )
    u   = C6*r^(1/3)
    up  = C6*g1
    upp = C6*g2
    uppp= C6*g3
    if n_terms >= 2
        # correction h(r) = (1/36) r^{-2}
        h0 = (1/36)*r^(-2)
        h1 = -(1/36)*(-2)*r^(-3)
        h2 = -(-1)*((1/36)*(-2)*(-3)*r^(-4))
        h3 = -(-1)*((1/36)*(-2)*(-3)*(-4)*r^(-5))
        u   += h0; up += h1; upp += h2; uppp += h3
    end
    return [u, up, upp, uppp]
end

# tritronquee reference value (correct sign): u ~ +C6 |x|^{1/3} for x<0
u_tri(x) = C6*abs(x)^(1/3) + (1/36)*abs(x)^(-2)

println("="^70)
println("STEP 11 — sign resolution: tritronquee is POSITIVE on x<0")
println("="^70)

# verify the correct-sign seed has small ODE residual
println("\n### ODE residual of the correct-sign asymptotic seed")
for x in (-40.0,-20.0,-10.0,-5.0)
    y = tri_seed(x; n_terms=2)
    # need u'''' ; reuse f: f returns y' = (u',u'',u''',u'''')
    fr = f(x, y)
    # the 4th component of fr is u'''' from the ODE; compare to actual u''''
    # actual u'''' by differentiating the seed once more numerically (central)
    dx=1e-4
    u4_num = (tri_seed(x+dx;n_terms=2)[4]-tri_seed(x-dx;n_terms=2)[4])/(2dx)
    println("  x=$x  u=", round(y[1],sigdigits=6),
            "  ODE-u''''=", round(fr[4],sigdigits=6),
            "  seed-u''''=", round(u4_num,sigdigits=6),
            "  mismatch=", round(abs(fr[4]-u4_num),sigdigits=3))
end

# --- plain 2+2 (u,u') BVP with the CORRECT-SIGN seed --------------------
println("\n### 2+2 (u,u') BVP, correct-sign seed at endpoints + interior guess")
function run(za,zb,N; nt=2, tol=1e-9, maxiter=120)
    Ba=zeros(4,4); Bb=zeros(4,4)
    Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
    ya=tri_seed(za;n_terms=nt); yb=tri_seed(zb;n_terms=nt)
    g=[ya[1],ya[2],yb[1],yb[2]]
    try
        sol=vector_bvp_solve(f,za,zb,Ba,Bb,g; N=N,tol=tol,maxiter=maxiter,
                             jacobian=Jf, initial_guess=z->tri_seed(z;n_terms=nt))
        xs=range(za,zb,length=25)
        e=maximum(abs(real(sol(x)[1])-u_tri(x)) for x in xs)
        umin=minimum(real(sol(x)[1]) for x in xs)
        mono=all(diff([real(sol(x)[1]) for x in xs]).<0)  # u decreasing as x increases
        println("  [$za,$zb] N=$N iters=$(sol.iterations) ODEres=",
                round(sol.residual_inf,sigdigits=3),
                " KKGerr=", round(e,sigdigits=4),
                " umin=", round(umin,digits=3),
                " monotone=", mono, "  ",
                e<5e-3 ? "*** TRITRONQUEE ***" : "wrong")
        return sol
    catch ex
        println("  [$za,$zb] N=$N threw: ", first(sprint(showerror,ex),120))
        return nothing
    end
end

for (za,zb,N) in [(-20.0,-2.0,60),(-20.0,-2.0,80),(-20.0,-2.0,128),
                  (-20.0,-2.0,200),(-40.0,-2.0,160),(-20.0,-10.0,80),
                  (-30.0,-5.0,128),(-50.0,-4.0,200),(-20.0,-1.0,128)]
    run(za,zb,N)
end

println("\n### Detailed dump: 2+2 [-20,-2] N=128, correct-sign")
let
    sol=run(-20.0,-2.0,128)
    if sol!==nothing
        for x in -20.0:2.0:-2.0
            un=real(sol(x)[1])
            println("  x=$x  u_num=", round(un,digits=6),
                    "  u_tri=", round(u_tri(x),digits=6),
                    "  |Δ|=", round(abs(un-u_tri(x)),sigdigits=3))
        end
    end
end
