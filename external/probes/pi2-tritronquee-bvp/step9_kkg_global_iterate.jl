# Step 9 — Approach A, part 1: KKG's ACTUAL global initial iterate.
#
# Re-reading KKG 2015 §7.1 (tritronquee_coeff.tex:3173-3177) the prior probe
# missed the decisive sentence:
#
#   "The resulting system of algebraic equations is solved using Newton's
#    method with the initial iterate  u = -6^{1/3} x/(1+x^2)^{1/3},  a smooth
#    function which satisfies for large x the asymptotic conditions, or a
#    linear interpolation between the boundary data."
#
# Every prior probe step seeded Newton with `pI2_tritronquee_ic` -- the
# OUTER-END asymptotic expansion u ~ -cbrt(6)|x|^{1/3}, accurate only near
# the large-|x| endpoint. Across the whole interval that vector seed is a
# poor guess and Newton slid to the +branch.
#
# KKG's iterate u = -6^{1/3} x/(1+x^2)^{1/3} is a GLOBALLY tritronquee-
# leaning smooth function. On the negative real axis the real branch the
# probe targets is u = -cbrt(6)|x|^{1/3} (negative). The sign-correct global
# iterate for that branch is
#
#       u_init(x) = -cbrt(6) * |x| / (1+x^2)^{1/3}
#
# which is negative for ALL x<0 and ~ -cbrt(6)|x|^{1/3} as |x|->infty.
# Its companion-form derivatives are taken analytically below.
#
# This step: 2+2 (u,u') BVP, seed the WHOLE interior with this global
# iterate (not the asymptotic seed), plain Newton. Test across N/segment.

using PadeTaylor.PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP: vector_bvp_solve
using LinearAlgebra

t   = 0.0
f   = painleve_hierarchy(:I, 2; t = t)
Jf  = painleve_hierarchy_jacobian(:I, 2; t = t)
u_kkg(x)   = -cbrt(6.0) * abs(x)^(1//3)
seed(z; n) = pI2_tritronquee_ic(z; t = t, n_terms = n)

# KKG global iterate, sign-corrected for the real negative-axis branch:
#   u(x) = -c * |x| * (1+x^2)^{-1/3},  c = cbrt(6)
# For x<0, |x| = -x, so u(x) = c * x * (1+x^2)^{-1/3}.
# Provide (u,u',u'',u''') analytically via that closed form.
const C6 = cbrt(6.0)
function kkg_iterate(x)
    # work with g(x) = x*(1+x^2)^{-1/3}; u = c*g  (x<0 so |x|=-x, u=-c|x|w=+c x w)
    p = (1 + x^2)
    w  = p^(-1//3)
    # derivatives of w = p^{-1/3}, p=1+x^2, p'=2x
    wp  = (-1//3) * p^(-4//3) * 2x
    wpp = (-1//3)*((-4//3)*p^(-7//3)*(2x)^2 + p^(-4//3)*2)
    wppp = (-1//3)*((-4//3)*((-7//3)*p^(-10//3)*(2x)^3 + p^(-7//3)*2*(2x)*2)
                    + (-4//3)*p^(-7//3)*2*(2x))
    g   = x*w
    gp  = w + x*wp
    gpp = 2*wp + x*wpp
    gppp = 3*wpp + x*wppp
    return C6 .* [g, gp, gpp, gppp]
end

println("="^70)
println("STEP 9 — Approach A pt1: KKG global initial iterate u=-c|x|(1+x^2)^{-1/3}")
println("="^70)

# sanity: iterate value vs KKG asymptotic
println("\n### KKG global iterate vs asymptotic, sample points")
for x in (-40.0,-20.0,-10.0,-4.0,-2.0)
    yi = kkg_iterate(x)
    println("  x=$x  u_init=", round(yi[1],sigdigits=6),
            "  u_kkg=", round(u_kkg(x),sigdigits=6),
            "  u'_init=", round(yi[2],sigdigits=5))
end

# --- 2+2 (u,u') BVP, seeded with the GLOBAL iterate everywhere -----------
println("\n### 2+2 (u,u') BVP, interior guess = KKG global iterate")
function run22(za,zb,N; nt=2, gmode=:global, tol=1e-9, maxiter=200)
    Ba=zeros(4,4); Bb=zeros(4,4)
    Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
    ya=seed(za;n=nt); yb=seed(zb;n=nt)
    g=[ya[1],ya[2],yb[1],yb[2]]
    guess = gmode === :global ? kkg_iterate : (z->seed(z;n=nt))
    try
        sol = vector_bvp_solve(f,za,zb,Ba,Bb,g; N=N,tol=tol,maxiter=maxiter,
                               jacobian=Jf, initial_guess=guess)
        e = maximum(abs(real(sol(x)[1])-u_kkg(x)) for x in range(za,zb,length=25))
        umin = minimum(real(sol(x)[1]) for x in range(za,zb,length=25))
        umax = maximum(real(sol(x)[1]) for x in range(za,zb,length=25))
        println("  [$za,$zb] N=$N gmode=$gmode iters=$(sol.iterations) ",
                "ODEres=", round(sol.residual_inf,sigdigits=3),
                " KKGerr=", round(e,sigdigits=4),
                " u in [", round(umin,digits=3),",",round(umax,digits=3),"]  ",
                e<5e-3 ? "*** TRITRONQUEE ***" : "wrong")
        return (sol,e)
    catch ex
        println("  [$za,$zb] N=$N gmode=$gmode threw: ",
                first(sprint(showerror,ex),120))
        return (nothing,Inf)
    end
end

for (za,zb,N) in [(-20.0,-2.0,80),(-20.0,-2.0,128),(-20.0,-2.0,200),
                  (-40.0,-2.0,160),(-20.0,-10.0,80),(-30.0,-5.0,128),
                  (-20.0,-2.0,64),(-50.0,-4.0,200)]
    run22(za,zb,N; gmode=:global)
end

println("\n### Control: same cases seeded with the OLD asymptotic seed")
for (za,zb,N) in [(-20.0,-2.0,80),(-20.0,-2.0,128)]
    run22(za,zb,N; gmode=:asym)
end

# detailed dump of one global-iterate case
println("\n### Detailed: 2+2 [-20,-2] N=128, global iterate")
let
    sol,e = run22(-20.0,-2.0,128; gmode=:global)
    if sol !== nothing
        for x in -20.0:2.0:-2.0
            un=real(sol(x)[1])
            println("  x=$x  u_num=", round(un,digits=5),
                    "  u_KKG=", round(u_kkg(x),digits=5),
                    "  |Δ|=", round(abs(un-u_kkg(x)),sigdigits=3))
        end
    end
end
