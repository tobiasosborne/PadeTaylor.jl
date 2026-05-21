# Step 10 — Approach A (asymptotic collar) and Approach B (deviation var).
#
# Step 9 confirmed (again): a 2+2 endpoint-value BVP genuinely ADMITS the
# +branch as a solution with the SAME 4 endpoint values; no initial guess
# can fix this because the +branch is a true solution of that BVP. The BC
# must constrain the INTERIOR.
#
# Approach A -- asymptotic collar. KKG's tau-method replaces the equations
# at j=0,1,N_c-1,N_c by boundary conditions. The deeper idea: pin the
# solution to the truncated asymptotic series u = Y + sum c_n Y^{-n} on a
# COLLAR of K nodes at each end. With K>=2 nodes per end the pinned band
# fixes u AND its slope on a finite stretch -- the +branch (u ~ +cbrt6|x|)
# cannot match a band that says u ~ -cbrt6|x|. We over-determine: keep all
# (N+1)*4 collocation rows, append 4*K_a + 4*K_b collar rows pinning the
# full companion state at the collar nodes, solve least-squares per Newton
# step with collar and collocation rows COMPARABLY scaled (collar weight
# chosen so the band is enforced but Newton still converges).
#
# Approach B -- deviation variable. u = u_a + v, u_a = -cbrt6*(-x)^{1/3}
# (the tritronquee leading order, real & negative for x<0). The companion
# RHS for v is g(x,Y_v) = f(x, Y_a+Y_v) - Y_a' where Y_a is the companion
# vector of u_a. Newton from v=0 (=> u=u_a everywhere) starts AT the
# tritronquee leading order across the whole interval. Test plain 2+2 in v.

using PadeTaylor.PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP: VectorBVP
using LinearAlgebra

t   = 0.0
f   = painleve_hierarchy(:I, 2; t = t)
Jf  = painleve_hierarchy_jacobian(:I, 2; t = t)
u_kkg(x)   = -cbrt(6.0) * abs(x)^(1//3)
seed(z; n) = pI2_tritronquee_ic(z; t = t, n_terms = n)
_D1 = getfield(VectorBVP, :_chebyshev_D1)
const C6 = cbrt(6.0)

# ---- tritronquee leading-order u_a = -cbrt6*(-x)^{1/3}, x<0 -------------
# r=-x>0; u_a = -C6*r^{1/3}; d/dx = -d/dr.
function ua_state(x)
    r = -x
    u   = -C6 * r^(1//3)
    up  =  (1//3)  * C6 * r^(-2//3)        # du/dx = -du/dr
    upp =  (2//9)  * C6 * r^(-5//3)
    uppp=  (10//27)* C6 * r^(-8//3)
    upppp= (80//81)* C6 * r^(-11//3)
    return (u,up,upp,uppp,upppp)
end
ua_vec(x) = collect(ua_state(x)[1:4])

function baryval(nt,vals,ts)
    np1=length(nt)
    for j in 1:np1; abs(ts-nt[j])<1e-13 && return vals[j]; end
    num=zero(eltype(vals));den=0.0
    for j in 1:np1
        sgn=iseven(j-1) ? 1.0 : -1.0
        half=(j==1||j==np1) ? 0.5 : 1.0
        wj=sgn*half;dl=ts-nt[j]
        num+=wj*vals[j]/dl;den+=wj/dl
    end
    num/den
end

println("="^70)
println("STEP 10 — Approach A (collar) + Approach B (deviation variable)")
println("="^70)

# =====================================================================
# Approach A — asymptotic collar (over-determined weighted least squares)
# =====================================================================
# collar: at the first Ka nodes (near z_b end, t=+1) and last Kb nodes
# (near z_a end, t=-1) pin the full 4-vector to the asymptotic seed.
function collar_bvp(za,zb,N; Ka=4, Kb=4, w=1.0, nt=2, tol=1e-9, maxiter=300)
    d=4; np1=N+1
    ntv=Float64[cos(j*pi/N) for j in 0:N]
    hs=(za+zb)/2; s=(zb-za)/2
    nz=Float64[hs+s*tt for tt in ntv]
    D1=_D1(ntv,Float64,N); Dop=kron(D1,Matrix{Float64}(I,d,d))
    # collar node indices (1-based): nodes 1..Ka (t near +1 => z near zb)
    # and nodes np1-Kb+1..np1 (t near -1 => z near za)
    collar = vcat(1:Ka, (np1-Kb+1):np1)
    # collar target = asymptotic seed at those nodes
    Y = reduce(vcat,(seed(z;n=nt) for z in nz))
    local stepn=Inf
    nC = length(collar)
    for k in 1:maxiter
        R=Dop*Y
        for j in 1:np1
            rng=(j-1)*d+1:j*d
            R[rng].-=s.*f(nz[j],Y[rng])
        end
        Rcol=Float64[]
        for c in collar
            rng=(c-1)*d+1:c*d
            append!(Rcol, w.*(Y[rng].-seed(nz[c];n=nt)))
        end
        Rfull=vcat(R,Rcol)
        J=copy(Dop)
        for j in 1:np1
            rng=(j-1)*d+1:j*d
            J[rng,rng].-=s.*Jf(nz[j],Y[rng])
        end
        Jcol=zeros(d*nC, np1*d)
        for (ci,c) in enumerate(collar)
            for a in 1:d
                Jcol[(ci-1)*d+a, (c-1)*d+a]=w
            end
        end
        Jfull=vcat(J,Jcol)
        dY=Jfull\Rfull
        stepn=maximum(abs,dY)
        Y.-=dY
        stepn<=tol && return (Y,nz,ntv,k,stepn,true)
    end
    (Y,nz,ntv,maxiter,stepn,false)
end

println("\n### Approach A: collar-matched over-determined BVP")
for (za,zb,N,Ka,Kb,w) in [
        (-20.0,-2.0,128,4,4,1.0),(-20.0,-2.0,128,4,4,10.0),
        (-20.0,-2.0,128,8,8,1.0),(-20.0,-2.0,128,8,8,30.0),
        (-20.0,-2.0,200,12,12,1.0),(-20.0,-2.0,200,16,16,1.0),
        (-40.0,-2.0,200,12,12,1.0),(-20.0,-10.0,128,8,8,1.0)]
    Y,nz,ntv,iters,stepn,conv=collar_bvp(za,zb,N;Ka=Ka,Kb=Kb,w=w,
                                         tol=1e-8,maxiter=400)
    hs=(za+zb)/2; sc=(zb-za)/2
    uvals=[Y[(j-1)*4+1] for j in 1:N+1]
    errs=[abs(baryval(ntv,uvals,(x-hs)/sc)-u_kkg(x)) for x in range(za,zb,length=25)]
    e=maximum(errs)
    umin=minimum(real(baryval(ntv,uvals,(x-hs)/sc)) for x in range(za,zb,length=25))
    umax=maximum(real(baryval(ntv,uvals,(x-hs)/sc)) for x in range(za,zb,length=25))
    println("  [$za,$zb] N=$N Ka=$Ka Kb=$Kb w=$w conv=$conv iters=$iters ",
            "KKGerr=", round(e,sigdigits=4),
            " u in [",round(umin,digits=2),",",round(umax,digits=2),"]  ",
            e<5e-3 ? "*** TRITRONQUEE ***" : "wrong")
end

# =====================================================================
# Approach B — deviation variable u = u_a + v
# =====================================================================
# RHS for v companion vector Yv:  Yv' = f(x, Ya+Yv) - Ya'
# where Ya = (u_a,u_a',u_a'',u_a''') and Ya' = (u_a',u_a'',u_a''',u_a'''').
function fv(x, Yv)
    ua = ua_state(x)
    Ya  = [ua[1],ua[2],ua[3],ua[4]]
    Yap = [ua[2],ua[3],ua[4],ua[5]]
    return f(x, Ya .+ Yv) .- Yap
end
# Jacobian of fv wrt Yv = Jf evaluated at Ya+Yv (Ya' term is const in Yv)
function Jfv(x, Yv)
    ua = ua_state(x)
    Ya = [ua[1],ua[2],ua[3],ua[4]]
    return Jf(x, Ya .+ Yv)
end

println("\n### Approach B: deviation variable, plain 2+2 (v,v') BVP, guess v=0")
function dev_bvp(za,zb,N; nt=2, tol=1e-9, maxiter=300, w_collar=0.0, Kc=0)
    d=4; np1=N+1
    ntv=Float64[cos(j*pi/N) for j in 0:N]
    hs=(za+zb)/2; s=(zb-za)/2
    nz=Float64[hs+s*tt for tt in ntv]
    D1=_D1(ntv,Float64,N); Dop=kron(D1,Matrix{Float64}(I,d,d))
    nodeN=N*d
    # deviation BC data: v = seed - u_a  (tiny at the asymptotic ends)
    va = seed(za;n=nt) .- ua_vec(za)
    vb = seed(zb;n=nt) .- ua_vec(zb)
    Ba=zeros(d,d);Bb=zeros(d,d)
    Ba[1,1]=1;Ba[2,2]=1;Bb[3,1]=1;Bb[4,2]=1
    g=[va[1],va[2],vb[1],vb[2]]
    collar = Kc>0 ? vcat(1:Kc,(np1-Kc+1):np1) : Int[]
    Y=zeros(np1*d)            # guess v=0  => u = u_a everywhere
    local stepn=Inf
    for k in 1:maxiter
        R=Dop*Y
        for j in 1:np1
            rng=(j-1)*d+1:j*d
            R[rng].-=s.*fv(nz[j],Y[rng])
        end
        n=nodeN
        R[1:d].=Ba*Y[n+1:n+d].+Bb*Y[1:d].-g
        Rcol=Float64[]
        for c in collar
            rng=(c-1)*d+1:c*d
            vtarget = seed(nz[c];n=nt).-ua_vec(nz[c])
            append!(Rcol, w_collar.*(Y[rng].-vtarget))
        end
        Rfull=vcat(R,Rcol)
        J=copy(Dop)
        for j in 1:np1
            rng=(j-1)*d+1:j*d
            J[rng,rng].-=s.*Jfv(nz[j],Y[rng])
        end
        Jbc=zeros(d,np1*d)
        Jbc[:,nodeN+1:nodeN+d].=Ba; Jbc[:,1:d].=Bb
        J[1:d,:].=Jbc
        Jcol=zeros(d*length(collar),np1*d)
        for (ci,c) in enumerate(collar), a in 1:d
            Jcol[(ci-1)*d+a,(c-1)*d+a]=w_collar
        end
        Jfull=vcat(J,Jcol)
        dY=Jfull\Rfull
        stepn=maximum(abs,dY)
        Y.-=dY
        stepn<=tol && return (Y,nz,ntv,k,stepn,true)
    end
    (Y,nz,ntv,maxiter,stepn,false)
end

for (za,zb,N) in [(-20.0,-2.0,80),(-20.0,-2.0,128),(-20.0,-2.0,200),
                  (-40.0,-2.0,160),(-20.0,-10.0,80),(-30.0,-5.0,128)]
    Y,nz,ntv,iters,stepn,conv=dev_bvp(za,zb,N;tol=1e-8,maxiter=400)
    hs=(za+zb)/2; sc=(zb-za)/2
    # u = u_a + v
    uvals=[Y[(j-1)*4+1]+ua_vec(nz[j])[1] for j in 1:N+1]
    errs=[abs(baryval(ntv,uvals,(x-hs)/sc)-u_kkg(x)) for x in range(za,zb,length=25)]
    e=maximum(errs)
    vmax=maximum(abs(Y[(j-1)*4+1]) for j in 1:N+1)
    println("  [$za,$zb] N=$N conv=$conv iters=$iters stepn=",
            round(stepn,sigdigits=3),
            " max|v|=", round(vmax,sigdigits=4),
            " KKGerr=", round(e,sigdigits=4), "  ",
            e<5e-3 ? "*** TRITRONQUEE ***" : "wrong")
end

println("\n### Approach B + collar: deviation var with collar pinning")
for (za,zb,N,Kc,w) in [(-20.0,-2.0,128,8,1.0),(-20.0,-2.0,200,12,1.0),
                       (-20.0,-2.0,128,4,1.0)]
    Y,nz,ntv,iters,stepn,conv=dev_bvp(za,zb,N;Kc=Kc,w_collar=w,
                                      tol=1e-8,maxiter=400)
    hs=(za+zb)/2; sc=(zb-za)/2
    uvals=[Y[(j-1)*4+1]+ua_vec(nz[j])[1] for j in 1:N+1]
    errs=[abs(baryval(ntv,uvals,(x-hs)/sc)-u_kkg(x)) for x in range(za,zb,length=25)]
    e=maximum(errs)
    println("  [$za,$zb] N=$N Kc=$Kc w=$w conv=$conv iters=$iters KKGerr=",
            round(e,sigdigits=4),"  ", e<5e-3 ? "*** TRITRONQUEE ***" : "wrong")
end
