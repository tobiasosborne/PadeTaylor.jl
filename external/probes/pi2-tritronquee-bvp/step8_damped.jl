# Step 8 — damped (Armijo) Newton + monitoring the branch sign-flip.
#
# Established: the wrong branch is u = +cbrt(6)|x|^{1/3} (the OTHER cube
# root); the tritronquee is u = -cbrt(6)|x|^{1/3}. Both solve the leading
# balance 40u^3 ~ -240x for x<0. Plain Newton from the negative seed flips
# to the positive branch on its FIRST FEW steps (Step 5d: no convergence
# at maxiter<=5, i.e. it is still mid-flight away from the seed).
#
# Hypothesis: the flip is a Newton OVERSHOOT -- the full Newton step is too
# large and jumps across the u=0 ridge into the +branch basin. A damped
# (Armijo line-search) Newton that only accepts steps reducing ||R|| should
# stay in the -branch basin. This step prototypes damped Newton on the
# square (u,u') 2+2 collocation system and monitors min(u) over the
# iteration -- if u stays negative throughout, damping is the recipe.

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

# square (u,u') 2+2 collocation residual + Jacobian, tau-method.
function build(za,zb,N; nt=2)
    d=4; np1=N+1
    ntv=Float64[cos(j*pi/N) for j in 0:N]
    hs=(za+zb)/2;s=(zb-za)/2
    nz=Float64[hs+s*tt for tt in ntv]
    D1=_D1(ntv,Float64,N);Dop=kron(D1,Matrix{Float64}(I,d,d))
    nodeN=N*d
    ya=seed(za;n=nt); yb=seed(zb;n=nt)
    Ba=zeros(d,d);Bb=zeros(d,d)
    Ba[1,1]=1;Ba[2,2]=1;Bb[3,1]=1;Bb[4,2]=1
    g=[ya[1],ya[2],yb[1],yb[2]]
    (;d,np1,ntv,nz,s,Dop,nodeN,Ba,Bb,g)
end
function resid(P,Y)
    R=P.Dop*Y
    for j in 1:P.np1
        rng=(j-1)*P.d+1:j*P.d
        R[rng].-=P.s.*f(P.nz[j],Y[rng])
    end
    n=P.nodeN
    R[1:P.d].=P.Ba*Y[n+1:n+P.d].+P.Bb*Y[1:P.d].-P.g
    R
end
function jac(P,Y)
    J=copy(P.Dop)
    for j in 1:P.np1
        rng=(j-1)*P.d+1:j*P.d
        J[rng,rng].-=P.s.*Jf(P.nz[j],Y[rng])
    end
    Jbc=zeros(P.d,P.np1*P.d)
    Jbc[:,P.nodeN+1:P.nodeN+P.d].=P.Ba
    Jbc[:,1:P.d].=P.Bb
    J[1:P.d,:].=Jbc
    J
end

# damped Newton with Armijo backtracking on ||R||_2
function damped_newton(P, Y0; tol=1e-9, maxiter=300, verbose=false)
    Y=copy(Y0)
    for k in 1:maxiter
        R=resid(P,Y); nR=norm(R)
        J=jac(P,Y)
        dY=J\R
        # backtrack
        lam=1.0; accepted=false
        local Ynew
        for bt in 0:25
            Ynew=Y .- lam.*dY
            if norm(resid(P,Ynew)) < (1 - 1e-4*lam)*nR
                accepted=true; break
            end
            lam/=2
        end
        if !accepted
            Ynew=Y .- lam.*dY  # take the smallest step anyway
        end
        umin=minimum(real(Ynew[i]) for i in 1:P.d:length(Ynew))
        umax=maximum(real(Ynew[i]) for i in 1:P.d:length(Ynew))
        verbose && k<=12 && println("    it$k lam=",round(lam,sigdigits=2),
            " ||R||=",round(norm(resid(P,Ynew)),sigdigits=3),
            " u in [",round(umin,digits=2),",",round(umax,digits=2),"]")
        step=maximum(abs,Y.-Ynew)
        Y=Ynew
        step<=tol && return (Y,k,true)
    end
    (Y,maxiter,false)
end

function baryval(nt,vals,ts)
    np1=length(nt)
    for j in 1:np1; abs(ts-nt[j])<1e-14 && return vals[j]; end
    num=0.0;den=0.0
    for j in 1:np1
        sgn=iseven(j-1) ? 1.0 : -1.0
        half=(j==1||j==np1) ? 0.5 : 1.0
        wj=sgn*half;dl=ts-nt[j]
        num+=wj*vals[j]/dl;den+=wj/dl
    end
    num/den
end

println("="^70)
println("STEP 8 — damped (Armijo) Newton, 2+2 (u,u') split")
println("="^70)

for (za,zb,N) in [(-20.0,-2.0,80),(-20.0,-2.0,128),(-20.0,-10.0,80),
                  (-20.0,-19.0,60),(-40.0,-2.0,160)]
    P=build(za,zb,N)
    Y0=reduce(vcat,(seed(z;n=2) for z in P.nz))
    Y,iters,conv=damped_newton(P,Y0; tol=1e-9, maxiter=400,
                               verbose=(za==-20.0&&zb==-2.0&&N==80))
    hs=(za+zb)/2;sc=(zb-za)/2
    uvals=[Y[(j-1)*4+1] for j in 1:N+1]
    errs=[abs(baryval(P.ntv,uvals,(x-hs)/sc)-u_kkg(x)) for x in range(za,zb,length=25)]
    e=maximum(errs)
    println("  [$za,$zb] N=$N conv=$conv iters=$iters KKGerr=",
            round(e,sigdigits=4),"  ", e<5e-3 ? "*** TRITRONQUEE ***" : "wrong")
end

# the verbose trace above shows whether u flips sign on step 1.
println("\n(see the verbose trace for [-20,-2] N=80 -- does u flip sign on step 1?)")
