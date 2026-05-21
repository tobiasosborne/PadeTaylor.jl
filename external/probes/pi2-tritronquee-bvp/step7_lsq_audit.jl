# Step 7 — audit the over-determined 4+4 solve: did it satisfy its BC?
#
# Step 6's 4+4 over-determined solve converged to u = +cbrt(6)|x|^{1/3}
# (POSITIVE) while g pinned u = -cbrt(6)|x|^{1/3} (NEGATIVE). KKG ground
# truth (tritronquee_coeff.tex:171, BMP): U_0(x,0) ~ -|6x|^{1/3} on the
# negative real axis -- the tritronquee IS negative there, ~ -4.93 at -20.
#
# So the 4+4 solve VIOLATED its BC. Diagnosis: in a plain dense
# least-squares  Jfull \ Rfull  with (N+1)*4 collocation rows and only 8
# BC rows, the BC rows are NUMERICALLY NEGLIGIBLE -- the residual minimum
# is dominated by the collocation block, the 8 BC equations contribute
# ~8/((N+1)*4) of the squared residual and are effectively ignored.
#
# This step:
#  (a) re-run the 4+4 solve, PRINT the BC residual explicitly;
#  (b) re-run with the 8 BC rows WEIGHTED UP (multiply BC block by a large
#      weight w) so least-squares actually enforces them -- the standard
#      fix for an over-determined collocation BVP;
#  (c) if weighting works -> that is the recipe.

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

# weighted over-determined collocation Newton.
# bcweight multiplies the 8 BC rows in BOTH the residual and Jacobian.
function wod_bvp(za, zb, g; N=80, tol=1e-9, maxiter=200, bcweight=1.0, nt=2)
    d=4; np1=N+1
    ntv = Float64[cos(j*pi/N) for j in 0:N]
    hs=(za+zb)/2; s=(zb-za)/2
    nz = Float64[hs+s*tt for tt in ntv]
    D1=_D1(ntv,Float64,N); Dop=kron(D1,Matrix{Float64}(I,d,d))
    nodeN=N*d
    Ba=Matrix{Float64}(I,d,d); Bb=Matrix{Float64}(I,d,d)
    ya=g[1:4]; yb=g[5:8]
    Y = reduce(vcat,(seed(z;n=nt) for z in nz))
    local stepn=Inf
    for k in 1:maxiter
        R=Dop*Y
        for j in 1:np1
            rng=(j-1)*d+1:j*d
            R[rng] .-= s.*f(nz[j],Y[rng])
        end
        Rbc_a = (Y[nodeN+1:nodeN+d] .- ya) .* bcweight
        Rbc_b = (Y[1:d] .- yb) .* bcweight
        Rfull = vcat(R, Rbc_a, Rbc_b)
        J=copy(Dop)
        for j in 1:np1
            rng=(j-1)*d+1:j*d
            J[rng,rng] .-= s.*Jf(nz[j],Y[rng])
        end
        Jbc=zeros(2d,np1*d)
        Jbc[1:d, nodeN+1:nodeN+d] .= bcweight*Matrix{Float64}(I,d,d)
        Jbc[d+1:2d, 1:d]          .= bcweight*Matrix{Float64}(I,d,d)
        Jfull=vcat(J,Jbc)
        dY = Jfull \ Rfull
        stepn=maximum(abs,dY)
        Y .-= dY
        stepn<=tol && return (Y,nz,ntv,k,stepn,true)
    end
    return (Y,nz,ntv,maxiter,stepn,false)
end

function baryval(nt,vals,ts)
    np1=length(nt)
    for j in 1:np1; abs(ts-nt[j])<1e-14 && return vals[j]; end
    num=0.0;den=0.0
    for j in 1:np1
        sgn=iseven(j-1) ? 1.0 : -1.0
        half=(j==1||j==np1) ? 0.5 : 1.0
        wj=sgn*half; dl=ts-nt[j]
        num+=wj*vals[j]/dl; den+=wj/dl
    end
    num/den
end

println("="^70)
println("STEP 7 — weighted over-determined 4+4 BC audit")
println("="^70)

za,zb=-20.0,-2.0
g=vcat(seed(za;n=2), seed(zb;n=2))
hs=(za+zb)/2; sc=(zb-za)/2

for bw in (1.0, 1e2, 1e4, 1e6, 1e8)
    Y,nz,ntv,iters,stepn,conv = wod_bvp(za,zb,g; N=120,tol=1e-8,
                                        maxiter=300,bcweight=bw)
    # BC residual (unweighted)
    YN = Y[end-3:end]; Y0 = Y[1:4]
    bcres = max(maximum(abs.(YN .- g[1:4])), maximum(abs.(Y0 .- g[5:8])))
    # KKG err
    uvals=[Y[(j-1)*4+1] for j in 1:121]
    errs=[abs(baryval(ntv,uvals,(x-hs)/sc)-u_kkg(x)) for x in range(za,zb,length=25)]
    e=maximum(errs)
    println("  bcweight=$bw  conv=$conv iters=$iters stepn=",
            round(stepn,sigdigits=3), "\n     BCres(unweighted)=",
            round(bcres,sigdigits=4), "  KKGerr=", round(e,sigdigits=4),
            "  ", e<5e-3 ? "*** TRITRONQUEE ***" : "wrong")
    println("     u(z_a) num=", round(Y[end-3],sigdigits=6),
            "  g pins u(z_a)=", round(g[1],sigdigits=6))
end

# also: weighted 4+0 (full state at outer end only) -- now well-posed?
println("\n### Weighted 4+0 (full state at z_a only)")
function wod_4plus0(za,zb,ya; N=120,tol=1e-8,maxiter=300,bcweight=1e4,nt=2)
    d=4; np1=N+1
    ntv=Float64[cos(j*pi/N) for j in 0:N]
    hs=(za+zb)/2;s=(zb-za)/2
    nz=Float64[hs+s*tt for tt in ntv]
    D1=_D1(ntv,Float64,N);Dop=kron(D1,Matrix{Float64}(I,d,d))
    nodeN=N*d
    Y=reduce(vcat,(seed(z;n=nt) for z in nz))
    local stepn=Inf
    for k in 1:maxiter
        R=Dop*Y
        for j in 1:np1
            rng=(j-1)*d+1:j*d
            R[rng].-=s.*f(nz[j],Y[rng])
        end
        Rbc=(Y[nodeN+1:nodeN+d].-ya).*bcweight
        Rfull=vcat(R,Rbc)
        J=copy(Dop)
        for j in 1:np1
            rng=(j-1)*d+1:j*d
            J[rng,rng].-=s.*Jf(nz[j],Y[rng])
        end
        Jbc=zeros(d,np1*d)
        Jbc[1:d,nodeN+1:nodeN+d].=bcweight*Matrix{Float64}(I,d,d)
        Jfull=vcat(J,Jbc)
        dY=Jfull\Rfull
        stepn=maximum(abs,dY)
        Y.-=dY
        stepn<=tol && return (Y,nz,ntv,k,stepn,true)
    end
    return (Y,nz,ntv,maxiter,stepn,false)
end
for bw in (1e2,1e4,1e6)
    ya=seed(za;n=2)
    Y,nz,ntv,iters,stepn,conv=wod_4plus0(za,zb,ya;N=120,bcweight=bw)
    uvals=[Y[(j-1)*4+1] for j in 1:121]
    errs=[abs(baryval(ntv,uvals,(x-hs)/sc)-u_kkg(x)) for x in range(za,zb,length=25)]
    e=maximum(errs)
    bcres=maximum(abs.(Y[end-3:end].-ya))
    println("  bcweight=$bw conv=$conv iters=$iters BCres=",round(bcres,sigdigits=4),
            " KKGerr=",round(e,sigdigits=4),"  ",
            e<5e-3 ? "*** TRITRONQUEE ***" : "wrong")
end
