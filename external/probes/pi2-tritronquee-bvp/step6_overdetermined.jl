# Step 6 — the missing ingredient: an over-determined (4+4) BC.
#
# Step 5(c) is decisive: even seeding the Newton guess with the KKG
# asymptotic at EVERY node, the (u,u') 2+2 BC slides to the wrong branch
# at N up to 256. The tritronquee is not a stable fixed point of the
# collocation-Newton map under ANY square 2+2 BC -- every 2+2 split (Step 4)
# converges to a wrong, pole-laden P_I^(2) solution that ALSO satisfies the
# 4 pinned components. The 3+1 / 4+0 splits give a singular Jacobian (the
# separatrix growing mode leaves the linearised collocation operator
# rank-deficient when too much asymptotic data is pinned at one end).
#
# KKG 2015 pinned the truncated asymptotic series at BOTH ends -> 4+4 = 8
# conditions for a d=4 system. That is an OVER-DETERMINED two-point BVP;
# the collocation system is then solved in LEAST SQUARES. The current
# vector_bvp_solve is square-only (exactly d BC rows replacing node-0's d
# collocation rows). This step PROTOTYPES the over-determined collocation
# solve by hand (probe-local, no src/ change) to confirm it selects the
# tritronquee -- the recipe the production fix should implement.
#
# Over-determined assembly: keep ALL (N+1)*d collocation rows, APPEND the
# 8 BC rows (B_a y(z_a) full + B_b y(z_b) full), solve the resulting
# ((N+1)*d + 8) x ((N+1)*d) least-squares system each Newton step.

using PadeTaylor.PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using LinearAlgebra
using PadeTaylor.VectorBVP: VectorBVP   # for _chebyshev_D1 access via getfield

t   = 0.0
f   = painleve_hierarchy(:I, 2; t = t)
Jf  = painleve_hierarchy_jacobian(:I, 2; t = t)
u_kkg(x)   = -cbrt(6.0) * abs(x)^(1//3)
seed(z; n) = pI2_tritronquee_ic(z; t = t, n_terms = n)

# reach the module-private D1 builder (probe only; production fix would
# expose it or share via the existing kron path)
_D1 = getfield(VectorBVP, :_chebyshev_D1)

# Hand-rolled over-determined Chebyshev-collocation Newton for d=4.
# BC: B_a*y(z_a) + B_b*y(z_b) = g, with B_a, B_b being w x 4 (w BC rows).
function overdet_bvp(f, Jf, za, zb, Ba, Bb, g; N=80, tol=1e-9, maxiter=120,
                     guess = nothing)
    d  = 4
    np1 = N+1
    nt = Float64[cos(j*pi/N) for j in 0:N]
    hs = (za+zb)/2; s = (zb-za)/2
    nz = Float64[hs + s*t for t in nt]
    D1 = _D1(nt, Float64, N)
    Dop = kron(D1, Matrix{Float64}(I,d,d))
    w  = length(g)
    # endpoint -> node: node 0 = t=+1 = z_b ; node N = t=-1 = z_a
    nodeN = N*d
    Jbc = zeros(w, np1*d)
    Jbc[:, nodeN+1:nodeN+d] .= Ba
    Jbc[:, 1:d]             .= Bb
    Y = guess === nothing ? zeros(np1*d) :
        reduce(vcat, (guess(z) for z in nz))
    local stepn = Inf
    for k in 1:maxiter
        # residual: collocation rows
        R = Dop*Y
        for j in 1:np1
            rng = (j-1)*d+1:j*d
            R[rng] .-= s .* f(nz[j], Y[rng])
        end
        # BC residual rows
        Rbc = Ba*Y[nodeN+1:nodeN+d] .+ Bb*Y[1:d] .- g
        Rfull = vcat(R, Rbc)
        # Jacobian
        J = copy(Dop)
        for j in 1:np1
            rng = (j-1)*d+1:j*d
            J[rng,rng] .-= s .* Jf(nz[j], Y[rng])
        end
        Jfull = vcat(J, Jbc)
        # least-squares Newton step
        dY = Jfull \ Rfull
        stepn = maximum(abs, dY)
        Y .-= dY
        stepn <= tol && return (Y, nz, nt, k, stepn, true)
    end
    return (Y, nz, nt, maxiter, stepn, false)
end

# barycentric eval (chebyshev-2) of component i
function baryval(nt, vals, tstar)
    np1 = length(nt)
    for j in 1:np1
        abs(tstar-nt[j]) < 1e-14 && return vals[j]
    end
    num = 0.0; den = 0.0
    for j in 1:np1
        sgn = iseven(j-1) ? 1.0 : -1.0
        half = (j==1||j==np1) ? 0.5 : 1.0
        wj = sgn*half
        delta = tstar - nt[j]
        num += wj*vals[j]/delta
        den += wj/delta
    end
    num/den
end

println("="^70)
println("STEP 6 — over-determined 4+4 BC: prototype the KKG recipe")
println("="^70)

function run_od(za, zb, N, nt_seed)
    ya = seed(za; n=nt_seed); yb = seed(zb; n=nt_seed)
    # 4+4: pin the FULL asymptotic state at BOTH ends
    Ba = vcat(Matrix{Float64}(I,4,4), zeros(4,4))    # 8x4
    Bb = vcat(zeros(4,4), Matrix{Float64}(I,4,4))    # 8x4
    g  = vcat(ya, yb)
    Y, nz, ntv, iters, stepn, conv = overdet_bvp(
        f, Jf, za, zb, Ba, Bb, g; N=N, tol=1e-8, maxiter=200,
        guess = z->seed(z; n=nt_seed))
    hs=(za+zb)/2; sc=(zb-za)/2
    # KKG error
    xs = range(za, zb; length=25)
    errs = Float64[]
    for x in xs
        tstar = (x-hs)/sc
        uvals = [Y[(j-1)*4+1] for j in 1:N+1]
        un = baryval(ntv, uvals, tstar)
        push!(errs, abs(real(un)-u_kkg(x)))
    end
    e = maximum(errs)
    println("  [4+4 OD] [$za,$zb] N=$N nt=$nt_seed conv=$conv iters=$iters ",
            "stepn=", round(stepn,sigdigits=3), " KKGerr=", round(e,sigdigits=4),
            "  ", e<5e-3 ? "*** TRITRONQUEE ***" : "wrong")
    return (Y, nz, ntv, e, conv)
end

for (za,zb,N,nts) in [(-20.0,-2.0,80,2),(-20.0,-2.0,128,2),
                      (-20.0,-2.0,200,2),(-40.0,-2.0,200,2),
                      (-20.0,-10.0,80,2)]
    run_od(za,zb,N,nts)
end

# detailed dump of the [-20,-2] N=200 case
println("\n### Detailed: 4+4 OD [-20,-2] N=200")
let
    Y, nz, ntv, e, conv = run_od(-20.0,-2.0,200,2)
    hs=-11.0; sc=9.0
    for x in -20.0:2.0:-2.0
        tstar=(x-hs)/sc
        uvals=[Y[(j-1)*4+1] for j in 1:201]
        un=baryval(ntv,uvals,tstar)
        println("  x=$x  u_num=", round(real(un),digits=5),
                "  u_KKG=", round(u_kkg(x),digits=5),
                "  |Δ|=", round(abs(real(un)-u_kkg(x)),sigdigits=3))
    end
end
