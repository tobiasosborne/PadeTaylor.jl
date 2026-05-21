# external/probes/pade-validity-radius/diag2.jl
#
# Second diagnostic.  Two questions the main probe leaves open:
#  (A) Is min|t*| (gate ii) dominated by spurious Froissart roots?  If a node's
#      smallest-|t*| Q-root is a numerical artefact, R_ii is meaningless.
#  (B) Does the oracle anchor get truncated early because the order-64 Taylor
#      sum DIVERGES (nearest pole) before the Pade actually fails?  If so R_iv
#      is a pessimistic lower bound and the honest radius is larger.

using Printf
using LinearAlgebra: norm, eigvals
using PadeTaylor
using PadeTaylor.PainleveHierarchy: painleve_hierarchy,
                                    painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP:         vector_bvp_solve
using PadeTaylor.VectorProblems:    VectorPadeTaylorProblem
using PadeTaylor.VectorPathNetwork: vector_path_network_solve
using PadeTaylor.VectorCoefficients: vector_taylor_coefficients

function poly_roots(c)
    cc = collect(ComplexF64, c)
    last = findlast(x -> abs(x) > 1e-14*(1+maximum(abs,cc)), cc)
    last === nothing && return ComplexF64[]
    cc = cc[1:last]; n = length(cc)-1; n < 1 && return ComplexF64[]
    a = cc ./ cc[end]; C = zeros(ComplexF64,n,n)
    for i in 1:n-1; C[i+1,i]=1.0; end
    for i in 1:n; C[i,n]=-a[i]; end
    eigvals(C)
end

f  = painleve_hierarchy(:I, 2; t = 0.0)
Jf = painleve_hierarchy_jacobian(:I, 2; t = 0.0)
CT = ComplexF64
Ba = zeros(CT,4,4); Bb = zeros(CT,4,4)
Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
sl = pI2_tritronquee_ic(-20.0+0im; t=0.0, n_terms=2)
sr = pI2_tritronquee_ic(-2.0+0im;  t=0.0, n_terms=2)
bvp = vector_bvp_solve(f, -20.0+0im, -2.0+0im, Ba, Bb,
        CT[sl[1],sl[2],sr[1],sr[2]];
        N=128, tol=1e-9, maxiter=40, jacobian=Jf,
        initial_guess = z -> pI2_tritronquee_ic(z; t=0.0, n_terms=2))
y_seed = ComplexF64.(bvp(ComplexF64(-3.0+0im)))
prob = VectorPadeTaylorProblem(f, y_seed, (ComplexF64(-3.0+0im),ComplexF64(8.0+0im)); order=24)
targets = ComplexF64[r*cis(θ) for r in (2.0,4.0,6.0,8.0) for θ in (-0.5,-0.25,0.0,0.25,0.5)]
walk = vector_path_network_solve(prob, targets; order=24, h=0.1)
N = length(walk.visited_z)
@printf("walk: %d nodes\n", N)

evalpoly_lh(c,t) = (s=zero(ComplexF64); for k in length(c):-1:1; s=s*t+c[k]; end; s)
function pade(k,z)
    z_v=walk.visited_z[k]; h_v=real(walk.visited_h[k]); t=(z-z_v)/h_v
    Q=walk.visited_denominator[k]; q=evalpoly_lh(Q,t)
    ComplexF64[evalpoly_lh(num,t)/q for num in walk.visited_numerators[k]]
end
rel(a,b)=norm(a.-b)/(norm(a)+norm(b))

# (A) For 6 representative nodes: list ALL Q-roots sorted by |t*|, with the
# residue magnitude as a Froissart indicator (a doublet root + tiny numerator
# residue = spurious).  We approximate "spuriousness" by how isolated the root
# is: a genuine pole root is not paired with a near-coincident numerator zero.
println("\n[A] Q-root structure (rescaled t-plane), sorted by |t*|:")
for k in (2, 50, 100, 200, 300, 389)
    ts = sort(poly_roots(walk.visited_denominator[k]); by=abs)
    h_v = real(walk.visited_h[k])
    @printf("  node %3d  deg(Q)=%2d  |y|=%.2e\n",
            k, length(walk.visited_denominator[k])-1, norm(walk.visited_y[k]))
    for (j,t) in enumerate(ts)
        j > 6 && break
        @printf("    |t*|=%7.3f  z-dist=%7.3f  arg=%+.2f\n",
                abs(t), abs(t)*h_v, angle(t))
    end
end

# (B) For the same nodes: along the +real direction, tabulate Pade-vs-oracle
# error AND the oracle's own convergence, to see which fails first.
println("\n[B] Pade vs order-64 Taylor oracle along +real ray:")
for k in (2, 200, 389)
    z_v=ComplexF64(walk.visited_z[k]); h_v=real(walk.visited_h[k])
    jets = vector_taylor_coefficients(f, z_v, ComplexF64.(walk.visited_y[k]), 64)
    resc = [ComplexF64[jets[i][q+1]*h_v^q for q in 0:length(jets[i])-1] for i in 1:4]
    tpart(ji,t,deg)=(s=zero(ComplexF64); for q in min(deg,length(ji)-1):-1:0; s=s*t+ji[q+1]; end; s)
    @printf("  node %d  (z=%.2f%+.2fi):\n", k, real(z_v), imag(z_v))
    @printf("    %-7s %-13s %-13s %-13s\n","r","Pade-vs-T64","T64-vs-T56","T64-vs-T48")
    for r in 0.0:0.1:1.2
        z=z_v+r*1.0; t=(z-z_v)/h_v
        p   = pade(k,z)
        t64 = ComplexF64[tpart(resc[i],t,64) for i in 1:4]
        t56 = ComplexF64[tpart(resc[i],t,56) for i in 1:4]
        t48 = ComplexF64[tpart(resc[i],t,48) for i in 1:4]
        @printf("    %-7.2f %-13.3e %-13.3e %-13.3e\n",
                r, rel(p,t64), rel(t64,t56), rel(t64,t48))
    end
end

# (C) Histogram: at tol 1e-8, count nodes where the oracle ray terminated
# because the Taylor sum DIVERGED (last_conv) vs an actual Pade crossing.
# We re-run the per-node oracle but tag the termination cause.
println("\n[C] R_iv termination cause @ tol 1e-8 (24-dir fan, min over dirs):")
function gate_iv_tagged(k, tol)
    z_v=ComplexF64(walk.visited_z[k]); h_v=real(walk.visited_h[k])
    jets = try vector_taylor_coefficients(f, z_v, ComplexF64.(walk.visited_y[k]), 64)
           catch; return (NaN,:fail) end
    resc = [ComplexF64[jets[i][q+1]*h_v^q for q in 0:length(jets[i])-1] for i in 1:4]
    tpart(ji,t,deg)=(s=zero(ComplexF64); for q in min(deg,length(ji)-1):-1:0; s=s*t+ji[q+1]; end; s)
    Rmin=Inf; cause=:none
    for d in 0:23
        dir=cis(2π*d/24); r=0.02; lastc=0.0; cross=Inf; thiscause=:diverge
        while r<=3.0
            z=z_v+r*dir; t=(z-z_v)/h_v
            t64=ComplexF64[tpart(resc[i],t,64) for i in 1:4]
            t56=ComplexF64[tpart(resc[i],t,56) for i in 1:4]
            if rel(t64,t56)>=1e-12; break; end
            lastc=r
            p=try pade(k,z) catch; cross=r; thiscause=:padeerr; break end
            if any(!isfinite,p); cross=r; thiscause=:padeerr; break; end
            if rel(p,t64)>tol; cross=r; thiscause=:cross; break; end
            r+=0.02
        end
        rr = isfinite(cross) ? cross : lastc
        if rr<Rmin; Rmin=rr; cause=(isfinite(cross) ? thiscause : :diverge); end
    end
    (Rmin,cause)
end
causes = Symbol[]
for k in 1:N
    _,c = gate_iv_tagged(k,1e-8)
    push!(causes, c)
end
ncross = count(==(:cross), causes)
ndiv   = count(==(:diverge), causes)
nerr   = count(c -> c in (:padeerr,:fail), causes)
@printf("  binding direction terminated by: Pade-crossing=%d  Taylor-diverged=%d  Pade-error=%d\n",
        ncross, ndiv, nerr)
println("  (if Taylor-diverged dominates, R_iv is a LOWER bound and true R_honest is larger)")
