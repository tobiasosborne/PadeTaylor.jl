# external/probes/pade-validity-radius/diag.jl
#
# Diagnostic: the main probe's overlap anchor pinned at exactly RAY_DR=0.02
# for every node.  Rule 2 -- investigate the root cause before concluding.
# This script dumps the raw disagreement profile between concrete node pairs.

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

const KKG_T = 0.0
f  = painleve_hierarchy(:I, 2; t = KKG_T)
Jf = painleve_hierarchy_jacobian(:I, 2; t = KKG_T)
CT = ComplexF64
Ba = zeros(CT,4,4); Bb = zeros(CT,4,4)
Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
sl = pI2_tritronquee_ic(-20.0+0im; t=KKG_T, n_terms=2)
sr = pI2_tritronquee_ic(-2.0+0im;  t=KKG_T, n_terms=2)
g  = CT[sl[1],sl[2],sr[1],sr[2]]
bvp = vector_bvp_solve(f, -20.0+0im, -2.0+0im, Ba, Bb, g;
        N=128, tol=1e-9, maxiter=40, jacobian=Jf,
        initial_guess = z -> pI2_tritronquee_ic(z; t=KKG_T, n_terms=2))
y_seed = ComplexF64.(bvp(ComplexF64(-3.0+0im)))
prob = VectorPadeTaylorProblem(f, y_seed, (ComplexF64(-3.0+0im), ComplexF64(8.0+0im)); order=24)
targets = ComplexF64[r*cis(θ) for r in (2.0,4.0,6.0,8.0) for θ in (-0.5,-0.25,0.0,0.25,0.5)]
walk = vector_path_network_solve(prob, targets; order=24, h=0.1)
N = length(walk.visited_z)
@printf("walk: %d nodes\n", N)

function evalpoly_lh(c, t)
    s = zero(promote_type(eltype(c), typeof(t)))
    @inbounds for k in length(c):-1:1
        s = s*t + c[k]
    end
    s
end
function node_pade_eval(k, z)
    z_v = walk.visited_z[k]; h_v = real(walk.visited_h[k])
    t = (z - z_v)/h_v
    Q = walk.visited_denominator[k]
    q = evalpoly_lh(Q, t)
    ComplexF64[evalpoly_lh(num,t)/q for num in walk.visited_numerators[k]]
end
rel(a,b) = norm(a.-b)/(norm(a)+norm(b))

# (1) Sanity: each node's Pade at its OWN centre should reproduce visited_y.
println("\n[1] node Pade at its own centre vs visited_y (should be ~0):")
for k in (1, 2, 50, 200, 389)
    yc = node_pade_eval(k, walk.visited_z[k])
    @printf("  node %3d: rel-err at centre = %.3e   |y| = %.3e\n",
            k, rel(yc, ComplexF64.(walk.visited_y[k])), norm(walk.visited_y[k]))
end

# (2) Pick a parent->child adjacent pair and dump the disagreement profile.
# Find a node k>=2 whose parent p sits ~h away.
function adjacent_pair()
    for k in 2:N
        p = walk.visited_parent[k]
        p == 0 && continue
        d = abs(walk.visited_z[k]-walk.visited_z[p])
        0.05 < d < 0.2 && return (p,k,d)
    end
    return (1,2,abs(walk.visited_z[2]-walk.visited_z[1]))
end
p,k,d = adjacent_pair()
@printf("\n[2] adjacent pair: parent %d -> child %d, |dz| = %.4f\n", p,k,d)
zp = walk.visited_z[p]; zk = walk.visited_z[k]
dir = (zk-zp)/abs(zk-zp)
@printf("  %-8s %-12s %-12s %-12s\n","r","|P_p|","|P_k|","rel-disagree")
for r in 0.0:0.02:0.30
    z = zp + r*dir
    yp = node_pade_eval(p, z)
    yk = node_pade_eval(k, z)
    @printf("  %-8.3f %-12.4e %-12.4e %-12.4e\n", r, norm(yp), norm(yk), rel(yp,yk))
end

# (3) The midpoint: both nodes evaluated at the geometric midpoint.
zm = (zp+zk)/2
@printf("\n[3] at midpoint z=(%.3f,%.3f): rel-disagree = %.4e\n",
        real(zm),imag(zm), rel(node_pade_eval(p,zm), node_pade_eval(k,zm)))

# (4) Oracle: order-64 Taylor jet at node p vs its own order-24 Pade.
println("\n[4] node p: order-24 Pade vs order-64 Taylor partial sum")
z_v = ComplexF64(walk.visited_z[p]); h_v = real(walk.visited_h[p])
jets = vector_taylor_coefficients(f, z_v, ComplexF64.(walk.visited_y[p]), 64)
rescaled = [ComplexF64[jets[i][q+1]*h_v^q for q in 0:length(jets[i])-1]
            for i in 1:length(jets)]
function tpart(ji, t, deg)
    s = zero(ComplexF64)
    @inbounds for q in min(deg,length(ji)-1):-1:0
        s = s*t + ji[q+1]
    end
    s
end
@printf("  %-8s %-14s %-14s %-14s\n","r","Pade vs T64","T64 vs T56","note")
for r in 0.0:0.05:0.6
    z = z_v + r*dir
    t = (z-z_v)/h_v
    pad = node_pade_eval(p, z)
    ref = ComplexF64[tpart(rescaled[i],t,64) for i in 1:4]
    low = ComplexF64[tpart(rescaled[i],t,56) for i in 1:4]
    conv = rel(ref,low)
    @printf("  %-8.3f %-14.4e %-14.4e %s\n", r, rel(pad,ref), conv,
            conv<1e-9 ? "T-converged" : "T-diverging")
end

# (5) How close are the parent and child Q-denominators?  If the canonical
# Pade is genuinely a local object the two should give nearly identical y on
# the overlap -- unless one node sits on a near-degenerate jet.
@printf("\n[5] |y| at the two node centres: node p=%.4e  node k=%.4e\n",
        norm(walk.visited_y[p]), norm(walk.visited_y[k]))
@printf("    deg(Q_p)=%d deg(Q_k)=%d\n",
        length(walk.visited_denominator[p])-1, length(walk.visited_denominator[k])-1)
