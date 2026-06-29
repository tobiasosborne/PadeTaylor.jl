# P2 at the ACTUAL figure scale [-50,50]^2 (panel e). Field-value jump map +
# two-seed path-dependence + pole dump (for the orientation grain-boundary), so
# we can see WHERE the field jump is large and whether it coincides with the
# pole grain boundary the user sees.
using PadeTaylor, Printf
const PS = PadeTaylor.PadeStepper
pI(z,u,up)=6u^2+z
const HALF=50.0; N=101
xs=range(-HALF,HALF;length=N); ys=range(-HALF,HALF;length=N)
grid=ComplexF64[complex(x,y) for y in ys for x in xs]
prob=PadeTaylorProblem(pI,(-0.1875,0.3049),(0.0,HALF);order=30)
SD="/tmp/claude-1000/-home-tobias-Projects-PadeTaylor-jl/e817ae45-f9a1-4b4a-818b-2d5e3e929f4b/scratchpad"
ev(P,t)= try PS._evaluate_pade(P,t) catch; ComplexF64(NaN,NaN) end

sol=path_network_solve(prob, grid; h=0.5, order=30, rng_seed=0, extrapolate=true)
VZ=sol.visited_z; VH=sol.visited_h; VP=sol.visited_pade
poles=extract_poles(sol)
@printf("nodes=%d poles=%d\n", length(VZ), length(poles))
open(SD*"/poles_e_full.txt","w") do io; for p in poles; @printf(io,"%.5f %.5f\n",real(p),imag(p)); end; end

function nearest2(z)
    d1=Inf;d2=Inf;i1=0;i2=0
    @inbounds for k in eachindex(VZ)
        d=abs2(z-VZ[k]); if d<d1; d2=d1;i2=i1;d1=d;i1=k elseif d<d2; d2=d;i2=k end
    end; i1,i2
end
function jmap()
    J=fill(-1.0,N*N); mx=0.0; mxd=0.0; ns=0; nsd=0; zmax=0.0+0im
    for (idx,z) in enumerate(grid)
        i1,i2=nearest2(z); (i1==0||i2==0)&&continue
        t1=(z-VZ[i1])/VH[i1]; t2=(z-VZ[i2])/VH[i2]
        u1=ev(VP[i1],t1); u2=ev(VP[i2],t2); (isfinite(u1)&&isfinite(u2))||continue
        jj=abs(u1-u2)/(abs(u1)+abs(u2)+1e-30); J[idx]=jj
        any(p->abs(p-z)<1.2,poles) && continue
        if jj>mx; mx=jj; zmax=z; end
        (jj>1e-2)&&(ns+=1)
        if abs(t1)<=1 && abs(t2)<=1; mxd=max(mxd,jj); (jj>1e-2)&&(nsd+=1); end
    end
    (J,mx,mxd,ns,nsd,zmax)
end
J,mx,mxd,ns,nsd,zmax=jmap()
@printf("\n== P2 jump map (e) [-50,50]^2, OFF-POLE ==\n")
@printf("max J = %.3e at z=%+.1f%+.1fi\n", mx, real(zmax),imag(zmax))
@printf("max J in-disc (H1 pure) = %.3e\n", mxd)
@printf("off-pole cells J>1e-2: total=%d in-disc=%d (of %d)\n", ns,nsd,N*N)
open(SD*"/jmap_full.txt","w") do io; @printf(io,"%d %f\n",N,HALF); for v in J; @printf(io,"%.6e\n",v); end; end

sol2=path_network_solve(prob, grid; h=0.5, order=30, rng_seed=42, extrapolate=true)
function seeddiff()
    dmax=0.0; nd=0; nt=0; zmx=0.0+0im
    for k in 1:N*N
        u0=sol.grid_u[k]; u1=sol2.grid_u[k]; (isfinite(u0)&&isfinite(u1))||continue
        nt+=1; rel=abs(u0-u1)/(abs(u0)+abs(u1)+1e-30); (rel>dmax)&&(dmax=rel;zmx=grid[k]); (rel>1e-2)&&(nd+=1)
    end; (dmax,nd,nt,zmx)
end
dmax,nd,nt,zmx=seeddiff()
@printf("\n== two-seed field diff (rng 0 vs 42) ==\n")
@printf("max rel = %.3e at z=%+.1f%+.1fi ; cells rel>1e-2 = %d/%d\n", dmax, real(zmx),imag(zmx), nd, nt)
