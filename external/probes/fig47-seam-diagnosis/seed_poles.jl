# Clincher: the grain-boundary seam is a RANDOM PATH artifact. Solve panel (e)
# with two different walk shuffles and extract poles from each. If the seam
# (orientation grain boundary) sits at a DIFFERENT place for each seed, it is
# not physical structure — it is path-dependence (the bug). A correct
# single-valued solver gives seed-independent poles.
using PadeTaylor, Printf
pI(z,u,up)=6u^2+z
const HALF=50.0; N=101
xs=range(-HALF,HALF;length=N); ys=range(-HALF,HALF;length=N)
grid=ComplexF64[complex(x,y) for y in ys for x in xs]
prob=PadeTaylorProblem(pI,(-0.1875,0.3049),(0.0,HALF);order=30)
SD="/tmp/claude-1000/-home-tobias-Projects-PadeTaylor-jl/e817ae45-f9a1-4b4a-818b-2d5e3e929f4b/scratchpad"
for (seed,fn) in ((0,"poles_s0.txt"),(42,"poles_s42.txt"))
    sol=path_network_solve(prob, grid; h=0.5, order=30, rng_seed=seed, extrapolate=true)
    ps=extract_poles(sol)
    open(SD*"/"*fn,"w") do io; for p in ps; @printf(io,"%.5f %.5f\n",real(p),imag(p)); end; end
    @printf("seed=%d nodes=%d poles=%d -> %s\n", seed, length(sol.visited_z), length(ps), fn)
end
