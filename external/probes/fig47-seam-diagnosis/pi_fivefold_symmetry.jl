# PI panel (a): u''=6u^2+z, u(0)=u'(0)=0  — the Z5-symmetric solution.
# Its pole field must have 5-fold rotational symmetry (z -> e^{2πi/5} z)
# AND conjugate symmetry (real IC).  Test both on the actual extracted poles.
using PadeTaylor
using Printf

pI(z,u,up) = 6u^2 + z
const HALF = 25.0
N = 31
xs = range(-HALF,HALF; length=N); ys = range(-HALF,HALF; length=N)
grid = ComplexF64[complex(x,y) for y in ys for x in xs]
prob = PadeTaylorProblem(pI,(0.0,0.0),(0.0,HALF); order=30)
t0=time()
sol = path_network_solve(prob, grid; h=0.5, order=30)
poles = extract_poles(sol)
@printf("plain path_network: %.1f s ; nodes=%d ; poles=%d\n", time()-t0, length(sol.visited_z), length(poles))

# dump coords for rendering
open("/tmp/claude-1000/-home-tobias-Projects-PadeTaylor-jl/e817ae45-f9a1-4b4a-818b-2d5e3e929f4b/scratchpad/pi_a_poles.txt","w") do io
    for p in poles; @printf(io,"%.6f %.6f\n", real(p), imag(p)); end
end

inwin = filter(p-> -HALF<=real(p)<=HALF && -HALF<=imag(p)<=HALF, poles)
tol = 0.6
cm(set,q)=count(x->abs(x-q)<tol,set)

# conjugate symmetry
let up=filter(p->imag(p)>1.0, inwin)
    s=count(p->cm(inwin,conj(p))>0, up)
    @printf("conjugate-symmetric: %d/%d = %.3f\n", s, length(up), s/max(length(up),1))
end
# 5-fold rotational symmetry: rotate by k*72°, fraction of poles whose image is also a pole
for k in 1:4
    rot = cis(2pi*k/5)
    s = count(p-> (r=p*rot; -HALF<=real(r)<=HALF && -HALF<=imag(r)<=HALF) ? cm(inwin, p*rot)>0 : true, inwin)
    # only count those whose rotated image stays in-window
    inw = filter(p->(r=p*rot; -HALF<=real(r)<=HALF && -HALF<=imag(r)<=HALF), inwin)
    s2  = count(p->cm(inwin,p*rot)>0, inw)
    @printf("rotation %d*72°=%3d°: %d/%d in-window-rotated poles map onto a pole = %.3f\n",
            k, round(Int,72k), s2, length(inw), s2/max(length(inw),1))
end
# angular histogram (12 bins) of in-window poles, |z|>5
let far=filter(p->abs(p)>5, inwin)
    h=zeros(Int,12)
    for p in far; b=mod(floor(Int, rad2deg(angle(p))/30),12)+1; h[b]+=1; end
    @printf("angular histogram (30° bins, |z|>5): %s\n", h)
end
