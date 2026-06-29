# Ground-truth pole-extraction test on the Weierstrass-℘ problem (FW §5.1.1).
# u'' = 6u^2, IC from FW (c1=-1, c2=2): exact poles on the equianharmonic
# (hexagonal) lattice  p(m,n) = 1 + 2ω(m + n e^{iπ/3}),  ω ≈ 1.3629.
# We solve a real-IC pole field and measure: recall vs the KNOWN lattice,
# precision, and conjugate (real-axis) symmetry of the extracted set.
using PadeTaylor
using Printf

f(z,u,up) = 6u^2
u0  = 1.071822516416917
up0 = 1.710337353176786

const HALF = 12.0
N = 25
xs = range(-HALF, HALF; length=N)
ys = range(-HALF, HALF; length=N)
grid = ComplexF64[complex(x,y) for y in ys for x in xs]

prob = PadeTaylorProblem(f, (u0,up0), (0.0, HALF); order=30)
t0 = time()
sol = path_network_solve(prob, grid; h=0.5, order=30)
poles = extract_poles(sol)
@printf("solve+extract: %.1f s ; visited nodes=%d ; poles extracted=%d\n",
        time()-t0, length(sol.visited_z), length(poles))

# --- known exact lattice -------------------------------------------------
ω = 1.3629079683730673                # Γ(1/3)^3 / (2^(13/6) π)
gen1 = 2ω + 0im
gen2 = 2ω*cis(pi/3)
exact = ComplexF64[]
for m in -12:12, n in -12:12
    p = 1 + m*gen1 + n*gen2
    if -HALF<=real(p)<=HALF && -HALF<=imag(p)<=HALF
        push!(exact, p)
    end
end
@printf("exact lattice poles inside window: %d\n", length(exact))

const tol = 0.3
countmatch(set, q) = count(x->abs(x-q)<tol, set)
recall_all() = count(e->countmatch(poles,e)>0, exact)
precision_all() = count(p->countmatch(exact,p)>0, poles)
let found=recall_all(), tp=precision_all()
    @printf("RECALL    = %d/%d = %.3f\n", found, length(exact), found/length(exact))
    @printf("PRECISION = %d/%d = %.3f  (=> %d extracted poles are NOT on the lattice)\n",
            tp, length(poles), tp/max(length(poles),1), length(poles)-tp)
end

# --- conjugate (real-axis) symmetry of the EXTRACTED set -----------------
let up_poles=filter(p->imag(p)>0.5, poles)
    sym = count(p->countmatch(poles,conj(p))>0, up_poles)
    @printf("CONJUGATE-SYMMETRIC: %d/%d upper-half poles have a lower-half mirror = %.3f\n",
            sym, length(up_poles), sym/max(length(up_poles),1))
end

# --- precision vs the FULL lattice (incl. out-of-window) -----------------
fulllat = ComplexF64[]
for m in -16:16, n in -16:16
    push!(fulllat, 1 + m*gen1 + n*gen2)
end
let
    dists = [minimum(abs(p-e) for e in fulllat) for p in poles]
    nspur = count(>=(tol), dists)
    @printf("\nvs FULL lattice (incl. out-of-window): %d/%d extracted are >=%.2f from ANY real pole => GENUINELY SPURIOUS\n",
            nspur, length(poles), tol)
    spur = poles[dists .>= tol]
    @printf("spurious-pole locations (z = x+iy), with |u| context:\n")
    for p in sort(spur, by=abs)[1:min(end,40)]
        @printf("   %+.2f %+.2fi   (|z|=%.1f, arg=%.0f°, dist-to-lattice=%.2f)\n",
                real(p), imag(p), abs(p), rad2deg(angle(p)), minimum(abs(p-e) for e in fulllat))
    end
end

# --- recall split by half-plane (is one half worse?) ---------------------
function recall_in(pred)
    ex = filter(pred, exact)
    (count(e->countmatch(poles,e)>0, ex), length(ex))
end
fu,nu = recall_in(p->imag(p)>1.0)
fl,nl = recall_in(p->imag(p)<-1.0)
fr,nr = recall_in(p->real(p)>1.0)
fL,nL = recall_in(p->real(p)<1.0)
@printf("recall upper-half = %d/%d = %.3f ;  lower-half = %d/%d = %.3f\n", fu,nu,fu/max(nu,1), fl,nl,fl/max(nl,1))
@printf("recall right(+x)  = %d/%d = %.3f ;  left(-x)    = %d/%d = %.3f\n", fr,nr,fr/max(nr,1), fL,nL,fL/max(nL,1))
