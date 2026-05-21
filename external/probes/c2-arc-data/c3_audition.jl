# external/probes/c2-arc-data/c3_audition.jl
#
# C3 audition (bead padetaylor-0ln.37.11).  Voter 1 (`surf_ray_eval`)
# bilinearly interpolates the ray-fan BVP solution from a polar (r,φ)
# raster onto the Cartesian grid.  Bilinear interpolation is NOT a
# harmonic reconstruction.  Question: at the production resolution, is
# the bilinear inter-ray interpolation error material — does it dominate
# the triple-method vote — and would a harmonic voter-1 be better while
# keeping the three voters independent?
#
# Method: measure the bilinear interpolation error of voter 1 against an
# INDEPENDENT exact reference — a DEDICATED ray BVP solved straight
# through each test point's own angle (no interpolation at all).  The
# midpoint of two fan rays is where bilinear interpolation is worst.
#
# Run:  julia --project=figures external/probes/c2-arc-data/c3_audition.jl

using PadeTaylor
using PadeTaylor.PainleveHierarchy: painleve_hierarchy,
                                    painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP: vector_bvp_solve

const T = 0.0
const D = 4
const RMIN = 2.0
const RMAX = 20.0
const RINNER = 1.05            # the C2 extended inner BC

f  = painleve_hierarchy(:I, 2; t = T)
Jf = painleve_hierarchy_jacobian(:I, 2; t = T)

function ray_bvp(φ; r_out = RMAX, r_in = RINNER, N = 96)
    CT = ComplexF64
    z_a = r_out * cis(φ); z_b = r_in * cis(φ)
    Ba = zeros(CT,D,D); Bb = zeros(CT,D,D)
    Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
    sa = pI2_tritronquee_ic(z_a; t=T, n_terms=2)
    sb = pI2_tritronquee_ic(z_b; t=T, n_terms=2)
    g = CT[sa[1], sa[2], sb[1], sb[2]]
    return vector_bvp_solve(f, z_a, z_b, Ba, Bb, g;
                            N=N, tol=1e-9, maxiter=40, jacobian=Jf,
                            initial_guess = z->pI2_tritronquee_ic(z;t=T,n_terms=2))
end

# Build the production fan (C2-extended: inner BC at RINNER, raster on
# [RMIN, RMAX] so every raster point is interior).
φ_lo = deg2rad(36.0 + 4.0)
φ_hi = deg2rad(360.0 - 36.0 - 4.0)
nφ   = max(2, round(Int, (φ_hi - φ_lo) / deg2rad(6.0)) + 1)
phis = collect(range(φ_lo, φ_hi; length = nφ))
nr   = 40
radii = collect(range(RMIN, RMAX; length = nr))
U = Matrix{ComplexF64}(undef, nr, nφ)
for (j, φ) in enumerate(phis)
    sol = ray_bvp(φ)
    for (i, r) in enumerate(radii)
        U[i, j] = sol(ComplexF64(r * cis(φ)))[1]
    end
end

# Bilinear polar interpolation — `surf_ray_eval`'s scheme verbatim.
function bilin(z)
    r = abs(z); φ = mod2pi(angle(z))
    (r < radii[1] || r > radii[end]) && return ComplexF64(NaN,NaN)
    (φ < phis[1]  || φ > phis[end])  && return ComplexF64(NaN,NaN)
    i = clamp(searchsortedlast(radii, r), 1, length(radii)-1)
    j = clamp(searchsortedlast(phis, φ), 1, length(phis)-1)
    tr = (r-radii[i])/(radii[i+1]-radii[i])
    tφ = (φ-phis[j])/(phis[j+1]-phis[j])
    return (1-tr)*(1-tφ)*U[i,j] + tr*(1-tφ)*U[i+1,j] +
           (1-tr)*tφ*U[i,j+1]   + tr*tφ*U[i+1,j+1]
end

println("="^72)
println("C3 audition — bilinear voter-1 interpolation error")
println("  production fan: $(nφ) rays @ 6°, $(nr) radial samples")
println("="^72)

# Worst-case: inter-ray midpoints (mid-φ between two fan rays, several
# radii).  An EXACT dedicated BVP straight through that angle is the
# independent oracle — no interpolation.
errs_mid = Float64[]
errs_node = Float64[]
for jr in 6:6:(nφ-6)
    φ_mid = (phis[jr] + phis[jr+1]) / 2
    sol_mid = ray_bvp(φ_mid)               # exact, dedicated ray
    for r in (3.0, 6.0, 10.0, 15.0, 19.0)
        z = ComplexF64(r * cis(φ_mid))
        u_exact = sol_mid(z)[1]
        u_bilin = bilin(z)
        push!(errs_mid, abs(u_bilin - u_exact))
    end
    # also at an on-ray angle (bilinear is exact in φ there — only the
    # radial interpolation error remains, and that is between raster
    # rows of one smooth barycentric ray solution).
    φ_on = phis[jr]
    sol_on = ray_bvp(φ_on)
    for r in (3.5, 6.5, 10.5)
        z = ComplexF64(r * cis(φ_on))
        push!(errs_node, abs(bilin(z) - sol_on(z)[1]))
    end
end
sort!(errs_mid); sort!(errs_node)
println("\n  inter-ray-midpoint bilinear error (the worst case):")
println("    median = ", round(errs_mid[cld(end,2)]; sigdigits=3),
        "   max = ", round(errs_mid[end]; sigdigits=3))
println("  on-ray radial-only bilinear error:")
println("    median = ", round(errs_node[cld(end,2)]; sigdigits=3),
        "   max = ", round(errs_node[end]; sigdigits=3))
println("\n  Reference scales:")
println("    bulk triple-method spread (worklog 057)  ≈ 5.5e-4")
println("    FEM voter (3) accuracy O(h²)             ≈ 5e-5 .. 5e-4")
println("    spectral voter (2) accuracy              ≈ 1e-8 (geometric)")

# ----------------------------------------------------------------------
# Audition the candidate voter-1 reconstructions.  Keep voter 1 a DIRECT
# ODE solve (independent of the two Laplace voters); only improve the
# polar-raster → Cartesian sampling.
#
#   bilinear        : current — 40-row raster + bilinear (r,φ).
#   exact-r-linφ    : store the ray VectorBVPSolutions; evaluate the two
#                     bracketing rays EXACTLY at the query radius
#                     (barycentric, ~1e-9), linearly interpolate in φ.
#   exact-r-cubicφ  : same, but a 4-ray cubic (Catmull-Rom) angular blend.
# ----------------------------------------------------------------------
println("\n" * "="^72)
println("C3 — voter-1 reconstruction audition")
println("="^72)

# Store the full ray BVP solutions (the fan already solves them).
sols = [ray_bvp(φ) for φ in phis]

# exact-r, linear-φ: two bracketing rays, each exact in r.
function v1_linphi(z)
    r = abs(z); φ = mod2pi(angle(z))
    (r < RMIN || r > RMAX) && return ComplexF64(NaN,NaN)
    (φ < phis[1] || φ > phis[end]) && return ComplexF64(NaN,NaN)
    j = clamp(searchsortedlast(phis, φ), 1, length(phis)-1)
    tφ = (φ-phis[j])/(phis[j+1]-phis[j])
    ua = sols[j](ComplexF64(r*cis(phis[j])))[1]
    ub = sols[j+1](ComplexF64(r*cis(phis[j+1])))[1]
    return (1-tφ)*ua + tφ*ub
end

# exact-r, cubic-φ: 4-ray Catmull-Rom in the angular direction.
function v1_cubicphi(z)
    r = abs(z); φ = mod2pi(angle(z))
    (r < RMIN || r > RMAX) && return ComplexF64(NaN,NaN)
    (φ < phis[1] || φ > phis[end]) && return ComplexF64(NaN,NaN)
    j = clamp(searchsortedlast(phis, φ), 1, length(phis)-1)
    j0 = clamp(j-1,1,length(phis)); j1 = j
    j2 = clamp(j+1,1,length(phis)); j3 = clamp(j+2,1,length(phis))
    tφ = (φ-phis[j1])/(phis[j2]-phis[j1])
    p0 = sols[j0](ComplexF64(r*cis(phis[j0])))[1]
    p1 = sols[j1](ComplexF64(r*cis(phis[j1])))[1]
    p2 = sols[j2](ComplexF64(r*cis(phis[j2])))[1]
    p3 = sols[j3](ComplexF64(r*cis(phis[j3])))[1]
    # Catmull-Rom basis.
    t = tφ
    return 0.5*((2p1) + (-p0+p2)*t + (2p0-5p1+4p2-p3)*t^2 +
                (-p0+3p1-3p2+p3)*t^3)
end

for (name, vfun) in (("bilinear      ", bilin),
                     ("exact-r linφ  ", v1_linphi),
                     ("exact-r cubicφ", v1_cubicphi))
    em = Float64[]; eo = Float64[]
    for jr in 6:6:(nφ-6)
        φ_mid = (phis[jr]+phis[jr+1])/2
        sol_mid = ray_bvp(φ_mid)
        for r in (3.0,6.0,10.0,15.0,19.0)
            z = ComplexF64(r*cis(φ_mid))
            push!(em, abs(vfun(z) - sol_mid(z)[1]))
        end
        sol_on = sols[jr]
        for r in (3.5,6.5,10.5)
            z = ComplexF64(r*cis(phis[jr]))
            push!(eo, abs(vfun(z) - sol_on(z)[1]))
        end
    end
    sort!(em); sort!(eo)
    println("  $name : inter-ray  median=", round(em[cld(end,2)];sigdigits=3),
            " max=", round(em[end];sigdigits=3),
            " | on-ray median=", round(eo[cld(end,2)];sigdigits=3),
            " max=", round(eo[end];sigdigits=3))
end
println("\n  A harmonic voter-1 (Laplace solve on the w-rect with ray-fan")
println("  edge data) ≡ voter 2's machinery — would destroy ADR-0024")
println("  triple-method independence.  The audition keeps voter 1 a")
println("  direct ODE solve and only fixes the raster→Cartesian sampling.")
