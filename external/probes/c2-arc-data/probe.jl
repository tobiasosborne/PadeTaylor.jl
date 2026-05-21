# external/probes/c2-arc-data/probe.jl
#
# C2 exploration probe (bead padetaylor-0ln.37.10).  Question: does the
# ray-fan BVP give an accurate inner-arc datum, or is its value at
# |x| = R_MIN merely the O(1e-2) asymptotic seed (because the BVP's 2+2
# split BC PINS (u,u') at the inner endpoint to that very seed)?
#
# If the inner endpoint of the BVP is pinned to the seed, then
# `fan.U[1,:]` IS the seed, and "source voters 2/3 arc data from the
# ray fan" buys nothing.  The genuine fix must place R_MIN in the BVP
# INTERIOR so the harmonic voters read the true ODE solution there.
#
# Run:  julia --project=figures external/probes/c2-arc-data/probe.jl

using PadeTaylor
using PadeTaylor.PainleveHierarchy: painleve_hierarchy,
                                    painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP: vector_bvp_solve

const T   = 0.0
const D   = 4
const RMIN = 2.0
const RMAX = 20.0

f  = painleve_hierarchy(:I, 2; t = T)
Jf = painleve_hierarchy_jacobian(:I, 2; t = T)

# A 2+2 split-BC ray BVP over [r_out, r_in], (u,u') pinned at both ends
# from the n_terms=2 asymptotic seed.  This is `surf_ray_bvp`'s recipe.
function ray_bvp(φ, r_out, r_in; N = 96)
    CT = ComplexF64
    z_a = r_out * cis(φ); z_b = r_in * cis(φ)
    Ba = zeros(CT, D, D); Bb = zeros(CT, D, D)
    Ba[1,1] = 1; Ba[2,2] = 1
    Bb[3,1] = 1; Bb[4,2] = 1
    sa = pI2_tritronquee_ic(z_a; t = T, n_terms = 2)
    sb = pI2_tritronquee_ic(z_b; t = T, n_terms = 2)
    g = CT[sa[1], sa[2], sb[1], sb[2]]
    return vector_bvp_solve(f, z_a, z_b, Ba, Bb, g;
                            N = N, tol = 1e-9, maxiter = 40, jacobian = Jf,
                            initial_guess = z -> pI2_tritronquee_ic(z; t = T,
                                                                    n_terms = 2))
end

println("="^72)
println("C2 probe — ray-fan BVP inner-arc datum vs the asymptotic seed")
println("="^72)

for φ in (deg2rad(180.0), deg2rad(120.0), deg2rad(240.0))
    println("\n--- ray φ = $(round(rad2deg(φ); digits=1))° ---")

    # The production BVP: [RMAX, RMIN], pinned at RMIN.
    prod = ray_bvp(φ, RMAX, RMIN)
    u_prod_inner = prod(ComplexF64(RMIN * cis(φ)))[1]

    # The seed at RMIN.
    u_seed_inner = pI2_tritronquee_ic(ComplexF64(RMIN * cis(φ)); t = T,
                                      n_terms = 2)[1]

    # Production BVP value at the inner endpoint vs the seed there.
    println("  |prod(RMIN) - seed(RMIN)|         = ",
            abs(u_prod_inner - u_seed_inner),
            "   (== 0  ⇒ inner endpoint is pinned to the seed)")

    # ORACLE-A: a BVP whose segment EXTENDS INWARD past RMIN, so RMIN is
    # an interior collocation point.  Its value at RMIN is the genuine
    # ODE solution interpolated between endpoints — not a pinned seed.
    # But its inner endpoint seed (smaller |x|) is even less accurate.
    for r_in in (1.5, 1.2, 1.05)
        try
            ext = ray_bvp(φ, RMAX, r_in)
            u_ext_inner = ext(ComplexF64(RMIN * cis(φ)))[1]
            println("  extend-inward to r=$r_in : u(RMIN) = ", u_ext_inner)
            println("      vs prod-endpoint        Δ = ", abs(u_ext_inner - u_prod_inner))
            println("      vs seed(RMIN)           Δ = ", abs(u_ext_inner - u_seed_inner))
        catch err
            println("  extend-inward to r=$r_in : BVP FAILED ($(typeof(err)))")
        end
    end

    # ORACLE-B: solve the SAME [RMAX,RMIN] segment but at much higher N
    # and with the inner BC RELAXED.  We cannot relax the BC count
    # (need 2d conditions), but we can test sensitivity: re-pin the
    # inner end with the seed at n_terms=1 vs n_terms=2 and see how much
    # the inner value moves.  Large movement ⇒ inner endpoint is
    # seed-dominated; small ⇒ ODE-constrained.
    CT = ComplexF64
    z_a = RMAX * cis(φ); z_b = RMIN * cis(φ)
    Ba = zeros(CT,D,D); Bb = zeros(CT,D,D)
    Ba[1,1]=1; Ba[2,2]=1; Bb[3,1]=1; Bb[4,2]=1
    sa  = pI2_tritronquee_ic(z_a; t=T, n_terms=2)
    sb2 = pI2_tritronquee_ic(z_b; t=T, n_terms=2)
    sb1 = pI2_tritronquee_ic(z_b; t=T, n_terms=1)
    g1 = CT[sa[1], sa[2], sb1[1], sb1[2]]
    bvp_n1 = vector_bvp_solve(f, z_a, z_b, Ba, Bb, g1; N=96, tol=1e-9,
                              maxiter=40, jacobian=Jf,
                              initial_guess = z->pI2_tritronquee_ic(z;t=T,n_terms=2))
    u_n1_inner = bvp_n1(z_b)[1]
    println("  inner BC n_terms=1 vs 2 : u(RMIN) moves Δ = ",
            abs(u_n1_inner - u_prod_inner),
            "  (seed n1↔n2 differ by ", abs(sb1[1] - sb2[1]), ")")

    # ORACLE-C: convergence of the extend-inward family.  Solve at high
    # N with the inner BC ever closer to |x|=1; the value at RMIN should
    # converge if RMIN-as-interior is a genuine ODE datum.
    refs = ComplexF64[]
    for (r_in, Nn) in ((1.5, 128), (1.2, 160), (1.05, 192))
        ext = ray_bvp(φ, RMAX, r_in; N = Nn)
        push!(refs, ext(ComplexF64(RMIN * cis(φ)))[1])
    end
    println("  extend-inward hi-N convergence : ",
            "Δ(1.5↔1.2)=", round(abs(refs[1]-refs[2]); sigdigits=3),
            "  Δ(1.2↔1.05)=", round(abs(refs[2]-refs[3]); sigdigits=3))
    println("  ⇒ best inner-arc oracle u(RMIN) ≈ ", refs[3])
    println("    seed(RMIN) error vs oracle    = ", abs(u_seed_inner - refs[3]))
end

println("\n" * "="^72)
println("Interpretation:")
println(" - If prod(RMIN)==seed(RMIN) exactly, the ray fan inner ARC datum")
println("   is identically the asymptotic seed — sourcing voters 2/3 from")
println("   `fan.U[1,:]` does NOT improve accuracy, only collapses spread.")
println(" - The extend-inward BVP places RMIN in the interior; its u(RMIN)")
println("   is the genuine ODE solution.  Compare its accuracy.")
println("="^72)

# ----------------------------------------------------------------------
# Robustness sweep: does the extend-inward BVP succeed on EVERY fan ray,
# and which inner-BC radius is the senior-grade choice?
# ----------------------------------------------------------------------
println("\n" * "="^72)
println("Robustness sweep — extend-inward BVP across the full fan")
println("="^72)
φ_lo = deg2rad(36.0 + 4.0)
φ_hi = deg2rad(360.0 - 36.0 - 4.0)
nφ   = max(2, round(Int, (φ_hi - φ_lo) / deg2rad(6.0)) + 1)
phis = collect(range(φ_lo, φ_hi; length = nφ))
for r_in in (1.05, 1.2, 1.5)
    nfail = 0; maxrange = 0.0
    for φ in phis
        try
            ext = ray_bvp(φ, RMAX, r_in; N = 96)
            u = ext(ComplexF64(RMIN * cis(φ)))[1]
            isfinite(u) || (nfail += 1)
        catch
            nfail += 1
        end
    end
    println("  r_in=$r_in : $(nfail)/$(length(phis)) ray BVPs failed")
end

# Accuracy of the extend-inward inner datum vs the converged oracle,
# averaged over the fan, for the chosen r_in.
for r_in in (1.05, 1.2)
    errs = Float64[]
    for φ in (deg2rad(60.0), deg2rad(120.0), deg2rad(180.0),
              deg2rad(240.0), deg2rad(300.0))
        oracle = ray_bvp(φ, RMAX, 1.02; N = 200)(ComplexF64(RMIN*cis(φ)))[1]
        ext    = ray_bvp(φ, RMAX, r_in; N = 96)(ComplexF64(RMIN*cis(φ)))[1]
        seed   = pI2_tritronquee_ic(ComplexF64(RMIN*cis(φ)); t=T, n_terms=2)[1]
        push!(errs, abs(ext - oracle))
        println("  φ=$(round(rad2deg(φ);digits=0))° r_in=$r_in : ",
                "ext err=", round(abs(ext-oracle); sigdigits=3),
                "  seed err=", round(abs(seed-oracle); sigdigits=3))
    end
    println("  r_in=$r_in : mean ext-vs-oracle inner-arc error = ",
            round(sum(errs)/length(errs); sigdigits=3))
end

# ----------------------------------------------------------------------
# Outer-arc check: the seed at RMAX=20 is already ~1e-6 accurate (VC-8).
# Does extending the BVP outward (RMAX interior) change u(RMAX)?
# ----------------------------------------------------------------------
println("\n" * "="^72)
println("Outer-arc — seed accuracy at RMAX and extend-outward effect")
println("="^72)
for φ in (deg2rad(120.0), deg2rad(180.0))
    prod_out = ray_bvp(φ, RMAX, 1.05; N=96)(ComplexF64(RMAX*cis(φ)))[1]
    ext_out  = ray_bvp(φ, 24.0, 1.05; N=110)(ComplexF64(RMAX*cis(φ)))[1]
    seed_out = pI2_tritronquee_ic(ComplexF64(RMAX*cis(φ)); t=T, n_terms=2)[1]
    println("  φ=$(round(rad2deg(φ);digits=0))° : ",
            "|prod_endpoint - seed| = ", round(abs(prod_out-seed_out); sigdigits=3),
            "  |extend-out - seed|  = ", round(abs(ext_out-seed_out); sigdigits=3),
            "  |extend-out - prod|  = ", round(abs(ext_out-prod_out); sigdigits=3))
end
