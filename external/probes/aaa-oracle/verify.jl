# =============================================================================
# verify.jl — cross-check the AAA oracle (aaa.jl) against the shared
# denominator of PadeTaylor.SharedPade.shared_denominator_pade.
#
# This is the V1d verification record — the *third independent oracle* for
# the shared-Q keystone, alongside Calgo 766 (calgo766-oracle/) and the
# block-Toeplitz determinant (shared-pade-determinant-oracle/).  Where those
# two oracles consume Taylor jets, AAA consumes *sampled function values*:
# the verification below samples the SP.1.1 and SP.1.2 test functions on a
# contour in their analyticity region, fits each component independently
# with AAA, and checks that the AAA pole locations land on the roots of the
# shared Q recovered by SharedPade.
#
# Per docs/v0p2_pillarA_hermite_pade_findings.md:354–373 (§6 "AAA as a
# benchmark"): AAA fits each component *separately*, so for the d=2 case
# SP.1.2 it produces d=2 INDEPENDENT pole sets.  The cross-validation
# evidence is the *coincidence* of those independent pole sets with each
# other and with the shared-Q roots: AAA cannot enforce a shared
# denominator, so when its per-component poles agree anyway, that agreement
# is independent confirmation the shared Q is correct.  AAA is a diagnostic
# pole cross-check, NOT a coefficient-level ground truth (§6, "limitation").
#
# This script does NOT add asserts to the project test suite — bead V1e
# wires aaa.jl into test/ separately.  It is a verification record.
#
# Run (one Julia process — CLAUDE.md Rule 7):
#   julia --project=../../.. external/probes/aaa-oracle/verify.jl
# =============================================================================

using LinearAlgebra
using Printf
using PadeTaylor: shared_denominator_pade

include(joinpath(@__DIR__, "aaa.jl"))
using .AAAOracle: aaa

# --- helpers -----------------------------------------------------------------

# Formal power-series long division: Taylor coefficients [c_0 … c_N] of P/Q.
# Mirrors test/shared_pade_test.jl `_ratio_jet` so the jets fed to SharedPade
# are bit-identical to the project test inputs.
function ratio_jet(P::Vector{Float64}, Q::Vector{Float64}, N::Int)
    c = zeros(Float64, N + 1)
    for k = 0:N
        s = k + 1 ≤ length(P) ? P[k + 1] : 0.0
        for j = 1:k
            qj = j + 1 ≤ length(Q) ? Q[j + 1] : 0.0
            s -= qj * c[k - j + 1]
        end
        c[k + 1] = s / Q[1]
    end
    return c
end

polyval(coef, z)    = sum(coef[k] * z^(k - 1) for k = 1:length(coef))
ratval(num, den, z) = polyval(num, z) / polyval(den, z)

# Roots of a degree-≤2 polynomial (low-to-high coeffs), sorted, complex.
# These are the SHARED POLES recovered by SharedPade — the reference set
# the AAA poles are cross-checked against.
function roots_of(q)
    a, b, c = q[3], q[2], q[1]
    disc = sqrt(complex(b^2 - 4a * c))
    sort([(-b + disc) / (2a), (-b - disc) / (2a)], by = z -> (real(z), imag(z)))
end

# Match each pole of `set` to its nearest neighbour in `ref` and return the
# worst pairing distance.  AAA returns its poles in an arbitrary order, so a
# nearest-neighbour pairing is the honest comparison.
function max_pairing_dist(set, ref)
    isempty(set) && return Inf
    worst = 0.0
    for p in ref
        d = minimum(abs(p - q) for q in set)
        worst = max(worst, d)
    end
    return worst
end

# =============================================================================
# Sampling contour.
#
# AAA fits from VALUES, so the sample points must avoid the poles (a sample
# at a pole is Inf — aaa.jl throws on that, Rule 1).  Both SP.1.1 and SP.1.2
# have their poles OUTSIDE the unit disk:
#   SP.1.1  Q = 1 - 0.5z + 0.3z²  → |roots| = 1/√0.3 ≈ 1.826
#   SP.1.2  Q = 1 - 0.4z + 0.2z²  → |roots| = 1/√0.2 ≈ 2.236
# so the unit circle |z| = 1 lies safely inside the analyticity region and
# the functions are smooth and bounded there.  We sample 200 equispaced
# points on |z| = 1; AAA's greedy selection then picks its own support
# points from this set.  (NST §3 needs m ≤ M/2; 200 points admits types up
# to (99,99), far more than the (2,2) approximants here.)
# =============================================================================

const N_SAMPLE = 200
const SAMPLE_Z = ComplexF64[exp(2im * pi * k / N_SAMPLE) for k = 0:(N_SAMPLE - 1)]

# =============================================================================
# Test cases — SP.1.1 and SP.1.2 verbatim from test/shared_pade_test.jl.
# =============================================================================

sp11 = (
    name = "SP.1.1 d=1 scalar rational",
    m    = 2,
    P    = [[1.0, 2.0]],
    Q    = [1.0, -0.5, 0.3],
)

sp12 = (
    name = "SP.1.2 d=2 shared denominator",
    m    = 2,
    P    = [[1.0, 0.3], [2.0, -0.7, 0.5]],
    Q    = [1.0, -0.4, 0.2],
)

# Pinned tolerances (documented in README.md).
const TOL_AAA_VS_AAA   = 1e-11   # per-component AAA pole sets agree this well
const TOL_AAA_VS_TRUTH = 1e-11   # AAA poles vs the exact known Q roots
const TOL_AAA_VS_SHARE = 1e-9    # AAA poles vs SharedPade's shared-Q roots

# =============================================================================
# Run the verification for one case.
# =============================================================================

function verify_case(case)
    println("=" ^ 74)
    println(case.name, "   (m = ", case.m, ", d = ", length(case.P), ")")
    println("-" ^ 74)
    m = case.m
    d = length(case.P)
    all_ok = true

    # --- (1) SharedPade: recover the shared Q from Taylor jets -------------
    jets = [ratio_jet(case.P[i], case.Q, 2m) for i = 1:d]
    _, den_shared = shared_denominator_pade(jets, m)
    roots_shared = roots_of(den_shared)
    roots_truth  = roots_of(case.Q)
    fmt_roots(rs) = join((@sprintf("%.6f%+.6fim", real(r), imag(r)) for r in rs), "  ")
    println("  SharedPade shared-Q roots:   ", fmt_roots(roots_shared))
    println("  exact known-Q roots:         ", fmt_roots(roots_truth))

    # --- (2) AAA: fit each component independently from sampled values ----
    # Per pillar A §6: AAA fits EACH component separately, so we obtain d
    # independent pole sets.  The cross-validation is their mutual + shared-Q
    # coincidence.
    aaa_pole_sets = Vector{Vector{ComplexF64}}(undef, d)
    for i = 1:d
        # Sample component i = P_i/Q on the unit-circle contour.
        Fvals = [ratval(case.P[i], case.Q, z) for z in SAMPLE_Z]
        res = aaa(SAMPLE_Z, Fvals; tol = 1e-13, cleanup = true)

        # AAA's type-(m-1,m-1) approximant of a genuine type-(≤2,2) rational
        # should converge with m = 3 support points and expose 2 poles.
        aaa_pole_sets[i] = res.poles
        derr_truth = max_pairing_dist(res.poles, roots_truth)
        @printf("  AAA component %d:  %d support pts, %d poles, |Δ poles vs truth| = %.3e\n",
                i, length(res.support_points), length(res.poles), derr_truth)
        all_ok &= (length(res.poles) == 2) && (derr_truth < TOL_AAA_VS_TRUTH)
    end

    # --- (3) per-component AAA pole sets must coincide with each other ----
    if d ≥ 2
        worst_cross = 0.0
        for i = 2:d
            worst_cross = max(worst_cross,
                              max_pairing_dist(aaa_pole_sets[i], aaa_pole_sets[1]))
        end
        @printf("  per-component AAA poles agree: |Δ| = %.3e   (tol %.0e)\n",
                worst_cross, TOL_AAA_VS_AAA)
        all_ok &= (worst_cross < TOL_AAA_VS_AAA)
    else
        println("  per-component AAA poles agree: (n/a — d = 1)")
    end

    # --- (4) AAA poles vs SharedPade's shared-Q roots — the cross-check ----
    worst_vs_shared = 0.0
    for i = 1:d
        worst_vs_shared = max(worst_vs_shared,
                              max_pairing_dist(aaa_pole_sets[i], roots_shared))
    end
    @printf("  AAA poles vs SharedPade Q:   |Δ| = %.3e   (tol %.0e)\n",
            worst_vs_shared, TOL_AAA_VS_SHARE)
    all_ok &= (worst_vs_shared < TOL_AAA_VS_SHARE)

    println("-" ^ 74)
    println(all_ok ? "  CASE PASSED — AAA poles coincide with the shared Q" :
                     "  CASE FAILED — investigate (CLAUDE.md Rule 2)")
    return all_ok
end

# --- failure-mode check: a sample on a pole must throw (Rule 1) --------------
function verify_throws()
    println("=" ^ 74)
    println("Failure-mode check: a non-finite sample value must throw (Rule 1)")
    println("-" ^ 74)
    # Sample SP.1.1's f on a contour passing through a pole: |root| ≈ 1.826,
    # so a circle of that radius hits the pole and produces an Inf value.
    Q = sp11.Q
    rpole = abs(roots_of(Q)[1])
    Zbad = ComplexF64[rpole * exp(2im * pi * k / 64) for k = 0:63]
    Fbad = [ratval(sp11.P[1], Q, z) for z in Zbad]
    # Force an exact Inf by placing one sample exactly on a real-axis root if
    # the contour does not already produce one numerically.
    if all(isfinite, Fbad)
        Fbad[1] = Inf + 0im
    end
    threw = false
    try
        aaa(Zbad, Fbad)
    catch e
        threw = true
        @printf("  threw as expected: %s\n",
                first(split(sprint(showerror, e), '\n')))
    end
    println(threw ? "  PASSED" : "  FAILED — aaa did not throw on an Inf sample")
    return threw
end

function main()
    println("\nAAA per-component pole cross-check oracle — verification (bead V1d)\n")
    ok = true
    ok &= verify_case(sp11)
    ok &= verify_case(sp12)
    ok &= verify_throws()
    println("=" ^ 74)
    if ok
        println("VERIFICATION PASSED: AAA — fitting each component independently")
        println("from SAMPLED VALUES on the unit circle — recovers, for every")
        println("component, the SAME pole pair, and those poles coincide with the")
        println("roots of SharedPade's shared denominator Q.  Independent of the")
        println("Taylor-jet route: orthogonal third oracle for the shared-Q keystone.")
    else
        println("VERIFICATION FAILED: see the per-case lines above.")
    end
    return ok
end

main() || exit(1)
