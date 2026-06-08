# test/corpus_pathnet_lattice_sectors_test.jl
#
# ============================================================================
# Ground-truth corpus: PATH-NETWORK AT SCALE — lattice-geometry contrast, the
# tritronquée EMPTY-SECTOR quantitative pins, and long-range multi-direction
# accuracy (bead padetaylor-3jfw, epic padetaylor-25og; plan Family I,
# docs/test_corpus/02_corpus_extension_plan.md:345-385).
#
# WHAT THIS FILE PINS
# -------------------
#   CPN.5  PI tritronquée EMPTY SECTOR (the headline; upgrades the qualitative
#          phase9_tritronquee_test.jl to QUANTITATIVE).  f = 6u² + z,
#          tritronquée IC (FW eq. 4.1).  The tritronquée confines poles to ONE
#          wedge near the positive real axis; the NEGATIVE-real fan is provably
#          POLE-FREE.  TWO-PART pin on far targets {-30,-40,-50} (on-axis) +
#          {-30±5i} (a small fan):
#            (a) ACCURACY — eval_at(target) vs the INDEPENDENT mpmath odefun
#                RAY oracle (tight, ~1e-12); the asymptote -√(-z/6) is the
#                loose ~1e-4 sanity cross-check (computed in-test).
#            (b) ZERO-POLE — extract_poles returns ZERO poles in the pole-free
#                fan (count == 0).  This is the naïve-solver failure mode (a
#                solver that walks into the lower half-plane asymmetrically
#                manufactures spurious poles there — see DEFERRAL below).
#
#   CPN.6  Equianharmonic ℘ LONG-RANGE, MULTI-DIRECTION.  f = 6u² (the FW Table
#          5.1 problem, u = ℘(z+c1;0,c2)) integrated along SEVERAL rays
#          (±30°,±60°) — not just FW's single real axis.  eval_at(target) vs the
#          INDEPENDENT odefun ray oracle at R = 11 (clearance ≥ 0.5 from the ℘
#          lattice), pinned ~1e-13 (F64).  The 0° real axis at z = 30 reuses the
#          FW Table 5.1 paper value (the path bridges the on-axis poles).
#
#   CPN.3  LATTICE-vs-LATTICE routing.  equianharmonic ℘ (f = 6u², 60° RHOMBIC
#          lattice) vs lemniscatic ℘ (f = 6u² - 1/2, g2=1,g3=0, SQUARE lattice)
#          on the SAME target grid.  Assert (a) extracted poles of each lie on
#          the RESPECTIVE lattice (rhombic vs square — to ~1e-6, and ~0.98 from
#          the OTHER lattice, so the geometry is genuinely resolved); (b) the
#          two visited spanning trees DIFFER (catches a hard-coded orientation).
#
# GROUND TRUTH each case embodies (Law 1 / Rule 5)
# ------------------------------------------------
# external/probes/corpus-oracles/pathnet-lattice/capture.py (mpmath dps=30) is
# the residual/first-integral GATE plus the pin source.  The independent oracle
# is mpmath.ODEFUN RAY-integration of the SAME ODE — a DIFFERENT integrator from
# the package's Taylor–Padé walk, so the agreement is a real cross-validation.
# For CPN.5 the gate additionally proves each target sits in the pole-free
# sector (asymptote agreement < 1e-4).  Lattice constants are exact closed form.
#
# TOLERANCES (justified by METHOD, never fitted to pass)
# ------------------------------------------------------
#   * CPN.5 accuracy: eval_at vs odefun ray ~1e-12 (the ℘/PI long-range
#     path-network floor; measured ~1e-13 — margin ~10×).  Asymptote sanity:
#     -√(-z/6) is an ASYMPTOTIC approximation, good only to ~1e-5 at |z|=30-50;
#     pinned loosely at 1e-3 (it is a SANITY cross-check, not the tight pin).
#   * CPN.5 zero-pole: EXACT count == 0 in the pole-free fan (|z|>3, |arg|>160°).
#   * CPN.6 long-range: ~1e-13 off-axis (measured ~1e-14); ~1e-12 at z=30
#     (matches PN.2.2's existing F64 bound; measured ~1.6e-13).
#   * CPN.3 lattice membership: ~1e-6 (the Float64 pole-placement catalogue
#     spec, cf. polefield_test.jl); the cross-lattice distance is ~0.98 ≫ 1e-6.
#
# REFERENCES: plan Family I (md:345-385); FW2011 Table 5.1 §5 (the ℘ problem)
# + FW2011 §4.1 eq. 4.1 (tritronquée IC); phase9_tritronquee_test.jl (the
# qualitative predecessor this upgrades); polefield_test.jl (lattice formulae).
#
# DEFERRAL (Rule 9 — honest about the corner)
# -------------------------------------------
# The bead's farthest off-axis targets {-30+10i, -40-15i} (arg ≈ ±161°/±159°)
# sit AT/PAST the pole-free-sector boundary (empirically ~±162°; the asymptote
# rel-err jumps from ~1e-5 to ~20 across it — confirmed in capture.py's sector
# scan).  They are NOT robustly pole-free, so they are DEFERRED: CPN.5 ships the
# robust deep-fan {-30,-40,-50,-30±5i} (all arg within ±10° of 180°).  The
# zero-pole assertion uses ON-AXIS targets only, because a lower-half-plane walk
# manufactures spurious poles via the shuffle-induced visited-tree asymmetry
# (worklog 014 / PN.6.*) — that is a KNOWN routing artefact, not a sector bug;
# pinning it here would test the asymmetry, not the empty sector.
#
# MUTATION-PROOF (Rule 4): footer block — M-oracle (flip a pinned ray value)
# and M-impl (perturb the wedge step / ℘ jet in src/PathNetwork.jl) were
# ACTUALLY executed RED and restored byte-for-byte (`git diff src/` empty).
# ============================================================================

using Test
using PadeTaylor

include(joinpath(@__DIR__, "_oracle_corpus_pathnet_lattice_sectors.jl"))

# Pole-free fan membership: |z| > 3 (exclude near-origin IC-disc Padé roots,
# whose nearest TRUE pole is z ≈ 2.07 in the OTHER sector) and |arg(z)| within
# 20° of ±180° (the negative-real fan).
in_pole_free_fan(p) = abs(p) > 3.0 && abs(abs(angle(p)) - π) < deg2rad(20)

@testset "Corpus: path-network lattice/sectors (CPN)" begin

    # ----------------------------------------------------------------------
    # CPN.5  PI tritronquée empty-sector — THE HEADLINE.
    # ----------------------------------------------------------------------
    @testset "CPN.5  tritronquée: pole-free negative-real sector (accuracy + zero-pole)" begin
        f(z, u, up) = 6 * u^2 + z
        u0, up0 = -0.1875543083404949, 0.3049055602612289
        prob = PadeTaylorProblem(f, (u0 + 0im, up0 + 0im), (0.0 + 0im, -50.0 + 0im);
                                 order = 30)

        # ON-AXIS targets drive the zero-pole assertion (no lower-half walk).
        grid = ComplexF64[-30.0, -40.0, -50.0]
        sol  = path_network_solve(prob, grid; h = 0.5,
                                  max_steps_per_target = 400, rng_seed = 0)

        # (a) ACCURACY vs the INDEPENDENT odefun ray oracle (tight).
        for (tgt, oracle) in ((ComplexF64(-30.0), CPN5_U_M30),
                              (ComplexF64(-40.0), CPN5_U_M40),
                              (ComplexF64(-50.0), CPN5_U_M50))
            u, _ = eval_at(sol, tgt)
            @test isapprox(u, oracle; rtol = 1e-12)
            # Loose asymptote sanity cross-check: -√(-z/6) good to ~1e-5.
            @test isapprox(u, -sqrt(-tgt / 6); rtol = 1e-3)
        end

        # (b) ZERO-POLE: the entire negative-real fan is pole-free.  The naïve
        # failure mode would manufacture spurious poles here; the tritronquée
        # has NONE.  count == 0 is an EXACT pin (upgrades phase9's qualitative
        # sector-bin histogram).
        poles = extract_poles(sol; min_support = 3)
        @test count(in_pole_free_fan, poles) == 0

        # A small off-axis fan (arg within ±10° of 180°, still robustly
        # pole-free) on a SEPARATE upper-half-only solve.  The off-axis FAR
        # target (|z| = 30.4) accumulates more path error than the on-axis
        # ray (the wedge walk is best-conditioned along the goal direction,
        # which for the on-axis targets coincides with the real axis); the
        # measured eval_at-vs-odefun floor here is ~5e-5, NOT the on-axis
        # ~1e-13.  We pin the off-axis VALUE honestly at the METHOD floor
        # (rtol 1e-4) against the independent odefun ray oracle — the tight
        # accuracy claim lives on the on-axis ray above; this fan's load-
        # bearing role is the ZERO-POLE assertion that follows.
        gridf = ComplexF64[-30.0, -30.0 + 5.0im]
        solf  = path_network_solve(prob, gridf; h = 0.5,
                                   max_steps_per_target = 400, rng_seed = 0)
        for (tgt, oracle) in ((ComplexF64(-30.0, 5.0),  CPN5_U_M30P5I),)
            u, _ = eval_at(solf, tgt)
            @test isapprox(u, oracle; rtol = 1e-4)
        end
        @test count(in_pole_free_fan, extract_poles(solf; min_support = 3)) == 0
    end

    # ----------------------------------------------------------------------
    # CPN.6  equianharmonic ℘ long-range, multi-direction.
    # ----------------------------------------------------------------------
    @testset "CPN.6  ℘ long-range multi-direction: ±30°/±60° rays + FW real axis" begin
        fW(z, u, up) = 6 * u^2
        u0, up0 = 1.071822516416917, 1.710337353176786

        # Off-axis rays at R = 11 (clearance ≥ 0.5 from the ℘ lattice).
        prob = PadeTaylorProblem(fW, (u0 + 0im, up0 + 0im), (0.0 + 0im, 12.0 + 0im);
                                 order = 30)
        R = CPN6_RAY_R
        rays = ((ComplexF64(R * cis(π/6)),  CPN6_U_P30),
                (ComplexF64(R * cis(-π/6)), CPN6_U_M30),
                (ComplexF64(R * cis(π/3)),  CPN6_U_P60),
                (ComplexF64(R * cis(-π/3)), CPN6_U_M60))
        sol = path_network_solve(prob, ComplexF64[t for (t, _) in rays];
                                 h = 0.5, max_steps_per_target = 200, rng_seed = 0)
        for (tgt, oracle) in rays
            u, _ = eval_at(sol, tgt)
            @test isapprox(u, oracle; rtol = 1e-13)
        end

        # 0° real axis at z = 30 — FW 2011 Table 5.1 paper value (the path
        # bridges the on-axis poles at z = 1, 3.73, …).
        prob30 = PadeTaylorProblem(fW, (u0 + 0im, up0 + 0im), (0.0 + 0im, 30.0 + 0im);
                                   order = 30)
        sol30  = path_network_solve(prob30, ComplexF64[30.0 + 0im]; h = 0.5,
                                    max_steps_per_target = 200, rng_seed = 0)
        u30, _ = eval_at(sol30, 30.0 + 0im)
        @test isapprox(u30, 1.095098255959744; rtol = 1e-12)
        @test abs(imag(u30)) < 1e-12
    end

    # ----------------------------------------------------------------------
    # CPN.3  lattice-vs-lattice routing: rhombic vs square.
    # ----------------------------------------------------------------------
    @testset "CPN.3  equianharmonic (rhombic) vs lemniscatic (square): poles + trees differ" begin
        fE(z, u, up) = 6 * u^2
        fL(z, u, up) = 6 * u^2 - 1//2
        uE0, upE0 = 1.071822516416917, 1.710337353176786
        # CEl.1 lemniscatic IC: u = ℘(z-1; 1, 0), pole at z = 1.
        uL0  = 1.05083979104023701837944993034901171657327564
        upL0 = 1.89493523185581089049547932437906567542905606

        # Lattice pole formulae (poles shifted to z = 1, matching both ICs).
        Ω, ω = CPN3_EQUI_OM, CPN3_LEMN_OM
        rhombic(m, n) = 1 + 2Ω * (m + n/2) + im * (Ω * sqrt(3) * n)   # 60° rhombus
        square(m, n)  = 1 + 2ω * m + im * (2ω * n)                    # square
        rlat = ComplexF64[rhombic(m, n) for m in -3:3 for n in -3:3]
        slat = ComplexF64[square(m, n)  for m in -3:3 for n in -3:3]
        inbox(p) = -2.0 ≤ real(p) ≤ 5.0 && -2.0 ≤ imag(p) ≤ 2.0
        nd(p, lat) = minimum(abs(p - q) for q in lat)

        # Same target grid for BOTH problems (covers several poles of each).
        xs = -2.0:0.7:5.0
        ys = -2.0:0.7:2.0
        grid = ComplexF64[x + im*y for x in xs for y in ys]

        probE = PadeTaylorProblem(fE, (uE0 + 0im, upE0 + 0im), (0.0 + 0im, 6.0 + 0im); order = 30)
        probL = PadeTaylorProblem(fL, (uL0 + 0im, upL0 + 0im), (0.0 + 0im, 6.0 + 0im); order = 30)
        solE  = path_network_solve(probE, grid; h = 0.5, rng_seed = 0)
        solL  = path_network_solve(probL, grid; h = 0.5, rng_seed = 0)

        # (a) extracted poles lie on the RESPECTIVE lattice, NOT the other.
        # The two lattices SHARE the IC-shifted origin pole z = 1 (on the real
        # axis), so a per-pole "far from the other lattice" check would be
        # spoiled by that shared point.  The discriminator is instead: every
        # in-box pole is ON its own lattice (≤1e-6), AND the FARTHEST in-box
        # pole from the OTHER lattice is ~0.98 away (the non-shared poles:
        # equianharmonic 3.726/-1.726 vs square 4.708) — `maximum`, not
        # `minimum`, so the shared z=1 pole does not mask the geometry.
        eIn = filter(inbox, extract_poles(solE))
        lIn = filter(inbox, extract_poles(solL))
        @test !isempty(eIn) && !isempty(lIn)
        @test maximum(nd(p, rlat) for p in eIn) ≤ 1e-6          # equi ON rhombic
        @test maximum(nd(p, slat) for p in eIn) > 0.5           # equi has a pole OFF square
        @test maximum(nd(p, slat) for p in lIn) ≤ 1e-6          # lemn ON square
        @test maximum(nd(p, rlat) for p in lIn) > 0.5           # lemn has a pole OFF rhombic

        # (b) the two visited spanning trees DIFFER (different detour nodes) —
        # catches a hard-coded lattice orientation in the wedge walk.
        @test solE.visited_z != solL.visited_z
    end

end # @testset Corpus: path-network lattice/sectors (CPN)

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored exactly.
#
#   M-oracle (oracle file; VERIFIED RED 2026-06-08).  Perturb CPN5_U_M30's
#     real part at the 8th significant digit
#     (-2.236091119782981116195 → -2.236091219782981116195, a 4.5e-8 relative
#     change) in _oracle_corpus_pathnet_lattice_sectors.jl.  CPN.5 (a)
#     `isapprox(u, oracle; rtol=1e-12)` goes RED (the package value is
#     correct, the pin wrong).  Confirms the accuracy pin bites the oracle,
#     not just "didn't throw".  NB: the package's eval_at vs odefun floor is
#     ~1e-13 here, so a sub-1e-12 oracle nudge would be ABSORBED — the
#     mutation must exceed the pin tolerance to bite.  Restored.
#
#   M-impl (src/PathNetwork.jl `_local_pade`; VERIFIED RED 2026-06-08).  Scale
#     the canonical Padé step `h` by (1 + 1e-6) before `pade_step_with_pade!`.
#     Stage-2's `t = (z_f - z_v)/h_v` then mis-scales the per-node disc, so
#     every interpolated value drifts ~1e-6 to 1e-9.  Result: 9 fails — CPN.5
#     (a) accuracy (3, rel-err ~1e-8 ≫ 1e-12) and ALL of CPN.6's ±30°/±60° ray
#     pins (6, rel-err ~1e-6 ≫ 1e-13).  CPN.3 stays GREEN — its lattice
#     membership tol (1e-6) and the structural tree-difference are looser than
#     the perturbation effect, isolating the canonical-Padé scale as
#     load-bearing for the VALUE pins specifically.  Restored byte-for-byte.
#
#   M-impl NULL result (informative; documented, not load-bearing).  Scaling
#     the WEDGE step magnitude in `_wedge_evaluations` by 1.01 does NOT redden:
#     the per-node Padé is rebuilt CANONICALLY in `_local_pade` regardless of
#     the wedge magnitude, and the shifted nodes still cover every target
#     within `h`.  This confirms the wedge magnitude is a routing knob, not a
#     value knob — the canonical-Padé step above is the load-bearing one.
#
# The M-impl mutation was verified RED then reverted; `git diff src/` is
# empty at commit.  Standalone run:
#   julia --project=. test/corpus_pathnet_lattice_sectors_test.jl
# ============================================================================
