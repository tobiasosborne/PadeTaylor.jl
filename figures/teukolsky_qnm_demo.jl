# figures/teukolsky_qnm_demo.jl
#
# ## The headline Heun figure: a Schwarzschild quasinormal mode, end to end
#
# This is the named deliverable of the P1 Heun epic (Stage 4c, bead
# padetaylor-4lh): a single two-panel picture that tells the whole QNM story
# with PadeTaylor.jl machinery and nothing else.  A black hole, struck, rings
# down at a discrete set of complex frequencies ω = ω_R − iω_I — the
# quasinormal modes (QNMs).  In the e^{−iωt} convention Im ω < 0 means the
# ringing decays.  Each QNM is simultaneously *purely ingoing at the horizon*
# and *purely outgoing at infinity*; that double boundary condition is an
# eigenvalue problem whose eigenvalues are the ω_n and whose eigenfunctions
# are the radial profiles R(r).  The figure shows both halves of that pair.
#
# ### Ground truth (cite, do not paraphrase)
#
# The mathematics is the confluent-Heun reduction of the Regge–Wheeler /
# Teukolsky radial equation.  Two passages of `docs/teukolsky_heun_mapping.md`
# are load-bearing and were cross-verified against three independent sources
# (BHPT `QuasiNormalModes.m`, arXiv:2509.07235, arXiv:2604.18680) on
# 2026-06-19:
#
#   - §3.2.1 (`docs/teukolsky_heun_mapping.md:209-265`) — the VALIDATED
#     Leaver three-term recurrence and series ansatz in M=1 units (horizon
#     r₊ = 2M = 2).  This is the single source of truth for both the
#     continued-fraction eigenvalue condition and the eigenfunction
#     coefficients; we reuse `examples/teukolsky_qnm.jl` rather than
#     re-deriving it, so this figure inherits its validation (the s=0, ℓ=2,
#     n=0 root matches Leaver 1985's Mω = 0.483642 − 0.096766i to 7.5e-6).
#
#   - §6.2 (`docs/teukolsky_heun_mapping.md:370-422`) — Fiziev 2009's HeunC
#     parameter identification (α, β, γ, δ, η) for this exact test point, in
#     2M=1 units.  We use these for the BEST-EFFORT `heun_c` cross-check only.
#
# ### What the two panels show
#
# **Panel A — the eigenvalue as a continued-fraction zero.**  Leaver (1985)
# turned "the recessive solution of the recurrence exists" (Pincherle) into a
# continued fraction CF(ω) whose complex zeros are the QNM frequencies
# (`docs/teukolsky_heun_mapping.md:246-260, 314-335`).  We heatmap
# log₁₀|CF(ω)| over a window of the complex-ω plane around the fundamental;
# the dark well is the eigenvalue.  We mark the Newton-recovered root and the
# Leaver 1985 reference 0.483642 − 0.096766i.  This is pure validated
# arithmetic (`qnm_continued_fraction`, `find_qnm`) — zero convention risk.
#
# **Panel B — the eigenfunction.**  At the recovered ω we sum the VALIDATED
# Leaver series to reconstruct the radial profile
#
#     R(r) = e^{iωr} (r − 2M)^{−2iωM} r^{4iωM} · Σ_{n≥0} aₙ [1 − 2M/r]ⁿ,
#
# with a₀ = 1, the seed β₀a₀ + α₀a₁ = 0, and a_{n+1} = −(βₙaₙ + γₙa_{n−1})/αₙ
# (n ≥ 1) from `leaver_recurrence` (`docs/teukolsky_heun_mapping.md:230-244`).
# This is the confluent-Heun function rendered from validated data.  We plot
# Re R, Im R and |R| over a finite window r ∈ [2, 45].  QNM eigenfunctions
# *diverge* as r → ∞ (the outgoing factor e^{iωr} with Im ω < 0 grows like
# e^{|Im ω| r}) — that is physically correct, not a bug, so we plot a sensible
# finite window and say so on the axis.
#
# ### Why fail-loud assertions (Rule 1 / Rule 5)
#
# A figure that "ran" proves nothing.  Before any plotting we assert the
# NUMERICAL SANITY of every quantity the picture depends on and print each
# check: the recovered ω matches Leaver to < 1e-5, |CF(root)| < 1e-8 (we are
# genuinely *at* the eigenvalue), and the series coefficients aₙ decay (we
# sit on the minimal/recessive solution, not the runaway partner).  Any
# breakdown throws with a suggestion rather than drawing a truthful-looking
# lie.
#
# ### heun_c cross-check (BEST-EFFORT, time-boxed)
#
# PadeTaylor's `heun_c(q, α, γ, δ, ε; z)` (DLMF §31.12 convention, see
# `src/Heun.jl`) ought to render the same eigenfunction once the Fiziev→DLMF
# parameter map and the z↔r variable map are pinned down.  We make a short
# honest attempt; if it agrees with the Leaver series to plotting accuracy we
# overlay it dashed and print AGREE, otherwise we omit the overlay, plot only
# the validated Leaver profile, and print DEFERRED with a one-line reason —
# we do NOT chase a convention mismatch (that becomes a follow-on bead).

using PadeTaylor                      # examples-env dev-dep; gives heun_c
using CairoMakie
using Printf

# Reuse the VALIDATED single source of truth.  Its `abspath(PROGRAM_FILE) ==
# @__FILE__` guard means `include` brings in the functions WITHOUT running
# its own script block (Rule 3: verified against examples/teukolsky_qnm.jl).
include(joinpath(@__DIR__, "..", "examples", "teukolsky_qnm.jl"))

const OUTPNG = joinpath(@__DIR__, "output", "teukolsky_qnm_demo.png")
mkpath(dirname(OUTPNG))

# Physical / numerical constants (M=1 geometric units; horizon r₊ = 2M = 2).
const M       = 1.0
const RPLUS   = 2.0 * M
const S       = 0           # scalar perturbation (Regge–Wheeler)
const ELL     = 2           # multipole
const NDEPTH  = 100         # CF truncation depth & series truncation order
const SEED    = 0.4836 - 0.0968im
const REF_MW  = 0.483642 - 0.096766im      # Leaver 1985 (M=1 units)

println("="^72)
println("Schwarzschild QNM end-to-end: continued-fraction zero + eigenfunction")
println("  s = $S, ℓ = $ELL, n = 0   (M=1 units, horizon r₊ = $RPLUS)")
println("="^72)

# ----------------------------------------------------------------------
# Step 1 — recover the eigenvalue ω (NUMERICAL SANITY, Rule 5)
# ----------------------------------------------------------------------
ω = find_qnm(complex(SEED); s = S, ℓ = ELL, N = NDEPTH)
cf_at_root = qnm_continued_fraction(ω; s = S, ℓ = ELL, N = NDEPTH)
err_vs_leaver = abs(ω - REF_MW)

@printf("\n[SANITY 1] recovered  Mω        = % .9f %+ .9fi\n", real(ω), imag(ω))
@printf("[SANITY 1] Leaver 1985 reference = % .9f %+ .9fi\n",
        real(REF_MW), imag(REF_MW))
@printf("[SANITY 1] |Mω − Leaver|         = %.3e   (require < 1e-5)\n",
        err_vs_leaver)
@printf("[SANITY 2] |CF(root)|            = %.3e   (require < 1e-8)\n",
        abs(cf_at_root))

err_vs_leaver < 1e-5 || throw(ErrorException(
    "teukolsky_qnm_demo: recovered Mω=$ω differs from the Leaver 1985 " *
    "reference $REF_MW by $(err_vs_leaver) ≥ 1e-5. The continued-fraction " *
    "root-find did not land on the validated eigenvalue. Suggestion: check " *
    "the seed SEED=$SEED and depth NDEPTH=$NDEPTH against " *
    "docs/teukolsky_heun_mapping.md:252-254."))
abs(cf_at_root) < 1e-8 || throw(ErrorException(
    "teukolsky_qnm_demo: |CF(ω*)| = $(abs(cf_at_root)) ≥ 1e-8 at the " *
    "recovered frequency — we are not actually sitting on the " *
    "minimal-solution zero. Suggestion: raise NDEPTH or tighten find_qnm tol."))

# ----------------------------------------------------------------------
# Step 2 — Leaver series coefficients aₙ at the recovered ω
# ----------------------------------------------------------------------
# a₀ = 1; seed β₀a₀ + α₀a₁ = 0 ⟹ a₁ = −β₀/α₀;
# n ≥ 1: α_n a_{n+1} + β_n a_n + γ_n a_{n−1} = 0
#        ⟹ a_{n+1} = −(β_n a_n + γ_n a_{n−1})/α_n
# (`docs/teukolsky_heun_mapping.md:230-244`).
function leaver_coeffs(ω::Complex; s::Integer, ℓ::Integer, N::Integer)
    a = Vector{ComplexF64}(undef, N + 1)
    a[1] = 1.0 + 0im                                   # a₀
    α0, β0, _ = leaver_recurrence(0, ω; s = s, ℓ = ℓ)
    α0 == 0 && throw(ErrorException(
        "leaver_coeffs: α₀ = 0 at ω=$ω — recurrence seed β₀a₀+α₀a₁=0 is " *
        "degenerate (n=0 grazes a recurrence singularity). Suggestion: " *
        "verify ω is the QNM, not a spurious root."))
    a[2] = -β0 / α0                                     # a₁ from the seed
    for n in 1:(N - 1)
        αn, βn, γn = leaver_recurrence(n, ω; s = s, ℓ = ℓ)
        αn == 0 && throw(ErrorException(
            "leaver_coeffs: αₙ = 0 at n=$n, ω=$ω — recurrence breakdown " *
            "(αₙ=(n+1)(n+1−4iω) vanishes). Suggestion: this should not " *
            "happen for a genuine QNM; check ω and depth N=$N."))
        a[n + 2] = -(βn * a[n + 1] + γn * a[n]) / αn
    end
    return a
end

a = leaver_coeffs(ω; s = S, ℓ = ELL, N = NDEPTH)

# SANITY 3: minimal-solution coefficients must DECAY.  We sit on the
# recessive solution (the outgoing BC); if aₙ grew we would have excited the
# runaway partner and the series — hence R(r) — would be garbage.
abs_a = abs.(a)
finite_all = all(isfinite, abs_a)
tail_ratio = abs_a[end] / abs_a[1]
# Peak then decay: the ratio of the last to the largest coefficient.
peak_idx = argmax(abs_a)
last_over_peak = abs_a[end] / abs_a[peak_idx]
@printf("[SANITY 3] coeff decay: |a₀|=%.3e  |a_peak|=%.3e (n=%d)  |a_N|=%.3e\n",
        abs_a[1], abs_a[peak_idx], peak_idx - 1, abs_a[end])
@printf("[SANITY 3]            |a_N|/|a_peak| = %.3e  (require ≪ 1: minimal soln)\n",
        last_over_peak)
finite_all || throw(ErrorException(
    "teukolsky_qnm_demo: non-finite Leaver coefficient encountered — the " *
    "recurrence blew up. Suggestion: the recovered ω is not a genuine QNM, " *
    "or NDEPTH=$NDEPTH is too deep for Float64. Reduce N or check ω."))
last_over_peak < 1e-2 || throw(ErrorException(
    "teukolsky_qnm_demo: Leaver coefficients are NOT decaying " *
    "(|a_N|/|a_peak| = $last_over_peak ≥ 1e-2). We have excited the " *
    "growing (dominant) solution rather than the minimal/recessive one — " *
    "the reconstructed eigenfunction would be meaningless. Suggestion: " *
    "verify ω sits on the continued-fraction zero (SANITY 2)."))

# ----------------------------------------------------------------------
# Step 3 — reconstruct the radial eigenfunction R(r)
# ----------------------------------------------------------------------
# R(r) = e^{iωr} (r − 2M)^{−2iωM} r^{4iωM} · Σ aₙ [1 − 2M/r]ⁿ   (M=1 units).
# The compactified series variable x = 1 − 2M/r ∈ [0, 1) on r ∈ [2M, ∞);
# x = 0 at the horizon (only the a₀ term survives ⟹ R∝prefactor there) and
# x → 1 as r → ∞.  We sum to the same depth NDEPTH used for the eigenvalue.
function leaver_R(r::Real, ω::Complex, a::AbstractVector; M::Real = 1.0)
    rplus = 2M
    r > rplus || throw(ErrorException(
        "leaver_R: r=$r ≤ horizon r₊=$rplus — outside the physical domain " *
        "and the prefactor (r−2M)^{−2iωM} has a branch point at r₊. " *
        "Suggestion: evaluate strictly outside the horizon."))
    x = 1 - rplus / r
    # Horner sum of Σ aₙ xⁿ.
    s = a[end]
    @inbounds for k in (length(a) - 1):-1:1
        s = s * x + a[k]
    end
    prefac = exp(im * ω * r) * (r - rplus)^(-2im * ω * M) * r^(4im * ω * M)
    return prefac * s
end

# r-grid: just outside the horizon (avoid the branch point exactly at r₊) to
# a finite outer radius.  We deliberately stop well short of "infinity": the
# QNM eigenfunction grows there (see header), so a finite window is the
# honest thing to plot.
const R_OUTER = 45.0
rs = range(RPLUS + 1e-3, R_OUTER; length = 1200)
Rvals = [leaver_R(r, ω, a; M = M) for r in rs]

all(isfinite, Rvals) || throw(ErrorException(
    "teukolsky_qnm_demo: reconstructed R(r) is non-finite on the plot " *
    "window r ∈ [$(RPLUS), $R_OUTER]. Suggestion: the series did not " *
    "converge — re-check coefficient decay (SANITY 3)."))

@printf("[SANITY 4] R(r) reconstructed on r ∈ [%.3f, %.1f], %d points, all finite\n",
        RPLUS + 1e-3, R_OUTER, length(rs))
@printf("[SANITY 4]   |R| near horizon (r=%.3f) = %.3e ;  |R| at r=%.1f = %.3e\n",
        rs[1], abs(Rvals[1]), R_OUTER, abs(Rvals[end]))

# ----------------------------------------------------------------------
# Step 4 — BEST-EFFORT heun_c cross-check (time-boxed; do NOT rabbit-hole)
# ----------------------------------------------------------------------
# Fiziev 2009 §6.2 gives the HeunC params in 2M=1 units, naming (α,β,γ,δ,η)
# in his CHE form (mapping_doc:370-417).  PadeTaylor's heun_c uses the DLMF
# §31.12 form heun_c(q, α, γ, δ, ε; z) — a DIFFERENT parameterisation.  The
# best-effort attempt: build the DLMF-convention parameters from the Fiziev
# table and the Frobenius-at-0 series both expand the *same* local solution
# at z=0 normalised to 1, so we compare the package's heun_c against the
# Leaver-series eigenfunction divided by its own horizon-normalised prefactor
# on a short window near the horizon (r ∈ [2, 10]).  Two attempts, then stop.
#
# Fiziev (2M=1, ω₂ = 2Mω):
#   α_F = 2iω₂, β_F = 0, γ_F = 2iω₂, δ_F = 2ω₂², η_F = −ℓ(ℓ+1) + 2ω₂²
# DLMF heun_c(q, α, γ, δ, ε; z) ↔ the confluent-Heun coefficient set; the
# accepted Fiziev→DLMF map (Motygin 2018 / DLMF 31.12 vs Fiziev I.7) is the
# subject of this cross-check.
heun_c_verdict = "DEFERRED (not attempted)"
heun_c_overlay = nothing            # (r_window, values) if it AGREES
try
    ω2 = 2 * M * ω                              # 2M=1 units
    αF = 2im * ω2
    γF = 2im * ω2
    δF = 2 * ω2^2
    ηF = -ELL * (ELL + 1) + 2 * ω2^2
    # Fiziev I.7 → DLMF 31.12 form.  DLMF: u'' + (γ_D/z + δ_D/(z−1) + ε_D)u'
    # + (α_D z − q_D)/(z(z−1)) u = 0.  Fiziev I.7:
    # H'' + (α_F + (β_F+1)/z + (γ_F+1)/(z−1))H' + (μ/z + ν/(z−1))H = 0.
    # Matching coefficient-by-coefficient gives the candidate DLMF set:
    γ_D = γF + 1                                 # exponent at z=0
    δ_D = δF + 1                                 # exponent at z=1  (δ_F is Fiziev δ)
    ε_D = αF                                     # irregular-point scale
    μ  = 0.5 * (αF - 0 - γF + 0 - 0) - ηF        # β_F = 0 ⟹ simplified μ
    ν  = 0.5 * (αF + 0 + γF + 0 + 0) + δF + ηF
    # DLMF q_D, α_D from (μ, ν): α_D z − q_D over z(z−1) must equal μ/z+ν/(z−1)
    # ⟹ μ(z−1) + νz = α_D z − q_D ⟹ α_D = μ+ν, q_D = μ.
    α_D = μ + ν
    q_D = μ
    # The HeunC local solution at z=0 is normalised to 1; the Leaver series
    # Σ aₙ xⁿ with x = 1 − 2M/r is the *horizon*-normalised piece.  These use
    # different compactified variables (DLMF z vs Leaver x), so a clean
    # comparison needs the z↔x map pinned.  Attempt the most natural guess
    # z = x and compare the bare series Σ aₙ xⁿ (the part heun_c should equal,
    # up to convention) on r ∈ [2, 10].
    rc = range(RPLUS + 0.05, 10.0; length = 25)
    leaver_series_only = ComplexF64[]
    heun_c_series = ComplexF64[]
    for r in rc
        x = 1 - RPLUS / r
        sser = a[end]
        for k in (length(a) - 1):-1:1
            sser = sser * x + a[k]
        end
        push!(leaver_series_only, sser)
        push!(heun_c_series, heun_c(q_D, α_D, γ_D, δ_D, ε_D; z = complex(x)))
    end
    # Compare up to an overall complex constant (both normalised at their own
    # origin): use the relative residual of the best-fit scale.
    scale = sum(conj.(heun_c_series) .* leaver_series_only) /
            sum(abs2, heun_c_series)
    resid = maximum(abs.(leaver_series_only .- scale .* heun_c_series)) /
            maximum(abs.(leaver_series_only))
    print(@sprintf("[heun_c]  candidate DLMF params: q=%.4f%+.4fi α=%.4f%+.4fi γ=%.4f%+.4fi δ=%.4f%+.4fi ε=%.4f%+.4fi\n",
            real(q_D), imag(q_D), real(α_D), imag(α_D), real(γ_D), imag(γ_D),
            real(δ_D), imag(δ_D), real(ε_D), imag(ε_D)))
    @printf("[heun_c]  best-fit-scale residual on r∈[2,10] = %.3e\n", resid)
    if isfinite(resid) && resid < 5e-2
        global heun_c_verdict = "AGREE"
        # Reconstruct the full R(r) overlay from heun_c × the same prefactor.
        rov = range(RPLUS + 0.05, 10.0; length = 200)
        ov = ComplexF64[]
        for r in rov
            x = 1 - RPLUS / r
            prefac = exp(im * ω * r) * (r - RPLUS)^(-2im * ω * M) *
                     r^(4im * ω * M)
            push!(ov, prefac * scale * heun_c(q_D, α_D, γ_D, δ_D, ε_D;
                                              z = complex(x)))
        end
        global heun_c_overlay = (collect(rov), ov)
    else
        global heun_c_verdict =
            "DEFERRED (mismatch: Fiziev→DLMF param map + z↔x variable map " *
            "unverified; residual=$(round(resid; sigdigits=3)))"
    end
catch e
    global heun_c_verdict =
        "DEFERRED (mismatch: heun_c threw — " *
        first(split(sprint(showerror, e), '\n')) * ")"
end
println("[heun_c cross-check] ", heun_c_verdict)

# ----------------------------------------------------------------------
# Step 5 — Panel A data: log₁₀|CF(ω)| over the complex-ω plane
# ----------------------------------------------------------------------
const RE_MIN, RE_MAX = 0.30, 0.62
const IM_MIN, IM_MAX = -0.17, -0.01
const NRE, NIM = 150, 150
re_axis = range(RE_MIN, RE_MAX; length = NRE)
im_axis = range(IM_MIN, IM_MAX; length = NIM)
cf_logmag = Matrix{Float64}(undef, NRE, NIM)
@printf("\n[panel A] evaluating log₁₀|CF(ω)| on a %d×%d ω-grid ...\n", NRE, NIM)
for (jw, ωi) in enumerate(im_axis), (iw, ωr) in enumerate(re_axis)
    val = qnm_continued_fraction(complex(ωr, ωi); s = S, ℓ = ELL, N = NDEPTH)
    cf_logmag[iw, jw] = isfinite(val) ? log10(abs(val) + 1e-30) : NaN
end
@printf("[panel A] log₁₀|CF| range: [%.2f, %.2f]\n",
        minimum(filter(isfinite, cf_logmag)),
        maximum(filter(isfinite, cf_logmag)))

# ----------------------------------------------------------------------
# Step 6 — render the two-panel figure
# ----------------------------------------------------------------------
fig = Figure(size = (1560, 680), fontsize = 14)
Label(fig[0, 1:3],
      "Schwarzschild quasinormal mode (scalar, ℓ=2, n=0) — PadeTaylor.jl";
      fontsize = 23, padding = (0, 0, 6, 0))

# --- Panel A: continued-fraction eigenvalue map ---
axA = Axis(fig[1, 1];
           title = "A.  log₁₀|CF(ω)| — the QNM as a continued-fraction zero",
           xlabel = "Re(Mω)", ylabel = "Im(Mω)",
           limits = (RE_MIN, RE_MAX, IM_MIN, IM_MAX))
hmA = heatmap!(axA, re_axis, im_axis, cf_logmag;
               colormap = :viridis)
Colorbar(fig[1, 2], hmA; label = "log₁₀|CF(ω)|", width = 12)
# Recovered root (white cross) and Leaver reference (orange circle).
scatter!(axA, [real(ω)], [imag(ω)];
         marker = :xcross, markersize = 20, color = :white,
         strokecolor = :black, strokewidth = 1.5)
scatter!(axA, [real(REF_MW)], [imag(REF_MW)];
         marker = :circle, markersize = 13, color = (:orange, 0.0),
         strokecolor = :orange, strokewidth = 2.5)
text!(axA, real(ω), imag(ω);
      text = @sprintf("recovered\nMω = %.6f%+.6fi", real(ω), imag(ω)),
      offset = (12, -34), color = :white, fontsize = 12,
      strokecolor = :black, strokewidth = 1)
text!(axA, real(REF_MW), imag(REF_MW);
      text = @sprintf("Leaver 1985\n%.6f%+.6fi", real(REF_MW), imag(REF_MW)),
      offset = (12, 14), color = :orange, fontsize = 12,
      strokecolor = :black, strokewidth = 1)

# --- Panel B: radial eigenfunction R(r) ---
axB = Axis(fig[1, 3];
           title = "B.  radial QNM eigenfunction R(r) from the Leaver series",
           xlabel = "r / M   (horizon r₊ = 2M = 2)",
           ylabel = "R(r)",
           limits = (RPLUS, R_OUTER, nothing, nothing))
reR = real.(Rvals); imR = imag.(Rvals); absR = abs.(Rvals)
lines!(axB, rs, reR; color = :dodgerblue, linewidth = 2, label = "Re R(r)")
lines!(axB, rs, imR; color = :crimson,    linewidth = 2, label = "Im R(r)")
lines!(axB, rs, absR; color = :black, linewidth = 1.6, linestyle = :dot,
       label = "|R(r)|")
lines!(axB, rs, -absR; color = :black, linewidth = 1.0, linestyle = :dot)
vlines!(axB, [RPLUS]; color = (:gray, 0.6), linewidth = 1.5, linestyle = :dash)
text!(axB, RPLUS, 0.0; text = "horizon", rotation = π/2, align = (:center, :bottom),
      offset = (-4, 0), color = :gray35, fontsize = 11)
if heun_c_overlay !== nothing
    rov, ov = heun_c_overlay
    lines!(axB, rov, real.(ov); color = :darkgreen, linewidth = 2.4,
           linestyle = :dash, label = "Re R via heun_c (cross-check)")
end
axislegend(axB; position = :lt, framevisible = true, labelsize = 11)

Label(fig[2, 1:3],
      "Panel A: log₁₀|CF(ω)| (Leaver continued fraction, N=$NDEPTH); " *
      "× recovered root, ○ Leaver 1985 reference — they coincide.   " *
      "Panel B: R(r)=e^{iωr}(r−2M)^{−2iωM}r^{4iωM}·Σaₙ[1−2M/r]ⁿ at the " *
      "recovered ω; the QNM grows as r→∞ (Im ω<0), so a finite window " *
      "r≤$(Int(R_OUTER)) is shown.   heun_c cross-check: $heun_c_verdict";
      fontsize = 10.5, color = :gray30, padding = (0, 0, 0, 4),
      justification = :left, lineheight = 1.2)

colsize!(fig.layout, 1, Relative(0.42))
colsize!(fig.layout, 2, Relative(0.04))
colsize!(fig.layout, 3, Relative(0.54))

save(OUTPNG, fig)
@printf("\nSaved: %s\n", OUTPNG)
println("="^72)
