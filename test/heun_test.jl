# test/heun_test.jl — Heun module golden-master test suite (Stage 5).
#
# Pins `heun_g` and `heun_c` against the wolframscript oracle at
# `external/probes/heun-oracle/oracles.txt` (42 records, 30-digit) plus
# the Motygin MATLAB cross-validation at
# `external/probes/heun-oracle/CROSSVAL.md`.
#
# Per CLAUDE.md Rules 1+4 and user mandate 2026-05-18:
#   - Each individual test asserts to a precise numerical tolerance.
#   - Each test enforces a 2-second wall-time gate.  Total file < 60 s.
#   - HeunC tests at z ≥ 1 deliberately use looser tolerance / off-axis
#     anchors per CROSSVAL.md verdict: "HeunC for z >= 1 — Motygin and
#     Mathematica are on different analytic sheets … cannot be pinned
#     by Motygin without connection-formula work."  We pin against
#     `HeunC[..., z+0.001i]` (Mathematica's upper-limit value) — the
#     same sheet our `enforce_real_axis_symmetry=true` walker takes.

using Test
using PadeTaylor: heun_g, heun_c

# ---------------------------------------------------------------------------
# Wall-time gate (CLAUDE.md user mandate 2026-05-18).
# Each assertion must complete in under WALL_TIME_GATE seconds; we wrap
# the value computation in @elapsed and assert.
# ---------------------------------------------------------------------------
const WALL_TIME_GATE = 2.0   # seconds, hard ceiling per timed call

# JIT warmup: trigger compilation for both code paths (Frobenius-only and
# path-network) BEFORE the timed tests start.  Without this, the first
# timed call includes ~2s of LLVM compilation and trips the gate.
let
    _ = heun_g(2, 1, 1, 2, 3, 4; z = 0.5 + 0im)        # walks
    _ = heun_g(2, 1, 1, 2, 3, 4; z = 0.01 + 0im)       # Frobenius-only
    _ = heun_c(1, 1, 2, 3, -1; z = 0.5 + 0im)
    _ = heun_c(1, 1, 2, 3, -1; z = 0.01 + 0im)
end

macro pinned(call, oracle, rtol)
    quote
        local v, dt
        dt = @elapsed (v = $(esc(call)))
        @test dt < WALL_TIME_GATE
        @test abs(v - $(esc(oracle))) / abs($(esc(oracle))) ≤ $(esc(rtol))
    end
end

# ---------------------------------------------------------------------------
# HG.1 — HeunG Frobenius regime (|z| ≤ 0.05 default → direct series)
# ---------------------------------------------------------------------------
@testset "HG.1 — HeunG Frobenius-only regime (|z| ≤ ε_start, no walk)" begin
    # At z=0 exactly, HeunG normalisation: HeunG(...; 0) = 1
    @test heun_g(2, 1, 1, 2, 3, 4; z = 0.0 + 0im) ≈ 1.0 atol = 1e-15

    # B1 degenerate identity: q=0, α=0 → HeunG ≡ 1 along positive real axis.
    # (This is a DLMF-derivable identity: from the recurrence with q=0, α=0
    # we have c_1 = 0 and P_j = 0 for j ≥ 2, so all higher c_j vanish.)
    for zt in (0.02+0im, 0.04+0im)
        @pinned heun_g(2, 0, 0, 2, 3, 4; z = zt) (1.0 + 0im) 1e-14
    end
end

# ---------------------------------------------------------------------------
# HG.2 — HeunG inside Frobenius disc via path-network (|z| < 1, walks)
# ---------------------------------------------------------------------------
@testset "HG.2 — HeunG inside disc via path-network walk" begin
    # Real-axis values (oracle pinned to 1e-13).
    @pinned heun_g(2, 1, 1, 2, 3, 4; z = 0.1 + 0im) 1.01698183974548593693413697618 + 0im 1e-12
    @pinned heun_g(2, 1, 1, 2, 3, 4; z = 0.5 + 0im) 1.09086834192356524567686750497 + 0im 1e-12
    @pinned heun_g(2, 1, 1, 2, 3, 4; z = 0.9 + 0im) 1.08501231904897540061995738851 + 0im 1e-12
end

# ---------------------------------------------------------------------------
# HG.3 — HeunG past z = a (= 2) on real axis (Wolfram principal sheet)
# ---------------------------------------------------------------------------
@testset "HG.3 — HeunG past z = a via Schwarz-symmetric walk" begin
    # z=3 is past both singular points (z=1 and z=a=2).  Walker uses
    # enforce_real_axis_symmetry; result matches Mathematica's
    # positive-imag principal-sheet value.
    @pinned heun_g(2, 1, 1, 2, 3, 4; z = 3.0 + 0im) 2.50048593259480864400191470232 + 1.13019151087573739762935104864im 1e-8

    # z=1.9 between cuts (sheet matches).
    @pinned heun_g(2, 1, 1, 2, 3, 4; z = 1.9 + 0im) 1.62938380388991529201995833848 - 0.0237636055330802896121718283173im 1e-8
end

# ---------------------------------------------------------------------------
# HG.4 — HeunG off real axis (Schwarz reflection / direct walk)
# ---------------------------------------------------------------------------
@testset "HG.4 — HeunG off real axis" begin
    # Upper half plane, inside disc.
    @pinned heun_g(2, 1, 1, 2, 3, 4; z = 0.5 + 0.5im)   1.08330447245477887174278307207 + 0.101277513394027204304026632048im 1e-12
    # Upper half, past z=1.
    @pinned heun_g(2, 1, 1, 2, 3, 4; z = 1.5 + 0.25im)  1.39692708192180109660303463890 + 0.0887329767721584520116820132696im 1e-8
    # Lower half (Schwarz mirror).
    @pinned heun_g(2, 1, 1, 2, 3, 4; z = 2.5 - 0.5im)   1.80699599100006028078189143898 - 0.598905987721527227756086898592im 1e-12
end

# ---------------------------------------------------------------------------
# HG.5 — HeunG special-case identities (DLMF 31 degenerate parameters)
# ---------------------------------------------------------------------------
@testset "HG.5 — HeunG special-case identities" begin
    # B1: q=0, α=0 → HeunG ≡ 1 everywhere (DLMF-verifiable via recurrence).
    for zt in (0.04+0im, 0.25+0im, 0.5+0im, 0.75+0im)
        @pinned heun_g(2, 0, 0, 2, 3, 4; z = zt) 1.0 + 0im 1e-12
    end
end

# ---------------------------------------------------------------------------
# HC.1 — HeunC Frobenius regime
# ---------------------------------------------------------------------------
@testset "HC.1 — HeunC Frobenius-only regime" begin
    @test heun_c(1, 1, 2, 3, -1; z = 0.0 + 0im) ≈ 1.0 atol = 1e-15
    # Recurrence sanity: c_1 = -q/γ = -1/2; expansion at z=0.02 gives
    # 1 + c_1·z + O(z²) ≈ 1 - 0.5·0.02 = 0.99.  We don't have a high-precision
    # oracle here, so just check the leading-order expansion holds.
    v = heun_c(1, 1, 2, 3, -1; z = 0.02 + 0im)
    @test abs(real(v) - (1 - 0.5*0.02)) < 1e-3   # leading-order match
end

# ---------------------------------------------------------------------------
# HC.2 — HeunC inside disc via path-network walk (2-source pinned!)
# Per CROSSVAL.md: z in (0,1) machine-precision agreement Mathematica ↔ Motygin.
# ---------------------------------------------------------------------------
@testset "HC.2 — HeunC inside disc (2-source pinned)" begin
    @pinned heun_c(1, 1, 2, 3, -1; z = 0.1 + 0im) 0.947195554219883278485 + 0im 1e-12
    @pinned heun_c(1, 1, 2, 3, -1; z = 0.5 + 0im) 0.612124644911960719121 + 0im 1e-12
    @pinned heun_c(1, 1, 2, 3, -1; z = 0.9 + 0im) (-3.94427764664311075396 + 0im) 1e-8
end

# ---------------------------------------------------------------------------
# HC.3 — HeunC off real axis
# ---------------------------------------------------------------------------
@testset "HC.3 — HeunC off real axis" begin
    @pinned heun_c(1, 1, 2, 3, -1; z = 0.5 + 0.5im)  0.883915756159538923864831581540 - 0.360286128606480829377541615098im 1e-12
    @pinned heun_c(1, 1, 2, 3, -1; z = 1.5 + 0.25im) 0.784466562289524476596119188175 - 0.355829161520796657924173715451im 1e-10
end

# ---------------------------------------------------------------------------
# HC.4 — HeunC past z=1 (upper-limit sheet via Mathematica probe)
# Pinned against HeunC[..., z+0.001i] per CROSSVAL.md sheet caveat.
# ---------------------------------------------------------------------------
@testset "HC.4 — HeunC past z=1 (upper-limit sheet)" begin
    # Pinned upper-limit oracle from wolframscript probe in /tmp/heunc_sheet_check.wl:
    #   HeunC[1,1,2,3,-1; 3+0.001i] = 0.6572021676876615 - 0.3193993171834178*I
    # Our enforce_real_axis_symmetry walker delivers this upper-limit value.
    @pinned heun_c(1, 1, 2, 3, -1; z = 3.0 + 0im) (0.6572021676876615 - 0.3193993171834178im) 1e-3
end

# ---------------------------------------------------------------------------
# HX.1 — Existence-condition throws (CLAUDE.md Rule 1)
# ---------------------------------------------------------------------------
@testset "HX.1 — Existence-condition throws" begin
    # γ = 0 (non-positive integer): HeunG existence condition fails.
    @test_throws ArgumentError heun_g(2, 1, 1, 2, 0, 4; z = 0.5 + 0im)
    @test_throws ArgumentError heun_g(2, 1, 1, 2, -1, 4; z = 0.5 + 0im)
    # a = 1: singular points collide.
    @test_throws ArgumentError heun_g(1, 1, 1, 2, 3, 4; z = 0.5 + 0im)
    # γ = 0 for HeunC.
    @test_throws ArgumentError heun_c(1, 1, 0, 3, -1; z = 0.5 + 0im)
end

# ---------------------------------------------------------------------------
# HX.2 — Wall-time aggregate gate
# ---------------------------------------------------------------------------
@testset "HX.2 — Aggregate wall-time discipline" begin
    # Confirm one fresh evaluation completes well under the gate after
    # any precompile/warmup is done.
    t = @elapsed heun_g(2, 1, 1, 2, 3, 4; z = 2.0 + 0im)
    @test t < WALL_TIME_GATE
    t = @elapsed heun_c(1, 1, 2, 3, -1; z = 0.3 + 0.2im)
    @test t < WALL_TIME_GATE
end
