# verify.jl -- triangulate Phase 6 Problems oracle sources, then emit
# test/_oracle_problems.jl.
#
# Two purposes (matching the two-section structure of capture.wl):
#
#   (A) Pole-bridge demo (worklog 004; new Phase-6 v1 acceptance).
#       Cross-checks closed-form ≡ NDSolve ≡ mpmath.odefun for
#       z ∈ {0.5, 0.95} (3-source); closed-form sole-source for
#       z ∈ {1.05, 1.4} with optional NDSolve restart at z=1.05 →
#       1.4 if it converged.  BF-256 round-trip at z=1.05.
#
#   (B) FW 2011 Table 5.1 reference values (deferred to v2).
#       Closed-form ≡ paper-pinned reference (2-source) for
#       z ∈ {30, 10⁴, 28.261}, plus BF-256 round-trip at z=30.
#
# Sources:
#   oracles_wolfram.txt — Mathematica closed-form WeierstrassP +
#                         NDSolve at WorkingPrecision=50.
#   oracles_python.txt  — mpmath.odefun at 40 dps; valid for z < 1
#                         only (pole at z=1).
#   FW 2011 paper       — reference values pinned at lines 301, 372.
#
# Convention (mirrors Phase-5 padestepper-oracle/verify.jl): the
# pinned test/_oracle_problems.jl carries the Wolfram content only.
# The mpmath cross-check is the ground-truth provenance, asserted
# here, and not duplicated as test-readable values.  Tests assert
# against the Wolfram-pinned closed-form.
#
# Usage (from project root):
#   julia --project=. external/probes/problems-oracle/verify.jl
# Run AFTER:
#   wolframscript -file external/probes/problems-oracle/capture.wl
#   python3 external/probes/problems-oracle/capture.py

using Printf

const HERE = @__DIR__
const PROJECT_ROOT = abspath(joinpath(HERE, "..", "..", ".."))
const O_WL = joinpath(HERE, "oracles_wolfram.txt")
const O_PY = joinpath(HERE, "oracles_python.txt")
const PINNED = joinpath(PROJECT_ROOT, "test", "_oracle_problems.jl")

module _OWL end
module _OPY end

Base.include_string(_OWL, read(O_WL, String), O_WL)
Base.include_string(_OPY, read(O_PY, String), O_PY)

# ─── Cross-source assertions ──────────────────────────────────────────

function check_agree(label, vals...; rtol)
    a = vals[1]
    for (i, v) in enumerate(vals[2:end])
        rel = abs(a - v) / max(abs(a), abs(v), 1e-300)
        rel < rtol || error("$label: source 1 = $a, source $(i+1) = $v, rel = $rel > rtol = $rtol")
    end
    @printf("%-55s  %.16g  (max rel %.2e)\n", label, a,
            maximum([abs(a - v)/max(abs(a), abs(v), 1e-300) for v in vals[2:end]]))
end

function check2(label, a, b; rtol)
    rel = abs(a - b) / max(abs(a), abs(b), 1e-300)
    rel < rtol || error("$label  $a  vs  $b  rel=$rel > rtol=$rtol")
    @printf("%-55s  %.16g  (rel %.2e)\n", label, a, rel)
end

# ===================================================================
# (A) Pole-bridge demo cross-checks (worklog 004)
# ===================================================================

println("(A) Pole-bridge demo cross-checks:")

# ---- z = 0.5: closed-form ≡ NDSolve ≡ mpmath.odefun (3-source) ----
check_agree("u  at z=0.5  (3-source: ℘ ≡ NDSolve ≡ mpmath)",
            _OWL.u_at_0_5,
            _OWL.u_at_0_5_via_NDSolve,
            _OPY.u_at_0_5_py;
            rtol = 1e-13)
check_agree("u' at z=0.5  (3-source: ℘ ≡ NDSolve ≡ mpmath)",
            _OWL.up_at_0_5,
            _OWL.up_at_0_5_via_NDSolve,
            _OPY.up_at_0_5_py;
            rtol = 1e-13)

# ---- z = 0.95 (near pole): 3-source at looser rtol ----
# Phase-5 5.1.4 used 1e-11 for z=0.95 in the closed-form-vs-mpmath
# check; we replicate the same near-pole tolerance bound here.
check_agree("u  at z=0.95 (3-source: ℘ ≡ NDSolve ≡ mpmath)",
            _OWL.u_at_0_95,
            _OWL.u_at_0_95_via_NDSolve,
            _OPY.u_at_0_95_py;
            rtol = 1e-11)
check_agree("u' at z=0.95 (3-source: ℘ ≡ NDSolve ≡ mpmath)",
            _OWL.up_at_0_95,
            _OWL.up_at_0_95_via_NDSolve,
            _OPY.up_at_0_95_py;
            rtol = 1e-11)

# ---- z = 1.05: closed-form sole source past pole; magnitude sanity ----
let
    u105 = _OWL.u_at_1_05
    @assert isfinite(u105) "u_at_1_05 must be finite; got $u105"
    @assert 100 < u105 < 1000 "u_at_1_05 should be O(1/(z-1)^2) ≈ 400; got $u105"
    @printf("z=1.05 closed-form sole source (magnitude sanity)  u=%.4f  (expected ~400 = 1/0.05^2)\n", u105)
end

# ---- z = 1.4: closed-form primary; if NDSolve restart present, 2-source ----
if isdefined(_OWL, :u_at_1_4_via_NDSolve_restart)
    check_agree("u  at z=1.4  (2-source: ℘ ≡ NDSolve restart from 1.05)",
                _OWL.u_at_1_4,
                _OWL.u_at_1_4_via_NDSolve_restart;
                rtol = 1e-10)
    check_agree("u' at z=1.4  (2-source: ℘ ≡ NDSolve restart from 1.05)",
                _OWL.up_at_1_4,
                _OWL.up_at_1_4_via_NDSolve_restart;
                rtol = 1e-10)
else
    @printf("z=1.4 closed-form sole source (NDSolve restart did not converge or absent)\n")
end

# ---- BigFloat-256 round-trip at z=1.05 ----
let
    setprecision(BigFloat, 256) do
        u_bf  = parse(BigFloat, _OWL.u_at_1_05_80dps_str)
        up_bf = parse(BigFloat, _OWL.up_at_1_05_80dps_str)
        @assert isapprox(Float64(u_bf),  _OWL.u_at_1_05;  rtol = 1e-14) "BF round-trip on u_at_1_05"
        @assert isapprox(Float64(up_bf), _OWL.up_at_1_05; rtol = 1e-14) "BF round-trip on up_at_1_05"
        @printf("z=1.05 80-dps BF-256 round-trip ≡ closed-form Float64\n")
    end
end

println()

# ===================================================================
# (B) FW Table 5.1 paper-pinned reference values (deferred to v2)
# ===================================================================

println("(B) FW Table 5.1 paper-ref cross-checks (deferred to v2):")

let
    check2("u(30)     vs FW ref (FW2011_*.md:301)",
           _OWL.u_at_30, _OWL.u_at_30_FW_ref; rtol = 1e-14)
    check2("u(10^4)   vs FW ref (FW2011_*.md:301)",
           _OWL.u_at_10000, _OWL.u_at_10000_FW_ref; rtol = 1e-14)
    check2("u(28.261) vs FW ref (FW2011_*.md:372)",
           _OWL.u_at_28261, _OWL.u_at_28261_FW_ref; rtol = 1e-14)
end

# BF-256 round-trip at z=30 (existing).
let
    setprecision(BigFloat, 256) do
        u_bf  = parse(BigFloat, _OWL.u_at_30_80dps_str)
        up_bf = parse(BigFloat, _OWL.up_at_30_80dps_str)
        @assert isapprox(Float64(u_bf),  _OWL.u_at_30;  rtol = 1e-14) "BF round-trip on u_at_30"
        @assert isapprox(Float64(up_bf), _OWL.up_at_30; rtol = 1e-14) "BF round-trip on up_at_30"
        @printf("z=30 80-dps BF-256 round-trip ≡ closed-form Float64\n")
    end
end

println()

# ─── Emit canonical pinned file ──────────────────────────────────────
#
# Convention (mirrors Phase-5 padestepper-oracle/verify.jl): include
# Wolfram content only.  The mpmath cross-check serves as ground-truth
# provenance verified above; tests consume the Wolfram-pinned values.
open(PINNED, "w") do io
    print(io, """
# test/_oracle_problems.jl -- Phase 6 oracle, pinned values.
#
# Auto-generated by external/probes/problems-oracle/verify.jl.
# Do NOT hand-edit.  Re-run capture.{wl,py} + verify.jl to update;
# verify.jl asserts cross-source agreement before writing.
#
# Sources (see verify.jl for the exact rtol bounds):
#
#   (A) Pole-bridge demo (worklog 004; Phase-6 v1 acceptance):
#       z ∈ {0.5, 0.95}     -- 3-source: closed-form ≡ NDSolve ≡
#                              mpmath.odefun (rtol 1e-13 / 1e-11).
#       z = 1.05            -- closed-form sole primary; magnitude
#                              sanity check; 80-dps BF-256 round-trip.
#       z = 1.4             -- closed-form primary; optional NDSolve
#                              restart 2-source check at rtol 1e-10.
#
#   (B) FW Table 5.1 (deferred to v2): closed-form ≡ paper-pinned at
#       z ∈ {30, 10⁴, 28.261}; rtol 1e-14.  Plus 80-dps BF-256 round-trip.

""")
    print(io, read(O_WL, String))
end

println("Wrote ", PINNED)
