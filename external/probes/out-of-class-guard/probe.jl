# external/probes/out-of-class-guard/probe.jl
#
# EMPIRICAL δ-SEPARATION PROBE for bead padetaylor-v1ub (out-of-class guard).
#
# Measures the two-order Padé convergence defect δ at every step for:
#   (a) the e^{1/z} essential-singularity APPROACH (z0=-1 → -0.02 at h=0.1);
#   (b) THREE legitimate pole-bridge cases drawn from the regression corpus.
#
# δ per step = |P(1) − P_lo(1)| / (|P(1)| + tiny), where P is the working-order
# [m,n] Padé of the SAME scaled jet the real step builds, and P_lo is the
# reduced-order [m-1,n-1] Padé of the jet truncated by two coefficients.
#
# Goal: CONFIRM the two populations SEPARATE — max in-class pole-bridge δ sits
# safely BELOW the e^{1/z} δ at the steps where the guard must fire. τ/K are
# DERIVED from this data, never fitted to pass.
#
# Run:  julia --project=. external/probes/out-of-class-guard/probe.jl

using PadeTaylor
import PadeTaylor: taylor_coefficients_2nd
import PadeTaylor.RobustPade: robust_pade
import PadeTaylor.PadeStepper: PadeStepperState, pade_step_with_pade!,
                               _rescale_by_powers, _evaluate_pade

const TINY = 1e-300

function _printline(step, z, δ, note)
    zr = round(z; digits = 5)
    δs = δ == 0 ? "0" : string(round(δ; sigdigits = 4))
    println("  step $step  z=$zr  δ=$δs  $note")
end

# Compute δ from the SAME scaled jet the step builds, WITHOUT advancing state.
function delta_at(state, f, order, h)
    coefs   = taylor_coefficients_2nd(f, state.z, state.u, state.up, order)
    scaled  = _rescale_by_powers(coefs, eltype(coefs)(h))
    m = n   = order ÷ 2
    P    = robust_pade(scaled, m, n)
    one_T = one(eltype(scaled))
    P1   = _evaluate_pade(P, one_T)
    # Reduced order: truncate the scaled jet by TWO coefficients, build [m-1,n-1].
    scaled_lo = scaled[1:(end - 2)]
    P_lo  = robust_pade(scaled_lo, m - 1, n - 1)
    P1_lo = _evaluate_pade(P_lo, one_T)
    δ = abs(P1 - P1_lo) / (abs(P1) + TINY)
    return δ, P1
end

# March a 2nd-order ODE with fixed signed h, printing δ at each step.
function march(label, f, u0, up0, z0, zend, h; order = 30)
    println("\n=== $label ===")
    println("  z0=$z0 → $zend  h=$h  order=$order")
    state = PadeStepperState{Float64}(z0, u0, up0)
    dir = sign(zend - z0)
    h_s = dir * abs(h)
    deltas = Float64[]
    step = 0
    while dir * (zend - state.z) > 0
        step += 1
        step > 1000 && break
        gap = zend - state.z
        h_step = dir * min(abs(h), abs(gap))
        local δ, P1
        try
            δ, P1 = delta_at(state, f, order, h_step)
        catch e
            println("  step $step  z=$(round(state.z; digits=4))  δ=THROW ($(typeof(e)))")
            break
        end
        push!(deltas, δ)
        # advance
        try
            pade_step_with_pade!(state, f, order, h_step)
        catch e
            _printline(step, state.z, δ, "ADVANCE-THROW")
            break
        end
        _printline(step, state.z, δ, "")
    end
    if !isempty(deltas)
        println("  max δ = $(maximum(deltas))   final δ = $(deltas[end])")
        # monotone-growth run length ending at the last step:
        run = 1
        for i in length(deltas):-1:2
            if deltas[i] > deltas[i-1]; run += 1; else; break; end
        end
        println("  trailing monotone-growth run length = $run")
    end
    return deltas
end

# ---------------------------------------------------------------------------
# (a) e^{1/z} essential-singularity APPROACH.
# ---------------------------------------------------------------------------
fess(z, u, up) = u * (1 + 2z) / z^4
const U0E  = exp(-1.0)
const UP0E = -exp(-1.0)
d_ess = march("e^{1/z} approach (OUT-OF-CLASS)", fess, U0E, UP0E, -1.0, -0.02, 0.1)

# ---------------------------------------------------------------------------
# (b1) ℘ double pole at z=1 — CPB.1 / problems_test 6.1.3.  Single seg h=1.5
#      brackets the pole; here we MARCH at the small h=0.1 to get per-step δ
#      trajectory THROUGH the pole region (more steps near the pole = harder).
# ---------------------------------------------------------------------------
fW(z, u, up) = 6 * u^2
const WP_U0  = 1.07182251641691742068161429481250292922453268
const WP_UP0 = 1.71033735317678620846433189149624950839231872
# March 0 → 0.95 (just before the pole at z=1) at h=0.1, then the headline
# single big step bridging it.
d_wp_march = march("℘ double pole z=1: march 0→0.95 h=0.1", fW, WP_U0, WP_UP0, 0.0, 0.95, 0.1)
d_wp_seg   = march("℘ double pole z=1: SINGLE seg 0→1.4 h=1.5 (bridges)", fW, WP_U0, WP_UP0, 0.0, 1.4, 1.5)

# ---------------------------------------------------------------------------
# (b2) tan simple pole at π/2 — CPB.2.  Single big step h=3 brackets π/2.
# ---------------------------------------------------------------------------
ftan(z, u, up) = 2 * u + 2 * u^3
d_tan = march("tan simple pole π/2: SINGLE seg 0→2.5 h=3 (bridges)", ftan, 0.0, 1.0, 0.0, 2.5, 3.0)
# also a fine march approaching the pole
d_tan_march = march("tan simple pole π/2: march 0→1.5 h=0.1", ftan, 0.0, 1.0, 0.0, 1.5, 0.1)

# ---------------------------------------------------------------------------
# (b3) coth real pole at z=0 — CPB.4.  u''=-2u+2u^3, u=coth(z).  Bridge from
#      z0=1 toward... coth has a pole at 0.  March from z0=1 down to 0.1 and a
#      single step bridging 0 (z0=-1 → 1, h=2, pole at z=0 interior).
# ---------------------------------------------------------------------------
fcoth(z, u, up) = -2 * u + 2 * u^3
const COTH_U0_at1  = 1.3130352854993313  # coth(1)
const COTH_UP0_at1 = 1 - 1.3130352854993313^2  # u' = 1 - u^2
d_coth_bridge = march("coth real pole z=0: SINGLE seg -1→1 h=2 (bridges)", fcoth,
                      cosh(-1.0)/sinh(-1.0), 1 - (cosh(-1.0)/sinh(-1.0))^2, -1.0, 1.0, 2.0)

# ---------------------------------------------------------------------------
# (b4) simple-pole-riccati u=1/(1-z), pole z=1 EXACT rational — CPB.3.
#      u''=2u^3.  IC at z=0: u=1, u'=1.  Single seg h=1.5 brackets z=1.
# ---------------------------------------------------------------------------
fric(z, u, up) = 2 * u^3
d_ric = march("1/(1-z) pole z=1 (exact rational): SINGLE seg 0→1.4 h=1.5 (bridges)",
              fric, 1.0, 1.0, 0.0, 1.4, 1.5)

# ---------------------------------------------------------------------------
# SUMMARY / separation analysis.
# ---------------------------------------------------------------------------
println("\n\n################ SEPARATION SUMMARY ################")
function summarize(label, d)
    isempty(d) && (println("  $label : (no steps)"); return)
    println("  $label")
    println("      n=$(length(d))  max δ=$(round(maximum(d); sigdigits=4))  ",
            "final δ=$(round(d[end]; sigdigits=4))")
end
println("OUT-OF-CLASS (must fire):")
summarize("e^{1/z} approach", d_ess)
println("\nIN-CLASS pole-bridge (must NOT fire):")
summarize("℘ march 0→0.95", d_wp_march)
summarize("℘ single-seg bridge", d_wp_seg)
summarize("tan single-seg bridge", d_tan)
summarize("tan march 0→1.5", d_tan_march)
summarize("coth single-seg bridge", d_coth_bridge)
summarize("1/(1-z) single-seg bridge", d_ric)

# Max in-class δ across ALL pole-bridge steps:
allinclass = vcat(d_wp_march, d_wp_seg, d_tan, d_tan_march, d_coth_bridge, d_ric)
println("\nMAX IN-CLASS δ (all pole-bridge steps) = ",
        isempty(allinclass) ? "n/a" : round(maximum(allinclass); sigdigits = 4))
if !isempty(d_ess)
    println("e^{1/z} δ tail (last 6 steps): ",
            [round(x; sigdigits = 3) for x in d_ess[max(1,end-5):end]])
end
println("###################################################")
