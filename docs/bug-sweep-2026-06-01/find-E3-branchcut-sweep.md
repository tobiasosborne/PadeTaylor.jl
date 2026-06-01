# Bug sweep E3 — branch-cut / principal-value errors across `src/`

Date: 2026-06-01
Auditor area: WHOLE `src/` — branch-cut / principal-value handling
(`sqrt`, `log`, fractional powers, `angle`, `atan`, `cis`, `exp`).
Read-only audit (no `julia` / `Pkg` / tests run, per Rule 7).

## Area

Every `src/*.jl` site that consumes a multi-valued elementary function
on complex data was located by grep and read in context:

- `(6z)^{1/3}` / `x^{1/3}` tritronquée ridge — `PainleveHierarchy.jl`
- `ζ = 2 log z`, `ζ = log z`, `η = log ζ` coordinate maps and their
  `exp` inverses — `CoordTransforms.jl`, `SheetTracker.jl`
- `z^{1/3}` PIII asymptotic IC — `IVPBVPHybrid.jl`
- pole-distance / direction `angle(...)` — `PathNetwork.jl`,
  `VectorPathNetwork.jl`, `SheetTracker.jl`, `BranchTracker.jl`
- the `adjoint`/`'`-vs-`.'` complex-transpose class (historical
  `padeapprox.m` bug) — `RobustPade.jl`, `SharedPade.jl`, `LinAlg.jl`
- fractional-power step formulas on real magnitudes — `StepControl.jl`,
  `VectorStepControl.jl`, `PadeStepper.jl`

## References checked

- `external/chebfun/padeapprox.m:106-138` — `[U,S,V]=svd(C,0)`,
  `b=V(:,n+1)`, `D=diag(abs(b)+sqrt(eps))`, `[Q,R]=qr((C*D)')`
  with the in-source comment **"until July 2018 there was an erroneous
  `.'` here"** (line 113).
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:368`
  — steepest-descent direction `θ = arg(-u(z₀)/u'(z₀))`: "If this
  direction falls inside the wedge, we accept it … If not, we choose
  the edge of the wedge closest to this steepest descent direction".
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:158`
  — five-direction wedge (straight at goal ±22.5°, ±45°), `h = 0.3`.
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:354`
  — §5.3 "select the lowest choice in terms of |u(z)|".
- `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:41`
  — PIII transform `z = e^{ζ/2}`, `u = e^{-ζ/2} w`.
- `.../FFW2017…md:45` — PV transform `z = e^ζ`, `u = w`.
- `.../FFW2017…md:103` — PIII sheet map: strip `-2π+4πs < Im ζ ≤ 2π+4πs`
  ↔ sheet `s`, `u = e^{-ζ/2} w`.
- `.../FFW2017…md:146,149,154` — PVI `z=e^ζ`, `ζ=e^η`, branch points
  `η = log|2πk| + i·arg(2πik)`.
- `.../FFW2017…md:157` — `Re η < log 2π` is branch-point-free.
- `.../FFW2017…md:222` — tronquée PIII pole-free region `-3π/4 < arg z
  < 9π/4` ↔ `-3π/2 < Im ζ < 9π/2`; `u ~ z^{1/3} = (e^{ζ/2})^{1/3} =
  e^{ζ/6}`; angular width `3π`, "any of the branches of `z^{1/3}`".
- `.../FFW2017…md:191` — tronquée seed obtained "far out on the positive
  real axis".

## Findings

### [MEDIUM] `:steepest_descent` wedge selection uses an un-wrapped angular distance

**Location:** `src/PathNetwork.jl:994-997` (function `_select_candidate`,
`:steepest_descent` branch).

```julia
θ_sd = abs(up_cur) > 0 ? angle(-u_cur / up_cur) : T(goal_dir)
offsets = (T(goal_dir) + T(θ) for θ in wedge_angles)
return argmin(abs(θ_sd - off) for off in offsets)
```

**Ground truth (cited):** `FW2011…md:368` — steepest-descent direction
`θ = arg(-u(z₀)/u'(z₀))`; "If this direction falls inside the wedge, we
accept it … If not, we choose the edge of the wedge closest to this
steepest descent direction." Closeness here is an *angular* relation on
the circle, so it must be taken mod 2π.

**Code behavior:** `θ_sd = angle(-u/u')` is wrapped to `(-π, π]`, but the
candidate `offsets = goal_dir + wedge_angle` are NOT wrapped:
`goal_dir = angle(target - z_cur) ∈ (-π, π]` (`PathNetwork.jl:531`) and
`wedge_angles ∈ [-π/4, π/4]` (`DEFAULT_WEDGE`, `PathNetwork.jl:153`), so
`off ∈ (-5π/4, 5π/4]`. The selection minimises the **linear** difference
`abs(θ_sd - off)` rather than the circular distance
`abs(rem(θ_sd - off, 2π, RoundNearest))`. When `θ_sd` and an `off`
straddle the `±π` cut, the linear difference is inflated by `2π` and the
true-nearest wedge can lose to a farther one.

Secondary fidelity gap (same lines): the reference accepts the
continuous `θ_sd` directly when it lies inside the 90° wedge and snaps to
the *wedge edge* otherwise; the impl always snaps to one of the five
discrete rays. This is an algorithmic deviation independent of the wrap.

**Mechanism (intermittent discontinuity):** The wrap only mis-fires when
`goal_dir ± wedge` crosses `±π` — i.e. when the current target lies in a
direction near the negative-real direction from `z_cur` — AND `θ_sd`
sits on the opposite side of the cut. As a Stage-1 walk wanders and the
target bearing rotates through `±π`, this condition switches on and off,
so an occasional step is steered into the wrong wedge ray. A single
wrong-direction step displaces that node (and everything chained off it),
injecting a localised, geometry-dependent kink into the computed field —
exactly an intermittent discontinuity. It is gated behind the non-default
`step_selection = :steepest_descent` (default is `:min_u`,
`PathNetwork.jl:380`), and the only test exercising it
(`test/pathnetwork_test.jl:294-320`, "PN.3.1") uses right-half-plane
targets (`1.4`, `1.2+0.4im`, `0.6+0.3im`) whose `goal_dir ≈ 0`, so the
`±π` wrap is never triggered and the bug is uncaught.

**Intermittent?** Yes — data/geometry dependent (only when bearing
crosses `±π`), and only under the opt-in `:steepest_descent` mode.

**Confidence:** 0.75 (the un-wrapped comparison is a demonstrable
mismatch against `FW2011…md:368`'s "closest" semantics on the circle;
confidence is not higher only because the affected mode is non-default
and the magnitude of the resulting kink was not numerically measured in
this read-only pass).

### [LOW] Public `pIII_asymptotic_ic` cube-root sheet trap behind a wrap-blind guard

**Location:** `src/IVPBVPHybrid.jl:200-207` (sector guard) and
`:215` (`s = z^(1/3)`).

**Ground truth (cited):** `FFW2017…md:222` — the tronquée PIII pole-free
region is `-3π/4 < arg z < 9π/4` with `u ~ z^{1/3}` over an angular width
of `3π`, and the seed is taken "far out on the positive real axis"
(`md:191`); the upper part `(π, 9π/4)` is a non-principal sheet of the
cube root. `FFW2017…md:103` makes the sheet index `s` explicit.

**Code behavior:** The guard tests `-3π/4 < angle(z) ≤ π` on Julia's
*principal* `angle` (`IVPBVPHybrid.jl:200-207`), then forms
`s = z^(1/3)` with Julia's *principal* cube root (`:215`). When the
hybrid driver computes a boundary point at `Im ζ = 9π/2 − ε`
(`IVPBVPHybrid.jl:487-490`), `z = exp(ζ/2)` wraps `arg` from `9π/4 − ε/2`
down to the principal `π/4 − ε/2`. That wrapped value PASSES the guard
(`π/4 < π`), so the guard cannot detect that the requested point is on a
non-principal cube-root sheet; `z^{1/3}` then returns the
`arg = π/12` branch instead of the physically required `arg = 3π/4`
branch (a `2π/3` rotation — a different sheet of the tronquée).

**Mechanism (intermittent discontinuity):** If a caller drives
`solve_pole_free_hybrid` with a sector whose `im_hi > 2π` while passing
the *default* `pIII_asymptotic_ic` as `asymptotic_ic_fn`, the upper-edge
seed silently lands on the wrong cube-root sheet, so the upper part of
the rendered sector jumps to a different branch than the lower part —
a discontinuity along the sheet seam. The shipped figure avoids this:
`figures/ffw2017_fig_5.jl` uses `IM_HI = 2π − 0.05`
(`ffw2017_fig_5.jl:174`) and a bespoke `tronquee_ic_sheet` that consumes
the *unfolded* `arg_z` directly (`ffw2017_fig_5.jl:213-240`), bypassing
the principal-`z^{1/3}` helper entirely. So this is a latent API trap,
not a live defect in shipped output.

**Intermittent?** Only if a future caller crosses the principal strip
with the default IC; dormant in current usage.

**Confidence:** 0.5 (mechanism is real and reference-confirmed, but it is
not reachable from any shipped path and is partly self-documented at
`IVPBVPHybrid.jl:196-207`).

### [LOW] `winding_delta` folds steps that subtend `|Δθ| ≥ π` (documented precondition, but reachable under adaptive/large steps)

**Location:** `src/SheetTracker.jl:283-292` (`winding_delta`), consumed by
`src/BranchTracker.jl:156-167` (`step_sheet_update`) in PathNetwork
`cross_branch = true` mode (`PathNetwork.jl:624-630`).

**Ground truth (cited):** `FFW2017…md:178` — on the pole-free sheet the
adaptive method "chose a few, relatively large steps"; sheet index is the
accumulated winding (`SheetTracker` docstring cites `md:187-189`).

**Code behavior:** `Δθ = angle(z_new−branch) − angle(z_old−branch)` is
normalised to `(-π, π]`. The docstring (`SheetTracker.jl:275-281`)
explicitly requires "no single path step has `|Δθ| ≥ π`" and justifies it
for `h = 0.5`, branch distance `> 1` (`|Δθ| ≤ arcsin(0.5) ≈ 30°`). Under
`:adaptive_ffw` with a `node_separation` `R(ζ)`, the accepted step can be
"relatively large" near the pole-free sheet, and the assumed branch
distance `> 1` is not enforced.

**Mechanism (intermittent discontinuity):** A single large step that both
crosses a branch cut and subtends `|Δθ| ≥ π` as seen from that branch
would be folded by the `±2π` normalisation; the per-step sheet counter
then advances by the `sign(Δθ)` of the folded (wrong-sign) value, so a
crossing is mis-counted and the post-crossing solution is read off the
wrong Riemann sheet — a discontinuity at the seam. Requires the
conjunction of large step + near branch + actual cut crossing, hence rare
and geometry-dependent.

**Intermittent?** Yes (only when the documented `|Δθ| < π` precondition is
violated by an adaptive/large step near a branch).

**Confidence:** 0.3 (the fold is real and the precondition is only
documented, not asserted; but I did not find a shipped configuration that
demonstrably violates it, and `cross_branch` is itself an opt-in mode).

## Areas verified correct

- **`RobustPade.jl:443` `adjoint(C * D)` is the CORRECT post-2018
  convention.** `padeapprox.m:113` uses `(C*D)'` (adjoint) with the
  in-source comment "until July 2018 there was an erroneous `.'` here";
  the Julia port matches the fixed MATLAB. The preceding
  `b_init = Vt[end,:]` (`:435`) is `conj(V[:,end])` vs MATLAB
  `V(:,n+1)`, but it feeds only `D = diag(|b|+√ε)` (`:442`), where the
  conjugation is annihilated by `abs`. Verified against
  `padeapprox.m:106-117`.
- **`CoordTransforms.jl:120-186` PIII/PV maps match FFW exactly.**
  `ζ = 2 log z`, `z = exp(ζ/2)`, `u = e^{-ζ/2} w` (`pIII_z_to_ζ`,
  `pIII_ζ_to_z`) reproduce `FFW2017…md:41,103`; the `Im ζ = 2 arg z`
  relation gives `arg z = -3π/4 ↔ Im ζ = -3π/2` and `9π/4 ↔ 9π/2`,
  matching `md:222`. PV `ζ = log z`, `z = exp(ζ)` matches `md:45`.
  Forward/inverse pairs are mutually consistent.
- **`SheetTracker.jl:236-261` PVI `η`/`ζ`/`z` composition matches FFW.**
  `ζ = log z`, `η = log ζ`, inverses `ζ = exp(η)`, `z = exp(exp(η))`
  reproduce `FFW2017…md:146-154`; the branch-point list and the
  `Re η < log 2π` free region (`md:149,157`) are cited verbatim.
- **`winding_delta` normalisation interval is correct.** `Δθ ≤ -π ⇒ +2π`,
  `Δθ > π ⇒ -2π` maps the raw difference (in `(-2π,2π)`) into `(-π,π]`,
  with the `Δθ = -π` boundary correctly sent to `+π`
  (`SheetTracker.jl:286-290`).
- **`BranchTracker.jl:112-120` `segment_crosses_cut` Cramer's-rule is
  correct.** Parametrising `z_cur + t·d = branch + s·u`, projecting on
  `conj(u)` and `conj(d)` and using `imag(u·conj(u)) = imag(d·conj(d)) =
  0` reproduces `t = imag(δ·conj(u))/imag(d·conj(u))` and
  `s = imag(δ·conj(d))/imag(d·conj(u))` exactly as coded (the sign works
  out via `imag(u·conj(d)) = -imag(d·conj(u))`). Half-open `0 < t < 1`,
  `s > 0` is a deliberate, consistent endpoint convention.
- **`PainleveHierarchy.jl:557-575` tritronquée sheet angle is continuous
  on the pole-free sector and discontinuous only on the pole wedge.**
  `θ = φ₀≤0 ? φ₀+4π : φ₀+2π` is continuous through `θ = 3π` as `arg`
  crosses the negative real axis (`φ₀ = π⁻ → θ = 3π⁻`; `φ₀ = -π⁺ →
  θ = 3π⁺`), i.e. the centre of the pole-free sector is smooth. The only
  jump is at `φ₀ = 0` (positive real axis), which is exactly the `V_0`
  pole wedge / cut per `md:222` ("any of the branches of `z^{1/3}`") and
  the module docstring `:77-94`. `Y = -∛6·|x0|^{1/3}·exp(iθ/3)` gives the
  ODE-consistent `+∛6|x|^{1/3}` on the negative real axis. The
  `n_terms=2` correction `∝ z^{-2k}` is an even power and single-valued —
  branch-independent, as documented `:496-499`. The negative-real legacy
  branch (`:539-551`) uses real `r = |x0|` powers — no complex cut.
- **`IVPBVPHybrid.jl:200-215` is correct WITHIN the principal strip.** For
  the shipped sector (`Im ζ ∈ (-3π/2+, 2π-)`, `figures/ffw2017_fig_5.jl:174-175`)
  the wrapped `arg z` stays in `(-3π/4, π]` and `z^{1/3}` is the correct
  branch; the `≤ π` boundary inclusion keeps `arg z = π` continuous from
  below. The `make_principal_sheet_ic` adapter (`ffw2017_fig_5.jl:248-254`)
  feeds principal `angle(z)`.
- **`StepControl.jl:162-198`, `VectorStepControl.jl`, `PadeStepper.jl:505`
  fractional powers are branch-safe.** Every `^(1/k)` / `^(1/(n+1))` acts
  on `abs(...)` magnitudes or real positive `(ε/aux)`, `(k·Tol/T_h)` — a
  positive-real base, so the principal real power is unambiguous. Verified
  the step formulas against the cited `TaylorIntegration.jl stepsize.jl`
  port comments and FFW `md:88-91`.
- **`LinAlg.jl:75-109` SVD dispatch is transpose-clean.** Returns
  `(F.U, F.S, F.Vt)` from `LinearAlgebra.svd` / `GenericLinearAlgebra.svd`
  uniformly for real and complex element types; no manual transpose that
  could confuse `'` vs `.'`.
- **`PainleveClosedForm.jl:234` `cbrt(4)` and `Painleve`/`CoordTransforms`
  `exp` calls** act on real positive or already-unfolded arguments — no
  cut crossing in solution arithmetic.
- **`VectorWedgeStep.jl:392-439` / `VectorPathNetwork.jl:808` wedge
  selection is wrap-free.** Directions enter only through
  `cos(θ)`, `sin(θ)` (`VectorWedgeStep.jl:402-403`), and selection scores
  on `_candidate_pole_disc` and `-|wedge_angle|` — no naive angle
  subtraction, so the `PathNetwork.jl:997` wrap defect does not recur in
  the vector path.
