# 087 — Hybrid driver evaluated the tronquée series on the wrong branch of z^{1/3} (bead w80i)

2026-09-29. Claim from a read-only scout; measured, confirmed, fixed here.

## The claim

`solve_pole_free_hybrid` (`src/IVPBVPHybrid.jl`) forms `z = exp(ζ/2)` and
calls `asymptotic_ic_fn(z)` at the two boundary anchors and, as the BVP
Newton initial guess, at every Chebyshev node of every slice.  FFW's sector
is `-3π/4 < arg z < 9π/4` (FFW md:222,
`references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:222`),
i.e. `-3π/2 < Im ζ < 9π/2`: 3π wide, so for `Im ζ > 2π` `exp` folds `z`
onto principal arguments and `z^(1/3)` lands on another branch.

## Measurements (all real runs)

Probe at FFW's own upper IC point `z₁ = 30e^{13πi/6}` (md:243), i.e.
`ζ₁ = 2 log 30 + 13πi/3`:

- `angle(exp(ζ₁/2)) = 0.5236` (= π/6); continuous arg `Im ζ₁/2 = 6.8068` (= 13π/6).
- principal `pIII_asymptotic_ic(z₁)` = `3.0653 + 0.5386i`, `|u − u_FFW| = 5.389`.
- continued series (`s = 30^{1/3} e^{13πi/18}`) = `−2.000738 + 2.376169i`,
  `|u − u_FFW| = 8.41e-6` (the dropped `a_3 z^{-7/3}` term).
- ratio principal/continued ≈ `−0.503 − 0.867i` ≈ `e^{−2πi/3}`.

What the pre-fix driver did at sector heights reaching `Im ζ > 2π`:

- `im_hi ∈ (2π, 5π/2]`: top anchor principal arg in `(−π, −3π/4)` →
  `pIII_asymptotic_ic` throws (e.g. `im_hi = 2π+0.3`: angle `−2.997`).
- `im_hi > 5π/2` (e.g. FFW's `13π/3 + 0.01`): top anchor ACCEPTED on the
  wrong branch (angle 0.519 vs continuous 6.803), wrong PFS walk launched;
  with `pIII_asymptotic_ic` the run then throws from the BVP initial guess
  (nodes in `2π < Im ζ < 5π/2`).  With a callable that has no sector check
  (e.g. the test/figure helper `z -> tronquee_ic_sheet_test(abs(z), angle(z))`,
  `test/ffw_fig_5_test.jl:170`) nothing throws: emulating that pre-fix
  path gives a hybrid solution with relative error vs. the continued series
  of `0.010, 0.075, 1.733, 1.742` at `Im ζ = 0, 3, 8, 13` (`Re ζ = 2 log 30 − 0.3`)
  — `√3 = |1 − e^{2πi/3}|`, i.e. a silently wrong sheet.
- Fixed driver, same sector: `0.00215, 0.00215, 0.00215, 0.00216` at the
  same points, 2.5 s run.

So the claim is TRUE: a throw for some sector heights, and a silently
wrong boundary IC / solution for others.

**The ODE-residual check suggested in the brief does not discriminate.**
The series is a Laurent series in `s = z^{1/3}` and `z = s³` on every
branch, so each branch is a formal solution of PIII (FFW md:222: "any of the
branches of z^{1/3}" gives a tronquée solution).  Measured in BR.3: residual
`< 1e-5` on both the principal and the continued branch at `z₁`, while the
two `u` values differ by `> 5`.  The discriminating oracle is FFW's
published `u(z₁)` (md:243) = analytic continuation.

## The fix

- `pIII_asymptotic_ic` gains `argz::Union{Nothing,Real} = nothing`
  (`src/IVPBVPHybrid.jl:217`).  Without it: unchanged (principal slice
  `(−3π/4, π]`, `s = z^(1/3)` verbatim, `:249`, `:278`).  With it: sector
  check `(−3π/4, 9π/4)`, consistency `|e^{i argz} − z/|z|| ≤ 1e-8`, and
  `s = |z|^{1/3} e^{i argz/3}` (`:278`).
- `_branch_arg(pp, ζ)` (`:637`): `Im ζ/2` for PIII, `Im ζ` for PV,
  `nothing` otherwise.
- `_call_asymptotic_on_branch(fn, z, θ, label)` (`:621`): plain `fn(z)`
  when θ is the principal arg; otherwise `fn(z; argz = θ)` if the callable
  accepts `argz` (`hasmethod`), else `ArgumentError` with suggestion.
- Used at both boundary anchors (`:545`) and in the BVP initial guess (`:751`).
- Docs: module chapter branch note, `pIII_asymptotic_ic` docstring,
  `solve_pole_free_hybrid` docstring, ADR-0014 addendum.

## Tests

New `test/ivp_bvp_hybrid_branch_test.jl` (19 assertions, BR.1–BR.5):
BR.1 argz branch vs FFW md:243 `u(z₁)`, `u'(z₁)`; BR.2 argz validation +
principal equivalence; BR.3 residual small on both branches (documents why
it is not an oracle); BR.4 driver helper; BR.5 full driver on
`Im ζ ∈ [−3π/2+0.05, 13π/3+0.01]` vs the continued series (`< 1e-2`,
measured 2.2e-3) and a plain callable throws.

Regression: `ivp_bvp_hybrid_test.jl` 65/65, `ffw_fig_5_test.jl` 23/23,
`corpus_bvp_hybrid_test.jl` 10/10.

## Mutation proof

- M1 `_call_asymptotic_on_branch` passes `argz = angle(z)`: BR.4 RED
  (`5.389 < 5e-5`), 14/15 pass.
- M2 `s = z^(1/3)` unconditionally: BR.1 ×2 + BR.4 RED (`5.389`, `0.0597`).
- M3 `_branch_arg` returns `nothing` for PIII: BR.4 ×3 + BR.5 ×1 RED
  (`0.0279 < 0.01`).
- M4 BVP initial guess passes `angle(z)`: BR.5 errors (ArgumentError from
  `pIII_asymptotic_ic`, principal arg `−2.630`).

All restored (`cmp` against backup) and re-run GREEN.

## Open

- The uniform 2.2e-3 relative gap of the hybrid vs the n_terms=2 series is
  the v1 hybrid accuracy (FFW Fig 5 `~10⁻²` corner, beads ykg/sn9a), not
  the branch.
- `src/IVPBVPHybrid.jl` was already far over the 200-LOC rule before this
  change; not split here.
- Bead `padetaylor-5k9` (fn(ζ) API) remains the broader redesign.
