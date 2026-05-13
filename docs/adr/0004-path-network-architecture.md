# ADR-0004 — PathNetwork architecture (Tier-2 path-network module)

**Status**: Proposed (2026-05-13) | **Bead**: `padetaylor-1jf`

## Context

`docs/worklog/004-phase-6-pivot.md` proved that fixed-h real-axis stepping
(Attempts A–D) cannot reach FW Table 5.1 accuracy.  The root cause: FW's
`7.62e-14` at z=30 requires the **5-direction complex-plane path-network**
(FW 2011 §3.1–§3.2, `FW2011_painleve_methodology_JCP230.md:155–166`).
`docs/figure_catalogue.md §6` lists 30+ Tier-2 figures requiring path-
network navigation — none reachable with the current `solve_pade` driver.

## Decision

**File**: `src/PathNetwork.jl` ≤200 LOC (CLAUDE.md Rule 6); split to
`src/PathNetworkStage2.jl` if needed.

**Public API**:
```julia
struct PathNetworkSolution{T}
    visited_z    :: Vector{Complex{T}}
    visited_u    :: Vector{Tuple{Complex{T}, Complex{T}}}
    visited_pade :: Vector{PadeApproximant{Complex{T}}}
    grid_z       :: Vector{Complex{T}}
    grid_u       :: Vector{Complex{T}}
    grid_up      :: Vector{Complex{T}}
end

path_network_solve(
    prob             :: PadeTaylorProblem,
    grid             :: AbstractVector{<:Complex};
    h                :: Real    = 0.5,
    order            :: Integer = 30,
    wedge_angles     :: AbstractVector{<:Real} = [-π/4,-π/8,0.0,π/8,π/4],
    step_selection   :: Symbol  = :min_u,
    step_size_policy :: Symbol  = :fixed,
    max_steps_per_target :: Integer = 1000,
    rtol             :: Real    = 1e-10
) -> PathNetworkSolution
```

**Composition**: calls `PadeStepper.pade_step_with_pade!` per Stage-1
step (same inner loop as `solve_pade`, promoted to `Complex{T}` state).
Does NOT use `step_jorba_zou` / `step_pade_root` — worklog-004 Attempt C
showed Jorba-Zou conflicts with FW's Padé-bridge paradigm.

**Stage 1** (`FW2011_painleve_methodology_JCP230.md:155–166`): initialise
`visited = {z_0}`, shuffle targets, for each target find nearest visited
node (Euclidean; lexicographic tiebreak for reproducibility), step with
5-direction wedge until within h of target, record each intermediate node
with its Padé.  `max_steps_per_target = 1000`; exceed → `IntegrationError`.

**Stage 2** (`FW2011_painleve_methodology_JCP230.md:166, 397`; RF 2014
`ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.md:161–165`): for
each fine-grid node find the nearest visited Stage-1 node with coverage
`|z_f - z_v| ≤ h`; evaluate stored Padé at `t = (z_f - z_v)/h`.
No new Taylor jets.  Gap → store NaN + diagnostic; no silent extrapolation.

**5-direction wedge** (`FW2011_painleve_methodology_JCP230.md:158–159`):
FW 2011 default `[-π/4,-π/8,0,π/8,π/4]` (±22.5°/±45°); RF 2014 users
override to `[-π/6,-π/12,0,π/12,π/6]` (±15°/±30°,
`ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.md:148–153`).
Parametrised; no internal conditional.

**Step selection**: `:min_u` (default) — select candidate with minimum
`|u|` heuristic (`FW2011_painleve_methodology_JCP230.md:159, 354`).
`:steepest_descent` — `θ = arg(-u/u')`, clip to wedge
(`FW2011_painleve_methodology_JCP230.md:362–368`).  Both implemented.

**Step-size policy**: `:fixed` (Tier-2 default, `h = 0.5`).
`:adaptive_ffw` (FFW 2017 §2.4, `FFW2017_painleve_riemann_surfaces_preprint.md:81–92`)
→ throws `ArgumentError` with Tier-4 bead reference; deferred.

**Failure handling** (`docs/unified_path_network_spec.md §7`): denominator
zero at endpoint → try remaining wedge directions → retry `h *= 0.8` →
throw with detail.  No vault; no Jorba-Zou shrinkage.

## Test plan (port-and-verify, CLAUDE.md Rule 4)

`test/pathnetwork_test.jl`, 5 testsets, all mutation-proven.

| id | oracle | mutation → RED |
|---|---|---|
| PN.1.1 | `℘`, 5×5 grid; all visited values match `WeierstrassP` to ≤1e-10 | flip wedge step sign → `max_steps_per_target` throw |
| PN.1.2 | PI tritronquée IC (FW 2011 eq. 4.1, line 215–222), 3×3 grid; all nodes finite | drop `:min_u` → pole-crossing failure on ≥1 node |
| PN.2.1 | Stage-2 11×11 fine grid from PN.1.1 tree; all ≤1e-9 rel; NaN=0 | replace `(z_f-z_v)/h` with `(z_f-z_v)` → O(10²) error |
| PN.2.2 | `℘` path-network to z=30; `u(30)` ≤1e-13 rel vs FW Table 5.1 (`FW2011_painleve_methodology_JCP230.md:385–391`) | disable h^k rescaling → error inflates past 1e-9 |
| PN.3.1 | `:steepest_descent` on PN.1.1; agrees with `:min_u` to ≤1e-10 | use `+u/up` in θ → path diverges |

## Figure-catalogue Tier-2 acceptance

When PN.1.1–PN.3.1 GREEN: FW 2011 Fig 3.1–3.3, 4.2–4.4, 4.7a–f, 4.8,
5.1, 5.2; FW 2014 Fig 3–9, 13–15, 17–24; FW 2015 Fig 2–6, 10, 12–15;
RF 2014 Fig 3, 6, 8–14.  Per-figure acceptance criteria in
`docs/figure_catalogue.md §§1–4`.

## Explicit deferrals

- **BVP composition (Tier 3)**: stub comment at `PathNetwork.jl` bottom:
  `# TIER-3 INTERFACE: bead padetaylor-804`.  No code.
- **CoordTransforms (Tier 4)**: wraps `PathNetwork` externally; zero
  change to internals.  FFW 2017 §2.1–§2.2.
- **SheetTracker (Tier 5)**: same wrapping discipline.  FFW 2017 §2.3.
- **Non-uniform Stage-1 nodes (Tier 4)**: `grid` param accepts any
  `AbstractVector{<:Complex}`; uniform generator only in Tier 2.
- **EdgeDetector (`padetaylor-c2p`)**: separate module; not in scope.

## Consequences

- `src/PathNetwork.jl` added; umbrella gains `include` + re-exports for
  `path_network_solve, PathNetworkSolution`.
- No existing `.jl` files modified beyond the umbrella.
- Tier-2 figure reproduction unblocked.  Tier-3 BVP interface planted.

## References

- `docs/unified_path_network_spec.md §1–§7` — full algorithm spec.
- `docs/figure_catalogue.md §6` — Tier-2 figure list.
- `docs/worklog/004-phase-6-pivot.md` — failure analysis.
- ADR-0001 — four-layer architecture extended at driver layer.
