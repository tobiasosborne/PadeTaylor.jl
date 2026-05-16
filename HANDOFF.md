# HANDOFF.md — for the next agent picking up PadeTaylor.jl

> **Read this top to bottom before doing anything else.** Then read
> `CLAUDE.md` (the discipline rules). Then `RESEARCH.md` §1.1, §2,
> §3.3, §5, §7. Then `DESIGN.md` §4. Then this file again.
>
> If you skim and skip the references, you will repeat mistakes the
> previous session already paid for. The frictions surfaced are
> recorded in `docs/worklog/001-stages-Z-1-2-handoff.md`.

## TL;DR — where we are

**Stage 0 (research) and Stage 1 (design) are complete.** Stage 2
(implementation) is in progress — Phases Z, 1–6 (original DESIGN.md
scope), plus the path-network and BVP modules (Phases 10–11, new
docs/design_section_6_path_network.md scope), are all shipped:

| Phase | Module | Status | Tests |
|---|---|---|---|
| Z | scaffold + 3 ADRs + CLAUDE.md | ✅ shipped | 12 placeholders |
| 1 | `LinAlg.pade_svd` | ✅ shipped, mutation-proven | 14 |
| 2 | `RobustPade.robust_pade` | ✅ shipped, mutation-proven | 35 |
| 3 | `Coefficients.taylor_coefficients_*` | ✅ shipped, mutation-proven | 125 |
| 4 | `StepControl.step_jorba_zou` + `step_pade_root` | ✅ shipped, mutation-proven | 7 |
| 5 | `PadeStepper.pade_step!` | ✅ shipped, mutation-proven | 16 |
| 6 | `Problems.PadeTaylorProblem`, `solve_pade` | ✅ shipped (Padé-vs-Taylor pole-bridge demo), mutation-proven; see worklog 005 | 11 |
| 7 | `CommonSolveAdapter` ext (`PadeTaylorAlg`) | ✅ shipped per worklog 009 | 34 |
| 8 | `PadeTaylorArblibExt` ext (SVD-only v1) | ✅ shipped per worklog 010 | 30 |
| **9** | Tier C: PI tritronquée pole-field qualitative | ✅ **shipped** — 4-of-5 pole-free sectors + conjugate symmetry + leading-pole magnitude at 25×25 [-4,4]² (worklog 012) | 12 |
| **10** | `PathNetwork.path_network_solve` | ✅ **shipped + FW Table 5.1 GREEN @ 2.13e-14 BF-256 (beats FW's 8.34e-14)**; see worklogs 006 + 008 | 33 |
| **11** | `BVP.bvp_solve` | ✅ **shipped (Tier-3 Chebyshev-Newton BVP)**, mutation-proven; see worklog 006 | 43 |
| **12 v1** | `Dispatcher.dispatch_solve` (1D chain) | ✅ **shipped (Tier-3 IVP↔BVP chain composition)**, mutation-proven; see worklog 007 | 45 |
| **12.5** | `EdgeDetector` — 5-point Laplacian pole-field classifier | ✅ **shipped** (worklog 011); FW 2011 §3.2.2 verbatim | 743 |
| **12 v2** | `LatticeDispatcher.lattice_dispatch_solve` — 2D-grid composition | ✅ **shipped** (worklog 013) — per-row BVP fill on smooth runs flanked by IVP cells (FW2011...md:190); FW Fig 4.1 quantitative pin deferred to follow-up `padetaylor-0c3` | 76 |
| 13 | `CoordTransforms` (FFW 2017 PIII/PV) | ✅ **shipped** (worklog 017) | Tier-4 |
| 14 | `SheetTracker` (FFW 2017 PVI) | ✅ **shipped** (worklog 018) | Tier-5 |
| 15 | `Painleve.PainleveProblem` — per-equation builder (ADR-0006) | ✅ **shipped** (worklog 026) | — |
| 15+ | `PoleField.extract_poles` — poles from a solved path-network | ✅ **shipped** (worklog 027) | — |
| 15++ | `EdgeGatedSolve` — IVP confined to the pole field (FW md:401) | ✅ **shipped** (worklog 028) | — |
| 15+3 | FW Fig 4.1 full BVP+IVP+BVP composition | ✅ **shipped** (worklog 029) | — |
| 15+4 | `PainleveSolution` wrapper + Makie ext (ADR-0007) | ✅ **shipped** (worklog 030) | — |
| 15+5 | `tritronquee` / `hastings_mcleod` named constructors (ADR-0008) | ✅ **shipped** (worklog 031) | — |
| **12.5 fix** | `EdgeDetector` h-aware default level (ADR-0009) | ✅ **shipped** (worklog 032); fixes `edge_gated_pole_field_solve` stall on any lattice finer than `h_grid ≈ 0.25` | +22 ED.5.* |
| 15+6 | `pii_rational`, `pii_airy`, `piv_entire` closed-form families (ADR-0010) | ✅ **shipped** (worklog 033); the second of the three ADR-0008 named-solution kinds — parametrised families with exact closed-form oracles | +56 CF.* |
| **A1** | Adaptive Padé step `:adaptive_ffw` (ADR-0011) — FFW 2017 §2.1.2 controller `q = (k·Tol/T(h))^(1/(n+1))` | ✅ **shipped** (worklog 034, bead `padetaylor-8ui`); foundation for all FFW 2017 figure reproduction | +43 AS.* |
| **A2** | Non-uniform Stage-1 nodes — `node_separation::Function` kwarg (ADR-0012) | ✅ **shipped** (worklog 035, bead `padetaylor-1a3`); composes with `:adaptive_ffw` (R seeds, controller may only shrink) | +295 NU.* |
| **A6** | `IVPBVPHybrid.solve_pole_free_hybrid` (ADR-0014) — FFW 2017 §3 PFS↔BVP coupling for pole-free sectors | ✅ **shipped** (worklog 039, bead `padetaylor-0co`); v1 corner: asymptotic-series helper hard-codes `a_3+ = 0` (bead `padetaylor-ykg` opens the lift) | +66 IB.* |
| **FFW Fig 6** | PV generic three-sheet — first FFW 2017 figure reproduced | ✅ **shipped** (worklog 036, bead `padetaylor-mbu`); per-sheet errors 3.3e-12 / 7.0e-8 / 9.8e-8 vs FFW 3e-9 / 7e-6 / 2e-5 — **beats FFW by 2-3 orders** | +24 FF6.* |
| **FFW Fig 1** | PIII three-sheet spiral pole field — FFW headline figure | ✅ **shipped** (worklog 037, bead `padetaylor-2pg`); 2664 visited nodes vs FFW's 2701 (-1.4%); sheet-0 conjugate-symmetry median `4.05e-15` vs FFW Exp-2 `1e-6` (9 orders better); Cartesian-resample rendering | +18 FF1.* |
| **FFW Fig 4** | PV tronquée three-sheet | ✅ **shipped** (worklog 038, bead `padetaylor-pt2`); pole-free sector visible on sheet 0; pole-density gradient 8× across sheets 1+2 | +43 FF4.* |
| **FFW Fig 5** | PIII tronquée + cond-number heatmap — uses A6 hybrid | ✅ **shipped** (worklog 040, bead `padetaylor-3ug`); cond-number pin `κ_r(z=30) = 157.30` matches FFW md:264 to ±1; v1 corner: sector error `~10⁻²` (4 orders looser than FFW Table 2), gated on `padetaylor-ykg` lift | +23 FF5.* |
| **A3** | η-plane PVI RHS (`pVI_eta_transformed_rhs` + `pVI_z_to_η`/`η_to_z` + `:transformed_eta` frame) | ✅ **shipped** (worklog 041, bead `padetaylor-riu`); FFW 2017 eq. 5 (md:154) verbatim; branch-point-free region `Re η < log(2π)`; closes worklog 018 deferral | +42 ET.* |
| **A4** | Constrained-wedge routing + per-branch sheet bookkeeping (ADR-0013, new `src/BranchTracker.jl`) | ✅ **shipped** (worklog 042, bead `padetaylor-bho`); both refuse mode (`branch_points` + `branch_cut_angles`) and cross mode (`cross_branch = true` + per-branch sheet counter); FFW 2017 §2.2.2 (md:163-189) | +105 BT.* + PNB.* |
| **A5** | Sheet-aware Stage-2 (`grid_sheet` kwarg + `eval_at_sheet` accessor) | ✅ **shipped** (worklog 043, bead `padetaylor-hed`); ADR-less additive extension to A4 substrate; unblocks B5 multi-sheet rendering | +30 SA.* |
| **FFW Fig 2** | PVI three-method reproduction (η + ζ-refuse + ζ-cross) — first B5 figure | ✅ **shipped** (worklog 044, bead `padetaylor-n93`); demonstrates full A1+A2+A3+A4+A5+A6 composition on a single tronquée P_VI problem | +19 FF2.* |
| **Stage-2 extrapolate** | `extrapolate=true` kwarg + new `eval_at` accessor (ADR-0015) — fixes white-gap artefact across all FFW figures | ✅ **shipped** (worklog 045, bead `padetaylor-afs`); aligns Stage-2 with FFW md:62 spec (no disc-radius cutoff); Fig 1/2/4/6 now render at **100% coverage** (was 0.9–98%) | +30 SX.* |
| **FFW Fig 3** | PVI phase portraits in ζ + z planes (sheet 0 + alt sheet) | ✅ **shipped** (worklog 046, bead `padetaylor-a1l`); reuses Fig 2's A4 cross-mode + A5 sheet-aware machinery on a wider ζ-window with HSV phase-portrait rendering + z-plane Cartesian resample (Fig 4 pattern); **zero new `src/` code** | +20 FF3.* |
| **FFW Fig 7** | Generic (non-tronquée) PVI in η + ζ + z planes — closes the 11-step FFW arc | ✅ **shipped** (worklog 047, bead `padetaylor-mgx`); first generic PVI figure; FFW md:281 `log10\|u\|` z-plane convention; first multi-branch-point routing (`branch_points = (0im, ±2πi)`) — confirmed `BranchTracker` tuple API works for the multi-point case; +41 FF7.* assertions | +41 FF7.* |

**2508 / 2508 tests passing** as of worklog 047 (`HEAD = 1331755`, all committed; push pending).  **v0.1.0 tagged** (`38a49ae`, `CHANGELOG.md` ships); **Documenter site** generated (`30b3298`).  Since v0.1.0: classical-Padé F64 default (worklogs 020+021); **thirteen FW 2011 figures reproduced** in a new `figures/` project (Fig 3.1–3.3, 4.1, 4.2–4.4, 4.7, 4.8, 5.1, 5.2; worklogs 022–029); **ADR-0006 + `src/Painleve.jl`** — the `PainleveProblem` per-equation builder layer (worklog 026); `PoleField` (027), `EdgeGatedSolve` (028), the FW Fig 4.1 full composition (029), the `PainleveSolution` wrapper + Makie extension (ADR-0007, worklog 030), the `tritronquee` / `hastings_mcleod` named constructors (ADR-0008, worklog 031), the `EdgeDetector` h-aware level fix that unblocked the N=641 hero render (ADR-0009, worklog 032), and the closed-form Painlevé families `pii_rational` / `pii_airy` / `piv_entire` (ADR-0010, worklog 033).  See the "Session 2026-05-15" entry below for the worklog-032–033 rundown; older sessions live in `HANDOFF_ARCHIVE.md`.

**Phase 6 shipped 2026-05-09 on the pivoted scope** — the v1
acceptance is a Padé-vs-Taylor pole-bridge demonstration (one stored
Padé bridges the lattice pole of `℘(z + c₁; 0, c₂)` at `z = 1`,
giving correct values at `z ∈ {0.5, 0.95, 1.05, 1.4}` while plain
truncation diverges past `z = 1`); see `docs/worklog/005-phase-6-
implementation.md`.  Pivot rationale + failure analysis in
`docs/worklog/004-phase-6-pivot.md`.

**Phases 10 + 11 shipped 2026-05-13** — Tier-2 path-network +
Tier-3 Chebyshev-Newton BVP solver, both GREEN + mutation-proven.
Three commits land (`910aab9` BVP ground truth; `0ada60f` Phase 10;
`cc7d8ca` Phase 11).  The session's deliverables include:

  - Phase 10 `PathNetwork.path_network_solve`: FW 2011 §3.1 5-direction
    wedge path-tree + Stage-2 fine-grid barycentric extrapolation.
    Generic in `T <: AbstractFloat` for `Float64` + `Complex{T}`.
    PN.1.1, PN.1.2, PN.2.1, PN.4.1 GREEN; PN.2.2 (FW Table 5.1
    quantitative ≤1e-13) + PN.3.1 (:steepest_descent test) deferred
    to follow-up bead `padetaylor-yt1`.
  - Phase 11 `BVP.bvp_solve`: Chebyshev spectral collocation per FW
    2011 §3.2 + Trefethen SMIM ch.6/13.  Generic in `T <: AbstractFloat`
    for `Float64` + `BigFloat-256` + `Complex{T}`.  Step-norm Newton
    convergence (eps^(3/4) default).  BV.1.1-BV.5.1 GREEN.
  - Full **figure-acceptance catalogue** at `docs/figure_catalogue.md`
    — 79 FW-family figures tiered T0-T5 with per-row acceptance.
  - **Unified path-network spec** at `docs/unified_path_network_spec.md`
    — 14 sections synthesised across all 5 FW-family papers.
  - **ADR-0004** at `docs/adr/0004-path-network-architecture.md`.
  - **BVP canonical recipe** at `references/bvp_recipe.md`
    — 8 sections covering Chebyshev nodes, D₁/D₂, Newton, barycentric.
  - **Octave oracle** at `external/probes/bvp-oracle/{capture.m, oracles.txt}`
    using DMSUITE chebdif/chebint; 47 pinned constants in
    `test/_oracle_bvp.jl`.
  - **Marker-converted primary refs**: Weideman-Reddy 2000 (DMSUITE),
    Trefethen SMIM 2000, Berrut-Trefethen 2004 (SIAM Review).  All
    paywalled, acquired via TIB VPN.
  - **DMSUITE source** at `external/DMSUITE/` (gitignored, MIT/public).

The session's worklog at `docs/worklog/006-phases-10-11-path-network-bvp.md`
covers orchestration pattern (Opus codes, Sonnet does grunt work),
algorithmic findings (step-norm vs residual-norm Newton convergence,
the Octave oracle's PI vs equianharmonic-℘ spec-drift catch), and
mutation-proof procedures for both modules.

## Project shape

- `references/` — load-bearing PDFs (FW 2011, GGT 2013, FW 2014/15,
  RF 2014, FFW 2017, Jorba-Zou 2005, Mezzarobba 2019), each with
  `marker_single`-converted markdown under `references/markdown/<name>/`.
- `external/chebfun/padeapprox.m` — the canonical Padé reference impl.
  **Do not** modify; it is the oracle source.
- `external/probes/` — empirical probes that informed the design:
  - `taylorseries-arb/` — proves `TaylorSeries.jl::Taylor1{Arb}` works
    at order 80 (5/5 tests); `probe.jl` + `SYNTHESIS.md` are the
    record. **Phase 3 builds on this empirical confirmation.**
  - `padeapprox-oracle/` — `capture.m` runs `padeapprox.m` under
    Octave, writes `oracles.txt` (copied to `test/_oracles.jl`).
- `docs/adr/` — ten accepted ADRs (0001–0010); read 0001–0004 before
  touching the core/path-network code.
- `docs/worklog/` — frictions / lessons learnt; 33 shards (001–033).
- `src/` — 20 modules (the 18 originally substantive across tiers 0–5
  + the Painlevé layer + `PainleveClosedForm.jl`, plus the umbrella
  `PadeTaylor.jl`).
- `test/` — `runtests.jl` orchestrates per-module test files.
- `RESEARCH.md` — Stage 0 deliverable (1165 lines). Read §1.1 (FW 2011
  algorithm), §2 (GGT 2013 + reweighting), §3.3 (TaylorSeries+Arb
  empirical), §5 (bigfloat-SVD landscape), §7 (resolutions of PRD
  open questions). Skim §3 (adjacent packages) and §4 (verified
  arithmetic).
- `DESIGN.md` — Stage 1 deliverable. **§4 is the granular execution
  plan; treat as the source of truth for what each phase must build.**
- `CLAUDE.md` — 13 numbered discipline rules. Re-read after any
  context compression.

## Hard rules — do not deviate

1. **Ground truth before code.** Before writing a line of Phase 3 or
   later, open `RESEARCH.md` and the relevant FW 2011 / GGT 2013
   markdown line range. Cite by `references/markdown/<file>.md:<lines>`
   in commits.
2. **TDD always.** Tests RED before impl. For port-and-verify
   modules (Phase 4 step formulas), capture an oracle from the source
   (Octave for `padeapprox.m`; TaylorIntegration.jl for Jorba-Zou) and
   assert against it.
3. **Mutation-prove every load-bearing test.** Before committing a
   GREEN phase: perturb the impl in a *meaningful* way, confirm
   tests RED, restore. Document the procedure in test comments.
   The naïve mutation may not bite — see worklog shard 001 for the
   `GenericLinearAlgebra` piracy discovery.
4. **`bd create` any friction**. Beads is the only tracker.
   `bd ready` at start; `bd close <ids>` at end. Never `bd edit`.
5. **Commit after each GREEN phase** with a message that:
   - lists what was added (file paths + LOC)
   - lists tests added + which oracle was used
   - records the mutation procedure that bit
   - cites the references by line number
6. **Pushing to remote is allowed any time** (updated 2026-05-14).
   Push freely once work is committed and the suite is GREEN; no need
   to wait for an explicit instruction. See CLAUDE.md "Session close".
7. **No parallel Julia agents.** One `julia` process at a time. The
   precompile cache is brittle to concurrent invocations.
8. **No author outreach.** Per project decision (CLAUDE.md Rule 12).

## Specific instructions for Phase 3 (`Coefficients`)

This is your immediate next task.

### Goal

Implement `PadeTaylor.Coefficients.taylor_coefficients_1st(f, z0, y0,
order)` and `taylor_coefficients_2nd(f, z0, y0, y1, order)` that
return the Taylor coefficients `c_0, …, c_order` of the local Taylor
expansion of the ODE solution at `z0`.

### Ground truth to read first

- `RESEARCH.md §1.1.1` — FW 2011 §2.1 method (b): "use one term,
  substitute, integrate, gain one correct coefficient per pass."
- `references/markdown/FW2011_painleve_methodology_JCP230/
  FW2011_painleve_methodology_JCP230.md:79-107` — the source.
  Specifically:
  - Lines 88-94 describe method (a) — symbolic differentiation. Don't
    use this for general `f`.
  - Lines 96-107 describe method (b) — operator-overloading. **Use this.**
- `RESEARCH.md §3.3` — empirical confirmation that
  `TaylorSeries.jl::Taylor1{T}` works for `T ∈ {Float64, BigFloat,
  Arb}`.
- `external/probes/taylorseries-arb/probe.jl` — concrete usage pattern.
- `external/TaylorIntegration.jl/src/integrator/jetcoeffs.jl` lines
  23-43 — the canonical coefficient-generation pattern in Julia for
  reference; we adapt the same pattern but expose it as a standalone
  function rather than as a step inside an integrator.
- `DESIGN.md §4 Phase 3` — the 5-test test plan.

### The algorithm (sketch — read FW 2011 §2.1.2 for full)

For first-order ODE `dy/dz = f(z, y)`:

1. Construct `z_taylor = z0 + Taylor1(T, order)` (a Taylor1 representing
   `z0 + h` symbolically, with `order` terms).
2. Bootstrap: start with `y_taylor = Taylor1([y0]; order)` (constant
   term only).
3. For `k = 1:order`:
   a. Compute `f_taylor = f(z_taylor, y_taylor)`. The result is a
      `Taylor1` whose `k`-th coefficient is now correctly computable
      from `y_taylor`'s first `k` coefficients (already known).
   b. Set `y_taylor[k] = f_taylor[k-1] / k` (the integration step:
      if `y' = f`, then `y_k = f_{k-1}/k`).
4. Return `y_taylor.coeffs` (the `Vector{T}` of length `order+1`).

For second-order ODE `d²u/dz² = f(z, u, u')`:

Similar bootstrap, but with two unknowns `u` and `up = u'` and the
relation `up_k = u_{k+1} · (k+1)` (formal differentiation of the
Taylor series). Update both per pass.

### Test plan (DESIGN.md §4 Phase 3)

| test | input | expected | shape |
|---|---|---|---|
| 3.1.1 | `f(z, y) = y`, `(z₀, y₀) = (0, 1)`, order 10 | `c[k] = 1/k!` (exp(z) Taylor) | spec-from-scratch |
| 3.1.2 | `f(z, y) = z² + y²`, `(z₀, y₀) = (0, 0)`, order 14 | matches `Taylor1` of the Bessel solution `t · J_{3/4}(t²/2) / J_{-1/4}(t²/2)` to ball-equality at order 14 | port-and-verify (FW 2011 §2.2.1) |
| 3.1.3 | `f(z, u, u') = 6u² + z` (PI), tritronquée IC, order 30 | finite, no NaN/Inf, leading coefficients agree with `WeierstrassP`-derived expansion to 1e-12 (sanity only — PI ≠ ℘ at large `z`) | spec-from-scratch |
| 3.1.4 | `T = BigFloat(precision=256)`, same as 3.1.1, order 60 | radii < `2⁻²⁰⁰` per `RESEARCH.md §3.3` empirical | port-and-verify |
| 3.1.5 | mutation-proof | replace one rule in TaylorSeries (e.g. negate `exp`) → 3.1.1 fails | RED |

For 3.1.2 the Bessel oracle:
- `J_{3/4}(t²/2)` at small `t` has Taylor expansion you can compute
  via `SpecialFunctions.besselj` over a `Taylor1`. Or use
  `Symbolics.jl` to derive the first 14 coefficients symbolically.
- Cleaner alternative: just integrate `dy/dz = z² + y²` symbolically
  by hand for the first ~6 coefficients (it's tractable: y_0=0,
  y_1=0, y_2=0, y_3=1/3, y_4=0, …).

### What `Coefficients.jl` should look like

Stub at `src/Coefficients.jl`:
```julia
module Coefficients
function taylor_coefficients_1st end
function taylor_coefficients_2nd end
end # module Coefficients
```

Replace with:
```julia
module Coefficients
using TaylorSeries: Taylor1
using ..PadeTaylor: ... (none needed at this layer)

export taylor_coefficients_1st, taylor_coefficients_2nd

function taylor_coefficients_1st(f, z0::T, y0::T, order::Int) where {T}
    # See module docstring for algorithm.
    z = z0 + Taylor1(T, order)
    y = Taylor1([y0], order)
    for k = 1:order
        f_t = f(z, y)
        y[k] = f_t[k-1] / k
    end
    return y.coeffs
end

function taylor_coefficients_2nd(f, z0::T, y0::T, y1::T, order::Int) where {T}
    # u'' = f(z, u, u'); two unknowns evolved jointly.
    # ...
end
end
```

**Watch the type stability.** `y0::T` should not be coerced to
`Float64` accidentally; use `T(y0)` explicitly when constructing
`Taylor1`. The empirical probe at `external/probes/taylorseries-arb/`
already validated this pattern works for `T = Arb`.

### What to commit

After Phase 3 GREEN + mutation-proven:
```
Phase 3 GREEN: Coefficients.taylor_coefficients_1st/_2nd; 5 tests + mutation-proof

src/Coefficients.jl (~80 LOC):
- ...

test/coefficients_test.jl:
- 3.1.1 ...
- 3.1.2 ...

Mutation-proven: ...

Test summary: 66/66 GREEN (61 + 5 new).
```

## Specific instructions for Phase 4 (`StepControl`)

### Status — Phase 4 SHIPPED 2026-05-09

Phase 4 is complete (commit recording StepControl, three-source oracle,
DESIGN.md correction). 192/192 tests GREEN. The original spec text
below contained an incorrect formula; the corrections (and the rest
of the paragraph) are kept here as a record of the spec drift caught
at implementation time. **Read `src/StepControl.jl`'s top docstring +
`docs/worklog/002-phase-4-spec-correction.md` before changing this
module.**

### Goal (as shipped)

Two step-size strategies in `src/StepControl.jl`:
- `step_jorba_zou(coefs, ε_abs; ε_rel = ε_abs)` — Jorba-Zou 2005
  **§3.3.1 eq. 11** in the fixed-order form (TI.jl-equivalent), NOT
  §3.2 eq. 3-8. The `0.7`-safety factor cited in the original spec
  has no source.
- `step_pade_root(P, z_current, target)` — FW 2011 §3.1 forward-
  projection-onto-direction heuristic.

### Ground truth (verified)

- `references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/
  JorbaZou2005_taylor_IVP_package_ExpMath14.md:613-645` — §3.3.1 first
  step-size control, eq. 11.
- `external/TaylorIntegration.jl/src/integrator/stepsize.jl:17-89` —
  ported verbatim (`stepsize` + `_second_stepsize`).
- `external/probes/stepcontrol-oracle/{capture.jl, capture.py,
  capture.wl, verify.jl}` — three-source pin: TI.jl ≡ mpmath ≡
  wolframscript at 47 decimal digits on the canonical case.

### Polynomials.jl roots with Arb element type — still open

`Polynomials.roots` is validated for `Float64` / `Complex{Float64}`
only.  Arb element type for the Padé-root step is deferred to v2; see
the deferred bead `padetaylor-…` (file at session close per Rule 9).

## Phase 5 — SHIPPED 2026-05-09

Phase 5 (`PadeStepper.pade_step!`) is complete. Three-source oracle
(closed-form ℘ ≡ Mathematica `NDSolve` ≡ `mpmath.odefun`) at
<1e-15 relative agreement; mutation-proven on both `u` and `u'`
evaluation paths.  See `docs/worklog/003-phase-5-padestepper.md`
for the two algorithmic choices that crystallised at impl time
(`h^k` rescaling necessary near poles; analytic differentiation
of `P_u` for `u'` beats re-Padé by ~40×).

The Phase-5 `PadeStepper` now also exports
`pade_step_with_pade!(state, f, order, h) -> (state, P_u)` — used
by `solve_pade` for per-segment Padé storage.

## Phase 6 — SHIPPED 2026-05-09

Phase 6 (`Problems.PadeTaylorProblem`, `solve_pade`,
`PadeTaylorSolution`, `taylor_eval`) is complete on the pivoted v1
acceptance.  218/218 tests GREEN.  Three-source oracle (closed-form
`WeierstrassP` ≡ Mathematica `NDSolve` ≡ `mpmath.odefun`); BF-256
fourth source via `WeierstrassP[..., 80]` for the high-precision
test.  Mutation-proven on both the `robust_pade` call and the dense-
interpolation formula.

**Headline empirical result.**  At `z = 1.05` (just past the lattice
pole at `z = 1`), the Padé conversion of the local Taylor jet beats
plain Taylor truncation by ≈ 9.86 orders of magnitude on identical
input coefficients (rel-err `3.45e-10` vs `2.50`).  Tests 6.1.3 +
6.1.5 are the compelling pair — same jet, two paths, ten orders of
magnitude apart.

**Algorithmic finding caught at impl time.**  Test 6.1.6's BF-256
acceptance was bumped from `order = 30` to `order = 40`: at
`order = 30` both Float64 and BF-256 saturate at the same `(15, 15)`
Padé approximation-error floor (~3.5e-10 at `t = 0.7`), so BF-256
adds no value over Float64.  Only at `order = 40` does Float64
degrade to ~1.6e-8 (SVD ill-conditioning of the 41-by-41 Toeplitz
block) while BF-256 cleanly converges to ~5e-15 — making BF-256
required, not optional.  Full sweep + rationale in
`docs/worklog/005-phase-6-implementation.md` §"Algorithmic finding".

**Deferred to v2** (bead `padetaylor-8cr`):
- FW 2011 Table 5.1 long-range integration via the FW §3.1
  path-network and/or a step-size selector that composes with
  Padé-bridge stepping.  Failure analysis of fixed-`h_max` +
  pole-vault + Jorba-Zou-adaptive in
  `docs/worklog/004-phase-6-pivot.md` §"What was tried, what failed".

## Specific instructions for Phase 7 (`CommonSolveAdapter` ext)

This is your immediate next task.  Phase 7 is an optional package
extension wiring `PadeTaylorProblem` / `solve_pade` into the
`CommonSolve.jl` `init`/`step!`/`solve!` interface, so PadeTaylor.jl
plays cleanly with the wider `SciML` solver ecosystem.

**Read first** (in order):
1. `DESIGN.md §4 Phase 7` — Phase-7 acceptance + interface sketch.
2. `docs/adr/0003-*.md` — the ADR governing extension boundaries.
3. `src/Problems.jl` — the public driver this layer wraps.
4. `docs/worklog/005-phase-6-implementation.md` — the most recent
   shipped phase; orchestration pattern + mutation-proof discipline.

**Suggested workflow** (mirror Phase 3's "read first → test plan →
implementation skeleton" pattern from this HANDOFF lines ~96-225):

  - Stub `ext/CommonSolveAdapter/CommonSolveAdapter.jl` (a Julia
    weakdep extension on `CommonSolve`); update `Project.toml`'s
    `[weakdeps]` and `[extensions]` tables.
  - Sketch the test plan: at minimum, `init(prob; h_max) → integrator`,
    `step!(integrator) → (z, y)`, `solve!(integrator) → solution`,
    each cross-checked against the equivalent direct `solve_pade`
    call on the same problem.
  - Mutation-proof: drop the integrator's `pade` accumulator; assert
    the cross-check fails.

Phase 7 is **optional** per `DESIGN.md §4`: skip it if Phase 8
(`PadeTaylorArblibExt`) is the higher priority for downstream use.
The Phase-4 friction bead on `Polynomials.roots` over `Arb` element
type is a prerequisite for Phase 8 but unrelated to Phase 7.

## Friction beads worth re-reading before starting

```
bd ready -n 30
```

Run `bd ready -n 30` for the live list.  As of worklog 047 there are
no open P0/P1 beads; the ready queue is all P2/P3.  Notable open beads
worth knowing about:

**FFW 2017 11-step arc**: ✅ **ALL 11/11 SHIPPED** — A1–A6 + B1–B4 + B5{Fig 2, 3, 7} all closed.  No remaining work in the FFW arc itself.

**FFW v1 corner lifts**:
- `padetaylor-ykg` **P3 [task]** — `pIII_asymptotic_ic` higher-order coefficients (`a_3, a_4, ...`).  Closes the FFW Fig 5 `~10⁻²` v1 gap.
- `padetaylor-i76` **P2 [feature]** — `BVP.jl` 3-arg-RHS standing artefact from A6.
- `padetaylor-5k9` **P3 [feature]** — sheet-aware `solve_pole_free_hybrid` API (lifts Fig 5 to full FFW sector).
- `padetaylor-zwh` **P3 [feature]** — FFW Fig 1 widened to full `Re ζ ∈ [-2, 8]` via Poisson-disk Stage-1 nodes.
- `padetaylor-1r0` **P3 [chore]** — FFW Fig 6 Cartesian-resample backport (visual polish).
- `padetaylor-qum` **P3 [chore]** — `PathNetwork.jl` LOC-cap split (now at 1010 LOC after A4 + A5 + ADR-0015 additions, vs Rule-6 cap of 200).

**Other open P2s (pre-session, still open)**:
- `padetaylor-0tj` **P2 [bug]** — `lattice_dispatch_solve` throws on plain-IVP-corrupted BCs in smooth sectors.
- `padetaylor-7zw` **P2** — BF-256 tritronquée pin: rerun FW Fig 4.1 step-(i) BVP at `BigFloat`-256.
- `padetaylor-bhw` **P2** — expand quantitative figure-pinning test coverage.
- `padetaylor-61j` **P2** — Willers 1974 acquisition.
- `padetaylor-8pi` **P2** — GenericLinearAlgebra `svd!` piracy friction.

**Other open P3s**:
- `padetaylor-8dg` **P3** — `bvp_fill_with_asymptotic`: 1D per-row BVP using the algebraic asymptotic `u → -√(-z/6)` as the outer BC for half-bounded smooth sectors.
- `padetaylor-fe9` **P3** — `smooth_fill_2d`: 2D elliptic Laplace BVP for analytic smooth-region fill.
- `padetaylor-p3l` **P3** — quantitative pole-count pin for FW Fig 4.3/4.4.
- `padetaylor-zzw` **P3** — performance profiling pass.

**Next work item**: FFW 11-step arc complete.  Pick from the remaining open P2/P3 beads above — most natural candidates are (a) **`padetaylor-ykg`** (lift the FFW Fig 5 `~10⁻²` v1 gap by adding `pIII_asymptotic_ic` higher-order coefficients), which unblocks both `padetaylor-i76` and `padetaylor-5k9`; (b) **`padetaylor-zwh`** (widen FFW Fig 1 to FFW's full `Re ζ ∈ [-2, 8]` via Poisson-disk Stage-1 nodes — would also help Fig 3 / Fig 7 / Fig 2 reduce their rng_seed fragility); or (c) **`padetaylor-qum`** (split `PathNetwork.jl` for Rule 6 LOC cap — now at ~1010 LOC).  Alternatively pick up the multi-branch-point routing surfaced by Fig 7 (worklog 047 §"multi-branch-point investigation") and convert the architecture-ready-but-unexercised `ζ = ±2πi` handling into a properly tested capability.

## Hard-won lessons (don't repeat these)

These come from the worklog shard `docs/worklog/001-stages-Z-1-2-handoff.md`.
Read that file in full before Phase 3.

Quick summary:
1. **`bd init` rejects dotted DB names.** `--prefix=padetaylor
   --database=padetaylor` is required for projects ending in `.jl`.
2. **`GenericLinearAlgebra.jl` pirates into `LinearAlgebra.svd!`.**
   Naïve `svd → LinearAlgebra.svd` mutation does not bite. Proper
   mutation for any "BigFloat dispatch" claim is `Matrix{Float64}(A)`
   downcast.
3. **Octave `i` ≠ Julia `im`** for the imaginary unit. Use
   `%.18eim` in `fprintf` format strings, not `%.18ei`.
4. **Mid-precision Padé table reduction is precision-dependent.**
   exp(z) (20,20) reduces to (7,7) at Float64 but stays at (20,20)
   at BigFloat-256. This is correct behaviour, not a bug.
5. **GGT 2013 §7 ill-posedness is real and visible empirically.**
   For log(1.2-z) (10,10), Octave-vs-Julia coefficients diverge at
   ~1e-3 *while the rational function values agree to 1e-15*. Test
   for function-value match, not per-coefficient match, in
   ill-posed-block cases.
6. **Test 1.1.1's mutation-proof procedure was wrong on first
   attempt.** Documented the *correct* mutation in the test comments.
   When you write Phase 3 mutation-proofs, verify the mutation
   actually bites before committing the procedure.
7. **Empirically sweep `order` before pinning `rtol` on near-pole
   tests.** The `(m, n)`-Padé approximation-error floor depends on
   `m + n` and the distance to the next singularity, NOT on
   arithmetic precision; BF-256 only helps once you push past
   Float64's SVD-conditioning limit (typically `order ≥ 40`).
   Phase-6 test 6.1.6 hit this: at `order = 30` Float64 and BF-256
   both saturated at `~3.5e-10` (the `(15, 15)` Padé floor at
   `t = 0.7`); only at `order = 40` did BF-256 buy two extra
   orders of accuracy.  See `docs/worklog/005-phase-6-implementation
   .md` §"Algorithmic finding".

8. **Step-norm Newton convergence, NOT residual-norm, for spectral
   BVPs.**  The discrete residual `‖R‖_∞` floors at `cond(D₂) ·
   eps(T) ≈ N² · eps(T)` — an irreducible "spectral truncation"
   floor.  A residual-norm Newton tolerance below that floor is
   structurally unachievable; convergence loops forever / errors.
   The right choice is step-norm `‖Δu‖_∞ ≤ tol` with default
   `eps(T)^(3/4)` (production-Newton standard, NLsolve.jl-style).
   Phase-11 hit this on first cut; default was `100·eps(T) ≈ 2.2e-14`
   while PI N=20 floors at ~1e-12.  See
   `docs/worklog/006-phases-10-11-path-network-bvp.md` §"Step-norm
   Newton convergence".

9. **Oracle subagents can catch spec drift you didn't write.**
   Phase-11's Octave oracle subagent caught that the existing
   Phase-6 problems-oracle's `u_at_0_5 = 4.0044…` was for
   `u'' = 6u²` (equianharmonic ℘) NOT `u'' = 6u² + z` (Painlevé I).
   Recomputed via mpmath at 40 dps, cross-validated to ~1e-10.
   Lesson: oracle files should state their ODE prominently;
   `test/_oracle_problems.jl` line 1 should add a header comment
   distinguishing PI from equianharmonic ℘.  Recorded in worklog 006
   §"Spec-drift catch by the Octave oracle subagent".

10. **`marker_single` lives in `~/Projects/<project>/.venv/bin/`,
    not on system PATH.**  Same for many other Python tools (mpmath
    via uv, etc.).  Scan venvs before declaring a tool absent.  Saved
    as user memory `user-python-tools-venv-pattern`.

11. **TIB VPN gives access to paywalled SIAM/ACM/Elsevier journals.**
    Phase-11 acquired Weideman-Reddy 2000 (ACM TOMS) + Berrut-Trefethen
    2004 (SIAM Review) directly from publishers.  When user activates
    the VPN, attempt direct publisher URLs before preprint fallback.
    Saved as reference memory `reference-tib-vpn-paywalls`.

12. **ADR-specified invariants need invariant tests, not approximation
    tests**.  PN.2.2 surfaced a Phase-10 bug where the visited-node
    canonical-Padé invariant (specified in ADR-0004 + worklog 006) was
    NOT enforced in code — only the FW Table 5.1 long-range test at
    z=30 caught it (~22% rel-err Float64 AND BF-256).  No test asserted
    the invariant directly (e.g., "evaluating visited[k]'s stored Padé
    at t=0 returns visited_u[k]"); such a test would have caught the
    bug immediately at Phase-10 GREEN.  See worklog 008 §"What changed"
    and §"Hard-won lessons".

13. **When numerical error is invariant to precision, suspect the
    algorithm**.  PN.2.2's bug showed identical 22% rel-err at Float64
    AND BF-256.  Float64-rel-err == BF-256-rel-err is a strong signal
    that the issue is structural, not roundoff/truncation.  A
    pre-emptive BF-256 probe is worth the ~50× wall-time cost when a
    result smells wrong.  See worklog 008 §"Frictions surfaced" F3.

## Last commit before this handoff

`HEAD = 1331755` ("B5 Fig 7 SHIPPED — generic PVI in η + ζ + z planes (bead padetaylor-mgx)").  Two commits this session, closing the FFW 11-step arc at 11/11: B5 Fig 3 PVI phase portraits (`55fd533`, worklog 046, bead `padetaylor-a1l`) → B5 Fig 7 generic PVI three-frame (`1331755`, worklog 047, bead `padetaylor-mgx`).  Test suite **2467 → 2508 GREEN** (+61 assertions this session; +20 FF3.* + +41 FF7.*).  Both commits land on `main`; push pending.

### Session 2026-05-16 (continuation 2) — FFW B5 closeout: Figs 3 + 7 (worklogs 046–047)

Orchestrated session: parent agent coordinates, opus subagents do plan + implementation per figure; one Julia process at a time per CLAUDE.md Rule 7.  Two commits, two worklogs, **+61 GREEN assertions**, **FFW 11-step arc complete (11/11)**.

  - **B5 Fig 3 — PVI phase portraits** (`55fd533`, worklog 046, bead `padetaylor-a1l` **closed**) — second of three B5 figures.  Three columns × two rows (sheet [0] top, alt sheet bottom): `|w(ζ)|` modulus + `arg w(ζ)` HSV phase + `arg u(z)` HSV phase on z-plane Cartesian-resample (Fig 4 pattern lifted).  Reuses Fig 2's PVI problem `(α,β,γ,δ)=(4,-4,8,-8)`, IC `u(10)=0.4295…` with a FRESH cross-mode A4 solve over Re ζ ∈ [-1.5, 2.5] (vs Fig 2's [0.05, 3.0]).  526 visited tree nodes; sheet distribution `[0] => 477, [1] => 49`; **100% Stage-2 coverage** on both panels.  Three mutations: M1 (corrupt α) 15 ERROR (walker destabilises), M3 (drop cross_branch) 15 ERROR, M4 (drop sheet bookkeeping) 5 RED + 1 ERROR.  Honest v1 corners: reduced z-window vs FFW `[1/100, 10]`, rng_seed=2 hard-coded for sheet [1] reproducibility, no explicit cut-line annotation.  **Zero new `src/` code** — pure composition.

  - **B5 Fig 7 — generic PVI in η + ζ + z planes** (`1331755`, worklog 047, bead `padetaylor-mgx` **closed**) — closes the FFW arc.  First GENERIC (non-tronquée) PVI figure.  Two solves: η-frame via A3 inside `Re η < log 2π` branch-point-free region; ζ-frame via A4 `cross_branch=true` + winding ring + A2 `node_separation = R(ζ) = max(0.05, (4 - Re ζ)/10)`.  z-plane column rendered as `log10|u(z)|` (FFW md:281: z-residues scale with |z₀|; log compresses dynamic range).  **First figure to exercise multi-branch-point routing**: `branch_points = (0im, +2πi, -2πi)` — confirmed `BranchTracker`'s tuple API supports the multi-point case (sheet counter is a length-3 `Vector{Int}` per node; ±2πi counters stay 0 in this window but architecture is ready).  Five mutations: M1 (α) 3 FAIL, M3 (cross_branch) 1 ERROR + cascade, M4 (sheet bookkeeping) 4 FAIL + 1 ERROR, M5 (log10 → bare abs in render kernel) 7 FAIL on the load-bearing FF7.1.8 invariant.  **M2 (drop node_separation) is HONEST NON-BITE** (Rule 9 documented): at this window/density `:adaptive_ffw` alone shrinks h sufficiently; the kwarg is figure-side ergonomic only.  rng_seed=1 picked from a 5-seed sweep (16 alt-sheet nodes).  Honest v1 corners: single principal-sheet η-strip (FFW shows 3); reduced z-window `|z| ∈ [0.45, 4.48]`; `ζ=±2πi` circumambulation architecture-ready but unexercised at this window.  41 GREEN assertions.  **Zero new `src/` code**; new `figures/_ffw2017_fig_7_helpers.jl` (~75 LOC + docstring) extracts the Makie-free render kernel so the test exercises the EXACT figure code path (M5 mutation bites the kernel, not a parallel copy).

**Open beads filed this session**: none new.  Both `padetaylor-a1l` and `padetaylor-mgx` **closed**.

**FFW 11-step arc final status**: ✅ **11/11 done** — A1 adaptive_ffw; A2 node_separation; A3 η-plane PVI; A4 constrained-wedge + sheet bookkeeping; A5 sheet-aware Stage-2; A6 IVP+BVP hybrid; B1 Fig 6; B2 Fig 1; B3 Fig 4; B4 Fig 5; B5 Figs 2 + 3 + 7.  Plus the user-flagged ADR-0015 `extrapolate=true` fix and the multi-branch-point capability surfaced by Fig 7.

**Hard-won lessons added this session**:

  42. **`node_separation` doesn't always bite when expected.** Plan predicted M2 would crash the walker on dense generic poles; reality is `:adaptive_ffw` alone handles it at our window/density (1536 visited nodes without `node_separation` vs 205 with).  `node_separation` becomes load-bearing only at wider windows or when the adaptive controller would need many rejection passes per step.  At Fig 7's scale it's figure-side ergonomic (faster solve, fewer nodes), not algorithmically required.  Documented in worklog 047 §"Mutation M2 honest non-bite."  Lesson: a v1 walker passing without an opt-in kwarg is not evidence the kwarg is unnecessary in general; it's evidence the v1 scope didn't reach the regime where the kwarg becomes load-bearing.

  43. **Multi-branch-point `branch_points` tuple has worked all along.** Fig 7 passed `(0im, +2πi, -2πi)` as a 3-tuple; `BranchTracker` correctly resolves a length-3 `Vector{Int}` sheet counter per node.  No `src/` change needed — the capability was there since A4 shipped (worklog 042).  Worth knowing for future PVI work where the `ζ = 2πi·k` lattice points become reachable from the walker (currently they don't for tronquée Figs 2/3 or for the limited-window generic Fig 7).

  44. **Render kernel extraction enables true rendering invariant tests.** Fig 7's FF7.1.8 (`log10|u|` rendering invariant) required the figure script and test to share the EXACT same code path so that an M5 mutation (`log10 → abs`) bites the actual code, not a parallel copy.  Solution: extract `figures/_ffw2017_fig_7_helpers.jl` as a Makie-free pure-Julia helper module included by both.  This is a precedent worth lifting for future figures with load-bearing rendering invariants.

### Session 2026-05-16 — A3/A4/A5 + B5 Fig 2 + extrapolate Stage-2 (worklogs 041–045)

Closes 7 of the remaining FFW-arc steps + a user-flagged rendering fix.  Five commits, two ADRs (0013, 0015), five worklogs (041–045), +227 GREEN assertions.

**FFW infrastructure shipped this session** (3 of 6 remaining A-steps):
  - **A3 η-plane PVI** (`acb7989`, worklog 041, bead `padetaylor-riu` **closed**) — popped previous session's WIP stash; reverted the M2 sign-flip marker at `src/SheetTracker.jl:248` (`vp = -z*ζ*up` → `vp = z*ζ*up`); re-verified all three mutations (M1: 4 RED, M2: 6 RED, M3: 2 RED) bite distinct subsets.  `pVI_eta_transformed_rhs` (FFW md:154 verbatim) + `pVI_z_to_η`/`pVI_η_to_z` coordinate helpers + `:transformed_eta` frame on `PainleveProblem(:VI)`.  +42 ET.* assertions.  ADR-less additive extension.
  - **A4 constrained-wedge routing + sheet bookkeeping** (`09217af`, worklog 042, ADR-0013, bead `padetaylor-bho` **closed**) — the largest remaining infra piece.  New module `src/BranchTracker.jl` (~158 LOC): `segment_crosses_cut` (2×2 parametric ray intersection via Cramer on imag parts), `any_cut_crossed`, `step_sheet_update` (per-step counter bump via `sign(winding_delta)`), `resolve_cut_angles` (scalar broadcast).  `path_network_solve` gains `branch_points` + `branch_cut_angles` + `cross_branch` + `initial_sheet` kwargs + new `visited_sheet` field on `PathNetworkSolution` (with backward-compat 9-arg constructor).  Refuse mode rejects any wedge candidate crossing a cut; cross mode allows crossings + bumps per-branch sheet counter.  6 mutations (B1, B2, B3 in BranchTracker; PB1, PB2, PB3 in PathNetwork wiring) each bite distinct subsets (1-52 RED).  +105 BT.* + PNB.* assertions.  ADR-0013 documents: half-line cut model, per-branch sheet tuple (FFW md:188), helper-module split (PathNetwork.jl now at 818 LOC; `padetaylor-qum` bead more pressing).
  - **A5 sheet-aware Stage-2** (`c685f7b`, worklog 043, bead `padetaylor-hed` **closed**) — ADR-less additive extension to A4.  `grid_sheet::Union{Nothing, Vector{Vector{Int}}}` kwarg on `path_network_solve` restricts Stage-2 lookup to matching-sheet visited nodes; fail-soft NaN when no match.  New public accessor `eval_at_sheet(sol, z, sheet) -> (u, up)` for post-hoc per-point queries.  New internal `_nearest_visited_on_sheet` helper.  3 mutations (S1, S2, S3) bite distinct subsets (4/6/2 RED).  +30 SA.* assertions.

**FFW figure shipped this session** (1 of 3 remaining B5 figures):
  - **B5 FFW Fig 2 — PVI three-method reproduction** (`1313297`, worklog 044, bead `padetaylor-n93` **closed**) — first end-to-end demo of the complete A1+A2+A3+A4+A5+A6 stack on a single problem.  Three columns of the tronquée P_VI solution `(α,β,γ,δ)=(4,-4,8,-8)` at IC `u(10)=0.4295…`: (1) η-plane via A3; (2) ζ-plane sheet 0 refuse mode via A4; (3a + 3b) ζ-plane two-sheet cross mode via A4+A5 with `eval_at_sheet`.  Sparse-Stage-1 + manual-Stage-2 pattern (lifted from Fig 1) avoids the walker-target/render-lattice trap.  19 assertions; mutations M1 (corrupt α): 2 RED, M3 (drop cross_branch): 1 ERROR, M4 (drop sheet bookkeeping): 1 RED.  Honest v1 corners: η-window single-strip (FFW shows two); alt-sheet coverage 0.9% in col 3b pre-fix.

**User-flagged fix** (followed up after Fig 2 ship):
  - **`extrapolate=true` Stage-2 kwarg** (`6dad724`, worklog 045, ADR-0015, bead `padetaylor-afs` **closed**) — root cause of the "lots of white gaps in all FFW 2017 figures" report: our `PathNetwork.jl:596` Stage-2 fail-soft (NaN when fine-grid cell > `visited_h` from nearest visited node, ADR-0004) diverged from FFW md:62's spec (evaluate Padé regardless of disc).  Opt-in `extrapolate::Bool = false` kwarg on `path_network_solve` + `eval_at_sheet`; new public `eval_at(sol, z; extrapolate=false)` sheet-blind accessor (the sibling of `eval_at_sheet`).  Default preserved per ADR-0015's Rule-1 contract.  3 mutations probed; X2 + X3 bite 6 + 1 RED.  X1 documented as untestable in normal `path_network_solve` geometry (walker visits every grid cell → Stage-2 NaN path is dead code).  +30 SX.* assertions.  Four FFW figure scripts updated (Fig 1, 2, 4, 6) — coverage went from 0.9-98% to **100% everywhere**.  Fig 5 unchanged (its NaN is genuine "outside BVP sector" from `solve_pole_free_hybrid`).
  
  **Friction noted in worklog 045**: initially spawned 4 parallel julia processes for the figure re-renders, violating CLAUDE.md Rule 7.  User flagged immediately; killed all four via TaskStop and re-ran sequentially.  Lesson: even with warm precompile cache, multiple concurrent `using PadeTaylor` calls can race; Rule 7 is the safe default.

**Open beads filed this session**:
  - `padetaylor-a1l` **P2 [feature]** — FFW Fig 3 (PVI phase portraits).  Blocked-by Fig 2 (now unblocked).
  - `padetaylor-mgx` **P2 [feature]** — FFW Fig 7 (generic PVI in η/ζ/z planes).  Independent of Fig 3.
  - `padetaylor-hed` (A5) **closed** in same commit.
  - `padetaylor-n93` (Fig 2) **closed** in same commit.
  - `padetaylor-afs` (extrapolate) **closed** in same commit.

**11-step FFW arc status**: A1-A6 + B1-B4 + B5(Fig 2) = **8 of 11 done**.  Remaining: B5 Fig 3 (gated on Fig 2 → ready) and B5 Fig 7 (independent).  After those, the 11-step plan is complete.

### Session 2026-05-15 (continuation) — FFW 2017 figure arc (worklogs 034–040)

A long orchestrated session: 11-step plan to reproduce all seven FFW 2017 figures.  **Phase A** (core infrastructure) shipped 3 of 6 planned steps; **Phase B** (figure reproduction) shipped 4 of 7 figures.  Seven commits, four ADRs (0011, 0012, 0014), seven worklogs (034–040), +512 GREEN assertions.

**Phase A shipped this session** (3 of 6):
  - **A1 adaptive Padé step** (`7c7cac7`, worklog 034, ADR-0011, bead `padetaylor-8ui` **closed**) — `:adaptive_ffw` step_policy on `path_network_solve` implementing FFW 2017 §2.1.2 (md:74-97) truncation-error controller `T(h) = |ε_{n+1}·h^{n+1}/a(h)|` + Newton rescale `q = (k·Tol/T(h))^(1/(n+1))`.  Three helpers in `src/PadeStepper.jl`; inline rescale loop in `path_network_solve` honouring FFW md:93's "scaled step length stored at the current point" + "re-run wedge, re-select min-|u|" semantics.  `:fixed` default byte-unchanged.  Mutation-proven M1-M4; +43 AS.* assertions.
  - **A2 non-uniform Stage-1 nodes** (`eb2a64b`, worklog 035, ADR-0012, bead `padetaylor-1a3` **closed**) — opt-in `node_separation::Function` kwarg on `path_network_solve`.  When supplied, each Stage-1 wedge step uses magnitude `R(ζ_current)` instead of fixed `h`; composes orthogonally with `:adaptive_ffw` (R seeds; controller may only shrink, not grow).  Backward-compat: `nothing` default reproduces the existing uniform tree byte-equally.  Mutation-proven M1-M4 (M3 off-by-one + M4 ignore-R both bite NU.1.7); +295 NU.* assertions.  Friction `padetaylor-qum` filed for `PathNetwork.jl` LOC-cap split (P3 chore, effective LOC 331 vs Rule 6's 200).
  - **A6 IVP+BVP hybrid driver** (`c41eaa2`, worklog 039, ADR-0014, bead `padetaylor-0co` **closed**) — new module `src/IVPBVPHybrid.jl` (~480 LOC) implementing FFW 2017 §3 (md:203-247) algorithm: run PFS along curved sector boundaries from asymptotic-IC anchors → harvest `u, u'` at discrete boundary points → call `bvp_solve` on the sector → glue.  New positional overload `bvp_solve(f, ∂f_∂u, ∂f_∂up, ...)` for 3-arg RHS required by P̃_III's `(w')²/w` term (BVP API extension additive only; 2-arg byte-unchanged; friction bead `padetaylor-i76` opens to clean up).  Asymptotic-series first coefficient `a_1 = -β/3` derived inline.  v1 corner documented in ADR §"Out of scope": at Float64 + n_terms=2 + N=10 the joint regime sits at FFW md:226's `~10⁻¹` (BVP self-consistency is `O(10⁻¹²)`; the asymptotic-IC error dominates).  Mutation-proven M1-M4; +66 IB.* assertions.

**Phase B shipped this session** (4 of 7 figures):
  - **B1 FFW Fig 6 PV generic three-sheet** (`803876e`, worklog 036, bead `padetaylor-mbu` **closed**) — first FFW 2017 figure reproduced.  Three panels × three sheets = nine.  Per-sheet errors 3.3e-12 / 7.0e-8 / 9.8e-8 vs FFW md:297's 3e-9 / 7e-6 / 2e-5 — **beats FFW by 2-3 orders of magnitude** at Tol=1e-10.  PoC validation that A1+A2 compose end-to-end.  Friction: scatter-on-polar-grid rendering produces visible spoke artifacts in the z-plane panels; backport bead `padetaylor-1r0` filed (P3 chore, lift Fig 1's Cartesian-resample helper).
  - **B2 FFW Fig 1 PIII three-sheet spiral** (`3afae14`, worklog 037, bead `padetaylor-2pg` **closed**) — FFW's headline figure.  Four panels (Stage-1 nodes, Stage-1 tree edges via `visited_parent`, ζ-plane `|w(ζ)|`, z-plane Riemann sheets with spiral pole pattern).  **Cartesian-resample bilinear-interpolation rendering** for the z-plane (deprecates B1's polar scatter — Panel C shows the FFW md:283 "oblique pole lines" beautifully).  2664 visited nodes vs FFW's 2701 (-1.4%); sheet-0 conjugate-symmetry median `4.05e-15` vs FFW Exp-2 `1e-6` (9 orders better).  Ran on reduced `Re ζ ∈ [-2, 5]` window; widening to FFW's `[-2, 8]` needs Poisson-disk node placement (Fornberg 2015 ref [8]) — friction bead `padetaylor-zwh` filed (P3 feature).
  - **B3 FFW Fig 4 PV tronquée three-sheet** (`2a37c40`, worklog 038, bead `padetaylor-pt2` **closed**) — three columns × three sheets.  Pole-free sector clearly visible on sheet 0 (smooth dark-blue plain in `|u|` panel, consistent with FFW md:209 asymptote `u → -1` for `-π < arg z < π`).  Pole-density gradient 8× across sheets 1+2.  Same four-class assertion template as B1/B2 (IC pin + sheet pins + sector qualitative + density gradient).
  - **B4 FFW Fig 5 PIII tronquée + cond-number heatmap** (`f05f13e`, worklog 040, bead `padetaylor-3ug` **closed**) — uses A6's hybrid driver.  Two panels: Panel A `|w(ζ)|` over principal-sheet sub-strip via `solve_pole_free_hybrid`, Panel B closed-form `log₁₀ κ_r` over full FFW sector.  Cond-number pin `κ_r(z=30) = 157.30` matches FFW md:264 to ±1.  v1 corner: per-strip error `~10⁻²` is 4 orders looser than FFW Table 2's `4e-6 / 3e-7 / 2e-8 / 3e-9 / 4e-8`; closing the gap is gated on bead `padetaylor-ykg` (asymptotic-series helper `a_3+` coefficients; the v1 helper hard-codes them to zero so n_terms ≥ 3 is currently a no-op).  M1/M3/M4 all bite; M2 (n_terms cap) is a v1 no-op as expected.

**Open beads filed this session**:
  - `padetaylor-qum` **P3 [chore]** — `PathNetwork.jl` split into helpers (Rule 6 LOC cap at 331 vs 200).
  - `padetaylor-zwh` **P3 [feature]** — FFW Fig 1 widened to full `Re ζ ∈ [-2, 8]` via Poisson-disk Stage-1 node placement.
  - `padetaylor-1r0` **P3 [chore]** — FFW Fig 6 Cartesian-resample backport (visual polish, deprecate spoke artifacts).
  - `padetaylor-i76` **P2 [feature]** — `BVP.jl` 3-arg-RHS standing artefact from A6.
  - `padetaylor-ykg` **P3 [task]** — `pIII_asymptotic_ic` higher-order coefficients via TaylorSeries.jl substitution; gates closing the FFW Fig 5 `~10⁻²` v1 gap.
  - `padetaylor-5k9` **P3 [feature]** — sheet-aware `solve_pole_free_hybrid` API taking `asymptotic_ic_fn(ζ)` instead of `(z)`, for full-FFW-sector multi-sheet rendering.

**What is NOT shipped from the 11-step plan** (4 of 11 steps remaining):
  - **A3 η-plane PVI RHS** (bead `padetaylor-riu` open) — needed for FFW Fig 2 column 1 and Fig 7 column 1.  **A3 was dispatched mid-session and stopped before completion**; the agent wrote `pVI_eta_transformed_rhs`, `pVI_z_to_η`, `pVI_η_to_z` (~150 LOC into `src/SheetTracker.jl`), wired `:transformed_eta` frame into `Painleve._build_VI`, drafted `test/eta_pvi_test.jl`, and was applying mutation M2 when stopped — leaving a deliberate `# M2 mutation: sign flip` comment at `src/SheetTracker.jl:248` (`vp = -z * ζ * up`).  **WIP is stashed** via `git stash push -u -m "A3 η-plane PVI WIP (stopped mid-mutation-test; has M2 marker)"`.  Next session: `git stash pop`, revert that M2 line back to `vp = +z * ζ * up`, finish mutation-proof, commit.
  - **A4 constrained-wedge routing** (bead `padetaylor-bho` open) — branch-cut-respecting wedge selector + sheet-counter increment on circumambulation.  Needed for FFW Figs 2/3/7.  Significant change to `PathNetwork.jl`; new ADR-0013 expected.
  - **A5 sheet-aware Stage-2 evaluation** (no bead yet) — `PathNetworkSolution(ζ; sheet=k)` accessor restricting barycentric pool to matching `visited_sheet`; cross-sheet queries fail-loud.
  - **B5 FFW Figs 2, 3, 7 (PVI multi-sheet)** (no beads yet) — render the three remaining PVI figures.  Requires A3 + A4 + A5.

### Session 2026-05-15 — EdgeDetector h-aware level fix + closed-form Painlevé families (worklogs 032–033)

Two beads close.  Six new files, two ADRs, +78 GREEN assertions (1630 → 1708), one new headline render at 4× the previous resolution.

  - **`EdgeDetector` h-aware default level** (`905762c` + `17c284c`, worklog 032, bead `padetaylor-f8l` **closed**) — `src/EdgeDetector.jl` + `src/EdgeGatedSolve.jl` + `src/LatticeDispatcher.jl` + ADR-0009.  The 5-point Laplacian residual of an analytic function scales as `|Δu| ∝ h²`, so in `log₁₀` the smooth-region floor AND the pole-annulus residual the flood-fill rides on BOTH shift by `2·log₁₀(h)` as the grid refines.  The fixed `level = 0.001` default (FW 2011 §3.2.2's contour at their single resolution) over-classified as smooth on fine grids: below `h_grid ≈ 0.25` the connecting annulus dropped under `10^level`, the flood-fill lost connectivity, and `edge_gated_pole_field_solve` stalled at the seed.  Fix: `level = :auto` resolving to `LEVEL0 + 2·log₁₀(min(h, H0)/H0)` with anchor `(H0, LEVEL0) = (0.25, 0.001)`; clamped at `H0` so `h ≥ H0` keeps FW's published `0.001` (preserves backward compat for existing test grids at `h ∈ {0.333, 0.5, 0.75}`).  Mutation M3 (drop the clamp) cascades into RED on the existing EG.1.* tests at `h_grid = 0.75` — proves the clamp is load-bearing, not aesthetic.  Threaded `:auto` through both wrapper modules.  **Headline render**: `examples/tritronquee_3d.jl` re-runs at N=641 / `h_grid = 0.0625` (the 4×-resolution README hero the bead blocked).  92 region-growing passes converged on a 75,646-cell pole field (18.41 % of the 410,881-cell grid — exactly the natural extent of the PI tritronquée's positive-real-axis wedge), 79-min wall.

  - **Closed-form Painlevé families** (`2b3c3cf`, worklog 033, bead `padetaylor-icf` **closed**) — `src/PainleveClosedForm.jl` + ADR-0010.  Three new constructors extending the ADR-0008 named-transcendent surface with the second of three taxonomy kinds (point-IC was kind 1 in worklog 031):
      - `pii_rational(n)` for `n ∈ {1, 2, 3}` — PII rational solutions `u_n` at `α = n` per FW 2014 eq. 6.
      - `pii_airy(n; θ = 0.0)` for `n ∈ {0, 1}` — PII Airy solutions `u_{n+½}` at `α = n + 1/2` per FW 2014 eqs. 8-9 (`θ` parameterises the one-parameter family per `n`).
      - `piv_entire(:minus_2z)` and `piv_entire(:minus_two_thirds_z)` — PIV entire solutions per RF 2014 md:91 at `(α, β) = (0, -2)` and `(0, -2/9)` respectively (RF 2014 cites the solutions but not the `(α, β)` pairs; derived algebraically inline in `_u_piv_entire`'s docstring).
    The IC is derived from the analytic closed form at `zspan[1]`; the same formula is the *oracle* for self-validating tests (exact `==` on the IC tuple + `≈` solver-vs-formula cross-check at downstream `z`).  Mutation-proven with four mutations (MC1 sign flip on `u_1`, MC2 wrong `(α, β)` for `piv_entire`, MC3 drop pole guard, MC4 Airy chain-rule sign flip) — each bites at least one CF.* assertion.  New runtime dep: `SpecialFunctions.jl` for the Airy `Ai/Bi`.  Deferred to follow-up (not yet a bead): `pii_airy(2)` (u_{5/2} hand-derivation), `pii_rational(n ≥ 4)` via Bäcklund recurrence, PIV Generalized Hermite/Okamoto rationals (Noumi-Yamada ref not in-tree).

  - **Three new beads filed** mid-session for the smooth-region BVP-fill arc (surfaced by a "can we BVP the smooth sectors?" smoke test that errored out on plain-IVP-corrupted BCs):
      - `padetaylor-0tj` **P2 [bug]** — `lattice_dispatch_solve` robustness on corrupted BCs.  Use edge-gated as the IVP source; fail-soft on BVP non-convergence with a `:bvp_fail` `region_tag`.
      - `padetaylor-8dg` **P3** — `bvp_fill_with_asymptotic`: 1D per-row variant using the algebraic asymptotic as outer BC.
      - `padetaylor-fe9` **P3** — `smooth_fill_2d`: 2D elliptic Laplace BVP for the analytic smooth-region fill (since `∇²u = 0` for analytic `u`, it's a single sparse LU per smooth component).  The principled alternative to `-8dg`.

### Session 2026-05-14 (continuation) — `PoleField`, `EdgeGatedSolve`, FW Fig 4.1/4.7/4.8, `PainleveSolution`, named constructors (worklogs 027–031)

Continues the figure-reproduction + Painlevé-infrastructure arc.  Five worklogs, six beads' worth of code, two ADRs (0007, 0008).  Test suite **1526 → 1630 GREEN**.

  - **`PoleField.extract_poles`** (`ace8db4`, worklog 027, bead `padetaylor-xvf`) — `src/PoleField.jl`: reads pole locations back out of a solved `PathNetworkSolution` / `PadeTaylorSolution` by rooting the stored Padé denominators, mapping roots to the `z`-plane, and clustering with a cross-node-support filter.  Verified vs the equianharmonic ℘ lattice to ~1e-8.
  - **`EdgeGatedSolve` + FW Fig 4.7/4.8** (`225b7f2`, worklog 028, beads `padetaylor-dmb` + `padetaylor-7na`) — `src/EdgeGatedSolve.jl`: region-growing pole-field solver (dilate → solve → classify → morphological-open → flood-fill → fixpoint).  Fixes the spurious-pole bloom of plain `path_network_solve` on solutions with large smooth sectors (the PI tritronquée's four empty sectors).  FW Fig 4.7 (6 panels) and 4.8 reproduced.
  - **FW Fig 4.1 full composition** (`2356b80`, worklog 029, bead `padetaylor-gky`) — `figures/fw2011_fig_4_1.jl`: FW's 3-step recipe end-to-end — (i) imaginary-axis `bvp_solve`, (ii) two `edge_gated_pole_field_solve` run-outs, (iii) per-row `bvp_solve` smooth-band fill.  All 119 step-(iii) BVPs converged, zero fallbacks.
  - **`PainleveSolution` wrapper + `:transformed` flip + Makie recipe** (`adb4425`, worklog 030, beads `padetaylor-p1b` + `padetaylor-26r` + `padetaylor-ylr`; `padetaylor-soi` closed as subsumed) — `src/PainleveSolution.jl` + ADR-0007: a self-describing solve-output wrapper with a uniform `z`-frame access surface; makes the `:transformed` equations (PIII/PV/PVI) work end-to-end through `path_network_solve`.  `PadeTaylorMakieExt` ships the `painleveplot` recipe.
  - **Named transcendent constructors** (`705674a`, worklog 031, bead `padetaylor-c86`; `padetaylor-icf` filed P3) — `src/PainleveNamed.jl` + ADR-0008: `tritronquee(:I)` and `hastings_mcleod()` with literature-pinned 16-digit ICs baked in.  Worklog 031's hard-won lesson: a literature subagent's job is to *sort*, not just *find* — the three-kinds taxonomy is what made the v1 scope defensible.

  **Doc round (this session)**: README.md rewritten as an engaging, progressively-explained, beginner-targeted document with embedded figures; `CHANGELOG.md` `[Unreleased]` section filled in; this HANDOFF refreshed for worklogs 027–031.

### Session 2026-05-14 — FW figures 3.1–5.2, ADR-0006, `PainleveProblem` layer, `9xf` fix (worklogs 022–026)

A long session.  Eight FW 2011 figures reproduced, one ADR written + accepted, two beads' worth of code shipped.

  - **FW 2011 figure reproduction** — new `figures/` project (own `Project.toml`, CairoMakie + dev-ed PadeTaylor; `figutil.jl` shared helper).  Shipped: Fig 3.1 + 3.2 (worklog 022), Fig 3.3 (worklog 023), Fig 4.2 + 4.3 + 4.4 (worklog 024), Fig 5.1 + 5.2 (worklog 025).  Output PNGs committed.  Fig 5.2 *meets* its quantitative catalogue acceptance (sweep min at `(order,h)=(30,0.40)`, min rel-err `9.1e-15`).
  - **`PathNetworkSolution.visited_parent`** added (worklog 022) — Stage-1 tree edges, needed for Fig 3.2; test `PN.5.1`, mutation-proven.
  - **ADR-0006** (`docs/adr/0006-painleve-problem-layer.md`, **Accepted**) — the `PainleveProblem` builder decision: a per-equation problem *builder*, explicitly **not** a `painleve()` value-function (infinitely many solutions ⇒ "which solution" stays an explicit caller choice, Rule 1).
  - **`padetaylor-9xf` fixed** (`a2f319c`, worklog 026) — `path_network_solve` / `dispatch_solve` / `lattice_dispatch_solve` carried their own `order::Integer = 30` kwarg and silently ignored `prob.order`.  Now default to `prob.order`.  Tests `PN.7.1` / `DP.1.3` / `LD.2.2`, all mutation-proven.  Caught while writing `figures/fw2011_fig_5_2.jl` (an `(order,h)` sweep that was accidentally all-order-30).
  - **`padetaylor-avl` shipped** (`6ce0a97`, worklog 026) — `src/Painleve.jl`: wrapper struct + 6 constructors + 3 new direct-frame RHS factories (PI/PII/PIV; PIII/PV/PVI wrap the existing `CoordTransforms`/`SheetTracker` factories) + `solve_pade`/`path_network_solve` forwarding methods.  `test/painleve_test.jl`, 47 assertions, 3 mutations all bite.

  **Beads filed this session**: `padetaylor-d26` (Fig 3.1/3.2 — closed), `padetaylor-ote` (Fig 3.3 — closed), `padetaylor-cvo` (Fig 4.2/4.3/4.4 — closed), `padetaylor-664` (Fig 5.1/5.2 — closed), `padetaylor-9xf` (closed), `padetaylor-avl` (closed), `padetaylor-p3l` **P3 open** (quantitative pole-count pin for Fig 4.3/4.4), `padetaylor-soi` **P3 open** (`:transformed` forwarding with point-value remapping).

  **Suggested next work**: continue FW figure reproduction — Fig 4.7/4.8 need a *pole-extraction* helper (roots of each visited node's Padé denominator) which is a real capability gap worth a bead; Fig 4.1 steps (ii)+(iii) and Fig 4.6 are T3 BVP-hybrid reproduction.  Or attack `padetaylor-soi` to make the `:transformed` PainleveProblem forwarding fully usable.

### Session 2026-05-13 (post worklog 020) — `padetaylor-txg` closed (worklog 021)

  - **`src/RobustPade.jl`** — new `classical_pade_diagonal(c, m)` implementing FW 2011 §5.1.4 eqs. (5.4)+(5.5) via `lu(... ; check=false)` + `\`; throws `SingularException` on `!issuccess(F)`.  New `method::Symbol` kwarg on `robust_pade` with element-type-driven defaults (`_default_pade_method`).  Off-diagonal `(m≠n)` and singular Toeplitz both auto-route to existing SVD path.
  - **`test/classical_pade_test.jl`** (new, ~250 LOC) — 9 testsets / 63 assertions CP.1.1-9; mutation-proven (4 mutations P1-P4, all bite, all restored).  Bite cascades: P1 (sign-flip) cascades to 60+ RED across Phase 5/6/9 + PathNetwork + LatticeDispatcher; P3 (revert F64 default) bites 8 RED including the new rtol tightenings (proves they're load-bearing); P2 (drop fallback) bites 5; P4 (drop off-diagonal route) bites 4.
  - **`test/robustpade_test.jl`** — tests 2.1.2-4 + 2.1.6 (SVD-specific behaviour: diagonal-hop reduction, Froissart suppression, noise-thresholded recovery) explicitly pass `method=:svd` to keep testing GGT 2013 Algorithm 2 specifically.
  - **`test/pathnetwork_test.jl`** — PN.2.2 z=30 F64 rtol `1e-9 → 1e-12` (1000×); PN.2.3 z=10⁴ F64 rtol `5e-5 → 5e-10` (100,000×).  Both hold under classical.
  - **`docs/adr/0005-classical-pade-default-at-float64.md`** (new, ~150 lines) — element-type dispatch decision + when each method is right.
  - **`docs/worklog/021-classical-pade-default.md`** (new) — implementation log + mutation-proof bite counts + frictions.
  - Test suite: 1314 → **1377 GREEN** (+63 CP.1.* assertions).  Wall ~1m45s at 8 threads.

### Beads filed / closed this session

  - `padetaylor-txg` **closed** — classical-Padé F64 default shipped.

### Open beads end-of-session (worklog 021)

One P2 with concrete scope:
  - `padetaylor-7zw` P2 — BF-256 tritronquée pin (parallel testset for `test/fw_fig_41_test.jl`).

Plus existing P2/P3 admin / figure-pinning / perf beads (`padetaylor-bhw`, `-bho`, `-8ui`, `-1a3`, `-gky`, `-8pi`, `-61j`, and six P3s).  No P0 / P1 open.

### Hard-won lessons added this session (worklog 021)

**39. Mutation cascades expose architectural load-bearing-ness.** Mutation P1 (sign-flip on classical's RHS) was expected to bite ~10 CP.1.* assertions; it actually bit 60+ across Phase 5/6/9 + PathNetwork + LatticeDispatcher.  The mutation count is a direct measure of how deep an algorithm sits in the dependency graph; classical Padé is consumed by every step in every IVP path-network walk.  Lesson: don't pre-judge mutation bite scope — run it and see.

**40. Worklog 020's wedge-walker number is NOT the full-PathNetwork number.** Worklog 020 probed via a custom 5-direction wedge walker (Stage 1 only); the full `path_network_solve` adds Stage 2 barycentric interpolation, which is an extra Padé eval per target.  At z=10⁴ F64 this adds ~2× to the rel-err (6.15e-11 → 1.4e-10).  Lesson: when a probe simulates a stage, the production code's full pipeline carries additional error.  Always re-measure under the production code path before pinning test bounds.

**41. Degenerate test inputs are easier to construct by accident than expected.** CP.1.7's first attempt used `c_k = (1/2)^k` (geometric series) intending to demonstrate classical's no-reduction property — but the resulting Toeplitz is rank-1 (every row is a scaled copy of row 1), so classical correctly throws SingularException.  The "well-conditioned" claim of a Taylor series is a property of its convergence radius, NOT of the Toeplitz it builds.  Lesson: when picking a non-degenerate test case, compute or estimate the Toeplitz's rank explicitly before pinning expected behaviour.


## Older session entries

Sessions before "Session 2026-05-13 (post worklog 020)" — including
the v0.1.0 release session, the worklog-019 z=10⁴ investigation, the
classical-Padé empirical probe (worklog 020), the SheetTracker /
CoordTransforms / FW Fig 4.1 step-(i) BVP / `enforce_real_axis_symmetry`
beads (worklogs 015–018), the four-P1 evening (worklogs 011–014), and
all pre-v0.1.0 history — live in [`HANDOFF_ARCHIVE.md`](HANDOFF_ARCHIVE.md).

Read the archive when you need: (a) the hard-won lesson behind a
specific Phase you're touching, (b) the bead that closed a particular
historical issue, (c) the context for why an ADR cites a session
you don't recognise.  For day-to-day work, this file (the recent three
sessions) is enough.
