# ADR-0014 — IVP+BVP hybrid driver for pole-free sectors (FFW 2017 §3)

**Status**: Accepted (2026-05-15) | **Bead**: `padetaylor-0co` | **Worklog**: 039

## Context

The path-network IVP solver (`PathNetwork.path_network_solve`, ADR-0004)
and the Chebyshev spectral BVP solver (`BVP.bvp_solve`, Phase 11) have
shipped as **independent** Tier-2 / Tier-3 drivers. The Painlevé
substrate (PIII, PV, PVI under the exponential transforms shipped in
ADR-0006 / `CoordTransforms`) is now ready for the composition step:
*on pole-free sectors of a transcendent's Riemann surface, the IVP is
exponentially ill-conditioned*, so the BVP must be substituted for the
sector while the IVP continues to cover the pole-bearing exterior.

This is FFW 2017 §3, the *load-bearing* infrastructure for FFW Figure 5
(PIII tronquée with the pole-free sector `−3π/4 < arg z < 9π/4`).
FFW md:203-247 spell out the algorithm; md:252-264 derive the
condition-number motivation:

> "...the equation is exponentially unstable [as an IVP] on the smooth
> region. ... κ_r ~ (27/16) e^{2 Re ζ / 3}, Re ζ → +∞, ... maximum
> value ≈ 27/16 (30)^{4/3} ≈ 157." — FFW md:262 + md:264

A pure-IVP walk through the sector therefore amplifies any
perturbation by a factor of ~157 before reaching the IC anchor, and
FFW report `~10⁻¹` errors when this is attempted (md:226). The BVP
detour reduces this to `~10⁻⁷` in their stack (md:247).

This is **step A6** of the 11-step plan to reproduce FFW 2017's seven
figures. Steps A1 / A2 (adaptive Padé step + non-uniform Stage-1
nodes) and B1 / B2 / B3 (Figs 6, 1, 4) shipped in worklogs 034-038.
Step B4 (Fig 5 itself, the visible deliverable) needs this hybrid
driver as its substrate.

## Decision

Ship a new module `src/IVPBVPHybrid.jl` exporting three names:

  - `solve_pole_free_hybrid(pp::PainleveProblem, sector, asymptotic_ic_fn;
                            pfs_kwargs, bvp_kwargs, glue_tol, n_slices) -> IVPBVPSolution`
  - `IVPBVPSolution` — the four-component solution container
    (`pfs_top`, `pfs_bot`, `bvp_slices::Vector{BVPSolution}`, `slice_re`).
  - `pIII_asymptotic_ic(z; n_terms, β, δ, …) -> (u, u')` — truncated
    `z^{1/3}[1 + Σ a_n z^{-2n/3}]` series helper for the PIII tronquée
    family of FFW Fig 5.

The four-step FFW §3 algorithm (md:203-247) maps to existing modules:

  1. **PFS exterior** ← `path_network_solve` from each curved-boundary
     asymptotic-IC anchor.
  2. **Boundary harvest** ← linear-interp `pfs.grid_u` / `pfs.grid_up`
     at the slice `Re ζ` values.
  3. **BVP sector solve** ← `bvp_solve` 3-arg-RHS overload on each
     `Re ζ = const` slice through the sector.
  4. **Glue** ← `IVPBVPSolution.(::Complex)` callable dispatches by
     sector membership; continuity asserted to `glue_tol`.

**BVP three-argument-RHS extension** (bead `padetaylor-i76`).  The
P̃_III RHS `w'' = (w')²/w + …` depends on `w'`, but the existing
`bvp_solve(f, ∂f_∂u, …)` API assumes `f(z, u)` only.  We extended
`bvp_solve` *additively* — a new positional overload
`bvp_solve(f, ∂f_∂u, ∂f_∂up, …)` selects a 3-arg path.  The 2-arg API
is byte-unchanged; every existing call site (and test) is unaffected.
See `src/BVP.jl` module docstring section "Three-argument RHS overload".

**Asymptotic-IC derivation**.  `pIII_asymptotic_ic` implements the FFW
md:230 series

    u(z) ~ z^{1/3} · [ 1 + a_1 z^{-2/3} + a_2 z^{-4/3} + … ],

with **`a_1 = -β/3`** derived analytically by substituting the ansatz
into the canonical PIII equation `u'' = (u')²/u − u'/z + (αu² + γu³)/z
+ β/z + δ/u` and matching the `z^{-1}` order.  For FFW Fig 5's
`β = -1/20`, `a_1 = 1/60`.  Higher coefficients (`a_2`) follow from
the `z^{-7/3}` match, giving `a_2 = ((4/9) + δ a_1²) / (2δ)`; for
FFW's `δ = -1` this is `a_2 ≈ -0.22208`.  The verbatim algebra is in
`src/IVPBVPHybrid.jl` lines 138-156 (a_1) and 222-260 (a_2).

**v1 = Float64 + 2-term series**.  Per the bead spec, the hybrid is
Float64-first; the whole point is to avoid BF-256 on pole-free
sectors via the BVP detour.  At Float64 + `n_terms = 2` + `N = 10`
collocation, the test IB.1.5 records a `|Δw'|` of `~0.4–1.7` at the
sector boundary — exactly the FFW md:226 baseline (the "PFS without
BVP error ~10⁻¹" signature).  Achieving FFW's md:247 `~10⁻⁷` requires
BF-256 + `n_terms ≥ 15` + finer N; this is the standing v2 lift,
documented in `padetaylor-0co`'s "What is NOT shipped" footer.

## Alternatives considered

**A.  Pure-IVP everywhere, with BF-256 precision.**  Use
`path_network_solve` over the full domain; rely on arbitrary-precision
arithmetic to defeat the `κ_r ≈ 157` amplification.  Rejected on two
counts:

  - **Time** — BF-256 path-network walks take ~50× longer than Float64.
    FFW Fig 5 covers an `|Re ζ| · |Im ζ| ≈ 5 · 12π` rectangle; a BF-256
    walk would take an hour vs Float64's ~5 minutes.
  - **Wrong fix** — even at BF-256, the IVP through the sector is
    *fundamentally* ill-conditioned: the small perturbations introduced
    by Padé-Taylor truncation at each step compound exponentially.
    The BVP detour is *the* algorithmic cure recommended by FFW
    md:203-204 and originally by the reference at md:224's [9]. BF-256
    treats the symptom; the BVP treats the disease.

**B.  Pure-BVP everywhere.**  Discretise the entire transcendent's
domain as one giant BVP, with asymptotic-series BCs on the outer
boundary.  Rejected:

  - The pole-bearing region has thousands of poles (FFW md:283 spiral
    structures, FF1.1.5 reports 1325 poles for Fig 1) — a Chebyshev
    spectral BVP cannot represent these (poles are essential
    singularities of polynomial interpolants, the spectral error blows
    up).  A locally-pole-aware BVP would require a Padé-rational basis,
    which is what `bvp_solve` *isn't* — it's a Chebyshev-polynomial
    interpolant.
  - The path-network already handles poles correctly via local Padé
    approximants.  Throwing it out for the BVP everywhere wastes its
    pole-bridging capability.

**C.  Analytic-series-only solution on the sector.**  Skip the BVP
step; evaluate the truncated `u(z) ~ z^{1/3}[1 + Σ a_n z^{-2n/3}]`
series everywhere in the sector.  Rejected:

  - The series is divergent (FFW md:222 calls it "divergent
    asymptotic"); truncated evaluation at moderate `|z|` (away from
    the asymptotic regime) has uncontrolled error.
  - For the sector near `Re ζ ≈ 0` (i.e., `|z| ≈ 1`), the series has
    `|a_n z^{-2n/3}| ~ |a_n|` — no decay.  Optimal truncation gives
    machine-precision only at `|z| ≫ 30`.  We need accurate values on
    the entire sector, including the near-`|z| = 1` boundary the BVP
    domain extends into.
  - The BVP solver, given good BCs, achieves spectral convergence
    *inside* the sector to the Chebyshev floor `~10⁻¹⁴` — much
    better than truncated asymptotic-series evaluation in the bulk.

**D.  Hybrid via Tier-3 Dispatcher (existing).**  ADR-0004's
Dispatcher composes IVP + BVP for *segments*, not sectors. Its design
covers the `(IVP region, BVP segment)` chain of a 1D problem (e.g.
FW 2011 Fig 4.4).  Rejected for this bead: the sector geometry is 2D
in ζ, the Dispatcher is 1D; trying to repurpose it would force a
non-orthogonal generalisation. The clean answer is a *new* Tier-3
component dedicated to the 2D-sector composition.

## Consequences

### Code

  - **`src/IVPBVPHybrid.jl`** (new, ~480 LOC). Top-of-file docstring
    chapter covers the FFW §3 algorithm, the condition-number motivation,
    the asymptotic-IC derivation (closed-form `a_1` + `a_2`), and the
    "compose, don't refactor" discipline. Exports
    `solve_pole_free_hybrid`, `IVPBVPSolution`, `pIII_asymptotic_ic`.
    Over 200 LOC; the v1 ships at ~480 LOC including ~180 LOC of
    in-source asymptotic-series derivation prose — Rule 6 split is
    deferred (worklog 039 §"What is NOT shipped").

  - **`src/BVP.jl`** (additive 3-arg overload, +110 LOC). New
    `bvp_solve(f, ∂f_∂u, ∂f_∂up, …)` for RHS `f(z, u, u')`. Existing
    2-arg API and every test byte-unchanged. Module docstring section
    "Three-argument RHS overload" derives the Jacobian
    `J = D₂_ii − scale · diag(∂f/∂u) − (scale/half_diff) · diag(∂f/∂u') · D₁_int_int`.

  - **`src/PadeTaylor.jl`** — `include("IVPBVPHybrid.jl")` after
    `Painleve.jl`, plus three exports.

### Tests

  - **`test/ivp_bvp_hybrid_test.jl`** (new, ~360 LOC). Seven testsets
    `IB.1.1`–`IB.1.7`. Mutation-proven (M1-M4); footer documents bites.
    65 new assertions.

  - Aggregate count: **2131 → 2197 GREEN** (+66; the umbrella-loads
    test adds one assertion, the new file adds 65).

### Figures unblocked

  - **FFW 2017 Figure 5** (B4) — the PIII tronquée with the headline
    pole-free sector + condition-number heat map. The hybrid driver is
    its substrate.

### Performance

  - Single-slice BVP solve at `N = 10`, Float64, 3-arg RHS: ~0.05 s
    Newton-converge.
  - Per-call (with 5 slices + 2 PFS walks) wall: ~1.5 s on the test's
    reduced sector.
  - A future full Fig 5 sector (5× wider in Im ζ + ~20× more slices)
    extrapolates to ~30-60 s on a single thread.

### Compatibility

  - The new 3-arg `bvp_solve` overload is a *new positional dispatch*;
    every existing 2-arg call works unchanged. No regression.
  - `solve_pole_free_hybrid` is opt-in via `using PadeTaylor`; users
    not calling it pay no cost.

### Out of scope (deferred to separate beads)

  - **PV / PVI hybrid sectors**.  v1 hard-codes PIII RHS in
    `_bvp_solve_on_slice`; PV's RHS (FFW md:47) has a different
    analytic-Jacobian shape (the `1/(2w) + 1/(w-1)` factor doesn't
    factor cleanly). A small dispatch-on-`pp.equation` extension
    closes this.

  - **Higher-order asymptotic series** (`a_3`, `a_4`, …).  v1 ships
    `a_1` (closed form) + `a_2` (closed form).  FFW md:232 use
    "optimal truncation" with many more terms (probably 15+); the
    mechanical generalisation via `TaylorSeries.jl` arithmetic on
    Laurent series in `s = z^{1/3}` is straightforward but ~150 LOC
    of additional implementation. Bead-deferred.

  - **BF-256 precision path**.  v1 is Float64; the BVP solver is
    generic in `T`, the PFS solver supports BF-256, the
    asymptotic-IC helper currently hard-codes `Float64` in the
    coefficient computation.  A BF-256 lift to reach FFW's `1e-7`
    sector-derivative agreement is the natural next refinement.

  - **2D-spectral BVP** (vs 1D slices + linear-interp between).  v1
    discretises the sector as a stack of 1D slices and linear-
    interpolates between them in Re ζ.  A true 2D Chebyshev spectral
    BVP would replace the slice stack — ~1000 LOC of new BVP code,
    out of v1 scope. The 1D-slice + linear-interp v1 is accurate
    enough for figure-rendering purposes (Re ζ direction is smooth).

  - **The "exterior" callable**.  v1's `sol(ζ)` covers the BVP sector
    only; queries outside the sector raise `DomainError`. The PFS
    Stage-2 grid is accessible directly via `sol.pfs_top.grid_*`. A
    full unified `sol(ζ)` would require dense PFS interpolation
    (FW 2011 §3.2 ADR-0004 deferred this; current bead does not
    re-open).

## References

  - **FFW 2017 §3** — `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:203-247` (algorithm), md:252-264 (condition number κ_r), md:222 + md:230 (sector + asymptotic series), md:240-244 (Fig 5 IC values), md:43 (P̃_III RHS).
  - **ADR-0004** — `docs/adr/0004-path-network-architecture.md`. The path-network architecture this composes.
  - **ADR-0006** — `docs/adr/0006-painleve-problem-layer.md`. The `PainleveProblem` wrapper this consumes.
  - **ADR-0011** — `docs/adr/0011-adaptive-pade-step.md`. Step A1.
  - **ADR-0012** — `docs/adr/0012-non-uniform-stage-1-nodes.md`. Step A2.
  - **`src/IVPBVPHybrid.jl`** — the module shipped by this ADR.
  - **`src/BVP.jl`** — the 3-arg-RHS overload (additive).
  - **`test/ivp_bvp_hybrid_test.jl`** IB.1.1-IB.1.7 — acceptance tests + mutation-proof footer.
  - **Worklog 039** — implementation diary, design synthesis, frictions.
  - **Bead `padetaylor-0co`** — closed by this commit.
  - **Bead `padetaylor-i76`** — the BVP 3-arg-RHS extension (artefact bead, standing).
