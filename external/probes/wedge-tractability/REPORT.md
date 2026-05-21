# Phase-A2 spike — wedge tractability probe

**Bead**: `padetaylor-0ln.37.2` (ADR-0025 Phase A2, Lever 2).
**Status**: complete. Throwaway exploration code; no `src/` or `figures/`
production file was modified.
**Run**: `julia --project=figures external/probes/wedge-tractability/probe.jl`

---

## 1. Mission

The headline P_I⁽²⁾ tritronquée figure renders the pole-rich wedge
(`|arg x| < 36°`, straddling the positive real axis) with a vector
path-network walk + Stage-2 fine-grid fill. The shipped V8b walk
(`figures/_kkg_pi2_helpers.jl`) uses **20 targets all within `|x| ≤ 8`**,
fixed step `h = 0.1`, Taylor order 24. The figure window is `|x| ≤ 20`.
ADR-0025 reframes bead `0ln.37` acceptance criterion 3 — "full `|x|≤20`
wedge coverage" — to "honest coverage to the numerically-supported
frontier, the frontier itself a measured and reported quantity".

This spike measures that frontier. "Honest" = no shared-Q Padé is ever
evaluated outside its **A1-verified disc of validity**. The A1 probe
(`external/probes/pade-validity-radius/REPORT.md`) resolved that disc to
gate **(v-b)**:

```
R_gate = min( s(tol)·h_JZ·h_v ,  0.5·h_v·min|t*| )
  h_JZ    = vector_step_jorba_zou(rescaled order-24 jet, tol)
  s(1e-8) = 0.36  (production default)
```

The four tasks: (2) extend the target fan toward `|x|=20`; (3) measure
honest **union** coverage of the wedge vs `|x|`; (4) probe the fixed-`h`
failure mode; (5) deliver a tractability verdict.

---

## 2. Methodology

### 2.1 Reproduction of the V8b substrate

The probe rebuilds the exact Stage-A BVP from `figures/_kkg_pi2_helpers.jl`
(2+2 `(u,u')` boundary condition on `[−20,−2]`, `N=128` Chebyshev,
3 Newton iterations, residual `1.5e−11`) and seeds the path-network walk
from `bvp(z=−3)`. Constants are copied verbatim (deliberately not
`include`d) so the spike is frozen.

### 2.2 Per-node honest disc — the A1 gate (v-b)

For each visited node `k`: rebuild the order-24 Taylor jet via the
production primitive `vector_taylor_coefficients(f, z_v, y_v, 24)` (the
same jet the canonical Padé was built from); rescale it to the
`t`-variable (`c_k·h_vᵏ`); feed it to `VectorStepControl.vector_step_jorba_zou`
with `eps_abs = 1e-8`; multiply by `h_v` for the z-plane radius; clamp by
`0.5·h_v·min|t*|` over the node's shared-Q roots. This is the A1 winner
verbatim — reuse of tested project machinery.

### 2.3 Union honest coverage

The wedge is the region `{|arg x| ≤ 36°}`. For each annulus `[r, r+2]` the
probe samples a deterministic `18×60` polar grid and counts the fraction
of sample points covered by **at least one** node's honest disc
(`|z − z_k| ≤ R_k`). The **honest frontier** is the largest `r` for which
the annulus `[r−2, r]` still reaches ≥ 50 % union coverage.

### 2.4 Fixed-h failure mode

`vector_path_network_solve` aborts the **whole** solve on a single throw
(no partial result). So a fan reaching past the walk's death radius
reports only "BLOCKED, 0 nodes" — useless for locating the frontier. The
probe therefore also runs a **manual wedge ray-walk**: it replicates the
production `_wedge_step` (the min-‖y‖ 5-direction `±22.5°/±45°` wedge,
`VectorPathNetwork.DEFAULT_WEDGE`) and marches step-by-step from the seed
toward a far wedge target, catching the throw and recording the `|x|` of
the last good node. Run at `h ∈ {0.10, 0.05, 0.03}`.

---

## 3. Results

### 3.1 The extended fan blocks past `|x| = 8` (Task 2)

A V8b-shaped fan (5 angles spanning `±0.5 rad ≈ ±28.6°`, radial spacing
`dr = 2`) at the production `h = 0.1`:

```
fan r_max= 8.0 (20 targets): OK      389 nodes,  honest-reach |z|=8.08
fan r_max=12.0 (30 targets): BLOCKED  0  nodes
fan r_max=16.0 (40 targets): BLOCKED  0  nodes
fan r_max=20.0 (50 targets): BLOCKED  0  nodes
```

Every fan reaching past `|x|=8` throws — `shared_denominator_pade: every
singular value is below τ` — i.e. a walk step landed on a pole and the
shared-Q jet degenerated. The whole solve aborts. **The fixed-`h=0.1`
production walk cannot extend the V8b fan past `|x|=8`.** This is exactly
ADR-0025 Lever-2 / C3: the hand-tuned `h=0.1` is region-specific.

### 3.2 The fixed-h failure mode is fan-density-driven, not a hard radius (Task 4)

The manual single-ray wedge walk tells a different story — a single ray
reaches `|x|=20` at **every** step size:

```
h      target            reach |x|   steps   cause
0.10   (20.0, 0.0)        20.09       300     reached
0.10   (18.4, ±7.8)       19.93       281     reached
0.05   (20.0, 0.0)        20.04       607     reached
0.03   (20.0, 0.0)        20.02      1019     reached
```

The min-‖y‖ wedge selection **does** successfully steer a single march
around the pole field all the way to `|x|=20`. The failure in §3.1 is not
"the wedge cannot reach `|x|=14`" — it is that a *multi-target fan* forces
many **inter-target bridging walks**, each of which crosses the pole field
on a chord, and a fixed-`h` chord step eventually lands on a pole.

The largest *completing* production fan, per step size (5-angle `dr=2`
fan, scanned downward from `r_max=20`):

```
h      max r_max   nodes   reach |z|   runtime
0.10   8.0          389     8.08        1.0s
0.05   20.0        2377    20.07        6.6s
0.03   20.0        4008    20.05       11.1s
```

A finer `h` **does** push the completing frontier out: `h=0.05` and
`h=0.03` complete the full `r_max=20` fan where `h=0.1` dies at `r_max=8`.
So an **adaptive / finer step (B2) reaches `|x|=20` geometrically.**

### 3.3 But honest *coverage* saturates near ~8 % — the decisive finding (Task 3)

Reaching `|x|=20` with nodes is not the same as honestly *covering* the
wedge. The A1 gate-(v-b) honest discs are small (`R_med ≈ 0.05–0.06`).
Union honest coverage of the wedge, V8b fan (`r_max=8`, `h=0.1`,
389 nodes):

```
annulus       honest coverage   nodes-in
[ 0.0, 2.0]        6.6%          24
[ 2.0, 4.0]        7.6%          45
[ 4.0, 6.0]        7.5%          92
[ 6.0, 8.0]        7.3%         161
[ 8.0,20.0]      ~0.0%           0
Honest frontier (>=50% coverage): |x| ~ 0.0
```

Coverage scaling — denser fans, finer `h`:

```
                                       nodes   cov(|x|<=20)   frontier
h=0.05  5ang dr2.0  ( 50 tgt)            2377      5.2%          0
h=0.05  9ang dr1.0  (171 tgt)         BLOCKED
h=0.03  5ang dr2.0  ( 50 tgt)            4008      5.6%          0
h=0.03  9ang dr1.0  (171 tgt)            6850      8.6%          0
h=0.03 13ang dr0.5  (481 tgt)         BLOCKED
```

The best completing dense walk — **6850 nodes**, 9-angle `dr=1.0` fan at
`h=0.03`, 20 s runtime — achieves only **8.6 %** honest coverage of the
`|x|≤20` wedge. Its per-annulus profile (`R_med = 0.0465`):

```
annulus       honest coverage   nodes-in   tile-need
[ 0.0, 2.0]        5.5%          88            618
[ 2.0, 4.0]       14.8%         311          1853
[ 4.0, 6.0]       14.4%         608          3089
[ 6.0, 8.0]       12.1%         757          4324
[ 8.0,10.0]        9.2%         853          5560
[10.0,12.0]        6.7%         805          6795
[12.0,14.0]        5.8%         803          8031
[14.0,16.0]        6.3%         807          9266
[16.0,18.0]        3.2%         817         10502
[18.0,20.0]        3.4%         780         11737
```

**No annulus, at any tested configuration, ever reaches 50 % honest
coverage. The honest frontier under a path-network-walk + per-node-disc
architecture is `|x| ~ 0` — there is no radius at which the wedge is
honestly *filled*.**

The denser fans BLOCK: adding targets to fill the 2D wedge multiplies the
inter-target bridging walks, each a fresh chance to step onto a pole. The
walk is structurally a **ray-filament network**, not a 2D tiling — and a
filament network of small discs cannot tile a 2D region.

### 3.4 Why coverage saturates — the tiling arithmetic

The `|x|≤20` wedge has area `≈ 0.628·rad × 20² = 251`. An A1 honest disc
has area `π·R_med² ≈ π·0.0465² ≈ 0.0068`. At a realistic 2D packing
efficiency `η ≈ 0.6`, tiling the wedge needs

```
N_tile ≈ 251 / (0.6 · 0.0068) ≈ 62,000 honest nodes.
```

The 6850-node walk has **~9× too few nodes**, and they sit on filaments
rather than a 2D lattice — hence the measured 8.6 %, consistent with the
arithmetic. The per-annulus `tile-need` column ranges `618` (inner) to
`11,700` (outer `[18,20]` annulus alone). This is the honest bottleneck:
the A1 gate makes each disc small *by design* (the truncation-limited
honest radius is genuinely `≈ 0.05`), so honest coverage demands a node
count a path-network walk does not and cannot produce.

---

## 4. Verdict

### (a) Honest frontier of the fixed-`h` walk

The fixed-`h=0.1` production fan walk **blocks past `|x|=8`** (§3.1) — its
geometric reach is `|x| ≈ 8`. But even within `|x|≤8` its **honest
coverage frontier is `|x| ~ 0`**: no annulus reaches 50 % honest coverage
(§3.3). The fixed-`h` walk does not honestly tile the wedge **anywhere**.

### (b) Frontier achievable with an adaptive walk (B2)

An adaptive / finer-`h` walk (B2) **does** fix the *geometric* reach:
`h=0.03` marches nodes out to `|x|=20` (§3.2). But B2 does **not** fix the
*coverage* frontier. The best B2-style walk measured (6850 nodes, `h=0.03`,
9-angle fan) still achieves only **8.6 %** honest coverage and a 50 %-
coverage frontier of `|x| ~ 0` (§3.3). **B2 is necessary for geometric
reach but not sufficient for honest coverage.** Denser fans, the obvious
next lever, BLOCK — so B2 cannot be pushed to a tiling node count by fan
density alone.

### (c) Node count to tile the wedge

To honestly tile the `|x|≤20` wedge with A1 gate-(v-b) discs
(`R_med ≈ 0.05`): **~62,000 nodes** (§3.4). The largest *completing*
walk reaches ~6,850. The gap is ~9×, and the walk architecture (ray
filaments + pole-hitting bridging segments) cannot close it.

### (d) Runtime tractability

Runtime is **not** the binding constraint. The 6850-node walk runs in
~20 s; a 62,000-node walk would extrapolate to a few minutes — tractable.
The binding constraint is that the path-network walk **cannot be made to
visit 62,000 nodes**: the dense fans needed to drive the node count up
BLOCK on pole hits in the bridging walks.

### (e) Is full `|x|≤20` honest coverage realistic?

**No.** Under the architecture ADR-0025 commits to — a Stage-1
path-network walk producing per-node shared-Q discs gated by A1 (v-b) —
full honest coverage of the `|x|≤20` wedge is **not realistic**, and
neither is honest coverage to *any* non-trivial frontier radius. The
honest 50 %-coverage frontier measured is `|x| ~ 0`. The figure honestly
rendered this way would show ~8 % filled wedge and ~92 % `NaN` gaps.

This is the honest finding ADR-0025 §Decision-2 anticipated ("whether
honest coverage *reaches* `|x|=20` is a genuine numerical unknown") — and
the answer is sharper than "frontier at `|x|=14`": the *coverage* frontier
is `|x|≈0` because the binding limit is the **A1 truncation-limited disc
size**, not the pole-field radius. The A1 report already flagged this
coupling ("B1 alone makes the figure *holier*, not bigger"); A2 confirms
that **B2 does not rescue it either** — denser walks block, and the disc
size is fixed by truncation error.

---

## 5. Recommendation for B3 (the production target fan) — and beyond

### 5.1 B3 as scoped cannot deliver `|x|≤20` honest coverage

Bead B3 ("principled extended target fan") as currently scoped — a wider
path-network target fan — **cannot** produce a honestly-filled `|x|≤20`
wedge. A wider/denser fan either (i) BLOCKS on pole hits (the 9/13-angle
dense fans), or (ii) completes but still tiles only ~8 % of the wedge.
The orchestrator must know this before B3 is built: **the path-network
walk + per-node A1 disc is the wrong tool for a *filled surface* over the
wedge.**

### 5.2 Recommended B3 target-fan design (the honest, deliverable scope)

If the figure's wedge panel is to be **honest**, B3 should produce the
densest walk that **completes**, and the figure should render exactly the
honest-covered fraction with explicit `NaN` gaps elsewhere:

- **Step size**: `h = 0.03` (fixed). It is the finest tested that still
  completes the full `r_max=20` fan; `h=0.05` completes the 5-angle fan
  but blocks the 9-angle one.
- **Fan**: the **9-angle (`±0.5 rad`) `dr=1.0` fan at `h=0.03` is the
  densest that completes** (6850 nodes, ~20 s, 8.6 % honest coverage).
  The 13-angle `dr=0.5` fan blocks; do **not** scope B3 to it.
- **Honest output**: render the ~8.6 % honestly-covered region; everywhere
  else is an honest `NaN` gap (ADR-0025 fail-loud). This is "as good as the
  path-network walk honestly manages" — and it is *thin*.

### 5.3 The architectural escalation B3 actually needs (raise to ADR-0025)

An 8.6 %-filled wedge panel is not a senior-grade headline figure. The
honest path to a *filled* wedge surface needs an architecture change
beyond B2/B3, and the spike's evidence points at three candidates the
orchestrator should weigh (each is a Phase-B re-scope, not a tweak):

1. **Higher Taylor order.** The honest disc is truncation-limited
   (A1 §3.3); `R_honest` grows with jet order. Going from order 24 to,
   say, 48–64 would enlarge each honest disc and cut the tile-need
   quadratically in `R`. This is the most direct lever — the A1 gate
   *already* shows order is what binds — and it keeps the existing
   architecture. **Recommended first escalation.** Cost: order-64 jets are
   ~7× the per-node compute of order 24 (A1's oracle ran them); the walk
   would slow but stay tractable, and far fewer nodes are needed.
2. **A genuine 2D fill**, not a ray-filament walk: a triangulated /
   grid-seeded set of *independent* local Padé expansions (each its own
   short IVP from a BVP-or-walk-provided nearby state), so coverage is a
   2D lattice of discs rather than filaments — and there are no
   pole-hitting *bridging* walks.
3. **Accept a thin honest wedge** and reframe the headline figure's wedge
   panel as a pole *scatter* (the V8b-shipped deliverable) over the
   honestly-covered fraction — explicitly *not* a filled surface. ADR-0025
   §Alternatives already rejected a pure scatter as "lower-ambition"; this
   spike is the evidence that the filled-surface ambition is, under the
   current architecture, **numerically unreachable**.

The honest one-line verdict: **the wedge cannot be honestly tiled to
`|x|=20` — nor to any non-trivial frontier — by a path-network walk with
order-24 A1-gated discs; the 50 %-coverage frontier is `|x|≈0`. B2
(adaptive `h`) is required for geometric *reach* but does not rescue
*coverage*. Full honest coverage needs a higher Taylor order (the
truncation-limited disc is the true bottleneck) or a 2D-fill architecture
— a Phase-B re-scope, recorded here for ADR-0025 Amendment 2.**

---

## 5b. Phase-A2b — the display-tolerance follow-up

### 5b.0 Why this section exists

§3–§5 measured honest coverage with the A1 gate (v-b) at **solver**
tolerances (`1e-6/1e-8/1e-10`); the headline 8.6 % used `tol = 1e-8`. But
the headline figure is a **display artifact** — a clamped Re/Im heatmap
plus a 3D surface. A display colour channel resolves only `~1e-2..1e-3`
relative accuracy: you cannot *see* better than that on a clamped
colormap. Rendering the wedge **surface** honest-to-`1e-3` is therefore
legitimate ("honest for display") — provided the pole **locations**, the
actual scientific content, stay separately extracted and validated to full
`1e-8+` accuracy. The A1 honest disc grows as `tol` relaxes
(`h_JZ = vector_step_jorba_zou(jet, ε)` grows with `ε`). This section
recomputes honest **union** coverage of the `|x|≤20` wedge at **display**
tolerances, and tests whether raising the Taylor order to 48 rescues it.

The A2b run reuses the §3.3 best-completing fan verbatim: the 9-angle
`dr=1.0` fan at `h=0.03`, **6850 nodes**, ~20 s. `s` is held at the
production default `0.36` across all display tolerances (A1 §7 calibrated
`s` only for `1e-6/1e-8/1e-10`; holding it fixed is the **conservative**
choice — the true display-tol `s` would, if anything, be slightly larger
as the over-rate band loosens, so the coverage below is a *lower bound*).

### 5b.1 Honest union coverage at display tolerances (Tasks 1 + 2)

Order-24 jets, A1 gate (v-b), 6850-node fan, `|x|≤20` wedge:

```
tol      R_med    R_p10/p90       honest cov(|x|<=20)   50%-frontier
1e-8     0.0465   0.0283/0.0647         8.6%   (baseline)   |x|~0
1e-2     0.0843   0.0510/0.1176        14.8%                |x|~0
1e-3     0.0764   0.0462/0.1065        13.5%                |x|~0
1e-4     0.0692   0.0419/0.0964        12.4%                |x|~0
```

Relaxing from the solver `1e-8` to the display `1e-3` lifts honest union
coverage from **8.6 % → 13.5 %** (a 1.6× gain) and the per-node honest disc
from `R_med 0.047 → 0.076` (1.6×). At the loosest display tol `1e-2` it
reaches 14.8 %. **The gain is real but modest** — `R` grows only
logarithmically-ish in `tol` (the Jorba–Zou radius scales like
`tol^(1/order)`), so three orders of magnitude of tolerance buys only a
~1.6× radius.

### 5b.2 Per-node honest disc `R_med` at each display tol (Task 2)

`R_med` = `0.0843 / 0.0764 / 0.0692` at tol = `1e-2 / 1e-3 / 1e-4`
(order 24); vs `0.0465` at the solver `1e-8`. Even the loosest display
disc (`R_med ≈ 0.084`) is small against the wedge: disc area
`π·0.084² ≈ 0.022` vs wedge area `251`.

### 5b.3 Node count to tile the wedge at tol = 1e-3 (Task 3)

With `R_med = 0.0764` (order 24, tol `1e-3`), honest disc area
`π·R_med² ≈ 0.0183`. At packing efficiency `η ≈ 0.6`:

```
N_tile  ≈  251.3 / (0.6 · 0.0183)  ≈  22,900 honest nodes.
```

The best *completing* fan has **6850 nodes — 3.3× short** (vs 9× short at
the solver `1e-8`). Display tolerance shrinks the gap from 9× to 3.3×, but
does **not** close it: the denser fans that would supply 22,900 nodes still
**BLOCK** on pole hits in the inter-target bridging walks (§3.3 — that
failure mode is geometric and tolerance-independent).

### 5b.4 Coverage-vs-radius at tol = 1e-3 — does the surface reach `|x|=20` ? (Task 4)

Order-24 jets, tol `1e-3`, per-annulus honest union coverage:

```
annulus       honest coverage   nodes-in
[ 0.0, 2.0]        9.4%           88
[ 2.0, 4.0]       25.6%          311      <- peak
[ 4.0, 6.0]       23.7%          608
[ 6.0, 8.0]       20.4%          757
[ 8.0,10.0]       16.9%          853
[10.0,12.0]       11.9%          805
[12.0,14.0]        9.8%          803
[14.0,16.0]        9.1%          807
[16.0,18.0]        6.0%          817
[18.0,20.0]        5.2%          780
```

The honest surface **does reach `|x|=20`** geometrically (nodes populate
every annulus, the outer `[18,20]` has 780 of them) — but only at **5.2 %
honest coverage** there. Coverage *decays monotonically* outward: the
inner annuli peak at ~25 %, the outer wedge is ~5 %. **No annulus, at any
display tolerance, reaches 50 %.** The 50 %-coverage frontier stays
`|x|≈0`.

### 5b.5 Raising the Taylor order 24 → 48 (Task 5)

**An order-48 *walk* BLOCKS.** Running the 6850-node fan at order 48 throws
`shared_denominator_pade: every singular value is below τ` — the shared-Q
Padé linear solve degenerates: with `m = order/2 = 24` numerator/denominator
degree, the order-48 jets produce an all-below-tolerance singular spectrum.
**An order-48 walk is not viable for this wedge under the shared-Q route.**
This is a new, independent obstruction: order 48 breaks the *walk*, not
just the coverage.

The A1 honest disc, however, is a property of the **jet**, not the shared-Q
Padé — `vector_step_jorba_zou` reads the order-`N` coefficient decay
directly. So we measured order 48's *disc* effect cleanly by rebuilding
**order-48 jets at the order-24 walk's (already honest) node states** and
re-gating — the truncation bound, decoupled from the failed order-48 Padé
construction:

```
ORDER-48 JETS on the 6850-node order-24 walk:
tol      R_med    R_p10/p90       honest cov(|x|<=20)   50%-frontier
1e-8     0.0746   0.0400/0.1096        13.3%                |x|~0
1e-2     0.0996   0.0528/0.1466        18.0%                |x|~0
1e-3     0.0949   0.0504/0.1397        16.9%                |x|~0
1e-4     0.0904   0.0481/0.1330        16.2%                |x|~0
```

Order 48 + display tol `1e-3` together: `R_med 0.0764 → 0.0949` (**1.24×**),
honest coverage `13.5 % → 16.9 %`. Per-annulus at order-48 jet / tol `1e-3`:
inner peak 31 %, outer `[18,20]` 6.5 % — still **monotone-decaying, still
no annulus past 50 %**.

The combined best case — order-48 jets, display tol `1e-2` — reaches
**18.0 %** honest union coverage. That is the ceiling this architecture
offers, and it is **a quarter of the 70 % viability bar**.

---

## 6. Files

- `probe.jl` — the full A2 + A2b spike: BVP + seed, extended fans, A1
  gate-(v-b) honest discs, union coverage, the manual wedge ray-walk
  failure-mode probe, coverage-scaling sweep, tiling estimate, and (§A2b)
  the display-tolerance recompute + the order-24→48 comparison.

## 7. Final verdict — display tolerance + order 48 (the follow-up)

### 7.1 The question

Phase-A2's headline (§3.3) said honest coverage **saturates at 8.6 %** —
but it used the *solver* tolerance `1e-8`. The follow-up asks: at **display
tolerance** (`1e-2..1e-3`, all a clamped colormap can resolve) and with the
Taylor order raised 24 → 48, is a **filled honest wedge surface viable**
(≥ 70 % coverage to `|x|=20`) — or does the architectural rescope (render
the wedge as an *extracted + validated pole-location field*, with an honest
surface only where discs cover) still stand?

### 7.2 The numbers

| configuration | `R_med` | honest cov `|x|≤20` | tile-need | 50 %-frontier |
|---|---|---|---|---|
| order 24, tol 1e-8 (A2 baseline) | 0.047 | **8.6 %** | ~62,000 | `|x|≈0` |
| order 24, tol 1e-3 (display)     | 0.076 | **13.5 %** | ~22,900 | `|x|≈0` |
| order 24, tol 1e-2 (display)     | 0.084 | **14.8 %** | ~18,700 | `|x|≈0` |
| order 48 jets, tol 1e-3          | 0.095 | **16.9 %** | ~14,800 | `|x|≈0` |
| order 48 jets, tol 1e-2 (best)   | 0.100 | **18.0 %** | ~13,400 | `|x|≈0` |

Every lever helps a little; none is decisive. Display tolerance: 8.6 → 13.5 %.
Order 48 jets on top: 13.5 → 16.9 %. Loosest display tol on top: → 18.0 %.
Across **three** orders of magnitude of tolerance **and** a doubled Taylor
order, honest union coverage moves **8.6 % → 18.0 %** — a 2.1× gain that
leaves the figure **82 % NaN**. No annulus ever crosses 50 %; the
50 %-coverage frontier stays `|x|≈0` in every configuration. The honest
surface *does* reach `|x|=20` geometrically, but at 5–6 % coverage there.

### 7.3 Two independent obstructions

1. **Coverage arithmetic.** Even at the most generous setting the honest
   disc is `R_med ≈ 0.10`; tiling the `area-251` wedge needs ~13,400 nodes,
   and the densest *completing* fan produces 6850. The denser fans that
   would supply the rest **BLOCK** on pole hits in inter-target bridging
   walks (§3.3) — a geometric failure mode, tolerance- and order-independent.

2. **The order-48 walk itself blocks.** Raising the order to 48 degenerates
   the shared-Q Padé linear solve (`m = 24` deg, all-below-τ singular
   spectrum). Order 48 cannot enlarge the *shipped figure's* discs at all
   unless the shared-Q construction is replaced — the order-48 numbers in
   §7.2 are a *what-if* on rebuilt jets, **not** a deliverable walk.

### 7.4 VERDICT

> **A filled honest wedge *surface* is NOT viable, even at display
> tolerance and even with order 48.** Best achievable honest union coverage
> of the `|x|≤20` wedge is **~18 %** (order-48 jets, tol `1e-2`) — well
> short of the 70 % bar, and no annulus anywhere reaches 50 %. Display
> tolerance and higher order *help* (8.6 % → 18 %, a 2.1× gain) but the
> binding constraints are architectural: (i) the completing fan tops out at
> ~6850 nodes while ~13,000–23,000 are needed, and the denser fans BLOCK;
> (ii) an order-48 *walk* degenerates the shared-Q Padé solve outright.
>
> **The architectural rescope STANDS.** The headline figure's wedge panel
> must be rendered as an **extracted-and-validated pole-location field**
> (the pole positions resolved separately to full `1e-8+` accuracy — the
> scientific content), with an **honest continuous surface only where the
> A1 discs genuinely cover** (~13–18 % at display tolerance) and explicit
> `NaN`/unrendered elsewhere. A filled honest wedge surface to `|x|=20`
> requires the Phase-B re-scope flagged in §5.3 — a genuine 2D-fill
> architecture (independent grid-seeded local expansions, no pole-hitting
> bridging walks) — not a tolerance or order tweak.

### 7.5 One-line summary for the ADR-0025 amendment

> Phase-A2b (display-tolerance follow-up): relaxing the A1 gate to display
> tolerance (`1e-3`) lifts honest union coverage of the `|x|≤20` P_I⁽²⁾
> tritronquée wedge from **8.6 % → 13.5 %**; adding order-48 jets reaches
> **16.9 %**, and the loosest display tol `1e-2` peaks at **18.0 %** — a
> 2.1× gain over the solver-tol baseline that still leaves the figure
> **~82 % NaN**, with **no annulus past 50 %** and the 50 %-frontier at
> `|x|≈0`. An order-48 *walk* additionally **BLOCKS** (the shared-Q Padé
> solve degenerates). **Verdict: a filled honest wedge *surface* is not
> viable at display tolerance or order 48; the architectural rescope
> stands** — render the wedge as an extracted-and-validated pole-location
> field (poles to full `1e-8+`), honest continuous surface only on the
> ~13–18 % the A1 discs cover, `NaN` elsewhere. A truly filled honest
> surface needs a 2D-fill architecture (§5.3), a Phase-B re-scope.
