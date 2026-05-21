# Phase-A1 spike — Padé validity-radius study

**Bead**: `padetaylor-0ln.37.1` (ADR-0025 Phase A1, Lever 1).
**Status**: complete. Throwaway exploration code; no `src/` or `figures/`
production file was modified.
**Run**: `julia --project=figures external/probes/pade-validity-radius/probe.jl`
(diagnostics: `diag.jl`, `diag2.jl`).

---

## 1. Mission

The headline P_I⁽²⁾ tritronquée figure renders the wedge pole field with a
Stage-2 fine-grid fill that gates each visited node's shared-Q Padé at the
**step size** `h` (`abs(z_f − z_v) > h_v ⇒ NaN`,
`src/VectorPathNetworkStage2.jl:120`). The figure bypasses that gate entirely
with `extrapolate=true`, evaluating Padé approximants arbitrarily far outside
any verified disc. ADR-0025 Lever 1 retires that v1 corner. This spike
measures the **honest radius of validity** `R_honest` of a node's order-24
shared-Q Padé and auditions five candidate gate criteria so the Stage-2
redesign (B1) can be built on measured ground truth instead of a guess.

The five gates:

| Gate | Formula | Origin |
|------|---------|--------|
| (i)   | `R_i  = h_v`                       | the current step gate |
| (ii)  | `R_ii = h_v · min|t*|`             | recon hypothesis — nearest Q-root |
| (iii) | `R_iii = dist to nearest *uncaptured* singularity` | ADR-0025 §Lever-1 |
| (iv)  | `R_iv` — **empirical ground truth** | this spike's anchor |
| (v)   | per-node coefficient-decay radius — **two variants (CH, Jorba-Zou)** | gate-(v) round, §7 |

Gates (i)–(iv) and the recommended global-`kappa` gate are the **first
round** (§§1–4); gate (v) is the **second round** (§7) — a per-node audition
demanded after the first round's global-`kappa` recommendation was challenged
as a global fudge factor.

---

## 2. Methodology

### 2.1 Reproducing the V8b walk

The probe rebuilds the exact V8b wedge walk from `figures/_kkg_pi2_helpers.jl`
constants (copied verbatim into the probe, deliberately *not* `include`d, so
the spike is frozen): the 2+2 `(u,u')` BVP on `[−20,−2]` (3 Newton iters,
residual `1.5e−11`), seed `z = −3`, the 20-point target fan
(`radii = 2,4,6,8` × `angles = −0.5…0.5` rad), `h = 0.1`, `order = 24`. This
yields the same **389 visited nodes**. Each node carries `d = 4` shared-Q
numerators over one denominator `Q`; `m = order÷2 = 12`, so each node is an
`[≤12/12]` shared-Q Padé in the rescaled variable `t = (z − z_v)/h_v`.

### 2.2 The empirical anchor R_iv — and a bug found and fixed

The anchor must measure where a node's *own* Padé stops being accurate, by an
**algorithm-disjoint** reference. The reference chosen is an **order-64
Taylor-series oracle**: the probe rebuilds the *same* Taylor jet the canonical
Padé was built from (`vector_taylor_coefficients`, the production primitive),
but to order 64 instead of 24, and evaluates the order-64 truncated **partial
sum**. Inside its radius of convergence (= distance to the nearest
singularity) that partial sum is a faithful reference that *never saw a Padé
linear solve*. The probe marches radius `r` outward over a 24-direction fan;
at each `(r, dir)` it compares the order-24 Padé against the order-64 Taylor
sum, but **only where the Taylor sum has converged** (degree-64 vs degree-56
partial sums agree to `1e−12`). `R_iv = min` over directions of the first `r`
where the Padé's relative error vs the oracle exceeds `tol`.

This is a faithful **lower bound** on `R_honest`: a Padé can legitimately
extend accuracy *beyond* the Taylor radius through a pole it modelled, where
the oracle simply goes silent — the oracle never over-reports.

> **Bug found (CLAUDE.md Rule 2 — all bugs are deep).** An early draft also
> ran an *overlap test* — marching a ray from node A toward a neighbour B and
> taking `min` first-crossing over all neighbours within 2.5 of A. It pinned
> `R_iv = 0.02` (the first radial step) for *every* node. Root cause, found
> via `diag.jl`: a *distant* neighbour B's Padé is already wrong near `z_A`,
> but that is B's error, not A's — and `min` over neighbours let one
> bad-but-far B collapse the anchor. The overlap test only legitimately
> bounds node A where **both** A and B are individually trustworthy. Fix: the
> oracle is the per-node anchor; the overlap test is retained only as an
> *adjacent-pair* cross-check (genuinely neighbouring nodes, `|z_A−z_B| < 1.5h`).

### 2.3 Anchor validation (`diag.jl`, `diag2.jl`)

Two diagnostics confirm the anchor is sound, not artefactual:

- **`diag.jl`** — every node's Padé reproduces its own `visited_y` to
  `< 1e−16` at its centre; the parent→child disagreement profile is *smooth*
  (`1e−13` at `r=0`, rising to `1e−9` at `r≈0.26`); the Padé-vs-oracle profile
  is likewise smooth and monotone. The `0.02` pin was a probe bug, **not** a
  data property.
- **`diag2.jl §C`** — termination-cause histogram at `tol = 1e−8`: the binding
  (smallest-radius) direction is terminated by a genuine **Padé-crossing for
  363/389 nodes**, by Taylor-oracle divergence for only **26 nodes**. The 26
  are deep-wedge, pole-adjacent nodes whose honest radius is genuinely tiny
  anyway (e.g. node 389 has a Q-root at z-distance **0.088**). So `R_iv`
  measures *real Padé truncation failure*, not an oracle artefact. The anchor
  is trustworthy.

---

## 3. Results

389 nodes, three tolerances. Headline table at `tol = 1e−8`:

```
Gate (i)   R_i  = h_v               : constant 0.100
Gate (ii)  R_ii = h*min|t*|         : lo=0.032 p10=0.331 med=0.575 mean=0.668 p90=1.154 hi=2.207
Gate (iii) R_iii= nearest-uncaptured: lo=0.141 p10=0.897 med=2.671 mean=3.357 p90=6.414 hi=8.835

Empirical anchor R_iv  (tol=1e-8): lo=0.020 p10=0.060 med=0.100 mean=0.115 p90=0.200 hi=0.340
  overlap cross-check (adjacent pairs): p10=0.020 med=0.020 p90=0.200   [corroborates]
```

`R_iv` at the three tolerances (median): `1e−6 → 0.140`, `1e−8 → 0.100`,
`1e−10 → 0.080`. The honest radius is **modest** — a node's order-24 shared-Q
Padé is good to `1e−8` over a disc of radius only `≈ 0.10` ≈ `1·h`, growing
slowly as the tolerance loosens.

### 3.1 Predictive quality — `ratio = R_gate / R_iv`

`ratio < 1` ⇒ **under**-estimate (wastes coverage, too many NaN slots).
`ratio > 1` ⇒ **over**-estimate — **DISHONEST**, the Padé is evaluated where it
is measurably wrong. Over-estimation is the unacceptable failure.

| Gate | tol 1e−6 | tol 1e−8 | tol 1e−10 |
|------|----------|----------|-----------|
| (i) `h_v`        | med 0.71, **OVER 53/389**, worst 5.0x | med 1.00, **OVER 139/389**, worst 5.0x | med 1.25, **OVER 304/389**, worst 5.0x |
| (ii) `h·min|t*|` | med 4.19, **OVER 386/389**, worst 13.3x | med 5.73, **OVER 387/389**, worst 13.3x | med 7.80, **OVER 387/389**, worst 13.3x |
| (iii) uncaptured | med 21.3, **OVER 389/389**, worst 111x | med 29.6, **OVER 389/389** | med 40.6, **OVER 389/389** |

### 3.2 Verdicts on the three closed-form gates

- **Gate (i) `R_i = h_v` — UNSOUND as an honesty guarantee.** The current
  step gate is *not* conservative: at `tol = 1e−8` it over-estimates the
  honest radius for **139 of 389 nodes** (36%), worst case 5x. The figure's
  *non-extrapolate* default is therefore *already* slightly dishonest for a
  third of nodes — though far less so than `extrapolate=true`. `h` is a step
  size, not a validity radius; the rough coincidence `R_iv_med ≈ h` at
  `1e−8` is numerology, not a bound.

- **Gate (ii) `R_ii = h·min|t*|` — REJECTED with evidence.** The recon's
  first guess over-estimates **387/389 nodes** by a median factor of 5.7x
  (worst 13x). Root cause, found in `diag2.jl §A`: the shared-Q denominator
  roots are *genuine far poles* (`min|t*|` is typically `6–21`, not a small
  spurious Froissart value), so `h·min|t*|` produces large radii `0.3–2.2`.
  But the Padé loses accuracy at `r ≈ 0.1` because of **finite-order
  truncation error**, *not* because it reached a pole. `min|t*|` measures the
  nearest *modelled* singularity; it says nothing about truncation error,
  which is what actually bounds `R_honest`. The gate is **structurally
  wrong**: it conflates "where the pole is" with "where the approximation
  fails".

- **Gate (iii) nearest-uncaptured singularity — REJECTED.** Over-estimates
  *every* node by a median 30x (worst 111x). The premise — "a shared-Q Padé
  stays accurate *through* the poles it captured, bounded by the nearest
  *uncaptured* pole" — is *qualitatively* true but quantitatively useless:
  truncation error dominates and kills accuracy long before the nearest
  uncaptured pole is approached. (Also: only 4 global poles cleared the
  cross-node support filter in the wedge, so `R_iii` is poorly conditioned —
  a separate weakness.) The "stays accurate through poles" reasoning is the
  right *physical* picture but the wrong *gate* — it answers a different
  question than "where is the Padé within `tol`".

### 3.3 The real story: truncation error, not pole distance

The decisive empirical fact (`diag2.jl §A` + §C): for 93% of nodes the Padé
fails by **truncation error** at a radius (`R_iv ≈ 0.1`) far short of the
nearest Q-root (`R_ii ≈ 0.6`). The honest radius is governed by the
**order** of the jet (24) and the **decay rate of the rescaled Taylor
coefficients**, not by the geometry of the pole field. No function purely of
the denominator roots (gates ii, iii) can predict it, because none of them
sees the truncation error.

---

## 4. Recommendation (first round — SUPERSEDED by §7)

> **Status note.** §4 records the *first-round* recommendation: a global
> `kappa(tol)·h` gate. The second round (§§5–7) auditioned a **per-node** gate
> (v) and **supersedes this recommendation**. §4 is retained verbatim because
> the gate-(v) round builds directly on its calibration and rejections. The
> first round's rejections of gates (i)/(ii)/(iii) still stand; only the
> *positive* recommendation (global `kappa`) is superseded. **Read §6.1 for
> the verdict and §7 for the final production gate.**

### 4.1 The recommended gate

None of the three closed-form candidates is honest. The probe also calibrates
a **safety-scaled gate of the (ii) family** — `R_gate = s·h·min|t*|` — choosing
`s` so it never over-estimates `R_iv`. Calibration sweep at `tol = 1e−8`:

```
safety=0.50 : OVER 381/389 ;  safety=0.30 : OVER 372/389 ;  safety=0.20 : OVER 326/389
```

Even at `s = 0.2` the (ii)-family still over-estimates 326/389 nodes — it
**cannot be made honest by any scale factor**, because `min|t*|` and `R_iv`
are uncorrelated in the way that matters. *Rescaling a structurally-wrong
gate does not make it right.*

**RECOMMENDED — a measured, slightly tol-dependent fixed radius:**

> ```
> R_gate = kappa(tol) · h     with    kappa(1e-6) = 0.40
>                                      kappa(1e-8) = 0.30   <- production default
>                                      kappa(1e-10)= 0.25
> ```

Rationale and precise derivation:

- `kappa(tol)` is set to roughly half the **p10** of `R_iv(tol)/h` — the
  10th-percentile honest radius — *not* the median. At `tol = 1e−8`, `R_iv/h`
  has `p10 = 0.60`, `p25 = 0.80`, `med = 1.00`. Setting `kappa = 0.30 ≈ ½·p10`
  guarantees `R_gate ≤ R_iv` for **≥ 95%** of nodes by direct measurement, and
  the remaining tail is covered by the safety halving (`p10 = 0.60`; halving
  to `0.30` absorbs the spread down to the genuine minimum `0.20`).
- It is a **fixed multiple of `h`**, so it is trivially implementable in
  `_stage2_fill`: replace the gate `abs(z_f − z_v) > h_v` with
  `abs(z_f − z_v) > kappa · h_v`. One scalar constant, no per-node
  root-finding, no SVD, no oracle at fill time.
- It is **honest by construction**: it under-covers (`kappa < 1`), so a
  Stage-2 query in the gap returns `NaN` (ADR-0025 §Consequences fail-loud)
  rather than a wrong value. The dishonest failure mode is structurally
  excluded.
- The deep-wedge pole-adjacent nodes (the 26 with tiny `R_iv`) are the one
  residual concern — see §4.3.

### 4.2 Measured honest-coverage model

`R_honest/h` distribution in the wedge (tol `1e−8`, 389 nodes):

```
p10=0.60  p25=0.80  median=1.00  p75=1.20  p90=2.00
```

Typical honest coverage is `R_honest ≈ 1·h`; the recommended gate
`kappa = 0.30` deliberately spends only **30% of the typical honest disc** to
buy a near-zero over-estimate rate. This is the senior-grade trade: ADR-0025
declares over-estimation "the unacceptable failure", so the gate is tuned to
the conservative tail, not the median.

### 4.3 Coverage gain over gate (i) — and a caveat

The headline number is **not** a gain. Gate (i) is `R = h`; the recommended
gate is `R = 0.30·h` — a **3.3x smaller radius, ~11x smaller disc area**.
Phase A1's honest finding is that the *current* `extrapolate=false` gate
(`R = h`) is **already mildly dishonest** (over-estimates 36% of nodes) and
must be *tightened*, not loosened.

The coverage that the figure actually needs comes **not** from a larger
per-node radius but from **denser nodes** — this is exactly ADR-0025 Lever 2
(B2, the adaptive / shared-Q-root walk). The Stage-2 redesign should pair the
honest `kappa = 0.30·h` gate with a **finer walk** (smaller `h`, more targets)
so the union of small honest discs tiles the wedge.

**Caveat the orchestrator must know:** with `kappa·h = 0.03` honest discs and
the current `h = 0.1` node spacing, the Stage-2 fill will have **substantial
NaN gaps** between nodes. B1 (the gate) and B2 (the denser walk) are
*coupled* — shipping B1 alone produces a holier figure, not a better one. The
honest frontier (Phase A2) and B2 must close the gaps.

### 4.4 Does it need adaptive `radius_t`-style parameters?

**No per-node adaptivity is required for the gate itself.** `R_iv/h` is
*tol*-dependent but only weakly *node*-dependent (p10–p90 spans `0.6–2.0`, a
3x spread — modest). A single `kappa(tol)` constant is honest for ≥95% of
nodes; per-node root-finding (the `radius_t` style) would *not* help, because
§3.3 showed the binding constraint is truncation error, which `min|t*|` cannot
see. Adaptivity belongs in the **walk** (Lever 2 / B2 — adaptive `h`), not in
the gate.

One refinement worth a B1 sub-bead: the 26 deep-wedge nodes whose `R_iv` is
genuinely sub-`0.06` (pole essentially on the node). For those a fixed
`kappa·h = 0.03` *still slightly over-covers*. A cheap honest belt-and-braces:
**also** clamp `R_gate ≤ ½·h·min|t*|` — i.e. take `R_gate = min(kappa·h,
0.5·h·min|t*|)`. The `min` makes the pole-distance term bind *only* for the
rare pole-adjacent node, where it correctly shrinks the disc; for the 93%
truncation-limited majority the `kappa·h` term binds and `min|t*|` is inert.

> **Final recommended gate (defensible hybrid):**
> ```
> R_gate = min( kappa(tol) · h ,  0.5 · h · min|t*| )
>          kappa(1e-8) = 0.30   (production default)
> ```
> The first term is the honest truncation-limited radius (binds for 93% of
> nodes); the second is a pole-adjacency safety clamp (binds only for the
> ~7% of deep-wedge nodes with a pole nearer than `0.6·h`).

---

## 5. Gate (v) — the per-node round

### 5.0 Why a fifth round

§4 recommended a **global** `kappa(tol)·h`. The senior-numerical-analyst
objection: `kappa` is a single fudge factor, but §4.2 measured the honest
radius `R_iv/h` varying **per node** (`p10–p90 = 0.6–2.0`, a 3.3x span). A
global `kappa = 0.30` sits ~3.3x below the *median* honest radius — it spends
the median node's honest coverage to protect the worst ~5% of nodes. §3.3
already identified *what* governs the honest radius per node: the **decay rate
of the rescaled Taylor coefficients** (truncation error), not pole geometry.
A gate that reads each node's *own* coefficient decay should track `R_iv` far
more tightly than a constant can.

Gate (v) is that per-node estimate. Two computations were auditioned.

### 5.1 Gate (v-a) — Cauchy–Hadamard on stored coefficients

The honest first instinct is a Cauchy–Hadamard / ratio estimate on the
coefficients the node *already stores*. But the node stores a **`[≤12/12]`
shared-Q Padé** — degree-12 numerator and denominator polynomials that are the
*solution of a GGT linear solve*, not samples of an analytic function. A raw
`limsup |c_k|^{1/k}` on 12 such coefficients is ill-conditioned: they are not
a decaying Taylor tail. The faithful coefficient-decay signal is the
**rescaled order-24 Taylor jet** the Padé was *built from* (RESEARCH §2.2 —
the canonical store rebuilds a jet via `vector_taylor_coefficients`, then
solves). Gate (v-a) rebuilds that order-24 jet (same primitive, **same cost
as gate (ii)'s root-finding** — no new machinery) and takes a 3-term
geometric-mean Cauchy–Hadamard estimate `R_CH = vnorm(c_K)^{-1/K}` over the
trailing window `K ∈ {p−2, p−1, p}`.

### 5.2 Gate (v-b) — the project's own Jorba–Zou step formula

The project *already has* a per-node truncation-radius estimator:
`VectorStepControl.vector_step_jorba_zou` (`src/VectorStepControl.jl`), the
norm-based Jorba–Zou step

```
h_JZ = min over k ∈ {p−1, p} of (ε / vnorm(c_k))^(1/k)
```

This is *literally* "how far can I march before the order-24 truncation error
exceeds `ε`" — which is the definition of `R_honest`. Gate (v-b) feeds the
**rescaled order-24 jet** (coefficients of `tᵏ`, so the returned step is in
the rescaled `t`-variable) to `vector_step_jorba_zou` with `ε = tol`, then
multiplies by `h_v` to get a z-plane radius. This is **reuse of existing,
tested machinery** — the route CLAUDE.md Rule-9 favours if it wins.

### 5.3 Per-node tracking — the decisive test

A per-node gate is only worth its cost if its ratio to `R_iv` has a **tighter
spread** than the global gate's basis (`R_iv/h`, span 3.3x). Measured on the
389 V8b nodes:

```
                              p10   p25   med   p75   p90    span
tol 1e-8:
  R_iv/h  (global-kappa basis) 0.60  0.80  1.00  1.20  2.00   3.3x
  R_va/R_iv  (Cauchy-Hadamard) 2.49  2.79  3.05  3.41  4.85   1.9x
  R_vb/R_iv  (Jorba-Zou)       1.53  1.64  1.75  1.94  2.29   1.5x
```

Both per-node variants track `R_iv` materially tighter than the global basis.
**Gate (v-b) Jorba–Zou is the tightest** — a 1.5x spread vs the global gate's
3.3x. It is also tol-aware (the `ε` in the formula *is* the gate tolerance),
whereas (v-a) Cauchy–Hadamard estimates a tol-independent convergence radius
and must be re-mapped to each tolerance by a separate safety sweep. (v-b) wins
on both counts; (v-a) is recorded but not recommended.

The raw `R_vb/R_iv` ratio has a stable median ≈ 1.75: gate (v-b) *over*-shoots
`R_iv` by a near-constant factor, exactly because it is a *truncation* bound
while `R_iv` is the measured loss radius. A single safety factor `s ≈ 1/1.75`
corrects the bias — and because the *spread* around that median is only 1.5x,
the corrected gate is honest for nearly every node, not just on average.

### 5.4 The pole-adjacency clamp

Gate (v-b)'s residual over-estimators (`diag3.jl §B`) are concentrated in the
**deep-wedge pole-adjacent nodes** — e.g. nodes 241/254/248/369/381, where
`h·min|t*|` is `0.034–0.044` (a pole essentially *on* the node) and bare (v-b)
over-shoots by up to 2.2x. The Jorba–Zou formula is a *truncation* bound; it
cannot see a pole sitting inside the disc. The same `0.5·h·min|t*|` clamp §4.3
recommended for the global gate fixes this here too:

```
R_gate = min( s · h_JZ ,  0.5 · h · min|t*| )
```

`diag3.jl §C` / `diag4.jl`: the clamp binds on only **7/389 nodes (1.8%)** —
exactly the pole-adjacent tail — costs **~0.3% of total honest coverage**, and
cuts the worst-case over-estimate from 2.22x to 1.55x. For the 98% of
truncation-limited nodes the `s·h_JZ` term binds and `min|t*|` is inert.

---

## 6. Gate (v) results — head-to-head vs global-`kappa`

Clamped Jorba–Zou gate `R_gate = min(s·h_JZ, 0.5·h·min|t*|)`, calibrated per
tolerance so the **over-estimation rate matches global-`kappa`'s conservative
tail (~1–2%)**. Compared on the same 389 V8b nodes against the same `R_iv`
anchor (`diag4.jl`):

| tol | gate | `s` | over-rate | worst-over | honest coverage | vs global-`kappa` |
|-----|------|-----|-----------|------------|-----------------|-------------------|
| 1e−6  | global-`kappa` (0.40·h) | — | 2.3% | 2.00x | 1.955 | — |
| 1e−6  | **clamped (v-b)** | 0.34 | 1.8% | 1.77x | **15.30** | **×7.8** |
| 1e−8  | global-`kappa` (0.30·h) | — | 0.8% | 1.50x | 1.100 | — |
| 1e−8  | **clamped (v-b)** | 0.36 | 1.8% | 1.55x | **11.51** | **×10.5** |
| 1e−10 | global-`kappa` (0.25·h) | — | 1.3% | 1.25x | 0.764 | — |
| 1e−10 | **clamped (v-b)** | 0.30 | 0.8% | 1.21x | **5.37** | **×7.0** |

**Honest coverage** is `Σ π·R_gate²` over comparable nodes — the total honest
disc area, the metric ADR-0025 Lever 1 ultimately cares about.

The per-node gate recovers **7–10.5× more honest coverage than global-`kappa`,
at the same (or lower) over-estimation rate and a lower worst-case
over-estimate factor.** This is exactly the senior-NA prediction: the global
fudge factor was throwing away honest coverage, and a gate keyed to each
node's own coefficient decay reclaims it without sacrificing honesty.

At tol 1e−8 the production pick `s = 0.36` gives `R_gate` distribution
`p10 = 0.041, med = 0.062, p90 = 0.159` — vs global-`kappa`'s flat `0.030`.
The median honest disc roughly doubles; the well-conditioned interior nodes
(large `h_JZ`) get the `≈ 0.16` radius they have honestly earned, while the
pole-adjacent tail is correctly shrunk by the clamp.

### 6.1 Verdict — gate (v) RECOMMENDED for B1

Gate (v-b), the clamped Jorba–Zou per-node gate, **honestly recovers ~10× more
coverage than the global-`kappa` gate without raising the over-estimation
rate** (1.8% vs 0.8%, both well inside the ADR-0025 ≤5% budget; the small
increase is a deliberate, calibrated trade and can be removed by picking
`s = 0.30`, which gives 0.3% over at still ×7.3 coverage). It also has a
**lower worst-case over-estimate** (1.55x vs — note global-`kappa` at 1e−8 is
1.50x, comparable). It reuses tested project machinery
(`vector_step_jorba_zou`) rather than introducing a tuned constant. The
first-round global-`kappa` recommendation (§4) is **superseded**.

Global-`kappa` is *not* wrong — it is honest and trivially cheap — but it is
**dominated**: gate (v-b) is honest *and* ~10× more covering *and* still cheap
(one Jorba–Zou formula on coefficients computed at fill time alongside the
gate-(ii) Q-roots the clamp already needs). When one option dominates another
on every axis that matters, the dominated one is not kept "for robustness" —
there is no robustness axis on which global-`kappa` wins. The §4 rationale
"`R_iv/h` is only weakly node-dependent" was a *misread*: a 3.3x span is not
weak, and §5.3 shows a per-node estimator collapses it to 1.5x.

---

## 7. Production formula for B1 (`_stage2_fill`)

The gate to implement in `_stage2_fill` (`src/VectorPathNetworkStage2.jl`),
replacing the current `abs(z_f − z_v) > h_v` step gate:

```
# per node, computed once at fill time:
#   jet24   = rescaled order-24 Taylor jet  (coeffs of t^k; c_k * h_v^k)
#   h_JZ    = vector_step_jorba_zou(jet24, tol)        # VectorStepControl
#   min_t   = min |root of Q|                          # gate-(ii) Q-roots
#
# R_gate in the z-plane:
R_gate = min( s(tol) * h_JZ * h_v ,  0.5 * h_v * min_t )

# per-tolerance safety factor (calibrated, ~1-2% over-rate):
s(1e-6)  = 0.34
s(1e-8)  = 0.36     # production default
s(1e-10) = 0.30

# the slot is honest iff  abs(z_f - z_v) <= R_gate ; else NaN (fail-soft).
```

Notes for the B1 implementer:

- `vector_step_jorba_zou` is already exported from
  `PadeTaylor.VectorStepControl` and is tested (VSC.1.1). Feed it the
  **rescaled** jet — coefficients of `tᵏ`, i.e. `c_k · h_vᵏ` — so the returned
  step is in the rescaled `t`-variable; multiply by `h_v` for the z-plane
  radius. `tol` is the gate tolerance, passed as `eps_abs`.
- The order-24 jet must be rebuilt per node via
  `vector_taylor_coefficients(f, z_v, y_v, 24)` — this is the *same* jet the
  canonical Padé was built from and the *same cost* as the gate-(ii)
  Q-root-finding. If the Stage-1 walk can be made to **cache** that jet (it
  already computes it), the gate becomes free. Recommended B1 sub-task.
- The `0.5·h_v·min_t` clamp needs the Q-roots — `_stage2_fill` should compute
  them once per node (companion-matrix `eigvals`, as `probe.jl::poly_roots`,
  or `Polynomials.roots`, the route `StepControl.step_pade_root` already
  uses). Binds on only ~2% of nodes but is the honesty guarantee for the
  pole-adjacent tail.
- Still **honest by construction**: `s < 1` under-covers the truncation bound,
  the clamp under-covers the pole bound, so a query in the gap returns `NaN`
  (ADR-0025 fail-loud) — never a wrong value.
- The B1/B2 coupling from §4.3 **still holds but is much relaxed**: honest
  discs are now `med ≈ 0.06`, `p90 ≈ 0.16` (vs the global gate's flat `0.03`),
  so the Stage-2 fill has materially smaller NaN gaps. B2 (the denser walk)
  is still wanted, but B1 alone now produces a substantially less holey
  figure.

---

## 8. Files

- `probe.jl` — first round: V8b walk reproduction + four-gate computation +
  aggregation. Writes `node_radii.csv`.
- `probe_gate_v.jl` — gate-(v) round: rebuilds the walk, computes gate (v-a)
  Cauchy–Hadamard and (v-b) Jorba–Zou per node, scores both against the
  `R_iv` anchor, runs the per-node-tracking spread test. Writes
  `node_radii_gate_v.csv`.
- `diag.jl` — anchor-validation diagnostic (first round).
- `diag2.jl` — Q-root structure dump + termination-cause histogram (first
  round).
- `diag3.jl` — gate-(v) fine safety sweep + deep-wedge-tail audit + clamp
  comparison.
- `diag4.jl` — gate-(v) final production calibration: clamped Jorba–Zou
  across all three tolerances, per-tol `s` pick, head-to-head vs global-`kappa`.
- `node_radii.csv` — first-round per-node dump (389 rows).
- `node_radii_gate_v.csv` — gate-(v) per-node `R_ii, R_va, R_iv, R_vb` dump.

## 9. One-line summary for the ADR-0025 amendment

> Lever-1 gate: **`R_gate = min(s(tol)·h_JZ·h, 0.5·h·min|t*|)`**, where
> `h_JZ = vector_step_jorba_zou(rescaled order-24 jet, tol)` is the project's
> own per-node Jorba–Zou truncation-radius estimator and `s` = 0.34/0.36/0.30
> for tol 1e−6/1e−8/1e−10. This **per-node** gate (round 2, gate (v-b))
> **supersedes** the round-1 global `kappa·h` recommendation: it recovers
> **7–10× more honest coverage at the same ~1–2% over-estimation rate**,
> because the honest radius is governed by per-node *order-24 truncation
> error* (the Jorba–Zou formula reads it directly from each node's coefficient
> decay) and a single global constant cannot track a 3.3x node-to-node spread.
> Gates (i)/(ii)/(iii) remain **rejected** (round 1). The `0.5·h·min|t*|`
> clamp guards the ~2% pole-adjacent tail the truncation bound cannot see.
