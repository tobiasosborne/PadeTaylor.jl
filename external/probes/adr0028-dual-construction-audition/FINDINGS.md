# ADR-0028 audition — FINDINGS

**Status:** prototype + audition complete (probe only; `src/` untouched, build
stays gated). Bead `padetaylor-0ql3`. Feeds the seven sign-off decisions in
`docs/adr/0028-sharedpade-dual-construction-pareto-dispatch.md`.

**What was built (`external/probes/adr0028-dual-construction-audition/`):**
`cells.jl` unifies cell (A) and the cell-(B) variants into one parameterised
block-Toeplitz solver (cell (A) asserted **bit-identical** to the shipped
`shared_denominator_pade`, Δ=0 — `smoke.jl`), plus the held-out-residual metric
and the Pareto/lexicographic selector. `battery.jl` builds jets via the
package's own `vector_taylor_coefficients` (rescaled `h^k` as `VectorStepper.jl:237`)
for harmonic, exp, ℘, Calogero–Moser, Noumi–Yamada A₄. `audition.jl` +
`addendum.jl` + `addendum2.jl` run F64 + BigFloat-256 + pole-crossing probes.

Three cell-(B) realisations were auditioned (the Mano–Tsuda scout flagged the
`d∤m` row-count as *source-silent* — it turned out decisive):

| cell | window | rows/block | denom degree (d∤m) | stack shape |
|---|---|---|---|---|
| **A** (shipped GGT) | `+2` | `m` | `m` | `(d·m)×(m+1)` over-determined |
| **B_ceil** (ADR's literal cell B) | `+1` | `⌈m/d⌉` | `m` (kept) | `(d·⌈m/d⌉)×(m+1)` **over-determined when d∤m** |
| **B_floor** | `+1` | `⌊m/d⌋` | `⌊m/d⌋·d` (lower) | `m_eff×(m_eff+1)` wide-square |
| **B_grow** | `+1` | `⌈m/d⌉` | `⌈m/d⌉·d` (≥ m) | `m_eff×(m_eff+1)` wide-square |

---

## Headline results

### 1. cell (B) recovers the worklog-067 accuracy cost — by 6–8 orders (F64)

True per-step error vs oracle (Float64), cell A vs the best square cell:

| system | kind | d | err_A | err_B (square) | recovery |
|---|---|---|---|---|---|
| harmonic | entire | 2 | 5.22e-9 | **1.4e-15** | ×3.6e6 |
| exp-pair | entire | 2 | 2.34e-9 | **2.2e-16** | ×1.1e7 |
| ℘ (regular step) | meromorphic | 2 | 3.93e-4 | **8.2e-13** | ×3.3e7 |
| Calogero–Moser | meromorphic | 4 | 9.66e-11 | **4.4e-16** | ×2e5 |
| Noumi–Yamada A₄ | meromorphic | 5 | 4.66e-14 | 5.6e-17 | ×840 |

The `err_A = 5.22e-9` on harmonic reproduces worklog-067's `~5e-9` exactly,
confirming the harness is faithful. **The square Mano–Tsuda cell recovers the
entire/high-d cost the ADR set out to recover — confirmed.**

### 2. 🔴 The ADR's literal cell B (`⌈m/d⌉`, keep-degree) is the WRONG realisation

When `d∤m`, keeping degree `m` with `⌈m/d⌉` rows makes the stack
**over-determined** (`d·⌈m/d⌉ ≥ m+1`) — re-introducing *exactly the
least-squares pathology cell (B) was meant to cure*:

| system (d∤m) | err_B_ceil | err_B_floor | err_B_grow |
|---|---|---|---|
| exp-pair (d=2,m=15) | **8.9e-7** | 2.2e-16 | 6.7e-16 |
| Calogero–Moser (d=4,m=15) | **2.5e-4** | 3.4e-13 | 4.4e-16 |
| Noumi–Yamada A₄ (d=5,m=12) | **1.0e-10** | 5.6e-17 | 1.7e-16 |

B_ceil is **×4e9 / ×7e8 / ×2e6 worse** than the genuinely-square variants. The
cell-B definition must be a **wide-square** system (`d·n = m_eff`, guaranteed
1-D null), not the over-determined `+1` analogue of A.

**Of the two square variants, `B_grow` (`⌈m/d⌉·d`, higher degree) is best:** it
matches B_floor on entire jets and beats it on meromorphic ones (CM 4.4e-16 vs
3.4e-13), because the extra degree captures more pole structure while staying
clean-square.

### 3. ✅🔴 A per-step dispatch IS justified — but the win-boundary is {pole × precision}, not entire-vs-meromorphic

Stress-testing "is B always best?" at **BigFloat-256** found the counterexample:

| system | precision | err_A | err_B_grow | winner |
|---|---|---|---|---|
| Calogero–Moser | **BF-256** | **2.3e-22** | 2.1e-19 | **A ×900** |
| harmonic | BF-256 | 1.6e-17 | 1.8e-30 | B ×1e13 |
| Noumi–Yamada A₄ | BF-256 | 4.7e-14 | 5.6e-17 | B ×840 |

Cell A wins **only** when a step has *genuine shared-pole structure* **and**
the precision is high enough that roundoff no longer masks A's clean-isolated-
null degree-`m` advantage (CM has real off-axis poles; A recovers them to
2.3e-22 where B_grow's ν=16 stalls at 2.1e-19). At F64 the same CM step favours
B (roundoff floor ~1e-13 swamps A's advantage). NY's step has *no* genuine pole
(regular stretch) → A does **not** improve at BF-256 (stuck at 4.7e-14) → B wins.

This **vindicates the ADR's core instinct** (sometimes A, sometimes B) but
**re-characterises the boundary**: it is `{genuine pole} × {precision}`, and the
A-win regime — *genuine poles at extended precision* — is precisely the FW 2011
Table 5.1 BF-256 long-range showcase. "Always B" is therefore **not** safe.

### 4. Pole-crossing alone does NOT favour A (transcendental tan, single pole)

| step | err_A=(m,m) | err_B=(m-1,m) |
|---|---|---|
| h=1.8 (PAST π/2), m=15 | 5.1e-11 | **3.8e-13** |
| h=2.0 (PAST π/2), m=15 | 1.1e-10 | **1.1e-12** |

Across a single transcendental pole at F64, the off-diagonal (m-1,m) ties or
**beats** the diagonal (m,m); at higher m the over-determined diagonal actually
degrades (Froissart-like). Multi-pole single steps (h=5,6) break *both* cells
(O(1) error) — not an A/B discriminator (FW always use small steps). So A's
advantage is *clean-null degree at high precision*, not pole-bridging per se.

---

## Mapped to the seven sign-off decisions

1. **Objective ordering / ε-bands — 🔴 the (R,g,K) Pareto is unreliable; ✅ a
   validated replacement exists (held-out POINT).** The ADR-as-written mis-ranks
   **6/8** test cases (`discriminator.jl`). The single-coefficient residual `R`
   is fooled (an over-determined LS fit nails one held-out coefficient while its
   function is inaccurate); the `ε_rel = 100·eps` band is *relative to R* (~1e-30,
   far below R's ~1e-16 noise floor) so on near-exact cells it selects on noise
   and defers to `g`, which is itself unreliable (item 3). **Recommendation: drop
   the (R,g,K) Pareto. Use the held-out-POINT discriminator** (`held_out_point`):
   evaluate each cell's `P_i(t*)/Q(t*)` at a point `t*` inside the Taylor radius
   and compare to the full-jet Taylor sum (a more-accurate local surrogate truth);
   pick the cell closer to it. It integrates ALL held-out coefficients into a
   function value, so the single-coefficient pathology vanishes. **Scored 8/8**
   across the F64 battery + BigFloat-256 — including correctly choosing **A** on
   the BigFloat CM genuine-pole case (where pure-R and the ADR pick B, wrongly)
   and **B** on BigFloat harmonic (where pure-R picks A, wrongly). **Open design
   item:** `t*` must be inside the radius (no in-step pole); for genuine
   pole-crossing vector steps the surrogate truth diverges — needs a radius check
   + fallback (e.g. default to A on a detected in-step pole, consistent with the
   F64-tan tie + BF-CM A-win). Not yet prototyped.

2. **Neighbouring-degree breadth — superseded by the cell-B correction.** The
   live question is not ±1 neighbours; it is *which square realisation*. **Use
   `B_grow` (`⌈m/d⌉·d` wide-square, degree `m_eff ≥ m`).** Amend ADR §1(B): the
   over-determined `⌈m/d⌉`-keep-degree `(m-1,m)` is wrong.

3. **Always-both vs conditional `g_A` gate — 🔴 conditional NOT supported.** The
   conditioning gap `g = σ_m/σ_{m+1}` does **not** separate the regimes: across
   the battery the "A-worse" and "A-ok" `g_A` ranges overlap entirely
   (`log10 g_A`: harmonic 1.5, exp 3.9, ℘ 1.7, CM 4.5, NY 1.7 — no threshold
   separates the A-wins case). Worse, `g` is **not consistently defined across
   cell shapes**: A is tall (`σ_{m+1}` = the null), B_grow is wide-square (its
   null is structural, not in the returned σ-vector) — so `g_A` and `g_B` measure
   different things and are not comparable. **Recommendation: build both always**
   (the cost is one extra small SVD, item 5) **or** gate on precision + genuine-
   pole detection, not `g`.

4. **Diagnostic surfacing — ✅ keep.** `pole_disagreement` between the A/B root
   sets tracks regime (entire 0.2–0.4, meromorphic 0.5–0.9); the opt-in 3-tuple
   `(numerators, denominator, ::SharedPadeDispatch)` is fine and breaks no
   `VectorStepper.jl:242,269` destructure. Add `chosen_cell`, `held_out_residual`,
   `pole_disagreement`, `reduced_degree`.

5. **Performance budget — ✅ acceptable, lean to always-both.** One extra
   `O(m_eff³)` SVD on a small (`m_eff×(m_eff+1)`) matrix per d≥2 step. Since the
   `g` gate (item 3) is unreliable, **always-both** is the safe default; revisit a
   cheaper gate only if profiling demands it.

6. **Held-out combine rule (max across d) — fine in principle, moot until the
   metric is fixed (item 1).**

7. **d=1 mechanism — ✅ short-circuit to A confirmed needed.** At d=1, cell B =
   `robust_pade(m-1,m)` ≠ cell A = `robust_pade(m,m)` (μ_A=6 vs μ_B=5 in
   `smoke.jl`); the bit-identity contract (SP.1.1 + ~49 tolerances) requires the
   unconditional `d==1 && return A` short-circuit. Keep exactly as ADR §4.

---

## Discriminator prototype (`discriminator.jl`) — the redesigned selector

Four selection rules scored against the ground-truth oracle winner (cell pair
A vs the corrected B_grow), F64 battery + BigFloat-256:

| rule | score | notes |
|---|---|---|
| ADR-as-written (R→g→K→A) | 6/8 | fails BF-harmonic + F64-harmonic (g/Rfloor deferral) |
| pure-R (argmin held-out residual, no g/floor) | 6/8 | fails BF-harmonic + F64-harmonic (single-coeff fooled) |
| **held-out-POINT** | **8/8** | correct on every case incl. BF-CM (→A) and BF-harmonic (→B) |
| hybrid (point primary, R tiebreak) | 8/8 | equivalent to point on this battery |

The held-out-point rule is the recommendation: it is the only tested rule that
correctly fires **A** in the genuine-pole×precision regime *and* **B**
everywhere else.

## Bottom line for the maintainer decision

- **Cell B must be the higher-degree wide-square Mano–Tsuda (`⌈m/d⌉·d`), not the
  ADR's over-determined `⌈m/d⌉`-keep-degree.** (Corrects the source-silent gap.)
- **The dispatch is justified** — a genuine A-win regime survives (genuine poles
  × extended precision, e.g. the CM BF-256 case central to FW Table 5.1), so
  "always B" is unsafe.
- **Replace the (R,g,K) Pareto with the held-out-POINT discriminator** (8/8 vs
  6/8). Drop the conditioning gap `g` entirely (item 3). Add the in-step-pole
  radius check + fallback as the one remaining design item.
- d=1 short-circuit, opt-in diagnostics, and the always-both cost are all fine.

## Untested regimes (before fully retiring/committing the design)

- The actual multi-step **path-network walk** (accumulating error, adaptive `h`),
  not just isolated single steps.
- The **49-tolerance contract**: switching the d≥2 path to B_grow lets the
  worklog-067 tolerances tighten back — must confirm no test asserts the loose
  values as a *lower* bound.
- B_grow **degenerate guards** (e.g. `d > m`, an identically-zero component jet,
  bead `padetaylor-0o9`); the prototype skips the ADR-0027 reduction loop for the
  square variants.
- A **redesigned discriminator** for the genuine A-win regime (item 1) — not yet
  prototyped.
