# ADR-0028 — Shared-Padé dual-construction + validated-Pareto dispatch (the numerator-degree axis)

**Status**: proposed (design; awaiting maintainer sign-off before build, per
worklog 066 correctness-first directive). **Amended 2026-06-04 by an empirical
audition (Amendment 1, below) — three load-bearing claims of the original design
are corrected; read the amendment before the original §1–§6.**
**Date**: 2026-06-02 (Amendment 1: 2026-06-04)
**Beads**: `padetaylor-flnr` (this ADR), depends on the shipped
`padetaylor-d3a` (C1, `+2` window) + `padetaylor-3p9c` (C2, graceful
reduction) — i.e. **on ADR-0027**. Re-scopes `padetaylor-unk`. Audition:
`padetaylor-0ql3`.

## Context

Worklog 066 settled (unanimous, sympy-verified) that the shipped
`shared_denominator_pade` is the **GGT diagonal `(m,m)`** cell: `+2` window
(top-left `c_{m+1}`, matching `z^{m+1}…z^{2m}`), `m` rows/block, degree-`m`
numerators. It is *meromorphic-correct* (worklog 067: ℘, Calogero–Moser, KKG
exact, no throws). But worklog 066 identified **two genuine forks the math
cannot decide** (`066:94-103`), tied by an **order-vs-conditioning trade-off**
(`066:105-110`):

- the **numerator-degree axis** — GGT diagonal `(m,m)` [order `O(z^{2m+1})`;
  the over-determined `(d·m)×(m+1)` stack is rank-deficient away from poles]
  vs the Mano–Tsuda square `(m−1,m)` [`+1` window, `n=m/d` rows/block,
  degree-`(m−1)` numerators; lower order `O(z^{m+n})` but a square `m×(m+1)`
  system, better-conditioned on entire/regular jets];
- the **`d≥2` row-count fork** (Mano–Tsuda's minimal `n=m/d` vs the code's
  `m`) — *the same fork*: choosing `(m−1,m)` with `n=m/d` rows resolves both.

The shipped `(m,m)` core pays an accuracy cost on entire/regular jets that
**grows with `d`** (worklog 067 table, `067:58-78`):

| system | d | pre-fix `+1` | correct `+2` | tolerance now |
|---|---|---|---|---|
| harmonic `[cos,−sin]` (entire) | 2 | ~1e-10 | ~5e-9 | 1e-8 (VS.1.2/1.3, VP.1.2/1.4) |
| exp pair (entire) | 2 | ~1e-11 | ~4e-8 | 1e-7 (VS.1.3) |
| harmonic BigFloat-256 | 2 | ~3.4e-21 | ~1.6e-17 | 1e-16 (VS.1.2, VP.1.2) |
| **Noumi–Yamada A₄** Σf=t (meromorphic, regular stretches) | **4** | ~2e-9 | **~6.4e-6** | 1e-5 (NYF.1.2) |

The NY A₄ row is the real target: a high-`d` meromorphic system loses ~3
orders on its *regular stretches* between poles, not just on entire toys.
This ADR turns both forks into a **per-step, data-driven choice**: build both
table cells, validate each against held-out Taylor data, pick one on a small
Pareto frontier with a deterministic tie-break. It **composes with, does not
replace, ADR-0027**: ADR-0027 is the *denominator-degree axis* ("reduce =
select a lower table cell" along `m`, `0027:89-96`); this ADR adds the
*numerator-degree axis* ((A) vs (B)). Both cells get the ADR-0027 reduction.

## Decision (the design)

### 1. The two constructions, precisely

Both build from the same `d` rescaled Taylor jets `jets[i] = [c⁽ⁱ⁾₀,…]`
(`VectorStepper.jl:237`); `c⁽ⁱ⁾_j = 0` outside the available range
(zero-padding, `SharedPade.jl:128`). Both reduce along the denominator-degree
axis exactly per ADR-0027 (isolated-1-D criterion + `z^λ` cancellation,
`0027:39-62`).

**(A) GGT diagonal `(m,m)` — the shipped cell.** Verbatim the current
`shared_denominator_pade` (`SharedPade.jl:173-314`):
- *window* `+2`: block `[r,c]=c_{m+r-c+1}`, top-left `c_{m+1}`
  (`SharedPade.jl:124-133`, `idx=m+rr-cc+2`);
- *rows*: `m`/block → `A_full` is `(d·m)×(m+1)`;
- *numerator degree* `m` (`_upper_block` is `(m+1)×(m+1)`, `:142-151`);
- *recovery*: `Q=b` = smallest-σ right vector of `A_full`, QR-reweighted
  (`:254-268`); `aᵢ = upper_i·b`;
- *matching residual* `O(z^{2m+1})` (window ends at `z^{2m}`);
- *conditioning*: over-determined; full column rank (no exact null) on
  entire/regular jets ⇒ smallest-σ right vector is a *least-squares* `Q`
  placing spurious poles (source of the cost above). Exact-pole recovery is
  sharp (worklog 067).

**(B) Mano–Tsuda square `(m−1,m)` — the new cell.** From the primary source
`hp_arXiv_final.tex`: the simultaneous-Padé system (eq:lspa `:1405-1428`)
stacks `L−1` blocks `A^i_m(n,m+1)` of the rectangular Toeplitz of eq:rTm
`:1060-1080` at subscript `j=m`, so **top-left `a^i_m = c_m` (the `+1`
window)**; each block has `n` rows; `m = n(L−1)` (`:1381`) ⇒ the stack is
**square `m×(m+1)`** (RHS brace `:1423` counts `m` rows); numerators
`P_i^{(0)}` are degree `≤ m−1` (`:411,:399`); matching order `O(w^{nL}) =
O(w^{m+n})` (`:417`). In our framing the `d` stacked components are
Mano–Tsuda's `L−1` ratio blocks, so:
- *window* `+1`: block `[r,c]=c⁽ⁱ⁾_{m+r-c}`, top-left `c⁽ⁱ⁾_m` (one lower
  than (A); implementer note: `idx=m+rr-cc+1`);
- *rows*: `n=⌈m/d⌉`/block → `A_full ≈ m×(m+1)` (square, or minimally
  over-determined when `d∤m`: round `n` up so `d·n ≥ m` — one extra row is
  harmless and keeps the best-conditioned realisation without dropping
  equations);
- *numerator degree* `m−1` (`upper` block `m×(m+1)`, one row shorter than (A));
- *recovery*: same SVD + QR-reweighting on the square `A_full`; `Q=b`;
  `aᵢ = upper_i·b` with the degree-`(m−1)` block;
- *matching residual* `O(z^{m+n+1})` (lower order than (A) when `d>1`);
- *conditioning*: square/minimally-rectangular ⇒ far better conditioned on
  entire/regular jets (worklog-065 F1: `+1`/square recovers the regular
  harmonic step ~1 order better, no spurious high-degree numerator over-fit).

At `d=1`, (B) → `n=m` rows, same *shape* as (A) but `+1` window + degree-
`(m−1)` numerator = scalar GGT `(m−1,m)` (off-diagonal). (A) at `d=1` is the
scalar `(m,m)` — the bit-identity target (§3 shows (A) wins at `d=1`).

### 2. Held-out-residual validation metric (genuine cross-validation)

Score each candidate `(Q,{Pᵢ})` (at its ADR-0027-reduced degree) against the
**next unused Taylor coefficient** — one its window did *not* see. The cells'
windows end at different orders: (A) matches through `z^{2m}` → first unused
`2m+1`; (B) through `z^{m+n}` → first unused `m+n+1`. With matching order `K`,
the candidate's first unconstrained matching-equation residual per component is

```
        | coeff_{K+1}( P_i − f_i·Q ) |
rᵢ  =  ─────────────────────────────────────────────────
        ‖Q‖ · (Σⱼ |c⁽ⁱ⁾_j|, j over the matched window) + ε
```

i.e. the magnitude of the first unmatched matching residual, normalised by
`‖Q‖` (the convolution weight) × the component's matched-window coefficient
mass (so a quiet component is not over-penalised), `ε = tol·‖c‖`. Combine
across the `d` components by **max** (worst component governs the shared `Q`,
per ADR-0019's shared-pole philosophy); report the RMS alongside for the
diagnostic.

**Fairness across cells of different order.** (A) probes `2m+1`, (B) probes
`m+n+1 < 2m+1` for `d>1`. To compare on one yardstick, take the **common
held-out order `Kc+1 := min(K_A,K_B)+1`** (the highest order *both* leave
unconstrained) for the *selection-decisive* residual — apples-to-apples,
neither candidate saw it. Use each cell's *own* first-unused order only for
the `K`-objective. Both are one extra convolution coefficient per component
per cell — cheap.

### 3. Pareto selection + deterministic lexicographic tie-break

**Candidate set** `{(A),(B)}`, each ADR-0027-reduced, optionally with a small
degree neighbourhood (sign-off #2). Three objectives:

1. **held-out residual `R`** — minimise (the common-order metric of §2);
2. **null-space isolation/conditioning `g`** — maximise. From the
   *already-computed* SVD `S` (`SharedPade.jl:211` returns `S,Vt`), with
   `σ₁≥…≥σ_{m+1}`: `g = σ_m/σ_{m+1}` (smallest *retained* σ over smallest σ,
   both from the small end). **Larger `g` is better** — a cleanly isolated
   1-D null space ⇒ well-determined `Q`; small `g` flags the least-squares /
   ill-posed regime (A) hits on entire jets;
3. **order `K`** — maximise (`K_A=2m > K_B=m+n` for `d>1`).

**Pareto frontier** over `(R↓, g↑, K↑)`: drop dominated candidates (another
is `≤R`, `≥g`, `≥K` with one strict). **Deterministic lexicographic
tie-break** among the non-dominated set (C7: same input → same cell; no
`rng`, no float-equality ambiguity), with explicit ε-bands so jitter cannot
flip the choice:

1. **primary — `R`**, smaller wins; tie iff `|R_X−R_Y| ≤ ε_rel·max(R_X,R_Y)`,
   `ε_rel = 100·eps(real(T))` (fixed type-scaled band, *not* an absolute
   float compare);
2. **secondary — `g`**, larger wins, same ε_rel band;
3. **tertiary — `K`**, larger wins (integer, exact);
4. **final — cell label**, prefer **(A)** over (B) (guarantees a total order).

*Justification.* Residual first: the held-out coefficient most directly
measures which cell generalises one step further — the quantity the stepper
cares about. Conditioning second: among near-equal fits the better-
conditioned `Q` gives more trustworthy poles and is less roundoff-sensitive
at extended precision. Order third: higher order is an asset only when
accuracy/conditioning do not already separate the cells. Cell-label last
ties the final knot to (A), protecting the `d=1` contract (§4) and
meromorphic exactness.

### 4. `d=1` consistency contract

**Contract**: for `d=1`, return *bit-identically* (modulo SVD sign) the
existing `robust_pade(jet,m,m;:svd)` result — the current (A) output, on
which `SP.1.1` and ~49 downstream tolerances depend (`SharedPade.jl:47-54`).

**Recommendation: short-circuit to (A) at `d=1`** (`d==1 && return (A)`),
rather than relying on selection: (i) makes the contract *unconditional*,
immune to any future objective-weight change; (ii) avoids building (B) and
running selection on the hot scalar path. Belt-and-braces: a test also
asserts the selection *would* pick (A) at `d=1` (it should — (A) is the
higher-order `(m,m)`, with the cell-label tie-break favouring (A) on any
residual tie). The short-circuit is load-bearing; the assertion guards
against the two diverging silently.

### 5. Free diagnostic — A/B disagreement flags an ambiguous step

When (A) and (B) disagree markedly the step is ill-conditioned — a free
byproduct of building both. Surface `residual_gap = R_A/R_B` and
`pole_disagreement` = matching distance between the root sets of `Q_A`,`Q_B`
(min-bipartite-matched, normalised by `|root|`).

**Surfacing** (mirroring ADR-0016's eager opt-in, `0016:82-94`): keep the
default return of `shared_denominator_pade` the **2-tuple
`(numerators,denominator)`** — `VectorStepper.jl:242,269` destructures
exactly that and must not break. Add `diagnostics::Bool = false`; when
`true`, return a **3-tuple** `(numerators, denominator, ::SharedPadeDispatch)`
whose struct carries `chosen_cell::Symbol` (`:diagonal`/`:square`),
`held_out_residual`, `isolation_gap`, `residual_gap`, `pole_disagreement`,
`reduced_degree`. On large `residual_gap`/`pole_disagreement`, also emit a
rate-limited `@warn` (`maxlog`). Default `false` preserves every current test
invariant byte-for-byte.

### 6. Performance — conditional build

**Recommendation: build (B) only when (A) is rank-deficient or poorly-
conditioned**, detected from (A)'s *already-computed* SVD: if (A)'s gap
`g_A = σ_m/σ_{m+1}` is below a threshold (e.g. (A) is full-column-rank with no
σ≤τ — the entire/regular least-squares regime, `SharedPade.jl:218-227`) **or**
(A) required ADR-0027 reduction, build (B) and select; else return (A)
directly (a genuine isolated null space ⇒ exact shared poles ⇒ (B) cannot
beat it). The meromorphic-clean-pole hot path stays at the current single-SVD
cost; the second SVD is paid only on the steps this ADR targets. The `d=1`
short-circuit (§4) skips all of it.

## Decisions requiring maintainer sign-off

Each fork the math cannot decide, with my recommendation. **Build is gated on
sign-off (Law 1).**

1. **Objective ordering/weights.** *Rec*: lexicographic `R→g→K→cell(A)` with
   ε_rel = `100·eps(real(T))` bands (§3). *Awaiting sign-off.*
2. **Neighbouring-degree breadth.** *Rec*: **exactly the two ADR-0027-reduced
   cells** for v1 (each cell's reduction already picks its honest degree; ±1
   multiplies cost for marginal gain). Add neighbours only if a probe shows a
   step where the reduced degree is the wrong cell. *Awaiting sign-off.*
3. **Always-both vs conditional.** *Rec*: **conditional** (§6). *Awaiting.*
4. **Diagnostic surfacing.** *Rec*: opt-in `diagnostics::Bool=false` →
   3-tuple + struct, plus rate-limited `@warn`; default unchanged (§5).
   *Awaiting.*
5. **Performance budget.** Second SVD (≈`m×(m+1)`) per *ill-conditioned* step
   acceptable? *Rec*: **yes** — the gate confines it to the targeted steps;
   one extra `O(m³)` SVD on a small matrix. *Awaiting.*
6. **Held-out combine rule.** *Rec*: **max** across `d` for the decisive
   residual (worst component governs `Q`); RMS in the diagnostic (§2).
   *Awaiting.*
7. **`d=1` mechanism.** *Rec*: **explicit short-circuit to (A)** + selection-
   would-pick-(A) test (§4). *Awaiting.*

## Consequences

- **Accuracy recovery (the goal).** On regular/entire and high-`d` stretches
  the square (B) cell recovers toward **~1e-11** (the `flnr` target), undoing
  the worklog-067 cost — most for high-`d` (NY A₄ ~6e-6 → ~1e-11), since (B)'s
  `n=m/d` rows are *less* over-determined exactly where (A)'s `(d·m)×(m+1)`
  stack hurts most.
- **Meromorphic exactness preserved.** On a genuine pole (A) has a sharp null
  space, large `g_A`, small `R_A` → wins; the conditional gate (§6) returns
  (A) without building (B). ℘/CM/KKG stay exact.
- **`d=1` bit-identity preserved** (short-circuit, §4).
- **Determinism preserved** (C7): same jets → same cell, by construction.
- **Cost**: one extra small SVD per ill-conditioned step; zero on `d=1` and
  clean-pole steps.

### Re-examining `padetaylor-unk`

`unk`'s NOTES flag that its ~4e-12 entire-system residual was "LIKELY a
symptom of the off-by-one (d3a) and/or the full-column-rank rank-break
(3p9c)," to re-examine after d3a lands. d3a/3p9c **have** landed (worklog
067), and the honest `(m,m)` residual is now *larger* (~5e-9 harmonic) — so
the ~4e-12 was the *incidental* accuracy of the buggy `+1`/`(m−1,m)` cell:
**`unk` was observing cell (B)'s accuracy through the off-by-one.** The honest
framing is that entire-system accuracy is *cell-dependent*, and this ADR's
dispatch is exactly the mechanism that recovers (B)'s tighter entire-jet
accuracy on purpose and validated. **Recommendation**: **re-scope `unk`** to
"entire-system accuracy recovered by ADR-0028's (B)-cell selection; close when
the dispatch lands and entire-system tolerances tighten back to (B)-level,"
and make it a **child of `flnr`**. Do not close yet — the recovery is this
ADR's deliverable.

### Alternatives considered, rejected

- **Always (A) [status quo].** The problem itself — pays the worklog-067 cost
  on every regular/high-`d` step.
- **Always (B).** Lower-order; loses (A)'s *sharp* exact-pole recovery on
  genuine meromorphic poles (would regress ℘/CM/KKG exactness).
- **Order-only selection.** Always picks (A) ⇒ identical to status quo.
- **Conditioning-only selection.** Would pick (B) even on clean poles where
  (A)'s higher order is genuinely better; the held-out residual is the arbiter.
- **A third per-step oracle (e.g. AAA).** Expensive, adds failure modes; the
  held-out Taylor coefficient is a free, structurally-aligned validator.

## Test plan (sketch; mutation-prove each)

- **Meromorphic exact recovery stays exact** (℘, CM, KKG): selection picks (A)
  on pole-bearing steps; eval error at current oracle values. *Mut*: force
  always-(B) ⇒ ℘ exact-recovery RED.
- **Entire/high-`d` recovered toward ~1e-11** (harmonic, exp, NY A₄):
  dispatched step beats the worklog-067 `(m,m)` numbers. *Mut*: disable the
  (B) build (gate always returns (A)) ⇒ those tolerances regress to the table
  values, RED.
- **`d=1` bit-identity** (SP.1.1) AND the selection-would-pick-(A) assertion.
  *Mut*: remove the short-circuit + bias tie-break to (B) ⇒ SP.1.1 RED.
- **Determinism**: dispatch twice on the same jets, and on a jet perturbed by
  `eps` → same `chosen_cell`. *Mut*: replace an ε-band with bare float `<` ⇒
  the eps-perturbed test flips the cell, RED.
- **Held-out metric correctness**: on a known shared-`Q` rational, the true
  cell's held-out residual is ~0 and beats a wrong-degree candidate. *Mut*:
  probe the *matched* coefficient instead of the next-unused ⇒ both ~0, metric
  loses discrimination, the high-`d` recovery test RED.
- **Diagnostic surfacing**: `diagnostics=true` → populated 3-tuple;
  `=false` → byte-identical 2-tuple. *Mut*: default `true` ⇒ `VectorStepper`
  destructure errors, downstream RED.

## References

- `docs/worklog/066-sharedpade-window-rootcause-resolved.md:94-110` (two
  forks; order-vs-conditioning; the dispatch motivation).
- `docs/worklog/067-sharedpade-c1-c2-minimal-core.md:58-89` (the shipped
  `(m,m)` core + accuracy-cost table this ADR must recover).
- `docs/adr/0027-sharedpade-graceful-reduction.md:39-62,89-96` (the
  denominator-degree reduction this composes on top of).
- `docs/adr/0019-shared-denominator-pade.md:69-77` (the construction both
  cells refine; SP.1.1 `d=1` oracle).
- `src/SharedPade.jl:124-151` ((m,m) blocks = cell (A)), `:204-313` (reduction
  loop, SVD `S,Vt`, QR-reweighting, `z^λ` cancellation), `:173` (the 2-tuple
  return contract).
- `src/RobustPade.jl:368` ((m,m) `:svd` `d=1` target), `:406-425` (scalar
  `ρ==n` reduction), `:459-489` (`_trim_and_normalise` / `z^λ`).
- `references/tex/hermite_pade/ManoTsuda2017_..._MathZ285/hp_arXiv_final.tex`:
  eq:rTm `:1060-1080` (Toeplitz block, top-left `a^i_m`), eq:lspa `:1405-1428`
  (square `m×(m+1)` stack, `n` rows/block), eq:spa `:414-417` (matching
  `O(w^{nL})`), deg cap `:411,:399` (numerator `≤ m−1`), `:1381` (`m=n(L−1)`).
- `docs/adr/0016-diagnostics-extension.md:82-94` (the `diagnose::Bool=false`
  eager-opt-in pattern §5 mirrors).
- `bd show padetaylor-flnr`, `bd show padetaylor-unk`.
- CLAUDE.md Laws 1 & 2, Rules 1, 4, 5, 9, 10.

---

## Amendment 1 (2026-06-04) — empirical audition findings

A prototype + audition of this design was run **before build** (bead
`padetaylor-0ql3`; probe `external/probes/adr0028-dual-construction-audition/`,
synthesis in its `FINDINGS.md`). The probe unifies cell (A) and three cell-(B)
realisations into one parameterised block-Toeplitz solver — cell (A) asserted
**bit-identical** to the shipped `shared_denominator_pade` (Δ=0, `smoke.jl`) —
and measures *true* per-step error vs oracle (closed-form ℘/Calogero–Moser,
conservation for Noumi–Yamada A₄) across Float64, BigFloat-256, and
pole-crossings. It **confirms the motivation but corrects three load-bearing
claims of §1–§6.** Where this amendment and the original text conflict, **this
amendment governs.**

### A1.1 — ✅ The accuracy recovery is real (the goal stands)

The square Mano–Tsuda cell recovers the worklog-067 cost by **6–8 orders** at
Float64: harmonic `5.22e-9 → 1.4e-15`, ℘ (regular step) `3.93e-4 → 8.2e-13`,
Calogero–Moser `9.66e-11 → 4.4e-16`. (`err_A = 5.22e-9` reproduces worklog
`067:65` `~5e-9` exactly — harness faithful.)

### A1.2 — 🔴 Cell (B) must be the higher-degree **wide-square** system, not §1(B)'s `⌈m/d⌉`-keep-degree

§1(B) and sign-off #2 specify `n = ⌈m/d⌉` rows **keeping denominator degree
`m`**. The audition shows this is the *wrong* realisation: when `d∤m` it makes
the stack **over-determined** (`d·⌈m/d⌉ ≥ m+1`), re-introducing the very
least-squares pathology cell (B) was meant to cure — `×4e9 / ×7e8 / ×2e6`
*worse* than the square variants on exp / CM / NY. (This is exactly the
"source-silent" `d∤m` gap the Mano–Tsuda scout flagged against §1(B):82-85 — it
is decisive, not minor.) **Corrected cell (B): the genuinely wide-square system
at degree `m_eff = ⌈m/d⌉·d ≥ m`** — `n = ⌈m/d⌉` rows, `d·n = m_eff`, an
`m_eff×(m_eff+1)` matrix with a guaranteed 1-D null, `+1` window, degree-`(m_eff−1)`
numerators (probe `cell_B_grow`). This both recovers the entire-jet accuracy and
beats the lower-degree `⌊m/d⌋` square variant on meromorphic steps (CM `4.4e-16`
vs `3.4e-13`).

### A1.3 — ✅🔴 The dispatch is justified, but the win-boundary is {genuine pole × precision}, not entire-vs-meromorphic

The original "(A) wins on genuine meromorphic poles, (B) on entire/regular" framing
(§ Consequences) is too coarse. Stress-testing at **BigFloat-256** found cell (A)
wins on Calogero–Moser by **×900** (`2.3e-22` vs the corrected B `2.1e-19`) — A's
clean-isolated-null degree-`m` advantage, masked by roundoff at Float64 (where B
wins ×280), emerges only at extended precision **and** only where a genuine shared
pole is present (NY's *regular-stretch* step does **not** flip at BF-256: A stays
`4.7e-14`, B `5.6e-17`). A single transcendental **pole-crossing** does *not*
favour A (F64 tan past π/2: off-diagonal ties/beats diagonal); multi-pole single
steps break *both* cells. So the A-win regime is narrow — *genuine shared pole +
extended precision* — but it **includes the FW 2011 Table 5.1 BF-256 long-range
showcase**, so "always B" is unsafe and the per-step dispatch is retained.

### A1.4 — 🔴 Replace the (R,g,K) Pareto (§2–§3) with a held-out-POINT discriminator; drop the conditioning gap

The §2–§3 selection mechanism is **not trustworthy as specified** (mis-ranks
**6/8** test cases):
- the **held-out residual `R`** (§2) is fooled — an over-determined LS fit nails a
  single held-out coefficient while its function is inaccurate;
- the **`ε_rel = 100·eps` band** (§3) is relative to `R` (~1e-30), far below `R`'s
  ~1e-16 noise floor, so near-exact cells select on noise;
- the **conditioning gap `g = σ_m/σ_{m+1}`** (§3 obj. 2, §6 gate) **does not
  separate the regimes** (the "A-worse" and "A-ok" `g_A` ranges overlap entirely)
  **and is not even consistently defined across the two cell shapes** — cell (A)
  is tall (its null is `σ_{m+1}`), the corrected cell (B) is wide-square (its null
  is structural, absent from the returned σ-vector) — so `g_A` and `g_B` measure
  different objects.

**Corrected mechanism — RESOLVED 2026-06-04 (Amendment 2, below):** select by the
**relative ODE defect** `‖ỹ'(t) − h·f(z₀+h·t, y(t))‖ / ‖h·f‖` (a reference-free,
precision- and pole-agnostic accuracy proxy; Enright/Higham defect control,
Corless–Kaya 2025). It is a **single uniform selector — 18/18** across pole-free,
pole-crossing, BigFloat, and BigFloat-crossing regimes (probe `pole_crossing.jl`,
`relative_defect`), needing **no `g`, no radius check, and no separate
pole-crossing fallback**. The held-out-point discriminator (8/8 on pole-free
steps) remains a valid alternative/cross-check but is only accidentally correct
across poles; the defect is principled everywhere. **Drop objective `g`
entirely.** See Amendment 2 for the full pole-crossing audition.

### A1.5 — unchanged

§4 (**d=1 short-circuit to (A)**) is *confirmed needed* — at d=1 cell (B) =
`robust_pade(m−1,m)` ≠ cell (A) = `robust_pade(m,m)` (`smoke.jl`: μ_A=6, μ_B=5).
§5 (opt-in `diagnostics` 3-tuple) stands. §6 perf: since the `g` gate is dropped,
**build both always** (one extra small `O(m_eff³)` SVD per d≥2 step) — the
unreliable conditional gate is not worth its risk.

### A1.6 — updated sign-off status

1. Objective ordering — **superseded**: use held-out-POINT, drop `g`/ε-floor.
2. Neighbour breadth — **superseded**: cell (B) = wide-square `⌈m/d⌉·d` (A1.2).
3. Always-both vs conditional — **decided: always-both** (the `g` gate fails, A1.4).
4. Diagnostic surfacing — **stands** (opt-in 3-tuple; add `chosen_cell`, point-errors).
5. Performance budget — **accepted**, always-both (A1.5).
6. Held-out combine rule — **moot/replaced** by the point metric (max over components).
7. d=1 mechanism — **stands** (short-circuit, A1.5).

### A1.7 — untested before build (carry-overs)

The audition measured **isolated single steps**. Before/with build, verify on:
the multi-step **path-network walk** (accumulating error, adaptive `h`); the
**~49-tolerance contract** (the corrected cell (B) lets the worklog-067 tolerances
tighten back — confirm none asserts the loose value as a *lower* bound); cell-(B)
**degenerate guards** (`d > m`, identically-zero component jet `padetaylor-0o9`;
the prototype skips the ADR-0027 reduction loop on the square variants); and a
**close-call probe** for the defect selector (all audition cases were decisive;
A1.4's in-step-pole fallback is now *resolved* — Amendment 2).

### A1 references

- `external/probes/adr0028-dual-construction-audition/FINDINGS.md` (synthesis),
  `cells.jl` (the unified cell builder + `held_out_point`), `audition.jl`,
  `addendum.jl`/`addendum2.jl` (BigFloat + multi-pole), `discriminator.jl`
  (the 8/8 selector validation), `smoke.jl` (cell-A bit-identity + d=1 oracles).
- `bd show padetaylor-0ql3` (audition bead, with the numeric summary).

---

## Amendment 2 (2026-06-04) — the selector is resolved (relative ODE defect; uniform)

A2 closes the one design item A1.4 left open (the in-step-pole fallback for the
held-out-point selector). A literature survey (defect/residual control — Enright;
Higham 1991; Corless–Kaya 2025 arXiv:2510.20117; AAA — Nakatsukasa–Sète–Trefethen
2018; singularity-from-jet — Domb–Sykes, Mercer–Roberts, Hunter–Guerrieri) plus a
four-discriminator audition (probe `pole_crossing.jl`, bead the A2 research bead)
settled it.

### A2.1 — the relative ODE defect is a single UNIFORM selector

For the candidate `ỹ_i(t) = P_i(t)/Q(t)` over a step (rescaled `t∈[0,1]`,
`z = z₀+h·t`), the **relative defect**
`score = max_t ‖ ỹ'(t) − h·f(z₀+h·t, ỹ(t)) ‖_∞ / (‖h·f‖_∞ + atol)` (sampled at
interior nodes off the roots of `Q`; `ỹ' = (P'Q − PQ')/Q²` exact) is a
reference-free accuracy proxy that needs **no surrogate truth**, so it stays valid
*across* a pole (the `‖h·f‖` normalisation cancels the blow-up). Scored against
the oracle across **four regimes** — pole-free, pole-crossing (tan-companion
`u''=2u+2u³` → `(tan z, sec²z)`, double pole at π/2, exact oracle), BigFloat-256,
and BigFloat pole-crossing — it is **18/18**, correctly firing **A** on the lone
genuine A-win (BF-CM regular step) and **B** everywhere else.

### A2.2 — cell A never wins a pole-CROSSING step

Even at BigFloat-256 and degree `m=28`, the over-determined diagonal A *degrades*
across a pole (`eA` grows `1e2→1e3` with `m`) while the square cell B stays clean
(to `4e-32`). A's only win is a *regular* high-precision step near (not across) a
pole. So the crossing regime is trivially B — which is why even the
theoretically-invalid held-out-point scores 18/18 there (every crossing is a
landslide). The relative defect is preferred because it is *principled* across
poles, not merely accidentally correct.

### A2.3 — revised recommendation + remaining work

- **Selector: the relative ODE defect (A2.1).** Drop the held-out-point + radius
  fallback of A1.4 (the defect needs neither); keep the held-out-point only as an
  optional independent cross-check on pole-free steps. Drop `g` (A1.4).
- **Resolved sign-off #1**: relative ODE defect; no `g`, no ε-floor, no fallback.
- **Remaining before build** (A1.7 minus the now-resolved fallback): the
  ~49-tolerance contract; cell-(B) degenerate guards (`d>m`, zero-component jet
  `padetaylor-0o9`, ADR-0027 reduction loop on the square variant); a multi-step
  path-network walk; and a **close-call probe** (all audition cases were decisive
  — a within-one-order A≈B step has not been stress-tested, though a misrank there
  costs ≈ nothing by definition).

### A2 references

- `external/probes/adr0028-dual-construction-audition/pole_crossing.jl` (the
  4-discriminator × 4-regime audition, 18/18), `cells.jl` (`relative_defect`,
  `q_defect`, `held_out_point`, `nearest_pole_modulus`), `FINDINGS.md`
  (§"Pole-crossing resolution").
- Defect/residual control: Enright (CS Toronto E11); Higham 1991 SIAM
  J.Sci.Comp. 12:991; Corless–Kaya 2025 arXiv:2510.20117. Singularity-from-jet:
  Domb–Sykes; Mercer–Roberts (arXiv:2406.04228 App. B); Hunter–Guerrieri 1980
  SIAM J.Appl.Math. 39(2):248. AAA: Nakatsukasa–Sète–Trefethen 2018
  arXiv:1612.00337. Repo precedent: `figures/_kkg_pi2_helpers.jl` `kkg_ode_residual`,
  `test/painleve_hierarchy_test.jl:154-186`.
