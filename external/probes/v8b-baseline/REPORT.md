# V8b Baseline Probe — Phase-A4 (bead padetaylor-0ln.37.4)

Throwaway exploration spike. BEFORE-picture for ADR-0025 re-resolution.

## 0. Reproduce the V8b wedge walk

Compute time: 10.0 s
march_ok    : true
message     : march OK — 389 visited nodes, 21 poles extracted
visited nodes : 389
poles extracted: 21

## Task 2 — Loop-closure ΔP_rel certificate

### 2a. Does `quality_diagnose` accept a `VectorPathNetworkSolution`?

`quality_diagnose(walk)` THREW — `MethodError`.
First line: `MethodError: no method matching quality_diagnose(::VectorPathNetworkSolution{Float64})`

Verdict: `quality_diagnose` does NOT accept the vector walk.

### 2b. Hand-computed equivalent loop-closure certificate

Delaunay edges     : 1155
tree edges         : 388
non-tree edges     : 767

non-tree edges evaluated (finite both ends): 767

### ΔP_rel distribution (vector-norm loop closure)

  median : 9.915e-04
  p90    : 9.955e-01
  p99    : 1.000e+00
  max    : 1.000e+00

### Category breakdown

  well_closed   : 40 (5.2%)  (ΔP_rel ≤ 1.0e-10)
  noisy         : 190 (24.8%)
  extrap_driven : 503 (65.6%)  (|t| > 1 at midpoint)
  depth_driven  : 34 (4.4%)  (in-disc loop-closure failure)

bad_centroid (mean midpoint of ΔP_rel > 1.0e-6 edges): 4.824916697515218 - 0.24191710076754872im
n bad edges : 537

### Worst 10 edges

   1. ΔP_rel=1.000e+00  cat=extrap_driven  extrap_max=51.60  tree_d=203  mid=2.087 + 1.969im
   2. ΔP_rel=1.000e+00  cat=extrap_driven  extrap_max=52.09  tree_d=204  mid=2.038 + 1.96im
   3. ΔP_rel=1.000e+00  cat=extrap_driven  extrap_max=52.58  tree_d=205  mid=1.989 + 1.951im
   4. ΔP_rel=1.000e+00  cat=extrap_driven  extrap_max=53.07  tree_d=206  mid=1.94 + 1.941im
   5. ΔP_rel=1.000e+00  cat=extrap_driven  extrap_max=53.57  tree_d=207  mid=1.985 + 1.963im
   6. ΔP_rel=1.000e+00  cat=extrap_driven  extrap_max=42.16  tree_d=108  mid=0.935 - 1.512im
   7. ΔP_rel=1.000e+00  cat=extrap_driven  extrap_max=3.75  tree_d=14  mid=5.559 + 2.909im
   8. ΔP_rel=1.000e+00  cat=extrap_driven  extrap_max=31.68  tree_d=134  mid=0.511 + 1.051im
   9. ΔP_rel=1.000e+00  cat=extrap_driven  extrap_max=31.19  tree_d=133  mid=0.56 + 1.059im
  10. ΔP_rel=1.000e+00  cat=extrap_driven  extrap_max=3.33  tree_d=13  mid=5.52 + 2.877im

## Task 3 — Conjugate-symmetry pole pairing

V_0 is real on the negative real axis ⇒ V_0(conj x) = conj(V_0(x))
⇒ wedge poles come in conjugate pairs.

Extracted poles (21):
   1. +4.22528 -2.10279im
   2. +4.25617 -3.09205im
   3. +5.05069 -2.60535im
   4. +5.20392 -3.52345im
   5. +5.21284 +2.05293im
   6. +5.33319 +2.93581im
   7. +5.43025 +3.21076im
   8. +5.44104 +2.39613im
   9. +5.80286 +0.71321im
  10. +5.86522 +2.93216im
  11. +5.96884 -2.92236im
  12. +6.03421 +2.40245im
  13. +6.22296 +3.42012im
  14. +6.25558 +1.40439im
  15. +6.62648 +0.66336im
  16. +7.05643 +3.81211im
  17. +7.16812 +0.03784im
  18. +7.43626 +0.82049im
  19. +7.92709 -1.88356im
  20. +8.05026 +1.38250im
  21. +8.22685 +0.56947im

### Pairing result

matched conjugate pairs : 6
unpaired / near-real    : 9

Per-pair residual |p_upper - conj(p_lower)|:
  pair  1: p_up=+5.4302+3.2108im  p_lo=+5.0507-2.6054im  resid=7.146e-01
  pair  2: p_up=+6.2230+3.4201im  p_lo=+5.2039-3.5235im  resid=1.024e+00
  pair  3: p_up=+5.3332+2.9358im  p_lo=+4.2562-3.0920im  resid=1.088e+00
  pair  4: p_up=+7.0564+3.8121im  p_lo=+5.9688-2.9224im  resid=1.405e+00
  pair  5: p_up=+5.8652+2.9322im  p_lo=+4.2253-2.1028im  resid=1.838e+00
  pair  6: p_up=+6.0342+2.4024im  p_lo=+7.9271-1.8836im  resid=1.963e+00

max residual    : 1.963e+00
median residual : 1.247e+00
mean residual   : 1.339e+00

Unpaired / near-real poles (Im x should be ~0):
  pole  4: +5.44104 +2.39613im   (|Im| = 2.396e+00)
  pole 18: +5.21284 +2.05293im   (|Im| = 2.053e+00)
  pole 14: +6.25558 +1.40439im   (|Im| = 1.404e+00)
  pole 12: +8.05026 +1.38250im   (|Im| = 1.382e+00)
  pole  7: +7.43626 +0.82049im   (|Im| = 8.205e-01)
  pole 11: +5.80286 +0.71321im   (|Im| = 7.132e-01)
  pole 10: +6.62648 +0.66336im   (|Im| = 6.634e-01)
  pole 15: +8.22685 +0.56947im   (|Im| = 5.695e-01)
  pole  6: +7.16812 +0.03784im   (|Im| = 3.784e-02)

## Task 4 — VC-4 dominant-balance A spot check

ADR-0025 A3: every P_I^(2) pole has u ~ A(x-x0)^-2 with A ∈ {-1,-3}.
Fit u(x) ≈ A(x-x0)^-2 + B(x-x0)^-1 + C on a 32-point ring.
u(x) is recovered from the nearest visited node's shared-Q approximant.

Spot-checked poles (nearest-origin 6):

  pole  8 @ +4.2253-2.1028im : A=-1.0000-0.0000im  |B|=9.908e-05  min(|A+1|,|A+3|)=1.034e-05
  pole 21 @ +4.2562-3.0920im : A=-1.0000-0.0000im  |B|=2.297e-04  min(|A+1|,|A+3|)=7.866e-06
  pole 18 @ +5.2128+2.0529im : A=-1.0000-0.0000im  |B|=2.951e-04  min(|A+1|,|A+3|)=1.330e-05
  pole  9 @ +5.0507-2.6054im : A=-1.0000-0.0000im  |B|=3.220e-05  min(|A+1|,|A+3|)=1.633e-06
  pole 11 @ +5.8029+0.7132im : A=-1.0000-0.0000im  |B|=1.594e-04  min(|A+1|,|A+3|)=9.348e-06
  pole  4 @ +5.4410+2.3961im : A=-0.0000-0.0000im  |B|=2.970e-09  min(|A+1|,|A+3|)=1.000e+00

NOTE: ring radius r=0.04 < walk step h=0.1, so most ring points
share ONE node's Padé — this is a coarse spot check, not the full
VC-4 acceptance test (which re-expands a dedicated jet per pole).

VC-4 verdict: 5 of the 6 spot-checked poles have A = -1 to <1.1e-5
and |B| < 3e-4 — the A=-1 generic family, exactly ADR-0025 A3. The
6th (pole 4 @ +5.44+2.40im) fits A ≈ 0 / |B| ≈ 3e-9 — NOT a genuine
double pole: a spurious / mis-located cluster. So the A3 result holds
on the genuine poles, and the VC-4 check already flags one of the 21
"poles" as spurious.

---

# Synthesis — the V8b baseline numbers and what they mean

## Finding 1 — `quality_diagnose` does NOT work on the vector walk

`quality_diagnose` is typed `quality_diagnose(sol::PathNetworkSolution)`
(scalar). The V8b walk produces `VectorPathNetworkSolution{Float64}`.
Calling it threw `MethodError: no method matching
quality_diagnose(::VectorPathNetworkSolution{Float64})`.

An adapter is required for bead D-VC7. The data-model gap, field by
field:

| concern              | scalar `PathNetworkSolution`        | vector `VectorPathNetworkSolution`              |
|----------------------|-------------------------------------|-------------------------------------------------|
| per-node approximant | `visited_pade::Vector{PadeApprox.}` | `visited_numerators` + `visited_denominator` (shared-Q) |
| evaluation           | `_evaluate_pade(pade, t)` → scalar  | `Pᵢ(t)/Q(t)` Horner over a *d*-vector           |
| Riemann sheets       | `visited_sheet` field present       | **no `visited_sheet` field** (companion P_I⁽²⁾ is meromorphic / single-sheeted) |
| state at a point     | scalar `u`                          | *d*-vector `y`                                  |

**Adapter D-VC7 needs (probe-and-report only — not implemented here):**

1. A new method `quality_diagnose(sol::VectorPathNetworkSolution; …)`
   in `ext/PadeTaylorDiagnosticsExt.jl` (the extension already `using`s
   `DelaunayTriangulation`; add `VectorPathNetworkSolution` to its
   imports). It is *additive* — no change to the scalar method, so v0.1
   tests stay bit-identical (ADR-0001).
2. Replace `_evaluate_pade(visited_pade[k], t)` with the shared-Q
   evaluation `Pᵢ(t)/Q(t)` — the `_eval_poly` Horner already in
   `VectorPathNetworkStage2.jl`. The result is a *d*-vector.
3. Generalise `ΔP_rel` to a vector norm:
   `ΔP_rel = ‖y_A(M) − y_B(M)‖ / (‖y_A(M)‖ + ‖y_B(M)‖ + ε)` — the
   2-norm over the *d* companion components (used in this probe's §2b).
4. Drop the sheet-0 mask: there is no `visited_sheet` field, and the
   companion system is single-sheeted, so "sheet 0 = all nodes". The
   `_sheet_mask` call must be skipped (not adapted) for the vector type.
5. Everything else of `ext/PadeTaylorDiagnosticsExt.jl` ports verbatim:
   Delaunay triangulation, tree-edge subtraction, non-tree edge set,
   midpoint `t = (M − z_X)/h_X` rescaling, `_build_depths` /
   `_tree_path_distance` LCA, `_categorise` thresholds, the
   `DiagnosticReport` aggregation. `EdgeReport`/`DiagnosticReport`
   structs are reusable as-is (they store `ΔP_rel`/`ΔP_abs` as
   `Float64`, agnostic to scalar-vs-vector origin).

This probe's §2b is a working reference implementation of that adapter
— D-VC7 can lift it directly into the extension.

## Finding 2 — loop-closure baseline (the §2b hand-computed certificate)

767 non-tree Delaunay edges on the 389-node V8b walk:

- `median(ΔP_rel) = 9.9e-4`, `p90 = 9.96e-1`, `p99 = 1.00`, `max = 1.00`.
- categories: well_closed 5.2 %, noisy 24.8 %, **extrap_driven 65.6 %**,
  depth_driven 4.4 %.
- `bad_centroid = 4.82 − 0.24im` (537 edges with `ΔP_rel > 1e-6`).

The distribution is **catastrophically worse** than the FFW Fig 1
scalar baseline (ADR-0016: median 2.6e-12, 6.3 % catastrophic). Here
**70 %** of loop closures are catastrophic (`ΔP_rel > 1e-6`), almost all
`extrap_driven` — the worst edges have `extrap_max` up to 53, i.e.
midpoints sit 50× outside the canonical Padé disc. This is the direct,
quantified signature of the ADR-0025 C2 corner (`h = 0.1` walk, no
true-radius gating): the `h = 0.1` step produces tiny canonical discs,
so every Delaunay edge longer than ~0.2 has its midpoint far outside
both endpoints' discs. The baseline loop-closure quality is essentially
"no consensus" — the figure the re-resolution must beat is, on this
metric, near-zero.

## Finding 3 — conjugate symmetry is BADLY violated (the headline defect)

`V₀` is real on the negative real axis ⇒ poles must come in conjugate
pairs. They do not:

- the target fan IS symmetric (5 angles `±0.5, ±0.25, 0`, 4 radii),
- but the 21 extracted poles split **15 upper-half / 6 lower-half**,
- forced nearest-neighbour conjugate pairing yields only 6 pairs with
  residual `median ≈ 1.25`, `max ≈ 1.96` — i.e. **no pole actually has
  a conjugate partner**; the pairing is meaningless.

The visited-node cloud is only mildly asymmetric (175 upper / 213
lower), so the walk reaches both half-planes — but pole *extraction*
(`extract_poles_shared_q`, `min_support = 2`) keeps 2.5× as many upper
poles. The lower half-plane is under-resolved: nodes there do not
accumulate the cross-node support to clear the `min_support` filter.

This is exactly the defect VC-5 (conjugate-symmetry pole pairing) is
designed to catch, and the V8b figure fails it outright. **Baseline
conjugate-symmetry residual: median ≈ 1.25 (target for re-resolution:
≲ pole-spacing, ideally ≪ 0.1).**

## Finding 4 — VC-4 already flags a spurious pole

5 of 6 spot-checked poles: `A = -1` to `<1.1e-5`, `|B| < 3e-4` — the
A=-1 generic family, ADR-0025 A3 confirmed. The 6th (pole @ +5.44+2.40im)
fits `A ≈ 0` — not a double pole at all. So **at least 1 of the 21 V8b
"poles" is spurious** by the VC-4a criterion. (The full VC-4 test
re-expands a dedicated order-N jet per pole; this spot check shares one
node's Padé per ring and is coarse — but the A=-1 signal is unambiguous
and the spurious pole is clear.)

## Baseline scorecard — what Phase F must beat

| metric                          | V8b baseline           | re-resolution target            |
|---------------------------------|------------------------|---------------------------------|
| pole count                      | 21 (≥1 spurious)       | conjugate-complete, no spurious |
| loop-closure median ΔP_rel      | 9.9e-4                 | ≪ 1e-6 (FFW scalar: 2.6e-12)    |
| loop-closure % catastrophic     | ~70 % (`>1e-6`)        | ≲ 6 % (FFW scalar baseline)     |
| extrap_driven edges             | 65.6 %, extrap_max≤53  | 0 % (B1 true-radius gate)       |
| conjugate-pair residual median  | ~1.25 (no real pairs)  | ≲ 0.1 (VC-5 pass)               |
| upper/lower pole split          | 15 / 6                 | balanced                        |
| VC-4 A=-1 on genuine poles      | confirmed (5/6)        | all poles VC-4a/b pass          |

