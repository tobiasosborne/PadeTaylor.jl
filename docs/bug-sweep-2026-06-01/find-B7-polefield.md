# Bug sweep B7 — PoleField + VectorPoleField

## Area

Pole extraction: rooting the local Padé denominator `Q`/`b`, mapping
roots from the node-local rescaled variable `t = (z - z_node)/h` back to
global `z`, far-root filtering, residue/Froissart filtering, and
cross-node clustering. Files audited:

- `src/PoleField.jl` (scalar, `_extract_poles_core` + two `extract_poles`
  adapters)
- `src/VectorPoleField.jl` (vector shared-`Q`, `extract_poles_shared_q`)

## References checked

- `external/chebfun/padeapprox.m:146,150,153-154` — MATLAB poly/roots
  convention (`polyval(a(end:-1:1),z)`, `roots(b(end:-1:1))`), residue
  estimate. This is the canonical reversed-order roots call the brief
  flagged.
- `src/RobustPade.jl:110-137` — `PadeApproximant` field convention:
  `a`, `b` stored **low-to-high**, `b[1] = 1`. This is the convention the
  pole extractor consumes.
- `src/PadeStepper.jl:36-48,160-162,300-327,342-350,367-409` — the Padé
  lives in the rescaled variable `t = h'/h` with `c̃_k = h^k c_k`;
  `P_u(t) ≈ u(z_node + h·t)`; `_evaluate_pade` evaluates `b` low-to-high
  via Horner. Establishes that a denominator root `t*` is a pole at
  `z = z_node + h·t*`.
- `src/PathNetwork.jl:29-32,92-93,206-210,500-619,680-691` — Stage-2
  evaluates `visited_pade[k]` at `t = (z_f - z_v)/h` (inverse of the pole
  map); `:adaptive_ffw` stores a per-node varying `visited_h[k] = h_step`.
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:147,281-318`
  — pole fields are composites of independent runs (cross-node-support
  rationale); Weierstrass-℘ test, second-order poles at `x = 1 + 2ωk`.
- `docs/adr/0019-shared-denominator-pade.md:8-39` — shared `Q`; every
  component's poles are the roots of one polynomial.
- `docs/adr/0026-vector-resilient-walk-dense-targets.md:410-455`
  (Amendments 6 & 7) — the **documented historical bug**: the old
  far-root filter `|t*| ≤ radius_t` is a "scale-fixing heresy" that
  **intermittently empties pole fields** when adaptive `h` shrinks; fixed
  in the vector path (S7) by `h_node·|t*| ≤ radius_t·h_max`, and the
  representative-selection facet by ordering on z-distance not `|t*|`.
- `test/ffw_fig_6_test.jl:85,113-120,189-197` — scalar `extract_poles`
  called on an `:adaptive_ffw` + `node_separation=R` walk where
  `R(ζ)=max(0.05,(4-Re ζ)/10)` varies `h` ~3.3× across the domain.

## Findings

### [HIGH] Scalar `PoleField` still uses the scale-fixing far-root filter `|t*| ≤ radius_t` that ADR-0026 root-caused as the cause of intermittently empty/sparse pole fields

**Location:** `src/PoleField.jl:144`
(`abs(t) ≤ radius || continue`), with the ordering key `c[2] = abs(t)`
at `:147,156`.

**Ground truth (cited):**
`docs/adr/0026-vector-resilient-walk-dense-targets.md:418-426` —
"`extract_poles_shared_q` filters shared-`Q` roots by `radius_t` in the
**rescaled variable `t = Δz/h`** — i.e. it keeps poles within a fixed
number of *steps*. A smaller `SAFETY` ⇒ smaller `h` ⇒ the same physical
pole maps to a larger `|t*|` ⇒ it falls outside the `radius_t` window ⇒
the `min_support ≥ 2` cross-node filter empties the field. **`radius_t`
is a scale-fixing heresy**." The fix (`:436-441,449-453`, Amendment 7
"S7 shipped") is `h_node·|t*| ≤ radius_t·h_max` and ordering the
greedy-cluster representative by z-plane distance, not `|t*|` — both
"mutation-proven." That fix landed only in `src/VectorPoleField.jl`
(`:287-297` far-root, `:320` sort by `z_dist`); the scalar `PoleField`
was not touched.

**Code behavior:** scalar `_extract_poles_core` keeps a denominator root
iff `abs(t) ≤ radius_t` (`:144`), where `t` is the root in the *per-node*
rescaled variable `t = (z - z_node)/h_node`. It also sorts candidates by
`abs(t)` (`:147` stores `RT(abs(t))`, `:156` sorts on it) and takes the
first-in-cluster as the cluster representative (`:159-167`). Both
quantities are per-node-step-relative, exactly the pattern Amendment 7
fixed in the vector module.

**Mechanism (intermittent discontinuity):** the scalar `extract_poles` is
exercised on adaptive walks with strongly varying per-node `h`
(`test/ffw_fig_6_test.jl:113-120` uses `:adaptive_ffw` +
`node_separation=R`, `R` ranging 0.5→0.15 across `Re ζ ∈ [-1,2.5]`, plus
further adaptive shrink toward `h_min`). A physical pole at fixed
z-distance `D` from a node is kept only if `D/h_node ≤ 5`. From a
coarse-`h` node (`h≈0.5`) it is kept out to `D≈2.5`; from a fine-`h` node
(`h≈0.15`) only out to `D≈0.75`. So whether a given node contributes a
root toward a pole's `min_support` vote depends on that node's local
adaptive `h`, which depends on the (RNG-shuffled, FW 2011 line 156) path
the walker happened to take to that region. A pole whose nearby nodes
happen to have shrunk `h` loses cross-node support and silently drops out
of the returned set; a tiny change in walk order or `adaptive_tol` flips
it back in. That is precisely the "intermittent discontinuity in the
rendered pole field" symptom, and ADR-0026 records it as *measured*
("lowering SAFETY 0.25→0.10 emptied four `src/` test fields"). The
ordering-by-`|t*|` facet adds a second intermittency: under varying `h`
the smallest-`|t*|` candidate is not the physically-closest node
(Amendment 7 measured a ℘-pole placed at 4e-5 from the min-`|t*|` node vs
6e-7 from the z-closest node), so the *reported location* of a pole jumps
between near-coincident cluster members as `h` varies.

**Intermittent?** Yes — data/order-dependent: the failure mode appears
only when per-node `h` varies (adaptive / non-uniform `node_separation`)
and depends on the RNG-shuffled walk geometry. For a uniform-`h` walk all
`h_node` are equal and the bug is dormant (which is why the default
`:fixed`-mode tests stay GREEN and the scalar bug was never caught).

**Confidence:** 0.78. The mismatch against the cited ADR is exact and the
ADR explicitly calls this the root cause of intermittent empty pole
fields; the residual uncertainty is whether the maintainer regards the
scalar path as deliberately frozen at v0.1 behaviour (the vector module's
docstring repeatedly frames the scalar `PoleField` as the "v0.1" object).
But ADR-0026's "scale heresies *must* be scale-covariant, so they are
fixed on principle" (`:289`) reads as a project-wide invariant, not a
vector-only one, and the scalar path is demonstrably run on adaptive
walks (`test/ffw_fig_6_test.jl`).

### [LOW] Scalar `PoleField` absolute `cluster_atol` default (0.1) is the "absolute cluster tolerance" the ADR calls a scale heresy that mis-merges/splits poles in a varying-scale field

**Location:** `src/PoleField.jl:160` (`abs(p - r) ≤ catol`), default
`cluster_atol = 1.0e-1` at `:176,207,246`.

**Ground truth (cited):**
`src/VectorPoleField.jl:125-155` (citing `ADR-0026 Amendment 3 §S4` and
project memory `scale-covariance-core-principle`): "An absolute cluster
tolerance is a *scale heresy* … in a varying-scale pole field it wrongly
*merges* distinct poles where the field is dense and *splits* one pole
where it is sparse." The vector module replaced it with
`CLUSTER_FRAC·min(h_i,h_j)`.

**Code behavior:** the scalar clusterer merges candidates within a fixed
absolute z-distance `catol` regardless of local pole spacing.

**Mechanism:** in a pole field whose spacing varies across the domain
(any transformed-plane figure, e.g. the `ζ = log z` PV field of
`ffw_fig_6_test.jl` where pole spacing compresses with `Re ζ`), a fixed
`catol` can merge two genuinely distinct poles in the dense region (one
pole vanishes from the output) or split one pole's near-coincident
cluster members in the sparse region (a phantom doubles). Which happens
depends on where the adaptive walk placed its nodes — intermittent.
Lower confidence than the far-root filter because the scalar default is
also user-overridable (the fig-6 test pins `0.15`) and the effect is
secondary to the far-root drop above.

**Intermittent?** Yes — scale/density dependent.

**Confidence:** 0.4. The "absolute tolerance is a heresy" judgement is the
project's own (cited above), and it is undeniably still present in the
scalar default; but unlike the far-root filter, ADR-0026 did not pin a
*measured* failure to the scalar `cluster_atol` specifically.

## Areas verified correct

- **Denominator coefficient ORDER — no reversed/reciprocal-pole bug.**
  `PoleField.jl:139` builds `D = Polynomial(P.b)` and `VectorPoleField.jl:285`
  builds `Polynomial(Q)`, both feeding `Polynomials.roots`. `P.b`/`Q` are
  stored **low-to-high** (`src/RobustPade.jl:115-118`; corroborated by the
  Horner sweep `den = Σ b[k] z^(k-1)` in `src/PadeStepper.jl:373`), and
  `Polynomials.Polynomial` takes coefficients low-to-high. So the
  polynomial rooted is identical to the one evaluated by the stepper, and
  its roots are the genuine denominator zeros. The MATLAB
  `roots(b(end:-1:1))` reversal (`padeapprox.m:150`) exists only because
  MATLAB's `roots` expects high-to-low; the Julia code correctly does
  **not** reverse, because `Polynomials` uses the opposite convention.
  Adding a reversal here would have produced reciprocal poles — it is
  correctly absent.

- **Rescaling root → global z: shift AND scale both present and correct.**
  `PoleField.jl:147` maps `z = z_ctr + h*t` and `VectorPoleField.jl:295`
  maps `z = z_node + h_node*t`. This is the exact inverse of the Padé's
  variable `t = (z - z_node)/h` (`src/PadeStepper.jl:160-162,280`;
  Stage-2 forward use `t = (z_f - z_v)/h` at `src/PathNetwork.jl:31,93,690`).
  The `h` factor is present (no missing scale) and the `+ z_ctr` shift is
  present. The vector module's `z_dist = abs(z_pole - z_node)` (`:296`) is
  algebraically `|h_node·t*|`, as its comment claims — verified.

- **`h_max` recovery for the vector far-root radius.**
  `VectorPoleField.jl:259-261` computes `h_max = maximum(abs, visited_h)`
  and `radius_z = radius_t·h_max`; ADR-0026 Amendment 7 (`:449-453`) and
  the module docstring (`:62-82`) confirm this is the intended
  scale-stable radius, and `VectorPathNetwork.jl:708` (`visited_h = C[C(h_max)]`)
  confirms the root node carries exactly `h_max`. Empty-solution guard
  (`zero(T)` placeholder, loop empty) is harmless.

- **Residue / Froissart filter (scalar).** `PoleField.jl:145` computes
  `res = N(t)/D'(t)` with `D' = derivative(D)` (Polynomials analytic
  derivative) and drops `|res| < min_residue` or non-finite. This matches
  the standard simple-pole residue `N(t*)/D'(t*)` and is consistent with
  `padeapprox.m:153-154`'s finite-difference residue estimate (different
  formula, same quantity). No sign or conjugate issue: all arithmetic is
  ordinary complex `Polynomial` evaluation, no adjoint/transpose anywhere
  in either module.

- **No conjugate-transpose / adjoint-vs-transpose hazard.** Neither module
  performs any matrix adjoint, `'`, `transpose`, or `conj`; the only
  complex operations are polynomial root-finding, `abs`, and affine maps.
  The `padeapprox.m:113` historical `.'`-vs-`'` bug class cannot arise
  here (that lives in the SVD/QR step, not in pole extraction).

- **No in-place buffer aliasing.** Both modules build fresh
  `candidates`, `reps`, `support`, `rep_h` vectors per call; `roots`
  returns a fresh array each node; `Set{Int}` mutation is per-cluster and
  isolated. No shared buffer is reused across nodes, so there is no
  aliasing-induced intermittency *within* the extractor (the intermittency
  in Finding B7-1 comes from upstream adaptive-`h` data, not aliasing).

- **Clustering greedy-order determinism (given fixed input).**
  `sort!(candidates; by = c->c[2])` plus `findfirst` is deterministic for
  a fixed candidate list; there is no RNG inside either extractor (the
  RNG lives in the upstream path-network walk). The intermittency is
  data-driven via `visited_h`, not a nondeterministic extractor.

- **Far-root `≤` boundary.** Both modules use `≤` (`PoleField.jl:144`,
  `VectorPoleField.jl:297`) and the support threshold uses `≥`
  (`:168`, `:339`). These boundary choices are internally consistent and
  match the docstrings; no off-by-one or strict-vs-nonstrict flip found.
