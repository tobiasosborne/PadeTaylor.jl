"""
    PadeTaylor.VectorPathNetworkStage2

The **Stage-2 fine-grid fill** of the vector path-network — the dense
evaluation step `VectorPathNetwork.vector_path_network_solve` delegates
to when the caller passes a `fine_grid` (bead `padetaylor-0ln.20`).

## Why this is a separate module

V7 (`src/VectorPathNetwork.jl`, bead `padetaylor-0ln.16`) shipped only
the Stage-1 *tree* construction — the wedge walk that builds a visited
tree of per-node shared-`Q` approximants.  The `A_4⁽¹⁾`/`P_I⁽²⁾`
pole-field *scatter* figures need only that: the pole locations are the
roots of each node's shared denominator `Q`, extracted by
`VectorPoleField`.

The headline `P_I⁽²⁾` tritronquée *surface* figure needs more — a dense
filled heatmap of `y(z)` over the pole-rich wedge.  That is Stage 2.
It is the direct vector lift of v0.1's
`PathNetwork.path_network_solve` Stage-2 step (`src/PathNetwork.jl:29-32`,
"Stage 2 — fine-grid extrapolation"; the fill loop at `:669-693`).

Stage 2 lives in its own module purely for the CLAUDE.md Rule 6 200
effective-LOC cap: `VectorPathNetwork.jl` was already ~145 effective
LOC after V7, and folding the Stage-2 fill + its Horner evaluator + the
expanded docstrings inline pushed it to the cap.  The split is along
the natural seam — Stage 1 builds the tree, Stage 2 reads it.

## No dependency on `VectorPathNetwork`

This module is deliberately free of any `VectorPathNetwork` import, so
it can be `include`d *before* `VectorPathNetwork` and the latter can
`using ..VectorPathNetworkStage2: _stage2_fill` at its own load time.
The one thing Stage 2 needs from Stage 1 — the nearest-visited-node
scan — is passed in as a *function argument* (`nearest_visited`) rather
than imported, so the dependency arrow points one way only.  It *does*
depend on the load-order-earlier sibling modules `VectorStepControl`
(the Jorba–Zou step formula) and `VectorCoefficients` (the jet
fallback) — both `include`d before this file in `src/PadeTaylor.jl`.

## The B1 true-radius validity gate (ADR-0025 Amendment 1)

**This is the substance of bead `padetaylor-0ln.37.5`.**  Through V0.1
the Stage-2 fill gated each visited node's validity at the *step size*
`h_v` (`abs(z_f − z_v) > h_v ⇒ NaN`).  The Phase-A1 spike
(`external/probes/pade-validity-radius/REPORT.md`) measured — against
an algorithm-disjoint order-64 Taylor oracle — that the `h_v` gate is
**not an honesty guarantee**: at `tol = 1e−8` it over-estimates the
honest radius for **36 % of nodes** (worst case 5×).  `h_v` is a step
*length*, not a validity *radius*; the rough coincidence `R_honest ≈ h`
is numerology.

ADR-0025 Lever 1 retires that corner.  The honest radius of a node's
order-`order` shared-`Q` Padé is governed by **finite-order truncation
error** — the decay rate of the rescaled Taylor coefficients — *not* by
the geometry of the pole field (the spike rejected both `h·min|t*|` and
the nearest-uncaptured-pole gate with evidence; they over-estimate 5.7×
and 30×).  The winning gate, **gate (v-b)**, reads each node's own
coefficient decay through the project's own per-node Jorba–Zou
truncation-radius estimator:

```
R_gate(node, tol) = min( s(tol) · h_JZ · h_v ,  0.5 · h_v · min|t*| )

  h_v     = real(visited_h[k])                  — canonical step magnitude
  h_JZ    = vector_step_jorba_zou(rescaled order-`order` jet, tol)
            — the truncation-limited step, returned in the rescaled
              t-variable (so it multiplies h_v to reach a z-plane radius)
  min|t*| = min |root of the node's shared-Q denominator|
            — the pole-adjacency clamp, in the same t-variable
  s(tol)  : 1e-6 → 0.34 ,  1e-8 → 0.36 (default) ,  1e-10 → 0.30
```

The **rescaled jet** is `c̃_k = c_k · h_vᵏ` — the coefficients of `tᵏ`
in the variable `t = (z − z_v)/h_v` the canonical Padé lives in, the
identical rescaling `VectorStepper` applies (FW 2011 §3.2).  Feeding the
rescaled jet to `vector_step_jorba_zou` returns `h_JZ` as a step in
*that* `t`-variable; `h_JZ · h_v` is the z-plane radius.  The
`0.5·h_v·min|t*|` term clamps the disc when a shared-`Q` root sits
inside it — it binds for only ~2 % of (deep-wedge, pole-adjacent) nodes;
for the ~98 % truncation-limited majority the `s·h_JZ·h_v` term binds
and `min|t*|` is inert.

The gate is **honest by construction**: `s < 1` under-covers the
truncation bound, the clamp under-covers the pole bound, so a Stage-2
query in the resulting gap returns `NaN` (fail-soft, Rule 1) — never an
out-of-disc, silently-wrong value.  The per-node spike head-to-head
(REPORT §6) showed gate (v-b) recovers **7–10.5× more honest disc area**
than a global `κ·h` constant at the same ~1–2 % over-estimation rate,
because a global constant cannot track the 3.3× node-to-node spread of
the honest radius.

### Where the order-`order` jet comes from

`vector_step_jorba_zou` needs the node's raw Taylor jet.  The Stage-1
walk now **caches** it: `_canonical_pade` already computed the jet to
build the shared-`Q` Padé, so `VectorPathNetworkSolution.visited_jets`
carries it at zero extra cost (ADR-0025 Amendment 1 §A1, "cache the
order-24 jet so the gate costs nothing extra").  When `visited_jets` is
empty — a solution built by a pre-B1 backward-compat constructor — the
gate **recomputes** the jet from the node `(z_v, y_v)` state via
`vector_taylor_coefficients`; this is the documented fallback, slower
but identical in result.

## The fill (mirrors v0.1 `PathNetwork`'s Stage-2 loop)

For each fine-grid point `z_f`:

  1. find the nearest visited node `idx_v` via the supplied
     `nearest_visited` scan (the same lexicographic-tiebreak
     nearest-node lookup Stage 1 uses);
  2. compute the node's true validity radius `R_gate` (the B1 gate
     above) — once per node, memoised across grid points that share a
     nearest node;
  3. if `!extrapolate` and `abs(z_f − z_v) > R_gate`, the grid point
     lies outside the verified disc of validity of its nearest local
     Padé — the slot gets a `d`-vector of `NaN + NaN·im` (fail-soft per
     ADR-0015 / Rule 1; **not** a throw — the visited tree is honest,
     the grid point is simply uncovered);
  4. otherwise rescale to `t = (z_f − z_v) / h_v` with the *real*
     canonical step `h_v = real(visited_h[idx_v])`, and evaluate
     component `i` as `Pᵢ(t) / Q(t)` — every component over the *same*
     shared `Q(t)`, the v0.2 keystone applied at evaluation time.

The `extrapolate` flag keeps its ADR-0015 escape-hatch semantics
verbatim: default `false` honours the fail-soft Rule-1 contract and the
B1 true-radius gate; `true` skips the gate entirely so a figure script
can request a gap-free panel, at the cost of evaluating Padé
approximants outside their verified disc.  ADR-0025 deletes the
*figure's* use of `extrapolate=true` (bead B3), but the kwarg itself
stays — it is a legitimate diagnostic escape hatch.

## Why nearest-node Voronoi, not a distance-weighted blend (B4 audition)

ADR-0025 bead B4 (`padetaylor-0ln.37.8`) auditioned this nearest-node
assignment against a smooth distance-weighted **blend** of *all* the
node discs that cover a given grid pixel.  The hypothesis (ADR-0025
Phase-B plan; the C4/C5 audition idiom) was that with the B1
true-radius gate adjacent nodes' discs overlap, so a blend would give
C⁰/C¹ continuity across the former Voronoi seams and the inter-disc
spread would double as a free loop-closure quality signal.

The audition (`external/probes/stage2-blend-audition/`) measured the
hypothesis on the *actual* `P_I⁽²⁾` tritronquée wedge underlay — the
2289-node B2 adaptive walk, a 380² grid fine enough (≈0.012 spacing) to
resolve the discs.  The hypothesis **failed on the evidence**:

  - **The discs barely overlap.**  The B1 gate is *truncation-limited*
    (ADR-0025 Amendment 1) — honest discs are small (median `R_gate ≈
    0.038`).  98.2 % of covered pixels lie in exactly **one** node's
    disc; only 1.8 % are multiply-covered; the maximum multiplicity is
    2.  ADR-0025 Amendment 1 already flagged this — "B1 alone makes the
    figure holier, not bigger" — and the underlay confirms it directly.

  - **Most Voronoi seams are not blendable.**  Of the covered
    adjacent-pixel pairs straddling a Voronoi boundary, only **14.9 %**
    lie inside a genuine node-pair overlap; the other 85 % are
    separated by `NaN` gaps — there is no C⁰ discontinuity there for a
    blend to smooth, only honest emptiness.

  - **Where discs do overlap the nodes already agree.**  The
    inter-disc evaluation spread at multiply-covered pixels has median
    `≈ 6·10⁻¹⁰` — machine-precision agreement.  There is no
    discontinuity to remove and the "spread as a quality signal" idea
    yields a uniformly `~10⁻¹⁰` map carrying no information.

  - **No visible seam to begin with.**  The C⁰ jump in `|u|` across a
    Voronoi boundary has median `0.153`, which is `0.73×` the
    *intrinsic* same-cell pixel-to-pixel field gradient (`0.21`): the
    boundary "jump" is *smaller* than the field's own local variation.
    The Voronoi and blend C⁰-jump statistics are identical to three
    significant figures.

A blend would therefore change 1.8 % of a ~5 % sparse underlay by
`~10⁻⁹` — a no-op dressed as continuity.  Per Rule 9 ("no
over-engineering; a well-evidenced retain is a legitimate audition
outcome") the audition **retains the nearest-node Voronoi fill**.  The
Voronoi assignment is already honest by construction: a pixel is
evaluated by node `idx_v` only when `abs(z_f − z_v) ≤ R_gate`, i.e.
strictly inside *that* node's own verified disc — never extrapolated;
an uncovered pixel stays `NaN`.  The blend's honesty contract (combine
only strictly-in-disc evaluations) is thus already met by the trivial
one-term case the evidence shows is overwhelmingly dominant.

## Fail-fast contract (Rule 1)

`_validity_radius` throws an `ArgumentError` with a `suggestion` /
`detail` for a degenerate node — an empty shared-`Q` denominator
(`Q` must carry at least the constant term) or a jet too short for the
Jorba–Zou formula (`< 2` coefficients).  A degenerate gate must never
silently return a wrong radius (Rule 1); a thrown error with context is
the honest signal that the Stage-1 walk produced a malformed node.

## References

  - `src/VectorPathNetwork.jl` — the Stage-1 walk that produces the
    visited tree (with the cached `visited_jets`) this module fills
    against; the `vector_path_network_solve` driver that delegates here.
  - `src/VectorStepControl.jl` — `vector_step_jorba_zou`, the per-node
    Jorba–Zou truncation-radius estimator the B1 gate is built on.
  - `src/VectorCoefficients.jl` — `vector_taylor_coefficients`, the jet
    primitive used for the no-cache fallback.
  - `src/PathNetwork.jl` — the v0.1 scalar Stage-2 fill (`:669-693`)
    this module lifts to the shared-`Q` vector case.
  - `src/VectorProblems.jl` — the `_eval_poly` Horner evaluator
    (`:339-345`) and the `Pᵢ(t)/Q(t)` dense-evaluation pattern.
  - `docs/adr/0015-extrapolate-stage-2.md` — the `extrapolate` policy.
  - `docs/adr/0019-shared-denominator-pade.md` — the shared `Q`.
  - `docs/adr/0025-headline-figure-re-resolution.md` Amendment 1 §A1 —
    the B1 true-radius gate (gate v-b); the production formula.
    Amendment 2 — the B4-rescoped honest `|u|` surface underlay.
  - `external/probes/pade-validity-radius/REPORT.md` §5–§7 — the A1
    spike: the head-to-head evidence and the production formula.
  - `external/probes/stage2-blend-audition/audition.jl` — the B4
    Voronoi-vs-blend audition: the measured evidence that the B1 discs
    barely overlap, so the nearest-node Voronoi fill is retained.
"""
module VectorPathNetworkStage2

using Polynomials:          Polynomial, roots
using ..VectorStepControl:  vector_step_jorba_zou
using ..VectorCoefficients: vector_taylor_coefficients

# -----------------------------------------------------------------------------
# B1 true-radius validity gate (ADR-0025 Amendment 1, gate v-b)
# -----------------------------------------------------------------------------

"""
    _safety_factor(tol::Real) -> Float64

The per-tolerance safety factor `s(tol)` of the B1 gate
(ADR-0025 Amendment 1; `external/probes/pade-validity-radius/REPORT.md`
§7).  Calibrated on the 389-node V8b wedge walk so the gate's
over-estimation rate matches the global gate's conservative ~1–2 % tail:

    s(1e-6) = 0.34 ,  s(1e-8) = 0.36 (default) ,  s(1e-10) = 0.30

The spike calibrated only three tolerances; an arbitrary `tol` is
mapped to the calibration point nearest in `log10(tol)`.  This is a
deliberate clamp, not an interpolation — `s` is a measured tail
statistic, not a smooth function, and extrapolating it would invent
calibration the spike never did.
"""
function _safety_factor(tol::Real)
    L = log10(float(tol))
    # Nearest calibration point in log10(tol): 1e-6, 1e-8, 1e-10.
    if L > -7.0
        return 0.34
    elseif L > -9.0
        return 0.36
    else
        return 0.30
    end
end

"""
    _validity_radius(jet, denominator, h_v, tol) -> real radius R_gate

The per-node true validity radius of the B1 gate (ADR-0025 Amendment 1,
gate v-b — see the module docstring for the full rationale):

```
R_gate = min( s(tol) · h_JZ · h_v ,  0.5 · h_v · min|t*| )
```

`jet` is the node's raw order-`order` vector Taylor jet (`d` coefficient
vectors, low-to-high, in the un-rescaled `h'`-variable — exactly
`VectorPathNetworkSolution.visited_jets[k]`).  `denominator` is the
node's shared-`Q` polynomial (low-to-high coefficients).  `h_v` is the
real canonical step magnitude.  `tol` is the truncation tolerance.

The jet is first **rescaled** to the `t`-variable (`c̃_k = c_k · h_vᵏ`)
so `vector_step_jorba_zou` returns `h_JZ` in that variable; `h_JZ · h_v`
is the z-plane radius.  `min|t*|` is the smallest-magnitude root of
`Q`; `0.5 · h_v · min|t*|` is the pole-adjacency clamp (binds only when
a pole sits inside the truncation disc).

Throws `ArgumentError` (Rule 1) if `denominator` is empty (a node with
no shared `Q` is malformed) or if every component jet has fewer than 2
coefficients (the Jorba–Zou formula needs the last two coefficients).
A jet whose coefficients are all zero, or a constant `Q`, are *not*
errors: a zero jet means the local solution is constant (the truncation
term is then `Inf` and the clamp — or, for a constant `Q`, nothing —
governs); a constant `Q` has no roots, so the clamp is simply absent
and the truncation term alone gives `R_gate` (the disc is bounded by
truncation error, never reaching a pole).
"""
function _validity_radius(jet::AbstractVector{<:AbstractVector{C}},
                          denominator::AbstractVector{C},
                          h_v::Real, tol::Real) where {C}
    R = float(real(C))

    isempty(denominator) && throw(ArgumentError(
        "VectorPathNetworkStage2._validity_radius: the node's shared-Q " *
        "denominator is empty — a visited node must carry at least the " *
        "constant term Q(0).  Detail: an empty Q is a malformed Stage-1 " *
        "node, not a degenerate-but-valid one.  Suggestion: check the " *
        "Stage-1 walk did not push a node with no canonical approximant."))
    isempty(jet) && throw(ArgumentError(
        "VectorPathNetworkStage2._validity_radius: the node's Taylor jet " *
        "is empty — no component series to estimate a truncation radius " *
        "from.  Suggestion: check the cached `visited_jets`, or the " *
        "no-cache fallback's `vector_taylor_coefficients` call."))
    minimum(length, jet) ≥ 2 || throw(ArgumentError(
        "VectorPathNetworkStage2._validity_radius: a component Taylor jet " *
        "has fewer than 2 coefficients; the Jorba–Zou truncation-radius " *
        "formula needs the last two Taylor coefficients.  Suggestion: " *
        "build the visited tree with `order ≥ 1`."))

    # --- truncation term: s(tol) · h_JZ · h_v ---------------------------------
    # Rescale each component jet to the t-variable: c̃_k = c_k · h_vᵏ.  This
    # is the FW 2011 §3.2 substitution the canonical shared-Q Padé lives in
    # (`VectorStepper._rescale_by_powers`); `vector_step_jorba_zou` then
    # returns h_JZ as a step in t, so h_JZ · h_v is the z-plane radius.
    rescaled = [_rescale_by_powers(j, C(h_v)) for j in jet]
    h_JZ = vector_step_jorba_zou(rescaled, R(tol))
    R_trunc = _safety_factor(tol) * h_JZ * R(h_v)

    # --- pole-adjacency clamp: 0.5 · h_v · min|t*| ----------------------------
    # min|t*| over the roots of the shared-Q denominator in the t-variable.
    # A constant Q (length 1) has no roots — the clamp is then absent and the
    # disc is purely truncation-limited (it never reaches a pole).
    R_gate = R_trunc
    if length(denominator) ≥ 2
        Q_roots = roots(Polynomial(collect(denominator)))
        if !isempty(Q_roots)
            min_t = minimum(abs, Q_roots)
            R_clamp = R(0.5) * R(h_v) * R(min_t)
            R_gate = min(R_gate, R_clamp)
        end
    end
    return R_gate
end

"""
    _rescale_by_powers(c::AbstractVector{C}, h::C) -> Vector{C}

Return `c̃[k+1] = hᵏ · c[k+1]` — the FW 2011 §3.2 variable substitution
`h' → h·t`.  A private copy of `VectorStepper._rescale_by_powers`,
duplicated (not imported) to keep this module free of any
`VectorStepper` dependency — the same module-boundary discipline that
keeps `_eval_poly` a local copy.
"""
function _rescale_by_powers(c::AbstractVector{C}, h::C) where {C}
    out   = Vector{C}(undef, length(c))
    h_pow = one(C)
    @inbounds for k in 1:length(c)
        out[k] = c[k] * h_pow
        h_pow  = h_pow * h
    end
    return out
end

# -----------------------------------------------------------------------------
# Stage-2 fine-grid fill
# -----------------------------------------------------------------------------

"""
    _stage2_fill(grid_z, visited_z, visited_h, visited_numerators,
                 visited_denominator, visited_jets, extrapolate, tol,
                 f, visited_y, order, ::Type{C}, nearest_visited)
        -> Vector{Vector{C}}

Evaluate the dense vector solution `y(z)` at every point of `grid_z`
against the nearest visited node's canonical shared-`Q` approximant,
gated by the B1 true-radius validity gate (ADR-0025 Amendment 1 — see
the module docstring).

`nearest_visited` is the Stage-1 nearest-visited-node scan
(`VectorPathNetwork._nearest_visited`), passed in so this module needs
no `VectorPathNetwork` import (see the module docstring).  It is called
as `nearest_visited(visited_z, z_f) -> Int`.

`visited_jets` is the cached per-node order-`order` Taylor jet; when it
is empty (a pre-B1 backward-compat solution) the gate recomputes each
node's jet from `(visited_z, visited_y)` via `vector_taylor_coefficients`
with the supplied RHS `f` and `order` — the documented fallback path.

`tol` is the B1 gate truncation tolerance.  Each node's `R_gate` is
computed once and memoised, since many grid points share a nearest node.

Returns one `d`-vector per grid point.  A grid point outside its
nearest visited node's true validity disc `R_gate` gets a `d`-vector of
`NaN + NaN·im` unless `extrapolate` is `true` (ADR-0015).
"""
function _stage2_fill(grid_z::Vector{C},
                      visited_z::Vector{C},
                      visited_h::Vector{C},
                      visited_numerators::Vector{Vector{Vector{C}}},
                      visited_denominator::Vector{Vector{C}},
                      visited_jets::Vector{Vector{Vector{C}}},
                      extrapolate::Bool,
                      tol::Real,
                      f,
                      visited_y::Vector{Vector{C}},
                      order::Integer,
                      ::Type{C},
                      nearest_visited) where {C}
    T      = real(C)
    nan_C  = C(T(NaN), T(NaN))
    grid_y = Vector{Vector{C}}(undef, length(grid_z))

    # Per-node R_gate memo: many grid points share a nearest node, and
    # the gate (a Jorba–Zou formula + a polynomial root-find) is the
    # expensive part of the fill.  `nothing` ⇒ not yet computed.
    has_cache = !isempty(visited_jets)
    n_nodes   = length(visited_z)
    radius    = Vector{Union{Nothing, T}}(nothing, n_nodes)

    for (i, z_f) in enumerate(grid_z)
        idx_v      = nearest_visited(visited_z, z_f)
        z_v        = visited_z[idx_v]
        h_v        = real(visited_h[idx_v])      # canonical step is real
        numerators = visited_numerators[idx_v]

        if !extrapolate
            R_gate = radius[idx_v]
            if R_gate === nothing
                # B1 gate: the node's true validity radius.  Use the
                # cached jet when present; otherwise recompute it (the
                # documented fallback for a pre-B1 solution).
                jet = has_cache ? visited_jets[idx_v] :
                      vector_taylor_coefficients(f, z_v, visited_y[idx_v],
                                                 Int(order))
                R_gate = _validity_radius(jet, visited_denominator[idx_v],
                                          h_v, tol)
                radius[idx_v] = R_gate
            end
            if abs(z_f - z_v) > R_gate
                # Grid point outside the node's verified disc — fail-soft.
                grid_y[i] = fill(nan_C, length(numerators))
                continue
            end
        end

        t   = (z_f - z_v) / h_v
        q_t = _eval_poly(visited_denominator[idx_v], t)
        grid_y[i] = C[_eval_poly(num, t) / q_t for num in numerators]
    end
    return grid_y
end

"""
    _eval_poly(c::AbstractVector{C}, t) -> C

Horner evaluation of the polynomial `Σ c[k+1] tᵏ` (coefficients
low-to-high) — the convention the shared-`Q` numerators and denominator
are stored in.  A private copy of `VectorProblems._eval_poly`
(`src/VectorProblems.jl:339-345`): copied rather than imported so this
module does not reach into a sibling module's private symbol
(CLAUDE.md — module-boundary discipline).
"""
function _eval_poly(c::AbstractVector{C}, t) where {C}
    s = zero(C)
    @inbounds for k in length(c):-1:1
        s = s * t + c[k]
    end
    return s
end

end # module VectorPathNetworkStage2
