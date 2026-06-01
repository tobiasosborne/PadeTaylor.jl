# Bug sweep B5 — CoordTransforms + transformed RHS (PIII/PV/PVI η,ζ maps)

## Area

The exponential coordinate maps `z ↔ ζ` (PIII, PV), the double-exponential
map `z ↔ η` (PVI), their inverses, the chain-rule Jacobian factors carried
through to the IC-conversion helpers, and the transformed-equation RHS
closures (`P̃_III`, `P̃_V`, ζ-plane `P̃_VI`, η-plane `P̃_VI`).  Files audited
in full:

  - `src/CoordTransforms.jl` (PIII / PV — both transforms, both inverses,
    both RHS factories)
  - `src/SheetTracker.jl` (PVI ζ-plane RHS, η-plane RHS, `z↔η` maps,
    winding-number primitives, `sheet_index`)

`src/BranchTracker.jl` (the direct downstream consumer of `winding_delta`)
was also read and its cut-crossing geometry verified as a courtesy, since it
is the place where the winding primitives feed intermittent path-dependent
sheet decisions.

Method: every algebraic relation was re-derived symbolically with SymPy and
diffed term-by-term against the FFW 2017 markdown equations.  No Julia/Pkg
was run (read-only mandate).

## References checked

  - `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:31`
    — PIII and PV equations in the z-plane.
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:35` — PVI equation.
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:41` — PIII transform
    `z = e^{ζ/2}`, `u(z) = e^{-ζ/2} w(ζ)`.
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:43` — P̃_III equation.
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:45` — PV transform
    `z = e^{ζ}`, `u(z) = w(ζ)`.
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:47` — P̃_V equation.
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:103` — PIII strip
    width: `-2π + 4πs < Im ζ ≤ 2π + 4πs` (4π per sheet).
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:141` — PVI `z = e^ζ`
    map; z=1 → ζ = 2πik.
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:144` — ζ-plane P̃_VI
    equation (eq. 3).
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:146` — second map
    `ζ = e^η`; PVI ζ↔η obtained by `w(ζ)=v(η)`, `ζ=e^η`.
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:149` — lattice
    `η = log|2πk| + i·arg(2πik)`, |k| ≥ 1 (eq. 4).
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:154` — η-plane P̃_VI
    equation (eq. 5), the `-1` chain-rule artefact in the `dv/dη` bracket.
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:157` — branch-point-
    free region `Re η < log 2π`.
  - `.../FFW2017_painleve_riemann_surfaces_preprint.md:180-189` — sheet
    parametrisation `(2k-1)π < θ_k ≤ (2k+1)π` (eq. 6).
  - `docs/worklog/017-coord-transforms-pIII-pV.md:84-88, 160-165, 191-197`
    — hand-derivation of PIII `wp = (zu + z²u')/2`; mutation-proof record.
  - `docs/worklog/041-eta-plane-pvi.md:37-42, 66-75` — M2 sign-flip
    restoration at `SheetTracker.jl:248`; η-plane `-1` artefact provenance.

## Findings

### [LOW] `sheet_index` banker's-rounding disagrees with FFW sheet convention at exact odd half-revolution boundaries

**Location**: `src/SheetTracker.jl:320`
(`sheet_index(total_winding::Real) = round(Int, total_winding / (2 * π))`)

**Ground truth (cited)**: FFW 2017 eq. 6,
`.../FFW2017_painleve_riemann_surfaces_preprint.md:187-189` gives the sheet
parametrisation `(2k-1)π < θ_k ≤ (2k+1)π` — i.e. the angular range for
sheet `k` is the half-open interval centred at `2kπ`, **closed at the upper
end** `(2k+1)π`.  So a winding angle exactly equal to `(2k+1)π` belongs to
sheet `k`.

**Code behavior**: `round(Int, x)` in Julia uses round-half-to-even
(banker's rounding) by default.  At an exact boundary `total = (2k+1)π`,
`total/2π = k + 0.5`, and banker's rounding sends `k+0.5` to the nearest
**even** integer.  For even `k` this gives `k` (agrees with FFW); for **odd
`k`** it gives `k+1` (disagrees — the point is assigned to sheet `k+1`
instead of the FFW-prescribed sheet `k`).  Verified: `θ=3π → round(1.5)=2`
(FFW sheet 1); `θ=7π → round(3.5)=4` (FFW sheet 3); whereas `θ=π →
round(0.5)=0` and `θ=5π → round(2.5)=2` agree.  The reference's half-open
`(...]` convention wants `floor((x - π)/2π) + 1`-style rounding (ties always
toward the lower-sheet's upper-closed boundary), not round-half-to-even.

**Mechanism**: This affects only the integer *sheet label* returned to a
caller, not the computed solution values `w`/`u` themselves, and only when
the accumulated winding lands **exactly** on `(2k+1)π` to floating-point
precision for odd `k`.  In a path-network walk that is a measure-zero event,
and the reference convention is itself a labelling choice at the branch cut
where the assignment is genuinely ambiguous.  It could in principle produce
an order-dependent (hence intermittent) off-by-one in a reported sheet index
if two independently-ordered Stage-1 walks accumulate winding that straddles
an odd `(2k+1)π` boundary from opposite sides — but it does not warp the
continued solution itself.  Not a plausible source of the maintainer's
intermittent *solution* discontinuities; reported for completeness.

**Intermittent?**: Yes (order/data-dependent, but measure-zero and label-only).

**Confidence**: 0.55 — the banker's-rounding-vs-FFW-convention mismatch at
exact odd boundaries is demonstrable, but its real-world reachability is
near nil and the impact is a sheet *label*, not solution geometry.

## Areas verified correct

All of the following were re-derived symbolically (SymPy) and diffed
term-by-term against the cited reference lines; each matches.

  - **PIII forward map `pIII_z_to_ζ`** (`CoordTransforms.jl:120-125`):
    `ζ = 2 log z`, `w = z u`, `wp = (z u + z² u')/2`.  SymPy on
    `w(ζ)=e^{ζ/2}u(z(ζ))` gives `dw/dζ = (z u + z² u')/2` (using
    `e^{ζ/2}=z`, `e^ζ=z²`).  Matches FFW transform at md:41 and the
    worklog-017 hand-derivation (worklog 017:160-165, 191-197), including
    the load-bearing `z²` (not `z`) factor on the `u'` term.

  - **PIII inverse `pIII_ζ_to_z`** (`CoordTransforms.jl:132-137`):
    `z=e^{ζ/2}`, `u=w/z`, `up=(2wp - w)/z²`.  SymPy `solve` of the forward
    relations yields exactly `{u: w/z, up: (2wp - w)/z²}`.  Round trip
    `z → ζ → z` is `exp(2 log z / 2) = z` (exact on principal branch).

  - **P̃_III RHS `pIII_transformed_rhs`** (`CoordTransforms.jl:150-156`):
    `wp²/w + (α w² + γ w³ + β e^ζ + δ e^{2ζ}/w)/4`.  Matches FFW eq. at
    md:43 term-by-term (the `1/4` divisor, all four parameter terms,
    `e^ζ`/`e^{2ζ}` powers).

  - **PV forward map `pV_z_to_ζ`** (`CoordTransforms.jl:169-174`):
    `ζ=log z`, `w=u`, `wp = z u'`.  SymPy on `w(ζ)=u(e^ζ)` gives
    `dw/dζ = e^ζ u' = z u'`.  Matches FFW transform at md:45.

  - **PV inverse `pV_ζ_to_z`** (`CoordTransforms.jl:181-186`):
    `z=e^ζ`, `u=w`, `up=wp/z`.  Round trip exact on principal branch.

  - **P̃_V RHS `pV_transformed_rhs`** (`CoordTransforms.jl:195-204`):
    `(1/(2w)+1/(w-1)) wp² + (w-1)²(α w + β/w) + γ e^ζ w +
    δ e^{2ζ} w(w+1)/(w-1)`.  Matches FFW eq. at md:47 term-by-term,
    including the `(w+1)/(w-1)` orientation (worklog-017 mutation N
    confirmed this denominator/numerator orientation is load-bearing).

  - **PVI ζ-plane RHS `pVI_transformed_rhs`** (`SheetTracker.jl:155-174`):
    the three terms `(1/2)(1/w+1/(w-1)+1/(w-e^ζ)) wp²`,
    `-(e^ζ/(e^ζ-1)+e^ζ/(w-e^ζ)) wp`, and
    `w(w-1)(w-e^ζ)/(e^ζ-1)² · (α + β e^ζ/w² + γ(e^ζ-1)/(w-1)² +
    δ e^ζ(e^ζ-1)/(w-e^ζ)²)` all match FFW eq. 3 at md:144 term-by-term,
    including every sign and every `(e^ζ-1)` vs `(e^ζ)` factor.

  - **PVI η-plane RHS `pVI_eta_transformed_rhs`** (`SheetTracker.jl:206-226`):
    independently re-derived the chain rule `ζ=e^η`, `v(η)=w(ζ)`:
    `d²v/dη² = F·e^{2η} + vp` with `w'(ζ)=e^{-η}vp`.  Term 1 → the
    `e^{2η}` cancels `(e^{-η})²` leaving `(1/2)(...)vp²` (code `vp^2/2`,
    correct, **no spurious e^{2η}**).  Term 2 → `-e^η(...)vp` and the
    leftover `+vp` combine to `-(e^η E/(E-1)+e^η E/(v-E) - 1)vp` — the
    `-1` artefact present and on the correct side (code line 216, matches
    FFW eq. 5 md:154).  Term 3 → `v(v-1)(v-E)·e^{2η}/(E-1)²·param`, where
    code's `eη^2 = exp(η)^2 = e^{2η}` is the correct factor (code line
    222).  `E = exp(eη) = e^{e^η}` (line 209) is the correct nested
    exponential.  Worklog-041:37-42 confirms M2/M1 mutations restored.

  - **PVI `z↔η` maps `pVI_z_to_η` / `pVI_η_to_z`**
    (`SheetTracker.jl:244-265`): forward `vp = z·ζ·up = z·log(z)·up`,
    inverse `up = vp/(z·ζ)`.  SymPy confirms both the `z→η→z` round trip
    and the `up` round trip are exact (residual 0) on the principal
    branch.  The Jacobian factor `z·ζ = z·log z` is present, on the
    correct side, and **not inverted** (forward multiplies, inverse
    divides) — this is the specific reciprocal-error class flagged in the
    assignment, and it is correct here.  Reuse of `pV_z_to_ζ` for the PVI
    ζ-frame (Painleve.jl:281) is the mathematically identical PV map
    (md:146), correct.

  - **`winding_delta` normalisation** (`SheetTracker.jl:283-292`): the
    single-correction `if Δθ ≤ -π … elseif Δθ > π` is sufficient because
    the difference of two `angle()` values (each in `(-π,π]`) lies in
    `(-2π,2π)`; verified the result lands in `(-π,π]` across edge cases
    (±π, ±(π+ε), ±1.9π).  Sign convention `arg(new-branch) -
    arg(old-branch)` is positive for counterclockwise (matches docstring
    and the FFW sheet-increment convention).  The `|Δθ| < π` per-step
    precondition is a documented caller contract (lines 275-281), not a
    bug.

  - **`accumulate_winding`** (`SheetTracker.jl:304-311`): cumulative sum
    with `out[1]=0`, no buffer aliasing (writes to a fresh `out`, reads
    `path[i-1]`/`path[i]` only).  No RNG/order dependence within the
    function.

  - **`segment_crosses_cut` Cramer's rule** (`BranchTracker.jl:112-121`,
    consumer — outside assigned files but verified): solved the 2×2 system
    `t·d - s·u = δ` symbolically; the code's `det = imag(d·conj u)`,
    `t = imag(δ·conj u)/det`, `s = imag(δ·conj d)/det` match the analytic
    solution exactly, including all `conj` placements and signs.  No
    transpose/adjoint confusion here.

  - **Sheet strip widths** (`CoordTransforms.jl:17,48,56`): PIII `4π·im·s`
    (matches md:103 strip `4π` per sheet); PV `2π·im·s` (matches `z=e^ζ`
    period `2πi`).  Correct.
