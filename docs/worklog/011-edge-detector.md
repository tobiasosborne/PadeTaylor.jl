# Worklog 011 — EdgeDetector: FW 2011 §3.2.2 pole-field classifier

**Date**: 2026-05-13
**Author**: Claude Opus (orchestrating); one Sonnet subagent (FW §3.2.2 spec dive)
**Scope**: Phase 12.5 / bead `padetaylor-c2p` — the 5-point Laplacian
pole-field edge detector per FW 2011 §3.2.2.  ~60-LOC module + 743
detector-specific assertions; clean mutation-proof on the two
load-bearing parameters (stencil centre coefficient, `1/h²` scaling).
Total suite 1148/1148 GREEN (404 prior + 743 new + 1 umbrella).

> Take-home: a tiny, scope-bounded diagnostic primitive — eq. (3.3)
> verbatim, no over-engineering, returns the raw residual matrix
> plus a thresholded BitMatrix.  Composes naturally with PathNetwork
> dense output and will be the partition oracle for Phase 12 v2
> (`padetaylor-k31`).

## Why now

`padetaylor-c2p` was originally filed alongside Phase 12 v1; worklog 007
noted that v2 (`padetaylor-k31`) would "consume the 5-point Laplacian
classifier formerly tracked as `padetaylor-c2p`".  Two reasons to ship
the detector as its own standalone deliverable rather than fold it into
v2:

  1. **Phase 9 (`padetaylor-kvi`) needs it independently** — the Tier-C
     qualitative PI tritronquée pole-field plot uses the detector to
     verify sectorial structure without requiring the full 2D
     dispatcher.  Shipping the detector first unblocks Phase 9.
  2. **The bead description has a slight misread** — it says
     "|∇²u|/h² > 0.001 threshold (FW's empirical level)".  FW's eq. (3.3)
     already includes the `/h²` inside `Δu`, and Fig. 3.3's level curve
     0.001 is on `log₁₀|Δu|`, not on bare `|Δu|`.  Catching this before
     v2 consumes the detector matters; documenting the correction in
     the module docstring is cleaner than a v2 follow-up patch.

## Ground-truth disambiguation (Sonnet research + Opus verification)

The Sonnet subagent (read-only Explore) read FW 2011 §3.2.2 lines
197-208 and reported the verbatim formulation.  One small disagreement
between Sonnet's interpretation and the bead's interpretation needed
resolving:

  - **Sonnet's read**: "Fig. 3.3 shows log₁₀|Δu| ... contour level 0.001"
    means the threshold is `log₁₀|Δu| > 0.001`, i.e., `|Δu| > 10^0.001
    ≈ 1.0023`.
  - **Bead description's read**: "|∇²u|/h² > 0.001" — the threshold is on
    bare `|Δu|` (since FW's `Δu` already has the `/h²` inside).

Opus opened the FW2011 markdown lines 197-208 directly and confirmed
Sonnet's reading: line 208 is verbatim "Fig. 3.3 shows log₁₀|Δu| ...
it is easy to select a contour level (here 0.001) that gives a suitable
pole field edge description."  The figure displays log₁₀|Δu| as a
surface; the contour at LEVEL=0.001 is on that displayed quantity.

Going with FW's reading.  The 1000× difference between the two
interpretations is large, but at the practical-classification level
both give plausible separation between "essentially zero" and
"non-zero" |Δu| — just at different sensitivity levels.  Made the
threshold a kwarg with default 0.001 (FW's value), so callers who
prefer the bead's reading can pass `level = -3` for the equivalent
linear threshold `|Δu| > 10⁻³`.

## What changed

`src/EdgeDetector.jl` (~60 LOC + ~115 lines of literate docstring):

  - `laplacian_residual(u_grid::AbstractMatrix{Complex{T}}, h::Real)`
    — eq. (3.3) verbatim.  Returns a same-shape complex matrix with
    interior cells populated and boundary cells `NaN + NaN·im`.  Uses
    `@inbounds` over the explicit 5-point neighbourhood loop; one
    multiplication by `inv(h²)` per cell.
  - `pole_field_mask(u_grid, h; level=0.001)` — convenience wrapper.
    Returns `BitMatrix` of the same shape; boundary cells `false`.
  - `pole_field_mask(Δu; level=0.001)` — variant for pre-computed
    residuals (so callers can threshold at multiple levels without
    re-computing the stencil).

`src/PadeTaylor.jl`: `include("EdgeDetector.jl")` + `using .EdgeDetector:
laplacian_residual, pole_field_mask` + `export ...`.

`test/runtests.jl`: `include("edge_detector_test.jl")` + extra
`isdefined` line in umbrella testset.

`test/edge_detector_test.jl` (7 testsets, 743 assertions):

  - **ED.1.1**: stencil exact-zero on a real harmonic quadratic `u(x,y)
    = x² − y²`.  The 5-point stencil is exact for any polynomial of
    total degree ≤ 3 (the truncation `O(h²)` term involves the 4th
    partial; absent here).  Asserts `|Δu| < 1e-12` on every interior
    cell of a 21×21 grid at h=0.1.  Boundary cells asserted NaN.
  - **ED.1.2**: same as 1.1 but with the complex-analytic `u(z) = z²`.
    Both real (`x²−y²`) and imaginary (`2xy`) parts are harmonic.
  - **ED.1.3**: pole separation — `u(z) = 1/(z − 0)` on a grid offset
    from `z = 0` (so the pole sits between lattice cells).  Asserts
    `|Δu| > 100` at the near-pole interior cell and `|Δu_corner| <
    |Δu_centre| / 10` at the farthest corner.
  - **ED.2.1**: `pole_field_mask` returns a `BitMatrix`, boundary
    cells `false`, near-pole centre cell `true`.
  - **ED.2.2**: separation property — `u = exp(z)` on a 41×41 grid
    gives a `count(mask) == 0` all-false bitmap; `u = 1/z` near the
    pole gives a non-trivial cluster (empirically ~188 cells flagged,
    bounded `[4, length÷2]`).
  - **ED.3.1**: fail-fast on grid `< 3×3` (no interior cells) and on
    `h ≤ 0`.
  - **ED.4.1**: pre-computed-residual variant agrees with one-step
    wrapper, and varying the level changes the mask monotonically.

Total: 743 assertions (most from the dense per-cell loops in 1.1
and 1.2).

## Mutation-proof (verified 2026-05-13)

**Mutation A** — change stencil centre coefficient from `-4` to `-3`.
Verified bite: 686 fails — ED.1.1 (324/365), ED.1.2 (360/361), ED.2.2
(2/3).  Even on the harmonic quadratic, the residual is no longer zero;
asserted `|Δu| < 1e-12` blows up to `|Δu| ≈ 100` (proportional to `u/h²`
at h=0.1).  ED.1.3's pole test still passes because the dominant
near-pole signal swamps the centre-coefficient bias.

**Mutation B** — drop the `1/h²` factor.  Verified bite: 1 fail —
ED.1.3 (`|Δu_centre| > 100` fails because the magnitude scales down
by `h² = 0.01`).  The harmonic tests (1.1, 1.2) still pass because
the residual is exact-zero with or without `/h²`.  The
`pole_field_mask` test (2.2) still passes because the near-pole cells
remain above the log₁₀ threshold even at the reduced magnitude — by
the margin the original mask had to spare.  This is a clean **targeted**
bite: ED.1.3 was specifically designed to catch the `/h²` scaling.

Both mutations restored before commit per CLAUDE.md Rule 4.

## Frictions surfaced

  - **F1. Bead description vs FW source disagreed on threshold interpretation**.
    Caught by reading FW2011...md:208 verbatim instead of trusting the
    bead text.  Recorded in the module docstring (§"What the bead
    description got slightly wrong") and as the threshold-units paragraph
    near the top.  Future bead authors: when transcribing an algorithm,
    quote the source paper line verbatim in the bead description.
  - **F2. ED.2.2 upper-bound assertion was too tight on first cut**.
    Asserted `count(mask_pole) ≤ 100`; empirically 188 cells flagged at
    h=0.05 on a 41×41 grid (the natural `|∇²(1/z)| > 1` disk extends to
    `|z| ≲ 0.35`).  Loosened to `< length(mask) ÷ 2` (838 of 1681) so
    the assertion still detects the "everything flagged" bug class
    while not being brittle to the specific cluster size.  Comment in
    the test documents the empirical count for future maintainers.

## Pointers

  - [`src/EdgeDetector.jl`](../../src/EdgeDetector.jl) — the module.
  - [`test/edge_detector_test.jl`](../../test/edge_detector_test.jl) —
    7 testsets + mutation-proof commentary at the bottom.
  - `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md`
    lines 197-208 (§3.2.2 + Fig. 3.3 caption), 395-410 (the manual-
    classification admission at line 401).
  - [`docs/figure_catalogue.md`](../figure_catalogue.md) §1 row Fig. 3.3
    — the catalogue's existing acceptance ("Edge-detector bitmap on
    20×20 lattice matches FW's level-0.001 contour to ≤2 lattice cells
    off").  Phase 9 (`padetaylor-kvi`) will pin this once it ships.

## Bead state

Closed in this session:
  - `padetaylor-c2p` — EdgeDetector v1 GREEN at this commit.

Still open from prior sessions (no change):
  - `padetaylor-rgp`, `padetaylor-k31`, `padetaylor-kvi`,
    `padetaylor-bvh`, `padetaylor-grc`, `padetaylor-61j`,
    `padetaylor-8pi`.

`padetaylor-k31` (Phase 12 v2 2D dispatcher) now has its detector
prerequisite shipped.  `padetaylor-kvi` (Phase 9 qualitative
pole-field plot) is the next-up item in this session.

## Hard-won lesson (for HANDOFF.md §"Hard-won")

**14. Trust the source paper over the bead description**.
`padetaylor-c2p`'s bead said the threshold is `|∇²u|/h² > 0.001`, but
FW 2011 line 208 says the threshold is on `log₁₀|Δu|`, which already
absorbs the `/h²` and gives a 1000× different effective sensitivity.
A 5-minute verbatim read of FW2011...md:197-208 caught this before any
code was written.  Future bead authors: include the source quote in
the bead description, not a paraphrase.  Future bead consumers: open
the citation before coding.
