# test/corpus_morphology_test.jl -- bead `padetaylor-rn1d`, epic
# `padetaylor-p1v0`.  Corpus bucket B9 (strong-regenerable morphology).
#
# WHAT THIS PINS — the discrete morphology internals that gate the
# edge-gated pole-field solve (FW 2011 md:401).  `EdgeGatedSolve`'s
# region-growing algorithm rests on three morphological primitives —
# dilation (grow the target ring), opening (= erode∘dilate, kill
# false-positive specks and thin bridges), flood-fill (keep only mask
# cells connected to the field) — plus `LatticeDispatcher._dilate_one`,
# the one-ring lift used to scatter derivative values onto `up_grid`.
# Capability map B9 records these internals as having ZERO direct
# tests; `edge_gated_solve_test.jl` / `lattice_dispatcher_test.jl`
# exercise them only transitively through full IVP solves.  A wrong
# structuring element, a flipped boundary convention, or a 4-vs-8
# connectivity mismatch between `_dilate` and `_dilate_one` would
# silently misplace pole-field cells (and the derivative values at
# diagonal frontier cells) with no failing test to catch it.  These
# fixtures are hand-computed exact BitMatrices — small enough to read
# by eye — cross-checked against BOTH the package's actual convention
# (read from src) AND the canonical Serra (1982) definitions.
#
# THE PACKAGE'S ACTUAL CONVENTIONS (read from src, not assumed):
#   * Structuring element: the 3×3 all-ones box (8-connected /
#     Chebyshev), i.e. every offset in {-1,0,1}² \ {(0,0)} plus the
#     centre.  `_dilate`/`_erode` iterate `for dj in -1:1, di in -1:1`
#     skipping `(0,0)` — the 8 neighbours; the centre is implicit
#     (a cell already true stays true on dilate; a true cell whose 8
#     neighbours are all true survives erosion). src/EdgeGatedSolve.jl
#     :140 (_dilate "Chebyshev-distance (8-connected)"), :163 (_erode).
#   * Boundary / out-of-grid: dilation treats out-of-grid as FALSE
#     (a missing neighbour cannot turn a cell on); erosion treats
#     out-of-grid as FALSE too — a window hanging off the lattice
#     ERODES THE CELL AWAY (src:173 `!(1 ≤ ii ≤ nx && 1 ≤ jj ≤ ny) ||
#     !prev[ii,jj]` → set false).  This is the standard zero-padding
#     convention, and it matches the _erode docstring verbatim
#     ("grid-edge cells (a window hanging off the lattice) erode away").
#   * Flood-fill connectivity: 8-connected (src:203 `for dj in -1:1,
#     di in -1:1`).
#   * `LatticeDispatcher._dilate_one` (src/LatticeDispatcher.jl:447):
#     8-connected, out-of-grid FALSE — CONSISTENT with `_dilate(·,1)`.
#
# Tests CMO.1 (dilation single pixel), CMO.2 (erosion edge-hang —
# pins the boundary convention), CMO.3 (erosion 5×5→single pixel),
# CMO.4 (opening idempotence + clean-block fixed point), CMO.5
# (flood-fill: two separated components), CMO.6 (cross-module:
# `_dilate` ≡ `_dilate_one` connectivity — the gap-#11 question).
# Mutation-proof procedure in the file footer.
#
# Catalogue records: docs/test_corpus/01_corpus_catalogue.md:1739-1818
# (morph-dilation-single-pixel, morph-erosion-edge-hang,
# morph-erosion-5x5, morph-opening-idempotent, morph-connected-
# components-8conn).  Catalogue oracles are 0-based row-major; the
# Julia BitMatrix is 1-based `mask[i,j]` with `i` the first index.
# Because the SE is the symmetric 3×3 box and connectivity is
# 8-symmetric, dilation/erosion outputs are index-convention-invariant;
# only the flood-fill membership labels shift by +1 per axis (0-based
# (r,c) ↦ 1-based (r+1,c+1)).
#
# References:
#   Serra, J. (1982) Image Analysis and Mathematical Morphology,
#     Academic Press, ch. 2 — dilation A⊕B = {q : (B̌)_q ∩ A ≠ ∅};
#     erosion A⊖B = {q : B_q ⊆ A}; opening γ_B(A) = (A⊖B)⊕B, which
#     is idempotent and anti-extensive.
#   Haralick & Shapiro (1992) Computer and Robot Vision Vol. I, ch. 2 —
#     8-connected connected-component labelling.

using Test
using PadeTaylor

const _EGS = PadeTaylor.EdgeGatedSolve
const _LD  = PadeTaylor.LatticeDispatcher

# Build a BitMatrix from a 0-based row-major nested-array catalogue
# layout: `rows[i][j]` is catalogue (row=i-1, col=j-1) ↦ Julia (i,j).
_bm(rows) = BitMatrix([rows[i][j] == 1 for i in 1:length(rows), j in 1:length(rows[1])])

# All `true` cells of a BitMatrix, as a sorted Set of 1-based (i,j).
_cells(M) = Set((i, j) for i in 1:size(M,1), j in 1:size(M,2) if M[i, j])

@testset "Corpus morphology (bead padetaylor-rn1d): B9 internals" begin

    @testset "CMO.1 dilation: lone pixel → full 3×3 SE footprint" begin
        # Serra dilation of a singleton by the 3×3 box = the box itself.
        # Input  [[0,0,0],[0,1,0],[0,0,0]]  ⊕ 3×3-box  =  all-ones 3×3.
        A = _bm([[0,0,0],[0,1,0],[0,0,0]])
        expected = _bm([[1,1,1],[1,1,1],[1,1,1]])
        @test _EGS._dilate(A, 1) == expected
        # The lone centre's 8 neighbours all switch on; nothing more to
        # switch (the 3×3 grid IS the footprint).  Standard Serra value.
    end

    @testset "CMO.2 erosion edge-hang: zero-padding kills the boundary" begin
        # An all-true 3×3 eroded by the 3×3 box.  Only the centre's
        # window fits entirely inside the grid; every boundary cell has
        # a window hanging off the lattice edge → out-of-grid = FALSE →
        # eroded away.  This fixture PINS the boundary convention.
        B = _bm([[1,1,1],[1,1,1],[1,1,1]])
        expected = _bm([[0,0,0],[0,1,0],[0,0,0]])     # code's convention
        @test _EGS._erode(B, 1) == expected
        # Standard Serra with zero-padding gives the SAME result, and it
        # matches the _erode docstring ("grid-edge cells … erode away").
        # An infinite-domain / clamped-boundary convention would instead
        # leave the whole 3×3 true — that is NOT what the code (or its
        # docstring) does.  Code ≡ docstring ≡ Serra zero-pad: no finding.
    end

    @testset "CMO.3 erosion: 5×5 with 3×3 centre → single pixel" begin
        # Input B: a 3×3 true block centred in a 5×5 zero field.
        #   [[0,0,0,0,0],
        #    [0,1,1,1,0],
        #    [0,1,1,1,0],
        #    [0,1,1,1,0],
        #    [0,0,0,0,0]]
        # Erosion by 3×3 box: only (3,3) (Julia) — catalogue (2,2) —
        # has an all-true 3×3 window (rows 2:4, cols 2:4 = the block).
        B = falses(5, 5); B[2:4, 2:4] .= true
        expected = falses(5, 5); expected[3, 3] = true
        @test _EGS._erode(B, 1) == BitMatrix(expected)
        # Every other block cell touches a border zero; every border
        # cell is zero already.  Serra erosion shrinks the block to a
        # point — confirmed.
    end

    @testset "CMO.4 opening: idempotent + clean-block fixed point" begin
        # Opening γ = dilate∘erode (src:186 `_open = _dilate(_erode)`).
        # The 3×3 block is a union of SE translates, so by Serra's
        # opening theorem it is a FIXED POINT: γ(B) = B.  And opening is
        # idempotent: γ(γ(X)) = γ(X) for any X.
        B = falses(5, 5); B[2:4, 2:4] .= true; B = BitMatrix(B)
        opened = _EGS._open(B, 1)
        @test opened == B                       # clean block is fixed
        @test _EGS._open(opened, 1) == opened   # idempotence on B
        # Idempotence on a NON-fixed input too (a speck the opening
        # removes), so the property is not vacuously true on fixed sets:
        X = falses(5, 5); X[2:4, 2:4] .= true; X[1, 1] = true; X = BitMatrix(X)
        gX = _EGS._open(X, 1)
        @test _EGS._open(gX, 1) == gX           # γ(γ(X)) = γ(X)
        @test gX == B                           # the lone speck is opened away
    end

    @testset "CMO.5 flood-fill: two separated components, 8-connected" begin
        # Catalogue CC (0-based row-major):
        #   [[1,1,0,0,0],
        #    [1,0,0,0,0],
        #    [0,0,0,0,1],
        #    [0,0,0,1,1],
        #    [0,0,0,0,0]]
        # Two L-shapes, non-adjacent under BOTH 4- and 8-connectivity
        # (min col-gap between the clusters is 2).  In 1-based Julia
        # coords: comp1 = {(1,1),(1,2),(2,1)}, comp2 = {(3,5),(4,4),(4,5)}.
        mask = _bm([[1,1,0,0,0],
                    [1,0,0,0,0],
                    [0,0,0,0,1],
                    [0,0,0,1,1],
                    [0,0,0,0,0]])

        # Seed comp1 from its corner (1,1); the flood-fill must reach
        # exactly comp1 and NOT leak across the gap into comp2.
        seed1 = falses(5, 5); seed1[1, 1] = true
        comp1 = _EGS._flood_fill(BitMatrix(seed1), mask)
        @test _cells(comp1) == Set([(1,1), (1,2), (2,1)])

        # Seed comp2 from (3,5); reaches exactly comp2.
        seed2 = falses(5, 5); seed2[3, 5] = true
        comp2 = _EGS._flood_fill(BitMatrix(seed2), mask)
        @test _cells(comp2) == Set([(3,5), (4,4), (4,5)])

        # The two are disjoint and together tile every mask cell ⇒
        # exactly two connected components, none stranded.
        @test isempty(intersect(_cells(comp1), _cells(comp2)))
        @test union(_cells(comp1), _cells(comp2)) == _cells(mask)

        # Independent count by relabelling, to pin n_components = 2.
        visited = falses(5, 5); ncomp = 0
        for i in 1:5, j in 1:5
            if mask[i, j] && !visited[i, j]
                ncomp += 1
                s = falses(5, 5); s[i, j] = true
                visited .|= _EGS._flood_fill(BitMatrix(s), mask)
            end
        end
        @test ncomp == 2

        # Diagonal frontier: seeding comp2 must NOT reach (1,2)/(2,1) of
        # comp1 — confirms 8-connectivity does not over-bridge the gap.
        @test !comp2[1, 2] && !comp2[2, 1]

        # 8-connectivity discriminator.  The two-L fixture above is also
        # 4-connected within each component, so on its own it does NOT
        # distinguish 4- from 8-connectivity.  This diagonal chain does:
        # (1,1) and (2,2) are 8-adjacent but NOT 4-adjacent, so an
        # 8-connected flood-fill reaches (2,2) from (1,1) — a 4-connected
        # one would strand it.  This pins the documented 8-connectivity.
        diag = falses(3, 3); diag[1, 1] = true; diag[2, 2] = true
        dseed = falses(3, 3); dseed[1, 1] = true
        dreach = _EGS._flood_fill(BitMatrix(dseed), BitMatrix(diag))
        @test dreach[2, 2]                       # 8-connected reaches diagonal
        @test _cells(dreach) == Set([(1,1), (2,2)])
    end

    @testset "CMO.6 cross-module: _dilate ≡ _dilate_one connectivity" begin
        # Gap-#11: a 4-vs-8 connectivity mismatch between
        # EdgeGatedSolve._dilate (the region-growth ring) and
        # LatticeDispatcher._dilate_one (the up_grid scatter ring) would
        # silently misplace derivative values at DIAGONAL frontier cells.
        # Assert they agree on a single-cell input — the case that
        # distinguishes 4- from 8-connectivity (a diagonal neighbour is
        # filled only under 8-connectivity).
        D = falses(5, 5); D[3, 3] = true; D = BitMatrix(D)
        diag_filled = _bm([[0,0,0,0,0],
                           [0,1,1,1,0],
                           [0,1,1,1,0],
                           [0,1,1,1,0],
                           [0,0,0,0,0]])
        @test _EGS._dilate(D, 1) == diag_filled        # 8-connected
        @test _LD._dilate_one(D) == diag_filled        # 8-connected
        @test _EGS._dilate(D, 1) == _LD._dilate_one(D) # ⇒ they AGREE

        # And on the asymmetric two-component fixture, for good measure:
        mask = _bm([[1,1,0,0,0],
                    [1,0,0,0,0],
                    [0,0,0,0,1],
                    [0,0,0,1,1],
                    [0,0,0,0,0]])
        @test _EGS._dilate(mask, 1) == _LD._dilate_one(mask)
        # FINDING: connectivity AGREES (both 8-connected, both
        # out-of-grid = FALSE).  Gap-#11 is NOT a latent bug here — the
        # diagonal-frontier derivative scatter is placed consistently.
    end
end

# ----------------------------------------------------------------------
# Mutation-proof procedure (CLAUDE.md Rule 4 — port-and-verify shape).
#
# Each mutation below was applied to src, this file run via
# `julia --project=. test/corpus_morphology_test.jl`, confirmed RED,
# then reverted byte-clean (`git diff src` empty afterwards).
#
#   M1 — flip a dilation neighbour offset.  In
#        `EdgeGatedSolve._dilate`, change the inner loop
#        `for dj in -1:1, di in -1:1` to `for dj in 0:1, di in 0:1`
#        (drop the negative offsets → a quadrant SE, not the full box).
#        CMO.1 (lone-pixel → all-ones 3×3) goes RED: the down/left
#        neighbours of the centre never switch on.  CMO.6 also REDs
#        (`_dilate` no longer matches `_dilate_one`).  Reverted.
#
#   M2 — flip the erosion boundary convention.  In
#        `EdgeGatedSolve._erode`, change the survival test
#        `if !(1 ≤ ii ≤ nx && 1 ≤ jj ≤ ny) || !prev[ii, jj]` to
#        `if (1 ≤ ii ≤ nx && 1 ≤ jj ≤ ny) && !prev[ii, jj]`
#        (treat out-of-grid as TRUE — the clamped/infinite-domain
#        convention).  CMO.2 (edge-hang) goes RED: the all-true 3×3 no
#        longer erodes its boundary, so the output is all-true, not the
#        single centre pixel.  (CMO.3's 5×5 block is interior — its
#        relevant windows never hang off the lattice — so CMO.3 stays
#        GREEN; the all-true-grid CMO.2 is the fixture that pins the
#        boundary, exactly as designed.)  Reverted — confirms the test
#        pins the zero-padding convention the docstring claims.
#
#   M3 — restrict flood-fill to 4-connectivity.  In
#        `EdgeGatedSolve._flood_fill`, add `(abs(di)+abs(dj)==2) &&
#        continue` inside the neighbour loop (skip diagonals).  CMO.5's
#        diagonal discriminator goes RED: (1,1) and (2,2) are only
#        8-adjacent, so a 4-connected flood-fill strands (2,2) —
#        `dreach[2,2]` is false and `_cells(dreach) == {(1,1)}`.
#        (The two-L fixture alone does NOT bite: each L is also
#        4-connected; that is why the explicit diagonal chain was
#        added.)  Reverted — confirms the test pins 8-connected
#        flood-fill.
#
#   (A fourth candidate — perturbing `_dilate_one` alone — makes CMO.6
#   RED by construction, which is exactly the gap-#11 mismatch this
#   test exists to catch; that is the intended detector, not a separate
#   mutation, so it is documented here rather than re-run.)
# ----------------------------------------------------------------------
