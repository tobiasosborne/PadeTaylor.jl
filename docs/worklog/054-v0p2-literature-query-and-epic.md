# Worklog 054 — v0.2 literature query, TeX-source fetch, and epic decomposition

**Date**: 2026-05-20
**Author**: Claude Opus 4.7 (orchestrator) + 4 Sonnet literature-query subagents + 1 Sonnet TeX-fetch subagent
**Bead**: `padetaylor-0ln` (EPIC: v0.2 — vector / multi-component Painlevé-type systems) — created this session, with child `padetaylor-0ln.1` (V0) covering this consolidation.
**Scope**: Stage 0+ second deep pass. The literature acquired in worklogs 052–053 had landed but never been *read*; this session queries it, fetches LaTeX source for every arXiv paper (a side quest replacing the slow `marker` step on this device), and decomposes v0.2 into a 20-bead epic. **No source code touched** — research + planning artefacts only.

> **Take-home**: Four pillar findings docs now carry the actionable
> algorithmic ground truth (the shared-Q stacked-Toeplitz SVD; the
> explicit A_n⁽¹⁾ system and its A_2⁽¹⁾⇒PIV reduction; the explicit
> P_I⁽²⁾ equation and its 4-vector form; the Calogero–Moser exact
> oracle). 48/48 arXiv TeX bundles fetched into `references/tex/`.
> The v0.2 epic `padetaylor-0ln` (+19 children, 29 dependency edges)
> is the executable plan.

## What shipped (artefacts)

### Four pillar findings docs

A read-only Sonnet subagent per load-bearing pillar, each writing a
path/page-cited findings doc (Law 1) and reading the PDFs directly:

  - `docs/v0p2_pillarA_hermite_pade_findings.md` — **shared-denominator
    robust Padé**. The algorithm: stack `d` component Toeplitz blocks
    into one `dm×(m+1)` matrix, take one SVD, the smallest right
    singular vector is the shared denominator `Q` (`‖Q‖₂=1`); the
    `padeapprox.m` QR-reweighting ports unchanged. `d=1` reduces
    exactly to GGT 2013 `:svd`. Load-bearing reference: Mano–Tsuda
    2017 §2.2 (block-Toeplitz determinant). Convergence underpinning:
    Fidalgo–López-Lagomasino–Medina 2013 Thm 1.3. Independent oracle:
    ACM Calgo 766 (Cabay–Jones–Labahn vector Padé, FORTRAN).
  - `docs/v0p2_pillarB_noumi_yamada_findings.md` — **Noumi–Yamada
    A_n⁽¹⁾**. Explicit system (even-`l`: `f_j' = f_j(Σf_odd − Σf_even)
    + α_j`; odd-`l` a quadratic variant), `Σf_j` and `Σα_j` constraints.
    A_2⁽¹⁾⇒scalar PIV reduction pinned (`y=−f_1/c`, `c=√(−3/2)`).
    Recommended n≥4 figure: A_4⁽¹⁾, `α_j=1/5`. Exact oracle: A_4⁽¹⁾
    "Type C" rational `f_j(t)=t/5`. W(A_n⁽¹⁾) self-test is algebraic
    (`s_i²=id`, `π^{n+1}=id`).
  - `docs/v0p2_pillarC_painleve_hierarchy_findings.md` — **PI
    hierarchy**. `P_I⁽²⁾`: `u_xxxx + 10u_x² + 20u·u_xx + 40(u³−6tu+6x)
    = 0` (KKG 2015 eq. 1.1 = Cosgrove F-V). 4-component first-order
    companion form transcribed. Tritronquée IC: seed from the KKG
    asymptotic series at `x₀=−20, t=0`. Acceptance figure: KKG 2015
    Figs 7.4/7.5.
  - `docs/v0p2_pillarEF_methodology_calogero_findings.md` — **adjacent
    methodology + Calogero–Moser**. PadeTaylor is *complementary* to
    Novokshenov 2014 (static Padé scatter, no stepper) and Klein–
    Stoilov 2018 (sector evaluation, no pole-field maps) — the PRD
    "no competing software" claim survives. CM smoke test: N=2
    rational KdV, exact pole trajectories `x_{1,2}(t)=±½√(1+4t²)`.

### `references/tex/` — 48 arXiv LaTeX source bundles (side quest)

`marker` PDF→markdown conversion is too slow on this device. A Sonnet
subagent fetched the arXiv e-print source (`arxiv.org/e-print/<ID>`)
for every arXiv-hosted paper in `references/`: **48/48 succeeded**,
62 `.tex` files, 40 MB, mirrored cluster structure under
`references/tex/<cluster>/<paper-slug>/`. v0.2 ground-truth reasoning
reads the `.tex` directly — no marker pass needed for any v0.2 paper.
Journal-only papers without an arXiv ID (FW2011, GGT2013, …) were
correctly skipped.

### v0.2 epic — `padetaylor-0ln` + 19 children

See `docs/v0p2_plan.md` for the full 10-phase plan and the
bead↔phase map. The epic decomposes into V0–V9 with 29 dependency
edges; V1 (shared-Q robust Padé) is the keystone that blocks the
vector lift.

## Decisions taken this session

  1. **Three oracles for shared-Q** (user, 2026-05-20: "More = better").
     V1's cross-validation is triple-redundant — AAA (V1d) + Calgo 766
     via `ccall` (V1b) + an independent Julia block-Toeplitz determinant
     construction (V1c) — all feeding the V1e mutation-proof suite.
  2. **Additive architecture.** The vector lift ships *new* modules
     (`SharedPade.jl`, `VectorCoefficients.jl`, `VectorStepper.jl`,
     `VectorProblems.jl`) alongside untouched v0.1 scalar modules.
     Rationale: backward compat is non-negotiable (v0.1 scalar tests
     must stay bit-identical), and the v0.1 Float64 default `:classical`
     must not be perturbed by routing through the SVD-based shared-Q.
     Also keeps every file under the Rule 6 200-LOC cap.
  3. **Garnier stays out of scope.** Per PRD, Garnier is multi-`z`
     (a PDE / multi-time problem) — a different project. The
     `references/garnier/` cluster is retained as context only.
  4. **Calogero–Moser smoke test is in.** Not a PRD acceptance item,
     but a cheap end-to-end vector validation with a genuine exact
     oracle — included as V4.

## Lessons

**42. Read-only literature subagents parallelise; two of five hit a
transient socket error.** Pillars B and E/F died mid-run with
"socket connection closed unexpectedly" — an infrastructure fault,
not a task failure. Re-dispatched fresh (now pointed at the faster
`references/tex/` sources) they completed in ~5 min each. Lesson:
a subagent API error is re-dispatchable; check whether the deliverable
file was written before assuming lost work.

**43. Fetch TeX source, not just PDFs, when marker is slow.** The
arXiv e-print endpoint serves the LaTeX source bundle for nearly every
paper (48/48 here, including pre-2007 slash-form IDs). Reading `.tex`
is faster and more exact for equations than either marker output or
PDF page reads — worth doing as a standing convention for new
reference acquisitions.
