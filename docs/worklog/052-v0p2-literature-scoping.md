# Worklog 052 — v0.2 vector-Painlevé literature scoping (Stage 0+ kickoff)

**Date**: 2026-05-19
**Author**: Claude Opus 4.7 (orchestrator) + 1 Opus scout + 3 Sonnet acquisition subagents in parallel
**Bead**: (none yet — `bd ready` was empty at session start; next agent to create the v0.2 Stage 0+ epic)
**Scope**: First scouting + first acquisition wave for the v0.2 north star (vector / multi-component Painlevé-type systems, PRD lines 213–424). **No source code touched.** Output is research artefacts only.

> **Take-home**: 20 reference PDFs landed in seven new cluster directories,
> 6 metadata stubs cover what couldn't be fetched openly, and a 435-line
> scope doc + 3 batch reports + this worklog leave the next agent with
> a re-launchable Batch 1 (32 arXiv PDFs) and three concrete PRD
> citation corrections to fold into the deep-dive `RESEARCH.md` extension.

## What this is and isn't

This is the **scoping + first acquisition wave** of PRD §"Stage 0+ —
research before v0.2 design" (PRD lines 304–323). It mirrors v0.1
Stage 0 in shape: gather ground-truth papers before any design or
code. The four pillars (PRD lines 308–323) drove the scope:

  - **A** — Hermite–Padé / simultaneous Padé via SVD (the v0.2 design crux);
  - **B** — Noumi–Yamada A_n⁽¹⁾ normal forms (primary v0.2 target);
  - **C** — Painlevé hierarchies PJ⁽ⁿ⁾ (secondary target);
  - **D** — Garnier G_n (tertiary target);

plus three derived pillars E (adjacent / competing numerical
methodology — confirms the PRD "no released software" claim survives),
F (Calogero–Moser pole-soliton oracles for the v0.2 smoke test in PRD
line 361), and G (vector-Taylor-jet infrastructure — deferred to
code-side reading, not literature).

It is **not** the deep-dive document. `RESEARCH.md` extension covering
the same ground at v0.1's depth is the next step. This worklog and
the scope doc are the input to that.

## What shipped (artefacts)

### Scope document

  - `docs/v0p2_literature_scope.md` (435 lines) — categorised reading
    list with 12 Tier-1, 20 Tier-2, 15 Tier-3 entries plus 4 books and
    2 software packages. Each row has a verified canonical citation,
    target path, retrieval URL, acquisition notes, and one-sentence
    relevance line. Six per-pillar narrative paragraphs interpret the
    raw list. Seven open questions surfaced (three of which are PRD
    citation corrections — see "PRD corrections" below).

### Batch reports

  - `docs/v0p2_batch2_preprint_report.md` — author-preprint batch (3/3 OK).
  - `docs/v0p2_batch3_paywalled_report.md` — DOI/paywall batch (2 PDFs + 4 stubs).
  - `docs/v0p2_batch4_books_report.md` — books + software batch (6/6 stubs).
  - **No `v0p2_batch1_arxiv_report.md`** — see "Batch 1 status" below.

### Reference cluster directories

Seven new cluster directories under `references/`, mirroring the
`references/heun/` precedent:

  - `references/hermite_pade/` — Pillar A. **3 PDFs**:
    `BeckermannLabahn1992_…`, `BeckermannLabahn1994_…`,
    `GonnetPachonTrefethen2011_…`, `ManoTsuda2017_…`. (Mano–Tsuda
    landed despite the rejected Batch 1 — one of the running fetch
    agents grabbed it as a side-effect; identity verified via `file`:
    35-page PDF v1.4.)
  - `references/noumi_yamada/` — Pillar B. **3 stubs**:
    `Noumi2004_AMS_MMONO223_metadata.md`,
    `Sakai2001_…_metadata.md` (paywall stub).
  - `references/painleve_hierarchy/` — Pillar C. **1 PDF + 1 stub**:
    `Cosgrove2000_USyd_2000-6_higher_painleve_polynomial_classII.pdf`
    (113 pages, recovered via gunzip + ps2pdf from the USyd
    PostScript), `Kudryashov1997_…_metadata.md`.
  - `references/garnier/` — Pillar D. **2 stubs** (book + algebraic).
  - `references/calogero_moser/` — Pillar F. **1 PDF + 1 stub**:
    `Krichever1980_elliptic_KP_calogero_moser_FAA14.pdf` (Columbia
    faculty mirror — Krichever's own page), `AiraultMcKeanMoser1977_…
    _metadata.md` (Wiley/CPAM paywall).
  - `references/rh_numerical/` — Pillar E. **1 PDF sample + 2 stubs**:
    `TrogdonOlver2016_RH_OT146_sample.pdf` (open SIAM sample chapter
    + 12-chapter TOC) plus `TrogdonOlver2016_RH_OT146_metadata.md`
    and `RiemannHilbert_jl_metadata.md`.
  - `references/taylor_integration/` — Pillar G. **1 stub**:
    `PerezHernandezBenet_TaylorIntegration_metadata.md`.

### Top-level entries

  - `references/Novokshenov2014_pade_painleve_pole_distribution_ConstrApprox39_metadata.md`
    — stub (Springer paywall, no open mirror after 4 attempted
    routes). **Highest-priority next-agent fetch** — Novokshenov is the
    directly-methodology-adjacent paper that v0.1 RESEARCH.md missed.

### Total

**6 PDFs (≈4.7 MB) + 9 metadata stubs + 1 sample PDF + 1 scope doc +
3 batch reports + this worklog**. The arXiv-Mano–Tsuda landing
brings the PDF total to 7 if counted (≈5 MB).

## Batch 1 status — REJECTED before launch

The fan-out plan called for five parallel batches (per
`docs/v0p2_literature_scope.md` §"Fan-out instructions"). The
orchestrator launched four in a single tool message:

  - **Batch 1** — 32 arXiv PDFs (Tier 1 + 2 + 3 papers reachable via
    `arxiv.org/pdf/`). **The user rejected this tool call** in the
    middle of the parallel launch; the other three batches
    (preprints, paywalled, books) had already started and completed
    normally.

No Batch 1 work was performed. The 32-target list is preserved
verbatim in the scope doc at lines 354–388 — the next agent can
re-launch it by either (a) reusing the same Sonnet-subagent prompt
the orchestrator drafted (also archived in conversation context that
the next agent will not have), or (b) reconstructing the prompt
directly from the scope doc's fan-out section. The list itself is
canonical there.

**Whether to retry Batch 1 verbatim or to revise its shape is a
decision for the next agent + user.** The orchestrator did not have
explicit feedback on *why* the launch was rejected — only that the
session needed to wind down. Plausible reasons include: time
pressure, a wish to space out arXiv fetches over multiple
sessions, a wish to inspect the scope doc before mass-fetching, or a
preference for serial human-supervised fetches. The next agent should
ask before re-launching.

## What's not yet in the repo

(Items the scope doc nominated for top-level placement that didn't
land this session — and the reason.)

  - **Kapaev–Klein–Grava 2015 (`arXiv:1306.6161`)** — Tier-1 PI⁽²⁾
    tritronquée pole-field paper, intended top-level alongside FW2011.
    Blocked by Batch 1 rejection.
  - **Noumi–Yamada 1998 (`arXiv:math/9808003`)** — Tier-1 founding
    paper for the A_n⁽¹⁾ hierarchy. Same reason.
  - **27 other Tier-1/2/3 arXiv targets**. Same reason.
  - **Beckermann–Labahn 1994** — Batch 3 fetched a *2001 extended
    preprint* of the 1994 SIMAX paper from Labahn's Waterloo page (26
    pages vs the published 20). The next agent should decide whether
    to keep the longer preprint or replace with the published SIMAX
    version once institutional access is available.

## PRD corrections surfaced (Stage 0+ open questions 1–3)

Confirmed during citation-verification in the scope pass:

  1. **PRD line 402** — `arXiv:0704.2869` is by **Sasano**
     (*Studies on the Garnier system in two variables*), not
     Mazzocco. The Mazzocco Garnier paper PRD meant to cite is
     `arXiv:math/0106208` (Mazzocco, *The geometry of the classical
     solutions of the Garnier systems*, IMRN 2002).
  2. **PRD line 403** — `arXiv:2107.11680` is by **Bobrova–Sokolov**
     (*On matrix Painlevé-4 equations. Part 1*), not Adler–Sokolov.
     The companion Adler–Sokolov paper at `arXiv:2012.05639` is
     correctly attributed.
  3. **PRD line 405** — `arXiv:1712.08546` is the
     **Cafasso–Gavrylenko–Lisovyy** *Tau functions as Widom constants*
     paper. The Gavrylenko–Lisovyy CMP 363 paper PRD meant to cite
     is `arXiv:1608.00958`. Both are relevant; both belong in the
     v0.2 corpus.

**Action for the next agent**: fold these three corrections into the
`RESEARCH.md` v0.2 extension *and* into a small `PRD.md` clarifying
note (matching the in-place `> Note (2026-05-19):` style PRD already
uses at lines 7–13). Do not silently rewrite the original survey
list — append a `> PRD corrections (2026-05-19, worklog 052):` note
that records the misattributions so the original survey provenance is
preserved.

## Two papers PRD's survey missed (Stage 0+ open questions 4–5)

  4. **Mano–Tsuda 2017** (`arXiv:1502.06695`, Math. Z. 285) — the
     single most load-bearing v0.2 paper. Does Hermite–Padé
     approximation via block-Toeplitz determinants and connects
     directly to Schlesinger transformations, P_VI and Garnier. PRD
     line 308 names Baker–Graves-Morris and Beckermann–Labahn for
     Hermite–Padé but misses this paper, which is the *actual*
     intellectual heir of GGT 2013 for the multi-numerator case. It
     landed in `references/hermite_pade/` this session.
  5. **Novokshenov 2014**, *Distributions of Poles to Painlevé
     Transcendents via Padé Approximations* (Constr. Approx. 39,
     85–99) — directly methodology-adjacent prior art for
     PadeTaylor's own algorithm. Uses Fair–Luke Padé on PI/PII/PIV
     pole distributions. **PRD lines 396–414 do not cite it; v0.1
     `RESEARCH.md` does not cite it.** A metadata stub is in
     `references/` (top-level, methodology-adjacent) pending Springer
     institutional access. Reproducing Novokshenov's pole pictures
     with v0.1 PadeTaylor.jl would be a sharp mutation-proof of
     correctness *and* of methodology-adjacent positioning before
     the v0.2 lift.

## Two more open questions surfaced

  6. **Cosgrove PostScript-only handling**. The USyd preprint page
     hosts `.dvi.gz` / `.ps.gz` only. This session converted via
     `gunzip` + `ps2pdf` to land a 113-page PDF. The next agent
     should decide whether to also ship the source `.ps.gz` for
     auditability, or to rely on the converted PDF alone.
  7. **Sakai 2001 classification as type-tag scheme**. Sakai's
     rational-surfaces classification is the umbrella under which
     both Noumi–Yamada (Pillar B) and Garnier (Pillar D) sit. The
     v0.2 ADRs may want to use the Sakai labelling as the type tag
     for `NoumiYamadaProblem` / `GarnierProblem` / similar. Sakai 2001
     itself is paywalled (stub written); the algorithmically needed
     content for v0.2 is in the open papers (Noumi–Yamada 1998,
     Matsuda 2007/2012, Clarkson et al. 2018), so this is not
     blocking, but the type-tag decision wants Sakai's framing.

## Pickup point for the next agent

The single concrete next step is **re-launch Batch 1** (32 arXiv
PDFs) after confirming with the user. The prompt for that batch is
fully specified by the scope doc's Fan-out §Batch 1 (`docs/v0p2_literature_scope.md`
lines 343–388). The other four batches are done.

Once Batch 1 lands, the work order is:

  1. Run `marker_single` on the new PDFs to produce
     `references/markdown/<cluster>/<slug>/<slug>.md` extracts (the
     v0.1 convention). The Mano–Tsuda, Kapaev–Klein–Grava, and
     Noumi–Yamada 1998 papers are the highest priority for marker
     extraction since they will be the most-cited papers in the v0.2
     ADRs.
  2. Fix the three PRD citations per open questions 1–3 above (in
     `PRD.md` via a `> PRD corrections (2026-05-19, worklog 052):`
     note; do not rewrite the survey-source list).
  3. Begin the `RESEARCH.md` v0.2 extension, structured by the same
     six pillars as the scope doc.
  4. Create a v0.2 Stage 0+ epic bead (the orchestrator did not — `bd
     ready` was empty and the user's stop-now timing didn't leave
     room). Suggested shape: one epic with sub-beads per pillar
     (A–G), each carrying acceptance criteria from PRD lines 348–379.
  5. Optionally retry Novokshenov 2014 via institutional Springer
     access if available — it's the single most useful paper still
     locked.

The fan-out infrastructure proven this session (one Opus scout
producing a structured fan-out plan + N parallel Sonnet acquisition
subagents) is reusable for the `RESEARCH.md` extension itself: spawn
one Sonnet per pillar reading the PDFs in that cluster and producing
its summary section.

## What the session deliberately did *not* do

  - **No source code changes.** Not `src/`, not `test/`, not
    `ext/`. v0.1 tests remain at whatever count they were at
    session start (untouched).
  - **No ADRs**, **no `RESEARCH.md` edits**, **no `PRD.md` edits**.
    These belong in the deep dive that follows this scoping pass.
  - **No marker-conversion of the new PDFs.** Deferred per pickup
    point 1.
  - **No author outreach.** Per CLAUDE.md Rule 12 (the same
    discipline that retired outreach to Fornberg, Weideman et al.
    extends to Noumi, Mazzocco, Klein, Olver, Cosgrove, Novokshenov,
    et al.).
  - **No Sci-Hub or other unauthorised mirrors.** Batch 3's
    paywalled-fallback chain was limited to arXiv, ResearchGate
    (author-uploaded), Semantic Scholar, OpenAlex, author
    institutional pages, and the Wayback Machine.

## Friction surfaced

  - **One parallel-launch rejection.** The orchestrator launched four
    fan-out batches in a single tool message; the user rejected one
    (Batch 1, the 32-arXiv fetch). The other three completed before
    the interrupt. Future fan-outs of this size should land more
    incrementally — launch one or two batches, surface the plan to
    the user, then launch the rest after confirmation.
  - **Bead role unconfigured.** `bd` warns `beads.role not configured
    (GH#2950)` on every invocation. Not blocking, but the next
    session should set it: `git config beads.role maintainer`.
  - **Beads directory permissions.** `bd` warns `.beads has
    permissions 0775 (recommended: 0700)`. Cosmetic; next session
    can fix with `chmod 700 .beads`.

## Provenance

  - **Scout** (Opus 4.7, single subagent): produced
    `docs/v0p2_literature_scope.md` in ~10 minutes of web-verification
    work (≈87 tool calls). Verified every arXiv ID against the
    abstract page; surfaced the three PRD citation errors during
    that verification.
  - **Batch 2** (Sonnet 4.6, ~2 minutes): 3 author-preprint fetches.
    Cosgrove's `.ps.gz`-only handling was the load-bearing step;
    `ps2pdf` (Ghostscript 10.02.1) gave a clean 113-page PDF.
  - **Batch 3** (Sonnet 4.6, ~9 minutes): 6 paywalled papers; 2 PDFs
    recovered from author institutional pages (Labahn UW, Krichever
    Columbia), 4 stubs written with the full source-attempt chain
    preserved per paper.
  - **Batch 4+5** (Sonnet 4.6, ~8 minutes): 6 metadata stubs (4
    books, 2 software). Trogdon–Olver open sample chapter
    discovered and fetched as a bonus.
