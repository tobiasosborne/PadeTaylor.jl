# Worklog 053 — v0.2 Stage 0+ second wave: Batch 1 retry + prior-art gap sweep

**Date**: 2026-05-19 (same calendar day as worklog 052; continuation session after a brief gap)
**Author**: Claude Opus 4.7 (orchestrator, serialised — one subagent at a time per user instruction) + 1 Sonnet acquisition agent + 1 Opus prior-art-sweep agent.
**Bead**: (still none — to be created when the v0.2 epic lands)
**Scope**: (a) **Resume Batch 1**, the 32-arXiv-PDF fetch rejected in worklog 052. (b) Run a **second-wave prior-art sweep** targeting the five gaps that surfaced during the user-orchestrator discussion after the first wave (Cuyt school, van Iseghem, Aptekarev / Stahl school, Klein-Bothner-Lisovyy-Joshi 2024–2026, GitHub code search).

> **Take-home**: Batch 1 lands cleanly (32/32 PDFs, 27.85 MB, zero
> failures). The prior-art sweep finds **no prior art that scoops v0.2's
> composition** ("shared-denominator Hermite–Padé + adaptive Taylor
> stepping + Painlevé / Garnier / Noumi–Yamada targeting"), but
> surfaces **two important qualifications**: Amore 2021 (*Physica D*)
> already uses the exact name "Padé–Taylor method" for a scalar van
> der Pol integrator, and Adler–Sokolov 2025 (*arXiv:2512.18828*,
> J. Geom. Phys.) is a brand-new vector-Painlevé *classification*
> paper that reinforces the v0.2 motivation without scooping it.

## What shipped (artefacts)

### Batch 1 — 32 arXiv PDFs (27.85 MB)

  - `docs/v0p2_batch1_arxiv_report.md` — fetch report.
  - 32 PDFs distributed across `references/{hermite_pade,
    noumi_yamada, painleve_hierarchy, garnier, rh_numerical}/` per the
    scope doc's batch-1 manifest. Headline targets (per worklog 052
    §"What's not yet in the repo"):
      - `references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf`
        — the v0.2 PI⁽²⁾ acceptance-test reference figure source
        (PRD line 374–375).
      - `references/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41.pdf`
        — the founding paper for A_n⁽¹⁾.
      - `references/hermite_pade/NakatsukasaSeteTrefethen2018_AAA_SISC40.pdf`
        — the Trefethen-school successor to GGT 2013.

  Fetched serially with 3-second arXiv-etiquette gaps. Total wall-clock
  ~229 s. All 32 verified: size > 50 KB, `%PDF-` magic, page count > 5.
  Pre-2007 slash-form arXiv IDs (`math/9808003`, `math/0106208`,
  `math/0212117`, `math-ph/9810007`, `q-alg/9708018`, `nlin/0306020`)
  all resolved at `https://arxiv.org/pdf/<full-ID>.pdf` with no
  fallback needed.

### Prior-art second-wave sweep — 6 PDFs + 5 stubs

  - `docs/v0p2_prior_art_sweep.md` — sweep report (~22 KB), one section
    per gap with searches, findings, and a verdict.
  - 4 Aptekarev-line convergence-theorem PDFs in `references/hermite_pade/`:
    `LopezLagomasino2019_intro_multiple_orthogonal_hermite_pade_…pdf`
    (open-access pedagogical entry), `FidalgoLagomasino_Medina_2013_HP_meromorphic_…pdf`,
    `LagomasinoMedina2014_HP_typeI_meromorphic_…pdf` (rational-modification
    case — closer to v0.2's actual setting than the pure Markov / Nikishin
    setup), and `IkonomovSuetin2026_convergence_rational_HP_…pdf` (most
    recent in the line).
  - **`references/rh_numerical/AdlerSokolov2025_vector_painleve_type_2512.18828.pdf`**
    — the 2025 *J. Geom. Phys.* vector Painlevé classification paper.
    Pure Lax-pair construction, no numerics, no code; **strengthens
    v0.2 motivation by independently confirming "the classification
    exists; no numerical software does."** This should be cited in
    PRD's §"v0.2 north star" §"Survey sources" as a 2025 reinforcement.
  - **`references/Amore2021_pade_taylor_van_der_pol_PhysicaD_arXiv2111.12198.pdf`**
    — Top-level placement (methodology-adjacent). **Name collision:**
    Amore 2021 explicitly calls his algorithm "the Padé–Taylor method
    (PTM)" for a scalar van der Pol integrator. Construction is Taylor
    jet + diagonal Padé + residual-step. **v0.1's docs should
    acknowledge Amore 2021 as parallel scalar prior art for the
    Padé–Taylor methodology name** — this is a v0.1 doc-debt item, not
    a v0.2 one, but it surfaces during this scoping.
  - 5 metadata stubs under `references/hermite_pade/`: Cabay–Jones–Labahn
    1997 (Calgo 766 — the FORTRAN vector-Padé oracle that v0.2's
    shared-Q port can mutation-prove against, exactly the way
    `padeapprox.m` was used for the scalar GGT 2013 port in v0.1),
    Cuyt 1982 (closest Cuyt title; verified non-ODE-stepping by
    inspection), Van Iseghem 1987 JCAM (algebraic foundation for
    shared-Q with one denominator across `d` components, recurrences
    of order `d + 1`), Aptekarev–Stahl 1992 (the load-bearing
    convergence-theorem chapter; paywalled book chapter — the
    citable open-access alternative is López Lagomasino 2019), and
    the Amore 2021 metadata companion.

## The five-gap verdict, condensed

  - **Gap 1 — Cuyt school.** Surveyed Cuyt 1982 / 1983 / 1984 / 1990,
    *Handbook of Continued Fractions* (2008), MPF software, and the
    Cuyt–Van der Cruyssen 1983 nonlinear-systems companion. Cuyt's
    *Padé in operator theory* programme is abstract operator iteration
    (Halley-method generalisation), not step-by-step Padé–Taylor ODE
    integration with movable-pole step control. Multivariate Padé is
    several-variable approximation, orthogonal to v0.2's
    single-z / vector-y setting. **Gap survives.**
  - **Gap 2 — van Iseghem / Brezinski Lille line.** Van Iseghem 1987
    is the *direct canonical* algebraic source for shared-denominator
    vector Padé — recurrences of order `d + 1`, vector QD-algorithm as
    a backend alternative to SVD-based GGT 2013. But no application
    to numerical ODE solving in any of the searched papers (1987 JCAM,
    2003 *Vector Stieltjes*, Brezinski–Van Iseghem 1994 Handbook
    chapter). **Gap survives**; v0.2 ADR for shared-Q must cite Van
    Iseghem 1987.
  - **Gap 3 — Aptekarev / Stahl school.** The convergence-theorem chain
    v0.2 needs is well-established: Aptekarev–Stahl 1992 chapter for
    the historical formulation, López Lagomasino 2019 as the
    open-access pedagogical re-proof, Fidalgo–Lagomasino–Medina 2013
    and Lagomasino–Medina 2014 for the rational-modification case
    (Nikishin system + finitely many added poles — the closest setup to
    "meromorphic ODE solution approximated by Hermite–Padé"), and
    Ikonomov–Suetin 2026 as the most recent. **None applies to
    numerical ODE solving.** Convergence guarantee for shared-Q
    Hermite–Padé of a system of meromorphic functions on compact
    subsets of the complement of the contour system is the citation
    v0.2 ADR needs. **Gap survives.**
  - **Gap 4 — Klein / Bothner / Lisovyy / Joshi 2024–2026.** All four
    groups remain active but no 2024–2026 paper is a vector / multi-
    component Painlevé numerical method, and **none releases code**.
    Klein's Painlevé-numerics line (Grava–Klein 2011, Kapaev–Klein–
    Grava 2015, Klein–Stoilov 2018) appears dormant; Klein's 2024–2026
    output is on fractional NLS / fractional Schrödinger. Bothner 2024–
    2026 is RH-asymptotics theoretical. Lisovyy 2024–2026 is tau-function
    / conformal-block theoretical. Joshi 2024–2026 is geometric /
    arithmetic / q-difference. **The surprise** was Adler–Sokolov 2025
    (`arXiv:2512.18828`) — not from any of the four groups asked
    about, but a 2025 *Vector systems of Painlevé type* classification
    via group reduction of vector NLS / mKdV / KdV, producing ODE
    systems with isomonodromic Lax pairs generalising PI, PII, P34, PIV.
    **Pure classification, no numerics, no code — confirms v0.2's
    target class is alive and unaddressed.** **Gap survives.**
  - **Gap 5 — GitHub code search.** Seven queries via `gh search code`
    and `gh search repos`. The most plausibly competing scalar
    Painlevé numerical code is `cool-japan/scirs`'s
    `scirs2-special/src/painleve` (RK45 PI/PII only, no Padé, no
    pole-field). The most plausibly competing *adjacent* Hermite–Padé
    code is `ore_algebra` (Mark Kauers / SageMath; scalar Hermite–Padé
    for *guessing* D-finite ODEs from data — opposite direction:
    recovering an ODE from terms, not solving an ODE forward through
    poles). The named `JuliaHolomorphic/Painleve.jl` is a 2019-created
    stub, no releases, no commits since 2021-01-04, 0 stars. **Gap
    survives — and the survey-supporting evidence for the PRD claim
    "no production-quality, registered, vector- or multi-component-
    Painlevé numerical package" is now hardened.**

## Honest framing implication for v0.2 documentation

The first-wave scoping (worklog 052) recommended cluster directories
and surfaced 3 PRD citation errors + 2 missed keystones. The
second-wave sweep adds two further refinements to the v0.2 framing:

  1. **The v0.1 docs (not v0.2) should acknowledge Amore 2021 as
     parallel scalar prior art** for the "Padé–Taylor method" name.
     Amore's construction is scalar, per-component, no FW / GGT
     robustification, no path networks, no Painlevé framing — so v0.1
     PadeTaylor.jl's contribution is unambiguous — but the name
     collision is real and v0.1's own provenance prose deserves a
     one-paragraph "the term has also been used by Amore 2021 for a
     scalar van der Pol integrator" note.

  2. **The v0.2 framing in PRD §"v0.2 north star" should be sharpened
     from "first general-purpose released software for any of these
     systems" to one of two more precise alternatives:**
       - "first packaged general-purpose released software, and first
         complex-plane pole-field maps of any higher-Painlevé-
         hierarchy member or any Noumi–Yamada A_n⁽¹⁾ solution beyond
         n = 2, 3"
       - "the first software to combine shared-denominator Hermite–Padé
         (Mahler / Coates / Beckermann–Labahn / Van Iseghem algebraic
         line; Aptekarev / Stahl convergence theory) with FW-2011
         path-network pole-bridging and applied to vector Painlevé-type
         systems"
     The first is sharper for the headline; the second is sharper for
     the technical ADRs. Both are *true*; the bare "first general
     method statement" claim from the post-first-wave discussion is
     **not** what should ship — Mano–Tsuda 2017, Beckermann–Labahn
     1992/1994, Van Iseghem 1987, and the Aptekarev / Stahl line each
     supply pieces of the analytic statement; v0.2's contribution is
     the *composition* and the *applied numerical method + pictures*.

## Total v0.2 corpus state (as of this commit)

| Cluster | PDFs | Stubs |
|---|---|---|
| top-level (methodology-adjacent / v0.1-relevant) | 2 (Trogdon–Olver sample chapter is in `rh_numerical/`; Amore 2021 PTM here) | 2 (Amore 2021 companion, Novokshenov 2014) |
| `references/hermite_pade/` | 9 | 5 |
| `references/noumi_yamada/` | 8 | 2 |
| `references/painleve_hierarchy/` | 8 | 1 |
| `references/garnier/` | 6 | 1 |
| `references/calogero_moser/` | 1 | 1 |
| `references/rh_numerical/` | 12 | 2 |
| `references/taylor_integration/` | 0 | 1 |
| **Total v0.2** | **46 PDFs** | **15 stubs** |

Plus 5 reports under `docs/v0p2_*.md` (scope, batch 1–4 fetch reports,
prior-art sweep) and 2 worklog entries (052 + this one).

## Pickup point for the next agent

The literature acquisition is **substantially complete**. The next
agent should focus on the *consumption* side:

  1. **Marker-convert priority papers** to
     `references/markdown/<cluster>/<slug>/<slug>.md`. Highest priority
     for marker conversion: Mano–Tsuda 2017, Kapaev–Klein–Grava 2015,
     Noumi–Yamada 1998, Cosgrove 2000-6, López Lagomasino 2019,
     Adler–Sokolov 2025, Amore 2021. (Worklog 052's pickup point 1.)
  2. **Fix the three PRD citations** identified in worklog 052 §"PRD
     corrections" (lines 402, 403, 405) via a `> PRD corrections
     (2026-05-19, worklog 052+053):` append-note that preserves the
     original survey provenance. **Also add Adler–Sokolov 2025
     reinforcement note** and Amore 2021 acknowledgement for v0.1.
  3. **Begin the `RESEARCH.md` v0.2 extension.** Structure by the same
     seven pillars as the scope doc. Suggested fan-out: one subagent
     per pillar, *one at a time* (the user's preference, see this
     worklog header). Each subagent reads the PDFs in its cluster
     directory and produces its `RESEARCH.md` section.
  4. **Create the v0.2 Stage 0+ epic bead** before writing
     `RESEARCH.md` (so deep-dive work is tracked from the start).
     Suggested structure: one epic with seven sub-beads (A–G) carrying
     PRD lines 348–379 acceptance criteria.

The headline analytic-underpinnings synthesis that surfaced during
the orchestrator-user discussion (the "shared-Q Hermite–Padé converges
to a Painlevé-property vector solution on pole-free compacts" theorem)
is **not in the literature in this exact form**; it is a corollary of
(Painlevé property) + (joint Laurent structure) + (Aptekarev / Stahl
convergence) that v0.2 should state explicitly in its ADRs. Writing
this corollary cleanly is genuine v0.2 mathematical work, half a
dozen pages of synthesis — Mano–Tsuda 2017 is the natural anchor and
the Aptekarev / Stahl line is the technical engine.

## What this session deliberately did *not* do

  - **No source code changes.** Not `src/`, not `test/`, not `ext/`.
    v0.1 tests remain at whatever count they were at session start
    (untouched).
  - **No marker conversion of the new PDFs.** Deferred per pickup
    point 1.
  - **No ADRs, no `RESEARCH.md` edits, no `PRD.md` edits.** Deferred
    per pickup points 2 + 3 + 4.
  - **No `MEMORY.md` files** (per CLAUDE.md / SessionStart hook —
    beads is the project memory).
  - **No author outreach.** Per CLAUDE.md Rule 12.
  - **No `bd` commands.** Bead creation deferred to next agent. `bd
    ready` showed empty at session start of 052; remains empty at end
    of 053. The `.beads has permissions 0775 (recommended: 0700)`
    warning and `beads.role not configured` warning from `bd`
    continue; both are non-blocking and noted in worklog 052.

## Friction surfaced (delta from worklog 052)

  - **Serialised subagents work fine and are cheaper to reason about.**
    The user constrained this session to one subagent active at a
    time after the worklog-052 4-parallel-launch was interrupted.
    With sequential subagents, each result lands cleanly with no
    coordination cost. This is a good default for literature work
    where the throughput is dominated by external network calls (arXiv
    fetches) rather than by orchestrator-side reasoning.
  - **The 32-vs-31 off-by-one in the Batch 1 brief**. The orchestrator
    told the subagent "31 net targets" because Mano–Tsuda was
    pre-fetched, but the list itself contained 32 entries (the
    orchestrator hadn't subtracted Mano–Tsuda). The subagent
    nevertheless correctly fetched the list and verified all 32 (`Fetched:
    32/31`). No harm done; Mano–Tsuda was not in the list because the
    orchestrator excluded it earlier. Trivial process bug, surfaced
    here for completeness.

## Provenance

  - **Batch 1 retry** (Sonnet 4.6, ~7 minutes wall-clock, ~14 tool
    calls): serial arXiv fetches with 3-second gaps. 32 PDFs, zero
    failures, total 27.85 MB.
  - **Prior-art second-wave sweep** (Opus 4.7, ~15 minutes wall-clock,
    ~99 tool calls): five-gap sweep + 6 PDF fetches + 5 stub writes
    + sweep-report drafting. The Opus-tier model was used here because
    the work required judgement on whether each candidate paper
    closes a gap — Sonnet would have over-claimed surface matches as
    closures. The Adler–Sokolov 2025 surprise was discovered via a
    Klein-author-page-adjacent search that returned the arXiv listing
    around late 2025; it had to be read carefully to confirm the
    "classification, no numerics" status before being filed.
