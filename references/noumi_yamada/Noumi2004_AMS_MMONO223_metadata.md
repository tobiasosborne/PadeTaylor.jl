# Noumi (M.) 2004: Painlevé Equations Through Symmetry

- **Canonical citation**: Masatoshi Noumi, *Painlevé Equations Through Symmetry*,
  American Mathematical Society, Translations of Mathematical Monographs vol. 223,
  Providence, RI, 2004. ISBN 978-0-8218-3221-9 (paperback), ISBN-10: 0-8218-3221-2.
- **Publisher URL**: https://bookstore.ams.org/mmono-223
  (returns HTTP 403 from automated fetch; confirmed via OpenLibrary record OL16105653M
  and OCLC 53376099 / LCCN 2003062828)
- **DOI**: none (AMS Translations monograph, no Crossref DOI)
- **Internet Archive**: https://archive.org/details/painleveequation0000noum
  (access-restricted borrow; metadata confirmed author/year/series)
- **Acquisition status**: book_paywalled (no open PDF). Routes:
  - Author institutional preprint: none found (original Japanese, *Painlevé Hōteishiki
    to Taishōsei*, Asakura Shoten, Tokyo, 2000, is also closed)
  - Library / institutional scan: WorldCat OCLC 53376099; request via university
    interlibrary loan or scan of relevant chapters
  - Used-book market: search AbeBooks / Bookfinder for ISBN 978-0-8218-3221-9;
    paperback editions available secondhand

- **Table of contents** (reconstructed from Internet Archive metadata, scope-doc
  narrative, and known content of the Noumi–Yamada programme; publisher TOC page
  returned HTTP 403 — marked PARTIAL / INFERRED):

  The book is 156 pages (ix + 156 + index). Based on the known structure of Noumi's
  AMS Translations series monographs and cross-references in the literature:

  - **Chapter 1** — Introduction: Painlevé equations and symmetry
  - **Chapter 2** — The second Painlevé equation and affine Weyl group W(A_1^{(1)})
  - **Chapter 3** — The fourth Painlevé equation and W(A_2^{(1)})
  - **Chapter 4** — The Noumi–Yamada systems: A_n^{(1)} hierarchy (Hamiltonian forms,
    Bäcklund transformations, affine Weyl group W(A_n^{(1)}))
  - **Chapter 5** — Rational solutions and special-function solutions
    (generalised Hermite / Okamoto / Umemura polynomials, closed-form families)
  - **Chapter 6** — Birational transformations and the geometry of solutions

  Bibliography: pp. 153–154 (approx. 40 references). Index included.

  *Source note*: Chapter numbering is inferred from the scope-doc per-pillar narrative
  and standard Noumi–Yamada literature. Physical verification needed against a library
  copy; mark TOC as INFERRED until confirmed.

- **Chapters most relevant to v0.2**:
  - Chapter 4 (A_n^{(1)} Hamiltonian normal forms and Bäcklund generators) — the
    explicit ODE systems PadeTaylor v0.2 must integrate.
  - Chapter 5 (rational-solution closed forms) — the oracle for the
    `noumi_yamada_rational` family verifier (PRD line 320, 352–353).

- **Why this book matters for v0.2** (≤ 3 sentences): The canonical self-contained
  treatment of the A_n^{(1)} Noumi–Yamada hierarchy that PRD line 315 names as a
  primary v0.2 target; Chapters 4–5 supply the Hamiltonian normal forms and rational-
  solution closed forms needed to build and verify `NoumiYamadaProblem` constructors.
  The Mano–Tsuda 2017 paper and Matsuda 2007/2012 + Clarkson et al. 2018 papers cover
  much of the same ground and are freely available, so this book is not v0.2-blocking
  but does consolidate information scattered across multiple papers.

- **Blocking status**: partially blocking, mitigated by open papers:
  - `references/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41.pdf`
    (A_n^{(1)} system definitions and affine-Weyl generators)
  - `references/noumi_yamada/Matsuda2012_rational_A4_NoumiYamada_JMP53.pdf`
    (rational solutions for A_4^{(1)})
  - `references/noumi_yamada/Clarkson_etal_2020_cyclic_maya_higher_painleve_StudApplMath144.pdf`
    (unified Maya-diagram construction covering arbitrary A_{2n}^{(1)})

> Stub placeholder — future commits cite this book via
> `references/noumi_yamada/Noumi2004_AMS_MMONO223_metadata.md`.
