# Iwasaki, Kimura, Shimomura, Yoshida 1991: From Gauss to Painlevé

- **Canonical citation**: Katsunori Iwasaki, Hironobu Kimura, Shun Shimomura,
  and Masaaki Yoshida, *From Gauss to Painlevé: A Modern Theory of Special
  Functions*, Vieweg, Aspects of Mathematics vol. E16, Braunschweig, 1991.
  ISBN 978-3-528-06355-9 (original hardback), ISBN-10: 3-528-06355-6.
  Reprint (Springer, 2012): ISBN 978-3-322-90163-7.
- **Publisher URL**: https://link.springer.com/book/10.1007/978-3-322-90163-7
  (redirects to Springer auth wall; body not fetchable without institutional login)
  Secondary: https://www.springer.com/book/9783528063498 (also redirects to auth wall)
- **DOI**: 10.1007/978-3-322-90163-7 (Springer reprint DOI; confirmed from redirect chain)
- **Pages**: approx. 347–360 pp. (sources vary between original and reprint)
- **Acquisition status**: book_paywalled (no open PDF). Routes:
  - Author institutional preprint: none found (1991 Vieweg monograph; no arXiv era)
  - Library / institutional scan: interlibrary loan via WorldCat (search ISBN
    978-3-528-06355-9 or the Springer reprint 978-3-322-90163-7)
  - Used-book market: AbeBooks listings confirmed, ~$35–$255 USD (original hardback
    and Springer softcover reprint both available); ISBN confirmed via AbeBooks scrape
    2026-05-19

- **Table of contents** (reconstructed from scope-doc narrative, standard references
  in the Garnier / isomonodromy literature, and the book's known subject matter;
  publisher page body not accessible — TOC marked INFERRED):

  The book covers the classical theory of special functions (Gauss hypergeometric
  through Fuchsian ODEs) and its modern continuation to the Painlevé equations and
  the Garnier system via isomonodromic deformations. Known structure:

  - **Introduction** — From Gauss's hypergeometric equation to the Painlevé equations
  - **Chapter 1** — Gauss's hypergeometric equation (classical theory; Fuchsian ODEs,
    monodromy, Riemann's P-function)
  - **Chapter 2** — Fuchsian equations and isomonodromic deformations (Schlesinger
    equations, Painlevé VI as isomonodromy condition)
  - **Chapter 3** — The Painlevé equations (classical derivation, Painlevé property,
    Bäcklund transformations for PI–PVI)
  - **Chapter 4** — The Garnier system (multi-variable generalisation; monodromy-
    preserving deformations of Fuchsian systems in two or more variables)
  - **Chapter 5** — Special solutions of the Garnier system (classical solutions,
    algebraic solutions, relation to hypergeometric functions)
  - **Appendices** — Monodromy data, connection formulae, supplementary tables

  *Source note*: Springer DOI page returned auth-redirect (303 See Other) with no
  readable body. Chapter structure above is inferred from the book's mathematical
  content as cited in the Garnier literature (Mazzocco, Sasano, Calligaris–Mazzocco
  papers all cite chapters 4–5 specifically). Physical verification required.

- **Chapters most relevant to v0.2**:
  - Chapter 4 (Garnier system) — foundational reference for the Garnier-tertiary
    target; supplies the Hamiltonian structure and monodromy-preserving-deformation
    formulation that v0.2 must encode in `GarnierProblem`.
  - Chapter 5 (special solutions of Garnier) — the classical-solution and algebraic-
    solution landscape that determines what oracle comparisons are available for v0.2
    Garnier acceptance tests.

- **Why this book matters for v0.2** (≤ 3 sentences): The canonical textbook
  reference for the Garnier system (multi-variable isomonodromic deformations), which
  is the tertiary v0.2 target listed in PRD §D; Chapters 4–5 define the Hamiltonian
  and classical-solution structures that the Mazzocco and Sasano papers cite as
  background. In practice, the Garnier target is non-blocking for v0.2-primary
  (Noumi–Yamada A_n^{(1)} and PI-hierarchy are higher priority), so the book is
  useful pre-reading rather than a computation-day dependency. Acquiring one chapter
  scan via interlibrary loan is sufficient for v0.2; the full book is desirable for v0.3.

- **Blocking status**: non-blocking for v0.2 (Garnier is tertiary). Useful for v0.3.
  Partially mitigated by:
  - `references/garnier/Mazzocco2002_geometry_garnier_classical_IMRN_math0106208.pdf`
    (Riemann–Hilbert framework for Garnier classical solutions)
  - `references/garnier/Sasano2007_garnier_two_variables_0704.2869.pdf`
    (Hamiltonian forms for 2-variable Garnier)
  - `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf`
    (explicitly treats the Garnier system in the Hermite–Padé context)

> Stub placeholder — future commits cite this book via
> `references/garnier/IwasakiKimuraShimomuraYoshida1991_from_gauss_to_painleve_metadata.md`.
