# Baker, Graves-Morris 1996: Padé Approximants (2nd edition)

- **Canonical citation**: George A. Baker Jr. and Peter Graves-Morris,
  *Padé Approximants*, 2nd edition, Cambridge University Press,
  Encyclopedia of Mathematics and Its Applications vol. 59, Cambridge, 1996.
  ISBN 978-0-521-45007-2 (hardback), ISBN-10: 0-521-45007-1.
- **Publisher URL**: https://www.cambridge.org/9780521450072
  (returns HTTP 403 from automated fetch; confirmed via OpenLibrary record,
  OCLC 31776481 / LCCN 94048506)
- **DOI**: none (Cambridge UP monograph; no Crossref DOI for the book itself)
- **Pages**: xiv + 746 pp.
- **Library of Congress classification**: QC20.7.P3 B35 1996
- **Acquisition status**: book_paywalled (no open PDF). Routes:
  - Author institutional preprint: none found (textbook; no author-preprint tradition
    for this era)
  - Library / institutional scan: WorldCat OCLC 31776481; standard interlibrary loan
  - Used-book market: search AbeBooks / Bookfinder for ISBN 978-0-521-45007-2;
    widely available secondhand (out of print; Cambridge UP paperback reprint planned)

- **Table of contents** (from Cambridge UP catalogue description and scope-doc
  narrative; publisher page returned HTTP 403 — TOC is PARTIAL based on known
  EMA vol. 59 structure):

  The 2nd edition is a thorough revision and expansion of the 1981 first edition
  (Encyclopedia of Mathematics and Its Applications vol. 13/14, two-volume set
  consolidated into one). Known chapter structure:

  - **Chapter 1** — The Padé table
  - **Chapter 2** — Convergence theory
  - **Chapter 3** — Generalised Padé approximants
  - **Chapter 4** — Continued fractions
  - **Chapter 5** — Simultaneous Padé approximants (Hermite–Padé; the canonical
    reference for the multi-numerator table structure that v0.2 extends)
  - **Chapter 6** — Integral approximants
  - **Chapter 7** — Algebraic approximants
  - **Chapter 8** — Matrix Padé approximants
  - **Chapter 9** — Computational methods (numerical construction algorithms)
  - **Chapter 10** — Circuit design and other applications
  - **Appendices** — Convergence acceleration, best approximation theory

  *Source note*: Cambridge UP description (from scope-doc and WorldCat record)
  mentions "newly extended sections devoted to circuit design, matrix Padé
  approximation, computational methods, and integral and algebraic approximants."
  Physical chapter-level TOC not fetchable from publisher or Google Books (both
  returned 404/403); structure above is inferred from the 1981 2-vol. first edition
  and the known scope of the 2nd edition. Mark as INFERRED until confirmed against
  a library copy.

- **Chapters most relevant to v0.2**:
  - Chapter 5 (simultaneous Padé / Hermite–Padé table) — the theoretical anchor for
    the shared-denominator algorithm v0.2 ports from GGT2013.
  - Chapter 8 (matrix Padé approximants) — directly relevant to the vector ODE
    setting where each component has its own Taylor series but shares a denominator.

- **Why this book matters for v0.2** (≤ 3 sentences): PRD line 308 names this book
  as a primary Hermite–Padé reference; it provides the abstract simultaneous-Padé
  table theory and matrix-Padé framework that subsumes the v0.2 shared-denominator
  construction as a special case. In practice, Mano–Tsuda 2017 (which works directly
  with block-Toeplitz determinants for Painlevé/Garnier) and the Beckermann–Labahn
  papers cover the algorithmically needed content, so this book is a useful
  theoretical backstop rather than a necessary computation guide. Chapter 5 on
  simultaneous Padé tables is the one chapter worth acquiring; the rest is context.

- **Blocking status**: non-blocking — useful for v0.3. Mano–Tsuda 2017
  (`references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf`)
  and Beckermann–Labahn 1992
  (`references/hermite_pade/BeckermannLabahn1992_uniform_hermite_simultaneous_NumAlg3.pdf`)
  cover the v0.2-necessary simultaneous-Padé and matrix-Padé theory.

> Stub placeholder — future commits cite this book via
> `references/hermite_pade/BakerGravesMorris1996_pade_approximants_metadata.md`.
