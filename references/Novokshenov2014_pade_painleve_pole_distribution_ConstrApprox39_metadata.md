# Novokshenov 2014: Distributions of Poles to Painlevé Transcendents via Padé Approximations

- **Canonical citation**: V. Yu. Novokshenov (2014), *Distributions of Poles to Painlevé Transcendents via Padé Approximations*, Constructive Approximation 39 (2014) 85–99.
- **DOI**: 10.1007/s00365-013-9190-6
- **Publisher URL**: https://link.springer.com/article/10.1007/s00365-013-9190-6
- **Acquisition status**: blocked_by_paywall
- **Acquisition routes tried**:
  - arXiv search (arxiv.org, multiple queries on title, author, year 2013/2014): no arXiv preprint found; the paper was not posted to arXiv.
  - ResearchGate (https://www.researchgate.net/publication/257406070): page exists and lists a PDF icon, but HTTP 403 on direct fetch — ResearchGate blocks automated access; no anonymous download URL available.
  - Semantic Scholar API (DOI lookup): status `CLOSED`, no openAccessPdf URL.
  - Author page (Institute of Mathematics, Ufa Scientific Center, Russian Academy of Sciences; mathnet.ru search for Novokshenov): mathnet.ru search URL returned 404; no accessible preprint found.
  - Wayback Machine: web.archive.org blocked by WebFetch; not accessible during this batch.
  - Google Scholar / general web: multiple queries returned only the Springer page and the ResearchGate listing (403); no freely downloadable PDF found after 4+ open-source attempts.
- **Why this paper matters for v0.2**: Flagged as the highest-priority find in `docs/v0p2_literature_scope.md` (Tier 1, row 12, Pillar E,A). Uses a Fair–Luke Padé algorithm to compute pole distributions of PI/PII/PIV transcendents in the complex plane — the *closest prior art* to PadeTaylor.jl's own algorithm. v0.2 docs should cite this as prior art; v0.2 should reproduce Novokshenov's pole-distribution figures as a mutation-proof cross-validation step. Scope doc notes: "Mutation-proof opportunity: reproduce Novokshenov's pole pictures with PadeTaylor and confirm agreement, then lift the same method to Painlevé hierarchy members."
- **Mitigation if unobtainable**: Scope doc (Tier 1 row 12) names no partial substitute; the Padé-for-Painlevé methodology is the scope doc's primary motivation for acquiring this paper. v0.2 work proceeds (algorithm is already in-tree); revisit when institutional Springer access is available. The related earlier paper (Novokshenov 2009, Theoretical and Mathematical Physics 159, DOI 10.1007/s11232-009-0073-8) covers PI and PII with Padé approximations and may partially substitute if this paper remains unavailable.

> This stub is a placeholder so that future commits can cite the paper by `references/Novokshenov2014_pade_painleve_pole_distribution_ConstrApprox39_metadata.md` even though the PDF is not in-tree.
