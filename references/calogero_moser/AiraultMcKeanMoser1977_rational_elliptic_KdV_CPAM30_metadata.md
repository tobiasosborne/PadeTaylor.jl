# Airault, McKean, Moser 1977: Rational and Elliptic Solutions of the KdV Equation and a Related Many-Body Problem

- **Canonical citation**: H. Airault, H. P. McKean, and J. Moser (1977), *Rational and elliptic solutions of the Korteweg–de Vries equation and a related many-body problem*, Communications on Pure and Applied Mathematics 30 (1977) 95–148.
- **DOI**: 10.1002/cpa.3160300106
- **Publisher URL**: https://onlinelibrary.wiley.com/doi/10.1002/cpa.3160300106
- **Acquisition status**: blocked_by_paywall
- **Acquisition routes tried**:
  - arXiv search (multiple queries by author, year, title fragments): no arXiv preprint exists; 1977 predates arXiv by 15 years and the paper was never retroactively posted.
  - ResearchGate (https://www.researchgate.net/publication/238225545 — Krichever 1980 page confirmed reachable in search, AMM page also listed): HTTP 403 on all ResearchGate fetch attempts; automated access blocked.
  - Semantic Scholar API (DOI lookup): status `CLOSED`, no openAccessPdf URL.
  - Author pages: NYU/Courant Institute archives (McKean, Moser); no public preprint server existed in 1977. AMS memorial for Moser (2000) confirms the paper but does not host a PDF. No accessible NYU or Courant preprint page found.
  - NASA ADS (https://ui.adsabs.harvard.edu/abs/1977CPAM...30...95A): record exists but links to Wiley paywalled version only.
  - Wayback Machine: web.archive.org blocked by WebFetch; not accessible during this batch.
- **Why this paper matters for v0.2**: Tier 2, row 31 of `docs/v0p2_literature_scope.md` (Pillar F — Calogero–Moser pole dynamics). The founding paper for the Calogero–Moser ↔ KdV pole-soliton correspondence; gives the closed-form rational and elliptic solutions that the v0.2 Calogero–Moser N-particle smoke test (PRD line 361) cross-validates against. Without the closed-form expressions from this paper, the smoke test lacks a ground-truth oracle.
- **Mitigation if unobtainable**: Scope doc notes Pillar F needs AMM 1977 + Krichever 1980 + Wilson 1993 for the CM smoke test. Krichever 1980 is now in-tree (`Krichever1980_elliptic_KP_calogero_moser_FAA14.pdf`). The rational-solutions part of AMM is also described in Adler–Moser (1978) *J. Math. Phys.* 19, which covers the same polynomial families, and in Satsuma–Ablowitz (1979); these do not require CPAM access. v0.2 smoke test can proceed with Krichever's elliptic case plus ADS abstract data for AMM rational case; revisit when institutional Wiley access is available.

> This stub is a placeholder so that future commits can cite the paper by `references/calogero_moser/AiraultMcKeanMoser1977_rational_elliptic_KdV_CPAM30_metadata.md` even though the PDF is not in-tree.
