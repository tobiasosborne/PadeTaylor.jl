# v0.2 Literature Acquisition — Batch 4+5 Report

Agent: Batch 4+5 subagent (Claude Sonnet 4.6), 2026-05-19.

## Stub status: 6/6 OK

| # | Stub | Path | Lines | Status | Summary |
|---|---|---|---|---|---|
| 1 | Noumi 2004 (book) | `references/noumi_yamada/Noumi2004_AMS_MMONO223_metadata.md` | 65 | OK | ISBN 978-0-8218-3221-9 confirmed via OpenLibrary OCLC 53376099; AMS publisher page returned HTTP 403; TOC inferred from scope-doc narrative and known Noumi–Yamada literature (marked INFERRED); Internet Archive record confirmed 156 pp. |
| 2 | Baker–Graves-Morris 1996 (book) | `references/hermite_pade/BakerGravesMorris1996_pade_approximants_metadata.md` | 71 | OK | ISBN 978-0-521-45007-2 confirmed via OpenLibrary (OCLC 31776481, 746 pp.); Cambridge UP page returned HTTP 403; chapter structure inferred from first edition and publisher description (marked INFERRED); blocking status: non-blocking. |
| 3 | Iwasaki–Kimura–Shimomura–Yoshida 1991 (book) | `references/garnier/IwasakiKimuraShimomuraYoshida1991_from_gauss_to_painleve_metadata.md` | 74 | OK | ISBN 978-3-528-06355-9 (original) and 978-3-322-90163-7 (Springer reprint) confirmed via AbeBooks scrape; Springer DOI page redirected to auth wall; chapter structure inferred from Garnier literature citations (marked INFERRED); blocking status: non-blocking for v0.2. |
| 4 | Trogdon–Olver 2016 (book) | `references/rh_numerical/TrogdonOlver2016_RH_OT146_metadata.md` | 96 | OK | ISBN 978-1-61197-419-5, DOI 10.1137/1.9781611974201 confirmed from SIAM epubs page (HTTP 200); complete 12-chapter TOC extracted from open sample PDF via pdftotext; sample PDF fetched and saved. |
| 5 | TaylorIntegration.jl | `references/taylor_integration/PerezHernandezBenet_TaylorIntegration_metadata.md` | 56 | OK | Canonical repo confirmed as `PerezHz/TaylorIntegration.jl` (scope doc had wrong org names); v0.18.4 current (2026-05-02); concept DOI 10.5281/zenodo.2562352; MIT licence; last commit 2026-05-17; actively maintained (139 stars). |
| 6 | RiemannHilbert.jl | `references/rh_numerical/RiemannHilbert_jl_metadata.md` | 70 | OK | Canonical repo confirmed as `JuliaHolomorphic/RiemannHilbert.jl` (scope doc listed wrong org `JuliaApproximation`); v0.1.0 only release (2019-12-09); last commit 2024-12-19; PII-only scope confirmed from README (Hastings–McLeod example only); dormant status confirmed. |

## Trogdon–Olver sample chapter

**Found and fetched.** The URL `https://www.math.uci.edu/~ttrogdon/publications/TrogdonOlver-Sample.pdf`
is live (701.7 KB PDF, confirmed HTTP 200 on 2026-05-19). It contains the book's
front matter, complete table of contents, preface, and full Chapter 1 (Classical
applications of Riemann–Hilbert problems, pp. 1–21). The file has been saved to:

```
references/rh_numerical/TrogdonOlver2016_RH_OT146_sample.pdf
```

The TOC in the stub (`TrogdonOlver2016_RH_OT146_metadata.md`) was extracted directly
from this PDF via `pdftotext` and is therefore primary-source verified (not inferred).

## Key findings and corrections

1. **TaylorIntegration.jl repo URL corrected.** The scope doc listed
   `PerezHernandez93/TaylorIntegration.jl` and `JuliaDiff/TaylorIntegration.jl` as
   candidate URLs — both return HTTP 404. The canonical URL is
   `https://github.com/PerezHz/TaylorIntegration.jl`, confirmed from the Zenodo record
   DOI 10.5281/zenodo.2562353 (v0.4.1).

2. **RiemannHilbert.jl org corrected.** The scope doc and PRD list
   `JuliaApproximation` as the GitHub org; the canonical org is `JuliaHolomorphic`.
   GitHub API query for `JuliaApproximation/RiemannHilbert.jl` returns no results;
   `JuliaHolomorphic/RiemannHilbert.jl` is the live repo.

3. **Zenodo concept DOI clarification.** The scope doc cites Zenodo DOI
   `10.5281/zenodo.2562353` (which is the version DOI for v0.4.1). The concept DOI
   covering all versions is `10.5281/zenodo.2562352` (visible in the README badge and
   confirmed via the Zenodo record). Both are valid; use the concept DOI in citable
   references.

4. **Three book TOCs are INFERRED.** Publisher pages for the AMS, Cambridge UP, and
   Springer/Vieweg books all returned HTTP 403 or auth redirects. Google Books also
   returned 404. The TOCs in stubs 1–3 are reconstructed from: (a) known structure
   of the respective series, (b) scope-doc per-pillar narrative, and (c) literature
   citations. Each stub is explicitly labelled INFERRED and should be verified against
   a library copy before being cited at chapter resolution.

5. **Trogdon–Olver TOC is primary-source verified** (pdftotext from sample PDF).

## File inventory

```
references/noumi_yamada/Noumi2004_AMS_MMONO223_metadata.md       (new)
references/hermite_pade/BakerGravesMorris1996_pade_approximants_metadata.md (new)
references/garnier/IwasakiKimuraShimomuraYoshida1991_from_gauss_to_painleve_metadata.md (new)
references/rh_numerical/TrogdonOlver2016_RH_OT146_metadata.md   (new)
references/rh_numerical/TrogdonOlver2016_RH_OT146_sample.pdf    (new, 702 KB)
references/taylor_integration/PerezHernandezBenet_TaylorIntegration_metadata.md (new)
references/rh_numerical/RiemannHilbert_jl_metadata.md           (new)
```
