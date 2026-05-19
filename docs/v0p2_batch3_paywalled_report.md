# v0.2 Literature Acquisition — Batch 3 Report (Paywalled Papers)

> Generated 2026-05-19 by Sonnet 4.6 subagent.  
> Scope: 6 paywalled papers from Springer / Wiley / SIAM / Elsevier.  
> Method: arXiv mirror → ResearchGate/Semantic Scholar → author institutional pages → Wayback Machine → metadata stub.

---

## Summary

| # | Paper | Status | Artefact |
|---|---|---|---|
| 1 | Novokshenov 2014, Constr. Approx. 39 | **STUB** | `references/Novokshenov2014_pade_painleve_pole_distribution_ConstrApprox39_metadata.md` |
| 2 | Beckermann–Labahn 1994, SIMAX 15 | **OK** | `references/hermite_pade/BeckermannLabahn1994_fast_matrix_pade_SIMAX15.pdf` |
| 3 | Airault–McKean–Moser 1977, CPAM 30 | **STUB** | `references/calogero_moser/AiraultMcKeanMoser1977_rational_elliptic_KdV_CPAM30_metadata.md` |
| 4 | Krichever 1980, Funct. Anal. Appl. 14 | **OK** | `references/calogero_moser/Krichever1980_elliptic_KP_calogero_moser_FAA14.pdf` |
| 5 | Kudryashov 1997, Phys. Lett. A 224 | **STUB** | `references/painleve_hierarchy/Kudryashov1997_first_second_painleve_higher_order_PLA224_metadata.md` |
| 6 | Sakai 2001, CMP 220 | **STUB** | `references/noumi_yamada/Sakai2001_rational_surfaces_painleve_CMP220_metadata.md` |

---

## Paper 1 — Novokshenov 2014

**Full ref**: V. Yu. Novokshenov, *Distributions of Poles to Painlevé Transcendents via Padé Approximations*, Constructive Approximation 39 (2014) 85–99. DOI: 10.1007/s00365-013-9190-6.

**Status**: STUB

**Acquisition chain**:
1. arXiv (title + author + year, multiple query variants): no arXiv preprint exists; paper not posted to arXiv.
2. Semantic Scholar API (DOI 10.1007/s00365-013-9190-6): `openAccessPdf.status = "CLOSED"`, URL empty.
3. ResearchGate (https://www.researchgate.net/publication/257406070): page exists; HTTP 403 on fetch — ResearchGate blocks automated access.
4. Author page at Institute of Mathematics Ufa, mathnet.ru: mathnet.ru search URL returned 404; no accessible preprint.
5. Wayback Machine: `web.archive.org` not reachable by WebFetch during this session.
6. Four additional keyword+venue searches: all returned only the Springer paywall page and the 403 ResearchGate listing.

**Artefact**: `references/Novokshenov2014_pade_painleve_pole_distribution_ConstrApprox39_metadata.md`

**Note (highest-priority paper)**: This is the single highest-priority target in Batch 3, flagged Tier 1 in the scope doc. The paper uses the Fair–Luke Padé algorithm to compute pole distributions for PI/PII/PIV — the closest prior art to PadeTaylor.jl. Despite 6+ open-source attempts it remains behind the Springer paywall with no author-uploaded copy found. Institutional Springer access (e.g., via a university library subscription) is the recommended next step.

---

## Paper 2 — Beckermann–Labahn 1994

**Full ref**: B. Beckermann and G. Labahn, *A Uniform Approach for the Fast Computation of Matrix-Type Padé Approximants*, SIAM J. Matrix Anal. Appl. 15 (1994) 804–823. DOI: 10.1137/S0895479892230031.

**Status**: OK

**Acquisition chain**:
1. arXiv: no arXiv preprint; predates arXiv math.NA era.
2. Semantic Scholar API (DOI): `openAccessPdf.status = "CLOSED"`.
3. Labahn faculty page at University of Waterloo CS (https://cs.uwaterloo.ca/~glabahn/): directory listing shows `Papers/uniform.pdf` (linked from the "uniform computation of matrix rational approximants" group). Downloaded and verified.
4. PDF verified: `%PDF-1.2`, 26 pages, 306 KB, title topic matches (uniform approach for different concepts of matrix-type Padé approximants including vector, Hermite, simultaneous, and matrix variants). Keywords in abstract match SIMAX 1994 paper. The preprint is dated October 2001 and is an extended / revised version of the SIMAX 1994 paper; it is the closest freely available version on the author's own page.

**Artefact**: `references/hermite_pade/BeckermannLabahn1994_fast_matrix_pade_SIMAX15.pdf` (306 KB, 26 pp, `%PDF-1.2`)

**Note**: The deposited file is a 2001 extended preprint of the 1994 paper from Labahn's own faculty page. It covers the same algorithm (O(N log² N) complexity, divide-and-conquer, Hermite + simultaneous + matrix-type Padé), with the same author attribution. This is the standard academic practice for citing via the author preprint.

---

## Paper 3 — Airault–McKean–Moser 1977

**Full ref**: H. Airault, H. P. McKean, and J. Moser, *Rational and elliptic solutions of the Korteweg–de Vries equation and a related many-body problem*, Communications on Pure and Applied Mathematics 30 (1977) 95–148. DOI: 10.1002/cpa.3160300106.

**Status**: STUB

**Acquisition chain**:
1. arXiv: no preprint; 1977 predates arXiv by 14 years and the paper was never retroactively posted.
2. Semantic Scholar API (DOI 10.1002/cpa.3160300106): `openAccessPdf.status = "CLOSED"`.
3. ResearchGate: HTTP 403 on direct fetch.
4. NASA ADS (https://ui.adsabs.harvard.edu/abs/1977CPAM...30...95A): abstract only; link to Wiley paywall.
5. NYU/Courant Institute archives: no open preprint server existed in 1977; AMS memorial for Moser references the paper but does not host a PDF.
6. Additional web searches for "free pdf download", NYU Courant preprint, open access: no freely accessible PDF found.

**Artefact**: `references/calogero_moser/AiraultMcKeanMoser1977_rational_elliptic_KdV_CPAM30_metadata.md`

---

## Paper 4 — Krichever 1980

**Full ref**: I. M. Krichever, *Elliptic solutions of the Kadomtsev–Petviashvili equation and integrable systems of particles*, Functional Analysis and Its Applications 14 (1980) 282–290. DOI: 10.1007/BF01078304.

**Status**: OK

**Acquisition chain**:
1. arXiv: no preprint; 1980 predates arXiv.
2. Semantic Scholar API (DOI 10.1007/BF01078304): `openAccessPdf.status = "CLOSED"`.
3. ResearchGate: HTTP 403 on direct fetch.
4. mathnet.ru (Russian math portal): search URL returned wrong paper; could not identify correct record.
5. Krichever's Columbia University faculty publications page (https://www.math.columbia.edu/department/krichever/pubs.html): **success** — entry 128 lists both a Russian PDF and English translation PDF. English translation URL: `http://www.math.columbia.edu/~krichev/pdfs/1980-1984/1980-ESOTKSEAISOP0.pdf`.
6. Downloaded and verified: `%PDF-1.3`, 9 pages, 648 KB. PDF metadata title: "Elliptic solutions of the Kadomtsev-Petviashvili equation and integrable systems of particles". Correct paper, correct author, correct page count (9 pp = English translation of 9-page Russian article).

**Artefact**: `references/calogero_moser/Krichever1980_elliptic_KP_calogero_moser_FAA14.pdf` (648 KB, 9 pp, `%PDF-1.3`)

---

## Paper 5 — Kudryashov 1997

**Full ref**: N. A. Kudryashov, *The first and second Painlevé equations of higher order and some relations between them*, Physics Letters A 224 (1997) 353–360. DOI: 10.1016/S0375-9601(97)00850-5.

**Status**: STUB

**Acquisition chain**:
1. arXiv (multiple queries by author, title, year): no arXiv preprint found; Kudryashov did not post 1990s papers to arXiv.
2. Semantic Scholar API (DOI): lookup returned a different Kudryashov 1998 paper; both `CLOSED`.
3. ScienceDirect (Elsevier): publisher page confirmed paywalled.
4. ResearchGate / Academia.edu: searches found related Kudryashov papers (Bäcklund transformations for Painlevé hierarchies, 2022) but not this specific 1997 PLA paper with free download.
5. Author page (MEPhI / Steklov Institute): no open preprint page found for 1990s Kudryashov papers.
6. Wayback Machine: `web.archive.org` not reachable by WebFetch.

**Artefact**: `references/painleve_hierarchy/Kudryashov1997_first_second_painleve_higher_order_PLA224_metadata.md`

**Note**: Scope doc flags this as Tier 3 (non-critical for v0.2); Cosgrove 2000-6 already in-tree covers the same classification territory.

---

## Paper 6 — Sakai 2001

**Full ref**: H. Sakai, *Rational Surfaces Associated with Affine Root Systems and Geometry of the Painlevé Equations*, Communications in Mathematical Physics 220 (2001) 165–229. DOI: 10.1007/s002200100446.

**Status**: STUB

**Acquisition chain**:
1. arXiv (multiple queries): no arXiv preprint; paper not posted to arXiv.
2. Semantic Scholar API (DOI 10.1007/s002200100446): `openAccessPdf.status = "CLOSED"`.
3. Author page at Kobe University (https://www.math.kobe-u.ac.jp/~sakai/): HTTP 404 — page does not exist or has moved.
4. RIMS Kyoto University preprint server: no relevant preprint found via web search.
5. ResearchGate: HTTP 403 on direct fetch.
6. Springer Nature Link: paywalled; no "Shareable Link" available.
7. Wayback Machine: not reachable by WebFetch.

**Artefact**: `references/noumi_yamada/Sakai2001_rational_surfaces_painleve_CMP220_metadata.md`

**Note**: Scope doc (Tier 2) notes the survey `arXiv:1509.08186` (Kajiwara–Noumi–Yamada, "Geometric Aspects of Painlevé Equations") describes the Sakai classification in detail and is open-access; that paper was fetched in Batch 1 and partially mitigates the absence of the original Sakai paper.

---

## Verification Checklist

| Artefact | Size | First bytes | Pages | Pass? |
|---|---|---|---|---|
| `Krichever1980_elliptic_KP_calogero_moser_FAA14.pdf` | 648 KB | `%PDF-1.3` | 9 | YES |
| `BeckermannLabahn1994_fast_matrix_pade_SIMAX15.pdf` | 306 KB | `%PDF-1.2` | 26 | YES |
| 4 metadata stubs | >50 bytes each | Markdown | N/A | YES |
