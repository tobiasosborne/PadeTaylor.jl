# v0.2 Literature Acquisition — Batch 2 Report

> Generated: 2026-05-19
> Agent: Sonnet 4.6 (Batch 2 — author preprint pages + format conversion)

---

## 1. Cosgrove 2000-6 — Higher-order Painlevé equations in polynomial class II (Bureau symbol P1)

**Status: OK**

**Landing page:** `https://ssh.maths.usyd.edu.au/u/ResearchReports/Nonlinear/Cos/2000-6.html`

The page offered four download options: `2000-6.dvi.gz` (142 kB), `2000-6.dvi` (427 kB), `2000-6.ps.gz` (264 kB), `2000-6.ps` (774 kB). No PDF was provided on the page.

**URL fetched:** `https://ssh.maths.usyd.edu.au/u/ResearchReports/Nonlinear/Cos/2000-6.ps.gz`

**Conversion pipeline:**
1. `curl -fL https://ssh.maths.usyd.edu.au/u/ResearchReports/Nonlinear/Cos/2000-6.ps.gz -o /tmp/Cosgrove2000-6.ps.gz`
   — Downloaded 265 kB; confirmed as `gzip compressed data, last modified: Mon Feb 21 22:03:59 2000`.
2. `gunzip -k /tmp/Cosgrove2000-6.ps.gz`
   — Decompressed to 774 kB PostScript; confirmed as `PostScript document text conforming DSC level 2.0`.
3. `ps2pdf /tmp/Cosgrove2000-6.ps /tmp/Cosgrove2000-6.pdf`
   — Converted using ps2pdf (Ghostscript 10.02.1 backend). No errors.
4. Copied to target path.

**File size:** 543 kB (555,850 bytes)

**Page count:** 113 pages (threshold: > 60)

**First bytes:** `%PDF-` confirmed by `file` (reported as `PDF document, version 1.4, 113 page(s)`)

**Target path:** `references/painleve_hierarchy/Cosgrove2000_USyd_2000-6_higher_painleve_polynomial_classII.pdf`

**Note:** The `.ps` was originally produced by `dvips(k) 5.85` (Copyright 1999 Radical Eye Software), consistent with a 2000 TeX/DVI source. The USyd page date is February 2000.

---

## 2. Beckermann–Labahn 1992 — A uniform approach for Hermite–Padé and simultaneous Padé approximants

**Status: OK**

**URL fetched:** `https://cs.uwaterloo.ca/~glabahn/Papers/tenerife.pdf`

**Conversion pipeline:** None required — file is already PDF (direct `curl -fL`).

**File size:** 211 kB (215,244 bytes)

**Page count:** 10 pages (threshold: > 5)

**First bytes:** `%PDF-` confirmed by `file` (reported as `PDF document, version 1.2, 6 page(s)`)

Note: `file` and `pdfinfo` gave slightly different counts (6 vs. 10); `pdfinfo` is authoritative and reports 10 pages.

**Target path:** `references/hermite_pade/BeckermannLabahn1992_uniform_hermite_simultaneous_NumAlg3.pdf`

**Metadata:** Creator: TeX; Producer: pdfTeX-0.13d; PDF version 1.2.

---

## 3. Gonnet–Pachón–Trefethen 2011 — Robust rational interpolation and least-squares

**Status: OK**

**URL fetched:** `https://people.maths.ox.ac.uk/trefethen/publication/PDF/2011_141.pdf`

**Conversion pipeline:** None required — file is already PDF (direct `curl -fL`).

**File size:** 1.6 MB (1,667,471 bytes)

**Page count:** 22 pages (threshold: > 5)

**First bytes:** `%PDF-` confirmed by `file` (reported as `PDF document, version 1.4, 22 page(s)`)

**Target path:** `references/hermite_pade/GonnetPachonTrefethen2011_robust_rational_interp_ETNA38.pdf`

**Metadata:** Creator: LaTeX with hyperref package; Producer: dvips + GPL Ghostscript 8.62; PDF version 1.4.

---

## Summary

| Paper | Status | File size | Pages | Path |
|---|---|---|---|---|
| Cosgrove 2000-6 | OK (PDF via ps2pdf) | 543 kB | 113 | `references/painleve_hierarchy/Cosgrove2000_USyd_2000-6_higher_painleve_polynomial_classII.pdf` |
| Beckermann–Labahn 1992 | OK (direct PDF) | 211 kB | 10 | `references/hermite_pade/BeckermannLabahn1992_uniform_hermite_simultaneous_NumAlg3.pdf` |
| Gonnet–Pachón–Trefethen 2011 | OK (direct PDF) | 1.6 MB | 22 | `references/hermite_pade/GonnetPachonTrefethen2011_robust_rational_interp_ETNA38.pdf` |

All three PDFs passed size (> 50 kB) and `%PDF-` header checks. Cosgrove passed the page-count > 60 threshold with 113 pages.
