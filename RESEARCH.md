# RESEARCH.md — Stage 0 deep research for PadeTaylor

> **Status: in progress.** Sections below are scaffolds; each will be
> filled with citations to primary literature in `references/`. Per the
> PRD's discipline ("ground truth via local PDFs"), every claim in this
> document carries a `references/<file>.pdf` citation, not a paraphrase
> from memory.

## 0. Scope of this document

Stage 0 of the project plan in `PRD.md`. Output of this document
informs `DESIGN.md` (Stage 1) and gates the move to implementation.
Resolves or sharpens every provisional claim in `PRD.md` against the
literature and the ecosystem. Ends with a **go/no-go recommendation**
for Stage 1.

## 1. Algorithmic core — the Fornberg–Weideman line

### 1.1 FW 2011 (J. Comput. Phys. 230, 5957–5973)

- Reference: `references/<TBD>.pdf`
- Algorithmic spec (§3–4):
- Step control:
- Pole-immune stepping mechanism:
- Boundaries / failure modes the paper documents:

### 1.2 FW 2014 (FoCM 14, 985–1016)

- Reference: `references/<TBD>.pdf`
- Refinements over FW 2011:
- PII survey results (figures we must reproduce qualitatively per PRD §"Definition of done"):

### 1.3 Fasondini–Fornberg–Weideman 2017 (J. Comput. Phys. 344, 36–50)

- Reference: `references/<TBD>.pdf`
- Riemann-surface extension — relevance to v1 vs v2:

### 1.4 Reeger–Fornberg PIV survey

- Reference(s): `references/<TBD>.pdf`
- Findings relevant to algorithm choices:

### 1.5 Willers 1974 — original IVP-with-poles

- Reference: `references/<TBD>.pdf`
- What FW generalises from Willers; what Willers' constraints were:

## 2. Robust Padé approximation

### 2.1 Gonnet–Güttel–Trefethen 2013 (SIAM Review 55)

- Reference: `references/<TBD>.pdf`
- The SVD-based algorithm (in detail; this is the canonical routine):
- Tolerance choices and their tradeoffs:
- Failure modes that GGT robustifies against (Froissart doublets, etc.):

### 2.2 Chebfun `padeapprox.m`

- Reference impl of GGT 2013. Source path / commit:
- Notes on deviations from the paper:

### 2.3 Hermite–Padé (Baker & Graves-Morris)

- Reference: `references/<TBD>.pdf`
- v1 vs v2 implications (shared denominator, branch-point detection):

## 3. Adjacent Taylor-IVP packages — design lessons

### 3.1 TIDES (Barrio et al.)

### 3.2 TaylorIntegration.jl (Pérez-Hernández, Benet)

- Real-axis only, no Padé. What does it get right / wrong?

### 3.3 TaylorSeries.jl (Benet, Sanders)

- **Critical question**: generic-type support. Does it work with `Arb`
  element type at order 60+ without breakage? Empirical test required
  before depending on it.

### 3.4 ATOMFT (Chang & Corliss) — historical reference

### 3.5 Jorba & Zou 2005 — design lessons for Taylor IVP packages

## 4. Adjacent verified / arb-prec ODE work

### 4.1 Mezzarobba — D-finite / linear ODEs in Arb (`ore_algebra`)

- Different scope but useful patterns for verified arithmetic in Arb.

### 4.2 COSY Infinity / Taylor models (Makino & Berz)

- Verified-bound thinking. Relevant if v2 pursues Arb-ball enclosures.

## 5. Ecosystem inventory

### 5.1 `@workbench/cas-core` (Scientist Workbench, TS)

- Taylor / FPS class: present? Absent?
- `Poly<T>` over `Field<T>` exists; can it host a "formal power series" use?
- Bigfloat SVD: **does not currently exist in the workbench substrate.**
  - `linalg-svd` is `Float64` numerical-tier (ADR-0014/0015/0016).
  - `@workbench/bigfloat` provides the `BigFloat` / `BigComplex` arb-prec
    scalars (per ADR-0020) but no matrix factorisation built on top.
  - Implication: porting a bigfloat SVD (Demmel-Kahan or one-sided
    Jacobi) is a sub-task of the TS impl, comparable in size to the
    Padé routine itself.

### 5.2 Julia ecosystem

- `Arblib.jl`: does `ArbMatrix` support SVD natively? Or fall back to
  generic `LinearAlgebra` SVD over `Arb` element type?
- `TaylorSeries.jl` with `Arb` at order 60+: hands-on test required.
- `Polynomials.jl`: roots of `Q` for step control; precision behaviour
  with `Arb` coefficients?

## 6. Source-code hunt

### 6.1 Paper appendices

- FW 2011 appendix code excerpts (if present):
- FW 2014 appendix:

### 6.2 Author pages and arXiv ancillary files

- Weideman Stellenbosch page (`appliedmaths.sun.ac.za/~weideman/`):
- Fornberg CU Boulder page:
- arXiv ancillary files for the FW papers:

### 6.3 Direct outreach (async; does not block Stage 1)

- Drafted emails to Weideman / Fornberg / Fasondini / Reeger:
- Status of replies:

## 7. Resolutions of `PRD.md` open questions

For each provisional claim or open question in `PRD.md`, this section
states the literature-backed resolution.

### 7.1 PRD §"Stage 1 — Design decisions to lock"

1. **Padé form**: scalar-per-component vs Hermite–Padé with shared
   denominator. **Recommendation**: TBD
2. **Step control**: coefficient-decay vs Padé-denominator-root distance
   vs both switchable. **Recommendation**: TBD
3. **Order strategy**: fixed (30–60) vs adaptive. **Recommendation**: TBD
4. **Path strategy in v1**: real-axis only with `step_complex` primitive
   vs 2D network from day one. **Recommendation**: TBD
5. **API surface**: standalone clean core + thin `CommonSolve` wrapper.
   **Recommendation**: TBD
6. **Verified arithmetic via Arb balls**: v1 / v2 / never. **Recommendation**: TBD
7. **Module boundaries** under ~200 LOC budget. **Recommendation**: TBD

### 7.2 PRD §"Open questions"

1. Hermite–Padé in v1 or v2: TBD
2. Step control default: TBD
3. Verified enclosures: TBD
4. Path-network consistency metric: TBD
5. Naming (`PadeTaylor.jl` vs alternatives): TBD
6. Licence: TBD

## 8. New questions surfaced during research

(Populated as research progresses.)

## 9. Recommendation

**Status**: TBD until §1–§7 are filled.

Possible outcomes:
- Proceed to Stage 1 (`DESIGN.md`).
- Identify specific blocking gaps; recommend further research before
  Stage 1.

## Appendix A — references inventory

| short name | full citation | local path | source |
|---|---|---|---|
| FW 2011 | Fornberg & Weideman, *A numerical methodology for the Painlevé equations*, J. Comput. Phys. 230 (2011) 5957–5973. DOI 10.1016/j.jcp.2011.04.007 | `references/FW2011_painleve_methodology_JCP230.pdf` | user-fetch (Cloudflare-blocked direct) |
| GGT 2013 | Gonnet, Güttel & Trefethen, *Robust Padé Approximation via SVD*, SIAM Review 55(1) (2013) 101–117. DOI 10.1137/110853236 | `references/GGT2013_robust_pade_via_SVD_SIREV55.pdf` | user-fetch |
| RF 2014 | Reeger & Fornberg, *Painlevé IV: A numerical study of the fundamental domain and beyond*, Physica D 280–281 (2014) 1–13. DOI 10.1016/j.physd.2014.04.006 | `references/ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.pdf` | user-fetch |
| FW 2015 | Fornberg & Weideman, *A computational overview of the solution space of the imaginary Painlevé II equation*, Physica D 309 (2015) 108–118. DOI 10.1016/j.physd.2015.07.008 | `references/FW2015_imaginary_PII_PhysicaD309.pdf` | user-fetch (bonus; complements PRD's FW 2014 FoCM) |
| FFW 2017 | Fasondini, Fornberg & Weideman, *Methods for the computation of the multivalued Painlevé transcendents on their Riemann surfaces*, J. Comput. Phys. 344 (2017) 36–50. Preprint version. | `references/FFW2017_painleve_riemann_surfaces_preprint.pdf` | user-fetch |

### Outstanding (still to acquire)

| short name | rationale | acquisition status |
|---|---|---|
| FW 2014 (FoCM) | PII survey; PRD-named. FW 2015 (Imag PII) covers similar ground but is not a substitute. Springer FoCM. | not yet attempted |
| Willers 1974 | Originating IVP-with-poles algorithm. Historical baseline. | not yet attempted |
| Chebfun `padeapprox.m` | Reference impl of GGT 2013. | not yet attempted (GitHub fetch should work) |
| TaylorIntegration.jl / TaylorSeries.jl source | Ecosystem audit (Phase 0.3). | clone via git, no auth needed |
| TIDES papers (Barrio et al.) | Mature Taylor-IVP design lessons. | not yet attempted |
| Jorba & Zou 2005 | Taylor IVP package paper, design lessons. | not yet attempted |
| Mezzarobba — `ore_algebra` / Arb papers | Verified-arithmetic patterns for arb-prec ODE. | not yet attempted |
| COSY Infinity / Berz Taylor models | Verified-bound thinking for v2. | low priority for Stage 0 |
| Baker & Graves-Morris (Hermite-Padé) | v2 / branch-point detection. | book; low priority for Stage 0 |

