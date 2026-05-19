# van Iseghem 1987 — Vector orthogonal relations, vector QD-algorithm

## Bibliographic record

- **Author:** Jeannette Van Iseghem
- **Title:** *Vector orthogonal relations. Vector QD-algorithm*
- **Venue:** Journal of Computational and Applied Mathematics
  **19** (1), 141–150 (July 1987)
- **DOI:** 10.1016/0377-0427(87)90182-8

## Companion references (not separately fetched)

- C. Brezinski, J. Van Iseghem, *Padé approximations*, in *Handbook
  of Numerical Analysis* vol. III, P.G. Ciarlet & J.L. Lions (eds.),
  Elsevier 1994 — chapter survey covering scalar, vector, and matrix
  Padé approximants.
- J. Van Iseghem, *Approximants de Padé Vectoriels*, PhD thesis,
  Université des Sciences et Techniques de Lille 1, 1987.
- J. Van Iseghem, *Vector Stieltjes continued fraction and vector
  QD algorithm*, Numerical Algorithms **33** (1–4), 485–498 (2003).

## What the paper does

Defines vector Padé approximants to a function `f : ℂ → ℂ^d`,
uniquely, by choosing numerator and denominator degrees (the same
denominator degree shared across all `d` components). The
denominators are vector-orthogonal polynomials `P_{s,r}` satisfying
recurrence relations of order `d + 1`. The paper proves a
generalised MNA algorithm and defines the vector QD-algorithm linking
adjacent diagonal sequences.

**This is the canonical 1987 source for the "shared-denominator" vector
Padé approximant that v0.2 proposes to use** — the algebraic object
v0.2 needs already exists in van Iseghem's setting, with a recurrence
algorithm and an orthogonality interpretation.

## What the paper does **not** do

It stays purely on the approximation-theory side: given a power series,
compute the vector Padé approximant. There is **no** application to
numerical ODE solving, **no** step control, **no** pole-bridging across
movable singularities of an analytic continuation, and **no** Painlevé /
isomonodromic-deformation framing.

Brezinski–Van Iseghem 1994 (the *Handbook of Numerical Analysis*
chapter) similarly stays inside Padé-approximation theory: scalar →
vector → matrix Padé approximants, convergence acceleration, and
biorthogonality, without lifting to ODE solving.

## Gap-2 closure verdict (van Iseghem / Brezinski school)

**Does not close v0.2's gap.** Van Iseghem provides the rigorously-defined
*shared-denominator vector Padé* object the v0.2 design proposes to
use, and the vector QD-algorithm is a candidate alternative backend to
the GGT2013-style SVD route, but the *composition with ODE stepping
and pole-bridging is absent from this body of work*. The "first
explicit vector Padé–Taylor for ODEs" claim survives the van Iseghem
sweep; v0.2 must cite van Iseghem 1987 (and Brezinski–van Iseghem 1994)
as the algebraic foundation it lifts.

## Acquisition status

`acquisition_status: blocked_by_paywall`

- Elsevier paywalled (https://www.sciencedirect.com/science/article/pii/0377042787901828).
- Not on author preprint pages (Brezinski's Lille page does not host
  van Iseghem's 1987 paper as a preprint).

For v0.2 ADR purposes, this metadata stub plus the survey treatment in
the López Lagomasino 2019 arXiv survey (now in this cluster) is
sufficient context.
