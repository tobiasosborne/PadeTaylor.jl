# Cabay, Jones & Labahn 1997 — Algorithm 766: VECTOR_PADE

## Bibliographic record

- **Authors:** S. Cabay, A. R. Jones, G. Labahn
- **Title:** *Algorithm 766: Experiments with a Weakly Stable Algorithm for
  Computing Padé–Hermite and Simultaneous Padé Approximants*
- **Venue:** ACM Transactions on Mathematical Software (TOMS) **23** (1),
  91–110 (March 1997)
- **DOI:** 10.1145/244768.244790
- **Algorithm number:** ACM Calgo 766
- **Code:** FORTRAN 90, distributed by ACM under `0766/` in the Calgo
  archive (https://calgo.acm.org/). Mirrored on GitHub by third parties
  at `ACM-TOMS/CALGO`, `Beliavsky/Alan-Miller-Fortran/VEC_PADE.F90`,
  `iajzenszmi/CodeCode/766.sh`, and `iajzenszmi/Decent-Code-1/766.sh`.

## Companion theory paper (not separately fetched; same authors)

- S. Cabay, A. R. Jones, G. Labahn, *A Fast, Reliable Algorithm for
  Calculating Padé–Hermite Forms* (precursor / theoretical companion;
  Cabay–Labahn ISSAC 1989) and Beckermann–Cabay–Labahn, *Fraction-Free
  Computation of Matrix Padé Systems* (ISSAC 1997). These are listed on
  Labahn's Waterloo publication page; not fetched here because the
  Calgo 766 paper alone is the load-bearing reference.

## Relevance to PadeTaylor v0.2

This is the closest existing prior art for the *core linear-algebra
substrate* of v0.2's shared-denominator design: a published, weakly
stable, FORTRAN-implemented algorithm for both Padé–Hermite and
simultaneous Padé approximants along the diagonal of the Padé table.

It is **not** an ODE solver:

- It computes the approximants `(P_1, …, P_k, Q)` (with `Q` the
  shared denominator) from a given vector of formal power series.
- It does **not** apply step control, residual-based stepping, pole
  detection on the analytic-continuation path, or pole-bridging.
- It does **not** target Painlevé / Garnier / Noumi–Yamada systems
  specifically; the test functions in §6 are generic.

Its role in v0.2 is exactly the role GGT 2013 played in v0.1: a
mutation-proof oracle for the rational-approximant output. The v0.2
shared-Q port should reproduce Calgo 766's `VECTOR_PADE` output (or
its Beckermann–Labahn-1994 successor at higher orders) to bit-level
agreement on common test inputs.

## Gap closure verdict

**Does not close v0.2's "first explicit vector Padé–Taylor for ODEs"
claim.** Calgo 766 is an approximation-theory routine; v0.2 is an ODE
integrator. The two compose: v0.2 *could* use a Calgo 766 backend if
its FORTRAN ABI were wrapped, but neither this paper nor the FORTRAN
package ever made that composition explicit, much less applied it to
Painlevé pole-field mapping. The combination remains the v0.2
contribution.

## Acquisition status

`acquisition_status: blocked_by_paywall`

- ACM Digital Library: paywalled (institutional access via
  https://dl.acm.org/doi/10.1145/244768.244790).
- Author preprint: not found on Labahn's Waterloo publication page
  (only the companion 1989/1992 papers are openly preprinted there).
- FORTRAN source is openly redistributable via the Calgo archive at
  https://calgo.acm.org/ (the journal paper accompanying it is not).

Fan-out action: if institutional library access is available, fetch
the PDF into this directory replacing this stub. Until then, the
stub plus the freely available FORTRAN source are sufficient to
mutation-prove against.
