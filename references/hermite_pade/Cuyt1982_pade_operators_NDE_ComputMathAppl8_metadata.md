# Cuyt 1982 — Padé approximants in operator theory for nonlinear DEs

## Bibliographic record

- **Author:** Annie A. M. Cuyt
- **Title:** *Padé-approximants in operator theory for the solution of
  nonlinear differential and integral equations*
- **Venue:** Computers & Mathematics with Applications **8** (6), 445–466
  (1982)
- **DOI:** 10.1016/0898-1221(82)90019-0
- **Related book:** A. Cuyt, *Padé Approximants for Operators: Theory
  and Applications*, Springer Lecture Notes in Mathematics **1065**,
  1984 (extended monograph treatment of the same programme).
- **Related companion:** A. Cuyt, P. Van der Cruyssen, *Abstract Padé
  approximants for the solution of a system of nonlinear equations*,
  Computers & Mathematics with Applications **9** (4), 617–624 (1983).

## What the paper does

Cuyt 1982 takes Padé-approximant ideas into the *abstract operator*
setting: given a nonlinear operator equation `F(x) = 0` (covering ODE
IVPs / BVPs, PDEs, and integral equations as instances), construct
rational-function variants of Newton / Chebyshev / Halley iteration in
the operator setting. The Halley scheme in particular is the headline
result and is shown to behave well in the neighbourhood of
singularities of the underlying problem.

## What the paper does **not** do

It does **not** present a step-by-step Padé–Taylor numerical ODE
integrator with movable-pole step control of the Fornberg–Weideman /
Jorba–Zou kind, and it does **not** propose Hermite–Padé / simultaneous
Padé with a shared denominator for vector-valued analytic ODE systems.
The "Padé approximants" appearing here are *operator-valued rational
functions* used inside fixed-point iteration; they do not approximate
the solution as an analytic function of `t` with polar singularities,
and there is no notion of bridging movable poles along an integration
path.

## Gap-1 closure verdict (Cuyt school)

**Does not close v0.2's gap.** Cuyt's operator-Padé programme stays in
the operator-iteration regime (Newton/Halley variants for nonlinear
operator equations); van der Cruyssen 1983 extends it to systems but
remains operator-iteration, not analytic-continuation step control.
Cuyt's later *multivariate Padé* work (Numer. Math. 1983, Numer. Math.
1984, BIT 1990, the 2008 *Handbook of Continued Fractions for Special
Functions* with Petersen, Verdonk, Waadeland, and Jones) targets
function-of-several-variables approximation — fully orthogonal to the
v0.2 *single-z, vector-y* setting where the Padé approximant is in one
variable and the *output* is vector-valued.

The composition "(vector Padé in single z, shared Q) + (Jorba–Zou /
FFW 2017 step control) applied to Noumi–Yamada / Painlevé hierarchies /
Garnier" is not in Cuyt's body of work. The "first explicit vector
Padé–Taylor for ODEs" claim survives the Cuyt sweep.

## Acquisition status

`acquisition_status: blocked_by_paywall`

- Elsevier paywalled (https://www.sciencedirect.com/science/article/pii/0898122182900190).
- Not on author preprint pages.

This stub is sufficient context for v0.2 ADRs; the underlying paper is
a context citation, not a load-bearing oracle.
