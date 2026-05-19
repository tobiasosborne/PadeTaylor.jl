# Amore 2021 — "Padé–Taylor method" for the van der Pol equation

## Bibliographic record

- **Author:** Paolo Amore
- **Title:** *Computing the solutions of the van der Pol equation to
  arbitrary precision*
- **Venue:** Physica D **435**, 133263 (2022) (published April 2022);
  preprint arXiv:2111.12198 (23 November 2021)
- **DOI:** 10.1016/j.physd.2022.133263
- **Affiliation:** Universidad de Colima

## Why this paper matters for PadeTaylor.jl

**This is a name-collision and methodology-adjacent paper that v0.1's
RESEARCH.md does not cite.** Amore 2021 names the same algorithm:
*Padé-Taylor method (PTM)*, abbreviated *VSVO-PT* (Variable-Stepsize
Variable-Order Padé-Taylor) (§2 of the paper). The algorithm shape is
also similar:

1. Compute a Taylor jet at the current point (eq. (3) in Amore;
   Jorba–Zou 2005 is cited as ref [4]).
2. Reroute the Taylor polynomial through a Padé approximant
   `P[X(t)]_{[m,n]}` (eq. (4)), to extend the radius of usefulness
   beyond the Taylor radius of convergence.
3. Pick step `τ_0` to keep the residual below tolerance `δ` (eq. (5)).
4. Iterate.

PadeTaylor.jl v0.1's RESEARCH.md and PRD describe a substantially
similar pipeline, citing FW 2011 / FW 2014 / FW 2015 / FFW 2017 / GGT
2013 / JorbaZou 2005 — Amore 2021 is missing from both. Fornberg /
Weideman / Fasondini do not appear in Amore's reference list either.
The two threads developed independently.

## Where Amore 2021 stops short of v0.2

Amore's paper is **scalar-only and PII-irrelevant**:

- Tested only on the van der Pol equation (one scalar second-order
  ODE, rewritten as a 2-dimensional first-order system).
- **No Hermite-Padé / shared-denominator** treatment — Amore applies
  per-component diagonal Padé `[n, n]` to the vector unknown of the
  van der Pol rewrite, which is the *per-component* design call PRD
  explicitly rejects.
- **No Painlevé context.** Cites no FW2011, no GGT2013, no
  Trefethen-school robust-Padé work — Amore's Padé core is naive
  Toeplitz, with spurious-pole avoidance implemented by *stepping
  back from the spurious pole* in §2, not by Froissart-doublet
  detection.
- **No pole-field maps** — the figures in §3 plot the real poles of
  the Padé approximant on the time axis (Fig. 3) to illustrate the
  spurious-pole hazard, not as a deliverable of the method.
- **No vector / multi-component / hierarchy / Garnier / Noumi–Yamada
  scope.**

## Gap survival verdict

**Amore 2021 does not close v0.2's gap.** It is the closest published
prior art on "Padé approximants + Taylor stepping = stiff scalar ODE
integrator with arbitrary precision," but the *vector lift with shared
denominator, applied to multi-component Painlevé-type systems, with
pole-field mapping* is absent.

It also has consequences for v0.1's framing: **the "Padé–Taylor"
name is shared with Amore 2021** for a method that overlaps in
algorithmic shape on scalar ODEs. The v0.1 README and PRD should
acknowledge Amore as parallel prior art for the scalar Padé–Taylor
method; v0.2 should also acknowledge it for the same reason. This is
**not** a v0.1 retroactive correctness issue (Amore handles only van
der Pol, FW2011 handles Painlevé pole fields — both are correct on
their own targets), but a citation hygiene issue.

## Acquisition status

`acquisition_status: open_access_arxiv`

- arXiv preprint: https://arxiv.org/abs/2111.12198 (free PDF).
- Published version: Elsevier Physica D 435, 133263 (2022) — paywalled,
  but the arXiv version is the canonical source.

The arXiv PDF was inspected at `/tmp/amore2021.pdf` during this sweep;
not retained inside `references/` because the paper is not Painlevé /
multi-component prior art *per se* — it is name-collision and scalar
Padé–Taylor parallel work. If the v0.1 docs decide to acknowledge it,
the fetch should land at:

```
references/Amore2021_pade_taylor_van_der_pol_PhysicaD_arXiv2111.12198.pdf
```

(top-level, alongside FW 2011 / GGT 2013, because it is
methodology-adjacent to v0.1, not v0.2-only).
