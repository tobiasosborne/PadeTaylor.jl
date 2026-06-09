# ADR-0029 — Certified ball oracles (Arb midpoint-radius) for the golden corpus

**Status**: **accepted — pilot built & GREEN 2026-06-09**. Pilot scope: ONE
oracle family (the equianharmonic Weierstrass ℘ of the headline `u'' = 6u²`
problem). Widening to the rest of the corpus is deferred (see Consequences).
**Date**: 2026-06-09
**Beads**: `padetaylor-krgy.5` (this ADR + `test/certified_oracle_test.jl`),
child of epic `padetaylor-krgy` (test-hardening sweep), Tier-2 §2.2 of
`docs/test_corpus/03_hardening_methodology.md:93-108`. Depends conceptually on
ADR-0002 (Arblib-has-no-SVD) and ADR-0003 (Arblib is a weakdep / test extra).

## Context

The corpus pins package output against high-precision **floats** — mpmath,
Mathematica `WeierstrassP`, and 80-digit BigFloat literals captured by
`external/probes/*/capture.{wl,py}` scripts. Those goldens are trustworthy
(cross-source agreement is gated before they are written) but they are not
**certified**: a golden could itself be subtly wrong — a transcription slip, a
convention mismatch, a capture-script bug — and no test in the suite would
notice, because every test compares the package against that same float.

Verified numerics close this gap. Arb **ball arithmetic** (midpoint ± radius)
computes a closed-form value together with a *rigorous* enclosure of the true
value: the radius is a machine-checked upper bound on the distance from the
midpoint to the truth. When that radius is provably tiny, the midpoint *is* a
proof-carrying golden.

## Decision

Adopt a **two-step certified-oracle assertion pattern**, piloted on the
Weierstrass ℘ family and recorded in `test/certified_oracle_test.jl`:

1. **Certify the oracle.** Compute the closed-form oracle (here ℘) in Arb
   *ball* arithmetic and assert `radius(ball) < ε_cert` (we use
   `ε_cert = 1e-40`; achieved radii are ~1e-74 at 256-bit). This is the
   certification — it proves the oracle midpoint equals the true value to
   ε_cert, upgrading "a high-precision float we trust" to "a provably-correct
   golden". It removes the "is our golden itself right?" risk.

2. **Pin the package.** Then pin the package's computed value to the ball's
   midpoint at the **package's own** tolerance:
   `abs(package − midpoint(ball)) < tol_pkg`, with `tol_pkg ≫ ε_cert`
   (e.g. `1e-13` at the Float64 Padé floor, `1e-6` past the bridged pole).
   The equivalent containment form — inflate the certified ball by `tol_pkg`
   via `Arblib.add_error!`, then `Arblib.contains(inflated, package)` — is
   asserted alongside so the two phrasings cannot silently diverge.

### Why this is NOT `contains(tight_ball, package)`

The naive route fails by construction. The package computes to ~1e-15
(Float64) / ~1e-30 (BF-256, here IC-limited to ~1e-15 by the 16-digit FW IC);
a 256-bit Arb ball has radius ~1e-74 — **far tighter** than the package error.
Containing the package value in the tight ball therefore FAILS: the package
value sits outside the certified ball. The certification value is the
*midpoint exactness* (step 1), with the package pinned to that midpoint at its
own, much looser, tolerance (step 2). The ~30–60 orders of magnitude between
ε_cert and tol_pkg is exactly the headroom the certification buys.

## The precision / width tradeoff (the crux of this ADR)

- **ε_cert** (oracle exactness) is bounded only by the Arb working precision.
  At 256-bit, ℘ enclosures are ~1e-74 wide — `ε_cert = 1e-40` is a comfortable,
  precision-independent bar. Raising precision shrinks the radius for free;
  this is the "width" knob and it is essentially as tight as we want.
- **tol_pkg** (package accuracy) is set by the *algorithm*, not the oracle:
  the order-30 (15,15) Padé floor, SVD conditioning, the published IC's digit
  count, post-pole error amplification. It is fixed by the package and must be
  honoured, not relaxed (CLAUDE.md "don't relax tolerances").
- The two live on opposite scales and must stay separated: `tol_pkg ≫ ε_cert`.
  If they ever collide, the oracle precision is too low — raise it, never lower
  `ε_cert` to paper over it.

## Why BALL, not INTERVAL

Inf-sup interval arithmetic suffers the **dependency problem**: repeated
occurrences of a variable are treated as independent, so enclosures
over-estimate (radius grows faster than the true error). Arb's midpoint-radius
form is tighter and auto-propagates a near-optimal radius — confirmed 3-0 in
the §2.2 research synthesis, and Arb's own randomized suite has even exposed
MPFR bugs. For a single special-function evaluation the difference is the
gap between a ~1e-74 ball and a uselessly wide interval.

## Why we certify the ORACLE, not the package's SVD output (ADR-0002)

Arblib ships **no SVD** (ADR-0002; verified by source inspection). The
package's Padé conversion routes `Matrix{Arb}` through `BigFloat` for the SVD,
which **discards the Arb radius** — so balls do not flow through the package's
numerical core. Rigorous error therefore cannot be propagated *through*
`solve_pade`. The correct, reconciled scope (bead NOTES, `wf_aefca141-353`) is
to certify the **closed-form oracle side** as a ball and pin the package value
(an ordinary high-precision float) to its midpoint. A fully Arb-rigorous SVD is
out of scope (ADR-0002 "Caveats", deferred).

## Pilot family + the Arblib call used

The headline problem `u'' = 6u²` has closed form `u(z) = ℘(z + c₁; g₂=0, g₃=c₂)`
with the FW 2011 IC `c₁ = -1, c₂ = 2` (an equianharmonic, g₂=0, lattice with a
real double pole at z = 1). Arblib exposes ℘ directly:

- `Arblib.elliptic_p!(res, z, τ; prec)` — ℘ on the lattice with periods (1, τ);
- `Arblib.elliptic_invariants!(G₂, G₃, τ; prec)` — the invariants of that lattice;
- `Arblib.elliptic_p_prime!` — ℘′, used for the ODE-invariant cross-check.

Arb's ℘ is parametrised by τ, giving *fixed* invariants. For g₂ = 0 the
fundamental-domain point is the equianharmonic τ = exp(iπ/3) (⇒ G₂ = 0,
G₃ ≈ 820.82). We rescale to the package's g₃ = 2 via the ℘ homogeneity
relation `℘(z; g₂, g₃) = λ⁻² ℘(z/λ; G₂, G₃)` with `λ = (G₃/g₃)^(1/6)` (real;
λ ≈ 2.726 is the real pole spacing). The whole construction runs in Acb ball
arithmetic, so the result carries a rigorous radius. The preferred ℘ pilot was
sufficient — the `robust_pade`-of-`exp` fallback was not needed.

## Consequences

- `test/certified_oracle_test.jl` (≤200 code LOC) ships four testsets: the IC
  golden certification (CO.1), package pins before (CO.2) and past (CO.3) the
  bridged pole, and a ball-residual check of the defining ODE invariant
  `(℘')² = 4℘³ − g₂℘ − g₃` (CO.4, guards the homogeneity rescale). Wired into
  `test/runtests.jl` after `ext_arblib_test.jl`.
- **No new dependency.** Reuses the existing `Arblib` test extra
  (`[extras]` / `[targets].test`); no ArbNumerics, no TaylorModels.
- Mutation-proofed BOTH directions (footer of the test): (a) perturb the pinned
  package value ⇒ the pin assertion goes RED while certify stays GREEN; (b)
  lower the working precision ⇒ the certify assertion goes RED while the pin
  stays GREEN — proving the two steps are orthogonal.
- **Deferred (Rule 9):** widening certified balls to the rest of the corpus
  (Jacobi cn/dn/sn via `Arblib` elliptic functions; Painlevé rational families;
  Padé-of-`exp`). Condition that forces the work: any new headline oracle whose
  float golden cannot be independently cross-sourced — certify it as a ball
  instead. Tracked under epic `padetaylor-krgy`.
- TaylorModels.jl as a "turnkey verified integrator" was REFUTED 0-3
  (`03_hardening_methodology.md:106`); not used here.

## Citations

- `docs/test_corpus/03_hardening_methodology.md:93-108` — Tier-2 §2.2 (the
  certified-oracle framing; ball-vs-interval; certify-oracle-not-SVD).
- `docs/adr/0002-bigfloat-svd-via-genericlinalg.md` — Arblib has no SVD; the
  radius is discarded through the BigFloat-SVD step.
- `docs/adr/0003-extensions-pattern.md` — Arblib is a weakdep / test extra.
- `test/problems_test.jl` 6.1.* and `test/_oracle_problems.jl` — the package
  values and float goldens this ADR certifies.
- Arblib `src/arbcalls/acb_elliptic.jl` — `elliptic_p!`,
  `elliptic_invariants!`, `elliptic_p_prime!`.
- DLMF 23.x — Weierstrass ℘ and the homogeneity relation 23.10.1.
