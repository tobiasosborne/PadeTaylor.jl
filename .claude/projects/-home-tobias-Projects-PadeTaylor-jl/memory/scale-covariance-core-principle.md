---
name: scale-covariance-core-principle
description: The solver must be scale-covariant — never hardwire an absolute length scale
metadata:
  type: feedback
---

The path-network solver must be **scale-covariant**. Never hardwire an
absolute length scale — not the step `h`, not a disc radius, not a pole
spacing, not a target spacing. Sparseness is scale-relative: zoom in and
any pole field is sparse, zoom out and any field is dense. `h` (and
every length the solver compares) must be *derived from the problem* —
specifically from the local pole spacing — not taken as a fixed input.

**Why:** User core principle, stated emphatically 2026-05-22. A fixed-`h`
"coverage ceiling" (e.g. the ~12–20% the headline-figure wedge probe
measured) is a **scale-mismatch artifact, not an architectural limit**.
FW 2011's method is general precisely because `h` tracks the local
scale; FW's fixed `h = 0.5` worked only because FW's P_I field was
roughly uniform in scale. A solver with a hardwired scale would imply FW
has no general method — an absurdity. FFW 2017 / ADR-0011 / ADR-0012
(adaptive step, non-uniform nodes, `CoordTransforms` log-maps) exist in
this very package for exactly this reason.

**How to apply:** Treat any absolute length compared to a threshold as a
red flag ("scale-fixing heresy"). The solver should compare *ratios*
(z-distances normalised by the local `h` or local pole spacing), never
raw lengths. When coverage/accuracy degrades, ask first "is `h` matched
to the local scale here?" before concluding anything architectural.
Reject "architectural ceiling" claims that have not been checked against
scale invariance. See [[schwarz-reflection-out-of-scope]].
