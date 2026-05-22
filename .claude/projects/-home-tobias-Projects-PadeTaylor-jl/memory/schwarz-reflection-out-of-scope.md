---
name: schwarz-reflection-out-of-scope
description: User is uncomfortable with Schwarz reflection — keep it out of scope
metadata:
  type: feedback
---

The user is uncomfortable with using Schwarz / conjugate-reflection symmetry
as a computational shortcut in the vector path-network solver. Do not propose
it for the headline-figure work (bead `padetaylor-0ln.40`) or scope it into
ADR-0026.

**Why:** Stated directly by the user 2026-05-22. FW 2011 itself does not use
Schwarz reflection to halve the work — it uses conjugate symmetry only as an
accuracy *diagnostic* (FW-md:303–310). So dropping it is also FW-faithful.

**How to apply:** When the vector solver needs more coverage or speed, reach
for resilience + denser targets + (if needed) genuine algorithmic fixes — not
symmetry halving. The scalar `PathNetwork._solve_with_schwarz_reflection`
exists but must not be ported to the vector stack.
