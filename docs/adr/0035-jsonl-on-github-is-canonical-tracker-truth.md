# ADR-0035 — The git-tracked `.beads/issues.jsonl` on GitHub is canonical tracker truth

**Status**: **Accepted** — 2026-09-07.  Maintainer decision.  Resolves the
durability leg of `padetaylor-05ti`; forcing conditions for revisiting are
carried by the deferred bead named in *Consequences*.
**Date**: 2026-09-07
**Worklogs**: none; the evidence is in this session's commits `d1e84d6`
(recovery), `c1dfdea` (merge of the stranded originals), `de06c15`
(pre-commit hook), `df5f20f` (the `created_at` regression and its repair).

## Context

`bd` (beads) stores its data in an embedded Dolt database under
`.beads/embeddeddolt/`, which is **git-ignored**.  The only artifact that
leaves a machine is the exported `.beads/issues.jsonl`.  There are therefore
two candidate homes for the truth, and until now the project had never chosen
between them — which is precisely how it lost three months of tracker state.

What actually happened, measured rather than recalled: the 2026-08-23 session
ran on a second checkout, filed ~23 beads, closed ~19, committed its **code**,
and never ran `bd export`.  `git reflog` on this machine jumps from
`d950914 HEAD@{2026-06-19}` straight to a `pull: Fast-forward` on 2026-09-07,
and no commit in `bed1e22..HEAD` touches `.beads/`.  By 2026-09-07 the two
checkouts disagreed by 15 issues and 10 statuses, **and neither store was a
superset of the other** (`c1dfdea`: the other machine's DB held 337 issues,
git held 338).  A bare `bd export` on either machine would have destroyed the
other's work, and `bd export` auto-stages, so `git diff` would have read empty
while it happened.

Two durability options were on the table.  A real off-machine Dolt remote
(`bd backup init` to DoltHub, or `bd dolt remote add`) would give the tracker
its own replicated history, but needs credentials this project does not have
and adds a second sync path to keep honest.  The alternative is to declare the
artifact that **already** crosses machines — the JSONL in the GitHub repo —
the single source of truth.

## Decision

**The git-tracked `.beads/issues.jsonl`, as pushed to GitHub, is canonical.**

The local `.beads/embeddeddolt/` database is a **derived working replica**, not
a store of record.  It may be rebuilt from the JSONL at any time, and where
the two disagree, **the JSONL wins**.

Mechanically:

- **Export on the way out.**  `scripts/git-hooks/pre-commit` (ADR-less, added
  in `de06c15`) re-exports and stages the JSONL on every commit, so the
  tracker can never ship staler than the code it travels with.  Fail-soft: a
  missing `bd`, a database locked by a concurrent agent, or a failed export
  warns and lets the commit through.  It skips mid-merge/rebase, where
  re-exporting the local DB over the JSONL would silently discard whichever
  side is not yet imported.
- **Import on the way in.**  `scripts/git-hooks/post-merge` runs `bd import`
  whenever a pull changes the JSONL, so the replica tracks canonical without
  anyone remembering to.
- **Install once per machine.**  `core.hooksPath` is per-clone configuration
  that no clone or pull carries: `./scripts/git-hooks/install.sh`.

## Two limitations of `bd import`, both measured this session

`bd import` has upsert semantics — "new issues are created and existing issues
are updated" (`bd import --help`).  Upsert is not reconciliation, and it fails
to converge in two specific ways.  Both are load-bearing for anyone treating
the JSONL as truth:

1. **`created_at` is preserved from the local row, not taken from the JSONL.**
   Import updates content fields (this session's first import corrected 10
   statuses, and the `c1dfdea` import brought in every recovered
   `--- ORIGINAL 2026-08-23 BEAD TEXT ---` notes block) but leaves the existing
   row's `created_at` alone.  Consequence, observed: 23 beads that had been
   re-created locally on 2026-09-07 as reconstructions kept that date even
   after importing the JSONL in which `c1dfdea` had restored their true
   2026-08-23 creation dates — and the next export wrote the wrong dates back
   over canonical.  Caught only because the new pre-commit hook staged a
   48-line diff on a commit that made no bead writes.
2. **Import never deletes.**  A bead removed upstream lingers in the local
   replica forever and will be re-exported, silently resurrecting it.

**Reconciliation procedure** when the replica diverges beyond upsert's reach
(and the repair actually used in `df5f20f`): identify the divergent ids by
diffing an export against the canonical JSONL, then
`bd delete --from-file <ids> --force` followed by `bd import <rows>`.  Import
honours `created_at` **on create**, so delete-then-import restores it.  Check
first that the affected beads carry no dependency links or comments; deleting
also discards that bead's local audit events, which are not part of the export
and are unrecoverable (40 were lost in `df5f20f`, no issue content depended on
them).

## Consequences

**Accepted, with eyes open.**  GitHub is now a single point of failure for the
tracker.  That is a deliberate v1 trade: it is the same single point of failure
the *code* already has, it needs no new credentials or sync path, and the loss
mode it replaces — two divergent local databases, neither a superset, silently
forking for three months — is strictly worse and actually happened.

**Good.**  One artifact to reason about.  Tracker state is reviewable in diffs,
travels with the commit that motivated it, and is restorable on any machine by
`git clone` + `bd import`.  The failure that cost 23 beads is now structurally
prevented rather than left to discipline.

**Bad, and known.**  No tracker history independent of git.  The replica can
still drift silently in the two ways listed above; the post-merge hook warns on
divergence but cannot repair it unattended.  A force-push or a lost GitHub
account takes the tracker with it.

**Forcing conditions to revisit** (deferred bead `padetaylor-jgt6`): a second
undetected divergence between replica and canonical; a need for tracker history
or branching independent of the code repo; a second GitHub outage blocking a
session's close; or the project acquiring a DoltHub account, at which point
`bd backup init` becomes cheap and this ADR should be superseded rather than
amended.
