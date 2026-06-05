# Corpus catalogue ERRATA

Corrections to scout-generated ground-truth recipes in
`01_corpus_catalogue.md` / `candidates.json` found during implementation
(epic `padetaylor-p1v0`). The raw scout output is kept verbatim for
provenance; the authoritative value is the one pinned in the test file.

These are **research-artifact recipe errors**, caught by independent
verification against canonical ground truth (Law 1) — *not* library bugs.

| Candidate | Catalogue value | Correct value | Notes |
|---|---|---|---|
| `jorba-zou-fallback-trigger-b6` (B6) | `(28!)^(1/27) = 12.3596…` | `(28!)^(1/28) = 11.298089984…` | Index/exponent mismatch in the scout's hand-derivation. The Jorba–Zou `_second_stepsize` criterion is `\|c[j]\|·h^j = 1` (same index `j` in coefficient and exponent). `src/StepControl.jl:194-197` is correct and faithful to TaylorIntegration.jl (confirmed in `docs/bug-sweep-2026-06-01/find-A3-stepcontrol.md:134-145`). Pinned the code-/TI.jl-faithful value in `test/corpus_taylor_jet_test.jl` (CTJ.3) with an anti-regression guard that the result must NOT equal the catalogue's `12.36…`. |
