# Worklog 067 — SharedPade C1+C2 minimal correct core (the (m,m) build)

**Date**: 2026-06-02
**Author**: Claude Opus 4.8 (1M)
**Beads**: `padetaylor-d3a` (C1, window), `padetaylor-3p9c` (C2, reduction);
new `padetaylor-l48l` (VPN.5.7 re-tune).
**Scope**: build the minimal correct shared-Padé core now that the math is
settled (worklog 066, ADR-0027) — the GGT diagonal `(m,m)` window + graceful
degree-reduction — with dispatch/Pareto deferred to ADR-0028.

> **Result**: shipped the `+2` window, graceful reduction, `z^λ` cancellation,
> and `Q=1` Taylor-fallback. The construction is now **correct for the real
> meromorphic workload** — the pole-bearing oracle suites pass exactly
> (Weierstrass ℘, Calogero–Moser, the P_I^(2)/KKG figure) with **no throws** —
> while **entire / regular jets** carry the honest method-set `(m,m)` accuracy
> (looser than the buggy `+1` incidentally was). 49 entire-system test
> tolerances were updated to that honest accuracy with documentation; the
> dispatch layer (ADR-0028) will recover the tighter accuracy.

## Code changes (`src/SharedPade.jl`)

1. **C1 window** `:117` `idx = m+rr-cc+1 → +2` (GGT diagonal `(m,m)`, top-left
   `c_{m+1}`); fixed the docstring arithmetic (`:27,:34,:105-116`,
   `c_{m+r-c} → c_{m+r-c+1}`). The jet-length guard stays `≥ m+1` (NOT raised
   to `2m+1` — `robust_pade` zero-pads short jets, so raising it would break
   the bit-identical `d=1` contract; verified `RobustPade.jl:378`).
2. **C2 reduction** — accept `ρ ≥ m_cur` (exact-null `ρ==m_cur` OR full-rank
   `ρ==m_cur+1` least-squares Q); reduce `m_cur = ρ` when rank-deficient. The
   research's `==` over-correction was reverted after the BigFloat empirics:
   at τ≈1e-74 an entire system is full column rank at every degree, so
   reduce-until-σ≤τ throws — accepting the legitimate least-squares Q (research
   construction verdict) is required.
3. **`z^λ` cancellation** replaces the `Q(0)≈0` throw — cancel the common
   leading factor from `b` and every numerator (mirrors `RobustPade.
   _trim_and_normalise:461-476`); the spurious origin pole–zero cancels.
4. **`Q=1` Taylor fallback** for `ρ==0` (the `+2` window reads the high-order
   coefficients, which vanish for a locally-regular / small-step jet) — reduce
   and return `Q=[1]` with Taylor-polynomial numerators (mirrors the scalar
   `n==0` branch `:407-410`), instead of the old eager Mode-2 throw. A
   genuinely-zero jet is still rejected by an early `‖c‖ ≤ tol` guard.
5. Removed the dead `n_near>1` guard; updated the module docstring failure-mode
   list (Rule 2 lockstep).
6. **Oracle 2** (`external/probes/shared-pade-determinant-oracle/oracle.jl:160`)
   bumped `+1 → +2` + docstring corrected — it mirrored the same bug and could
   never have caught it.

## Verification

- **`shared_pade_test.jl` 150/150**, including the new transcendental
  **SP.1.7** (the only window-independent `d=1` oracle; the determinant oracle
  shares the window). C1 mutation-proven (revert `+2`→`+1` ⇒ SP.1.7 RED).
- **C2 mutation-proven**: disabling the `z^λ` cancellation (reinstating the
  `Q(0)≈0` throw) ⇒ `vector_stepper` harmonic step throws (RED); restored.
- **Meromorphic oracles GREEN** (the real workload): `vector_pipeline_oracle`
  (℘), `calogero_moser`, `kkg_pi2_figure`, `vector_path_network`,
  `vector_path_network_stage2`. No throws anywhere.

## The accuracy cost (honest, documented — the dispatch motivation)

The correct `(m,m)` construction is *less* accurate than the buggy `(m−1,m)`
on **entire / regular jets** (no genuine pole → least-squares Q at a reduced
degree). The cost **grows with the component count `d`** (the stack is
`(d·m)×(m+1)`, more over-determined for larger `d`):

| system | d | pre-fix `+1` | correct `+2` | tolerance now |
|---|---|---|---|---|
| harmonic `[cos,−sin]` (entire) | 2 | ~1e-10 | ~5e-9 | 1e-8 (VS.1.2/1.3, VP.1.2/1.4) |
| exp pair (entire) | 2 | ~1e-11 | ~4e-8 | 1e-7 (VS.1.3) |
| harmonic BigFloat-256 | 2 | ~3.4e-21 | ~1.6e-17 | 1e-16 (VS.1.2, VP.1.2) |
| **Noumi–Yamada A₄** Σf=t (meromorphic, regular stretches) | **4** | ~2e-9 | **~6.4e-6** | 1e-5 (NYF.1.2) |

This matters: the cost reaches the **real high-`d` NY target's conservation
invariant** (~6e-6, ~3 orders looser), not just entire toys. It is *correct,
not wrong* (pole structure recovered exactly; bounded, stable), and acceptable
for the qualitative vector figures. **The dispatch layer (ADR-0028) — selecting
the square `(m−1,m)` cell on regular/off-pole stretches — is expected to recover
this, most for high-`d` systems.** Each updated tolerance cites ADR-0027 +
`padetaylor-unk`; none masks a wrong answer (Rule 5).

## Known follow-ups

- `padetaylor-l48l` — VPN.5.7 `skip_covered` benefit demonstration needs
  re-tuning for the correct (tighter) pole discs (assertion relaxed to
  correctness + comparability for now; the feature itself is unaffected).
- `diagnose_test.jl` does not run in this environment — the optional
  `DelaunayTriangulation` extension is not instantiated (pre-existing, unrelated
  to this change).
- **ADR-0028 (next)** — dual-construction + validated-Pareto dispatch to recover
  the regular/high-d accuracy.

## References

- `docs/adr/0027-sharedpade-graceful-reduction.md`; `docs/worklog/066-...md`.
- `src/SharedPade.jl`; `src/RobustPade.jl:378,407-410,461-476`.
- `external/probes/sharedpade-offbyone-confirm/{confirm,c2_reduction_study}.jl`.
