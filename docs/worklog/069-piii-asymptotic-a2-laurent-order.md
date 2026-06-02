# Worklog 069 — PIII asymptotic a₂ fixed (wrong Laurent order; a₂=0 at δ=−1) (C4)

**Date**: 2026-06-02
**Author**: Claude Opus 4.8 (1M) — orchestrated (Opus coding subagent), derivation
cross-verified independently by the orchestrator.
**Bead**: `padetaylor-qsj` (sweep C4, P1 bug, conf 0.95).
**Scope**: re-derive the tronquée-PIII asymptotic IC coefficient a₂ from the
canonical FFW equation, fix the shipped value + its three duplicates + the ADR
+ worklog, and replace the self-referential IB.1.1 test (Rule 5 violation) with
an FFW-published-IC oracle.

> **Result**: a₂ was matched at the **wrong order (s⁻⁷)**, giving the spurious
> `a₂ ≈ −0.22208`. The correct **s⁻⁵ balance** gives `a₂ = −a₁²(α+δ)/(2α−δ)`,
> whose `(α+δ)` factor vanishes under the tronquée validity constraint `δ=−α`
> ⇒ **a₂ = 0** for the Fig-5 family. The corrected n_terms=2 IC matches FFW's
> independently-published `u(z₂)` to **8.41e-6** (the old value was **7.4e-3**
> off — ~880× worse). `ivp_bvp_hybrid_test.jl` 65/65, `ffw_fig_5_test.jl` 23/23.

## Root cause (sweep C4)

`src/IVPBVPHybrid.jl:_pIII_asymptotic_coeffs` matched a₂ at the `s⁻⁷ (z^{−7/3})`
order via lengthy (and internally inconsistent) hand-algebra, shipping
`a₂ = ((4/9) + δ a₁²)/(2δ) ≈ −0.22208`. a₂ is fixed one order **earlier**, at
`s⁻⁵`. The exported `pIII_asymptotic_ic` is load-bearing at the PFS anchor
corners (`:491-496`), so the wrong a₂ injected a fixed ~7e-3 IC bias.

## Derivation (sympy, cross-checked by hand against FFW md:31/230/243)

Ansatz `u = z^{1/3}[1 + Σ aₙ z^{−2n/3}]` into the **FFW md:31** equation
`u'' = (u')²/u − u'/z + (αu²+β)/z + γu³ + δ/u` (γ=0 family), residual collected
in `w = z^{−1/3}`:

| order | residual coefficient | fixes |
|---|---|---|
| `z^{−1/3}` | `−(α+δ)` | **consistency ⇒ δ=−α** (FFW md:222: α=1=−δ) |
| `z^{−1}` | `−a₁(2α−δ) − β` | `a₁ = −β/(2α−δ) = −β/(3α)` (=1/60) |
| `z^{−5/3}` | `−a₁²(α+δ) − a₂(2α−δ)` | `a₂ = −a₁²(α+δ)/(2α−δ)` **= 0 at δ=−α** |
| `z^{−7/3}` | `a₁³δ − 2a₁a₂(α+δ) + (4/9)a₁ − a₃(2α−δ)` | a₃ = β(β²−4)/81 ≈ 2.47e-3 |

The bead's hypothesis `a₂ = β²(δ+1)/(δ−2)³` shares the load-bearing `(δ+1)`
numerator (so 0 at δ=−1) but its `(δ−2)³` denominator is **wrong** — the
correct denominator is `9(δ−2)` at α=1 (they agree only at δ=−1). We ship the
sympy-verified `a₂ = −a₁²(1+δ)/(2−δ)` (the α=1 specialisation, the helper's
validated domain).

## Fix (Law 2 — code + 3 duplicates + ADR + worklog)

- `src/IVPBVPHybrid.jl`: `a[2] = -(a[1]^2*(1+δ))/(2-δ)`; the s⁻⁷ comment block
  replaced with the s⁻⁵-balance derivation table above; **fail-loud `δ == -α`
  guard** added (Rule 1 — the ansatz is invalid otherwise) with a `suggestion`;
  docstring + a₃ note corrected.
- `figures/ffw2017_fig_5.jl:~224` + two comment copies — fixed.
- `test/ffw_fig_5_test.jl:~84` — a **third** duplicate (not in the bead spec)
  found and fixed; FF5.1.1 tol tightened **1e-2 → 5e-5** (toward accuracy,
  measured 8.41e-6).
- `docs/adr/0014-ivp-bvp-hybrid.md:~67-97` — transcription (FFW md:31 form),
  the s⁻⁷ wording, and the −0.22208 claim corrected.
- `docs/worklog/040-ffw2017-fig-5.md` — error-budget corrected (a₃≈2.47e-3 ⇒
  ~8e-6 IC floor, not the old O(1)→3e-4 over-estimate).

## Test (TDD, independent oracle)

New **IB.1.1** asserts `pIII_asymptotic_ic(30·e^{−2πi/3}; n_terms=2, β=−1/20,
δ=−1)` against FFW eq.6 `u(z₂)=2.384379236170−1.993845650158i` (atol 1e-5,
measured 8.41e-6) and `u'(z₂)` (atol 1e-6, measured 4.62e-7) — FFW's published
values, **not** recomputed from a₂. Structural cross-check: a₂=0 ⇒ n_terms=1 and
n_terms=2 give bit-identical u (which is exactly why the old `abs(u₂−u₁)>1e-3`
assertion was structurally impossible). **Mutation-proven**: re-introducing the
old a₂ ⇒ IB.1.1 RED (z₂ off by 7.4e-3). New IB.1.7 covers the δ≠−α throw.

## Follow-ups

- New bead (P3): thread α through `_pIII_asymptotic_coeffs` (latent α=1 hardcode;
  currently safe behind the α==1 + δ=−α validation).
- Historical diagnosis docs (worklog 039/064, bug-sweep-2026-06-01) still cite
  −0.22208 — left untouched as accurate *records of the bug*, not design docs.

## References

- `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/...md:31,222,230,243`.
- `src/IVPBVPHybrid.jl`; `docs/adr/0014-ivp-bvp-hybrid.md`; worklog 040, 065 (C4).
