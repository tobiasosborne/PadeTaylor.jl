# v0.2 plan — vector / multi-component Painlevé-type systems

> Stage 0+ deliverable. Consolidates the four pillar findings docs
> (`docs/v0p2_pillar{A,B,C,EF}_*.md`) into an executable 10-phase plan.
> The plan is tracked as beads epic `padetaylor-0ln` (+19 children).
> Supersedes nothing in the PRD — it *implements* the PRD §"v0.2 north
> star" acceptance checklist.

## Goal

Lift PadeTaylor from scalar analytic ODEs `y' = f(z, y)`, `y ∈ ℂ`, to
**vector** `y ∈ ℂ^d` — the first general-purpose numerical software
for multi-component Painlevé-type systems. Primary targets: the
Noumi–Yamada higher-order Painlevé systems `A_n⁽¹⁾` and the first
Painlevé hierarchy member `P_I⁽²⁾`.

## The keystone primitive — shared-denominator robust Padé

A vector meromorphic solution's components all blow up at the *same*
poles. The right representation is therefore `d` numerator polynomials
`P_1,…,P_d` over **one shared denominator** `Q` — a simultaneous
(type-II Hermite–Padé) approximant. The algorithm
(`docs/v0p2_pillarA_hermite_pade_findings.md`; Mano–Tsuda 2017 §2.2):

> Stack the `d` component Toeplitz blocks (each `m×(m+1)`, built from
> that component's Taylor coefficients exactly as in scalar GGT 2013)
> into one `dm×(m+1)` matrix. Take **one** SVD. The right singular
> vector of the smallest singular value is `Q` (`‖Q‖₂ = 1`). Recover
> each `P_i` by convolution. Port the `padeapprox.m` QR-reweighting
> unchanged. At `d = 1` this reduces *exactly* to GGT 2013 `:svd`.

Failure modes that must throw (Rule 1): jet too short; all singular
values below `tol·‖c‖`; `Q(0) ≈ 0`; non-isolated null space after
degree reduction.

## Architecture decision — additive

The vector lift is **additive**: new modules alongside untouched v0.1
scalar modules. Rationale — (a) v0.1 scalar tests must stay
bit-identical (backward compat is non-negotiable); (b) the v0.1
Float64 default `:classical` must not be perturbed by the SVD-based
shared-Q; (c) the Rule 6 ≤200-LOC cap. New modules:
`SharedPade.jl`, `VectorCoefficients.jl`, `VectorStepper.jl`,
`VectorProblems.jl`, plus a Noumi–Yamada and a Painlevé-hierarchy
problem layer.

## The 10 phases

| Phase | Bead | Deliverable | Depends on |
|---|---|---|---|
| V0 | `0ln.1` | Stage 0+ consolidation: findings docs ✓, `references/tex/` ✓, worklog 054, this plan | — |
| V1a | `0ln.2` | `SharedPade.jl` — stacked-Toeplitz SVD shared-Q core + QR-reweight | V0 |
| V1b | `0ln.3` | Oracle: ACM Calgo 766 (FORTRAN) via `ccall` | V1a |
| V1c | `0ln.4` | Oracle: independent Julia block-Toeplitz determinant | V1a |
| V1d | `0ln.5` | Oracle: AAA per-component pole cross-check | V1a |
| V1e | `0ln.6` | Shared-Q mutation-proof test suite (3 oracles) + ADR-0019 | V1a–d |
| V2 | `0ln.7` | `VectorCoefficients.jl` — `Vector{Taylor1{T}}` jet + ADR-0020 | V0 |
| V3a | `0ln.8` | `VectorStepper.jl` — vector `pade_step!` | V1e, V2 |
| V3b | `0ln.9` | `VectorProblems.jl` — vector `solve_pade` driver | V3a |
| V3c | `0ln.10` | Vector norm-based Jorba–Zou step control + ADR-0021 | V3a |
| V4 | `0ln.11` | Calogero–Moser N=2 exact-oracle smoke test | V3b, V3c |
| V5a | `0ln.12` | `NoumiYamadaProblem` builder + `A_n⁽¹⁾` RHS + ADR-0022 | V3b, V3c |
| V5b | `0ln.13` | `A_2⁽¹⁾ ⇒ PIV` backward-compat self-validation | V5a |
| V5c | `0ln.14` | `noumi_yamada_rational` oracle + W(A_n⁽¹⁾) self-test | V5a |
| V6 | `0ln.15` | `painleve_hierarchy(:I,m)` + `P_I⁽²⁾` 4-vector RHS | V3b, V3c |
| V7 | `0ln.16` | Vector path network + pole extraction from shared `Q` | V3b, V5a, V6 |
| V8a | `0ln.17` | Figure: `A_4⁽¹⁾` Noumi–Yamada pole field (PRD n≥4) | V5a, V7 |
| V8b | `0ln.18` | Figure: `P_I⁽²⁾` tritronquée vs Kapaev–Klein–Grava 2015 | V6, V7 |
| V9 | `0ln.19` | v0.2 docs: README/CHANGELOG/RESEARCH, ADR review | V8a, V8b |

`V1` is the critical path: nothing in the vector stack proceeds until
the shared-Q primitive is built and triple-oracle-validated.

## Ground-truth equations (transcribed for the implementer)

**Noumi–Yamada A_n⁽¹⁾, even `l = 2n`** (`docs/v0p2_pillarB_*`):
`f_j' = f_j(Σ_{r} f_{j+2r-1} − Σ_{r} f_{j+2r}) + α_j`, indices mod
`l+1`, with `Σf_j = t`, `Σα_j = 1`. Phase-space dimension `2n`.
A_2⁽¹⁾⇒PIV via `y = −f_1/c`, `c = √(−3/2)`.

**P_I⁽²⁾** (`docs/v0p2_pillarC_*`; KKG 2015 eq. 1.1):
`u_xxxx + 10 u_x² + 20 u u_xx + 40(u³ − 6tu + 6x) = 0`.
4-vector form `y = (u, u', u'', u''')`:
`y' = (y_2, y_3, y_4, −10y_2² − 20y_1y_3 − 40(y_1³ − 6t y_1 + 6x))`.

**Calogero–Moser N=2 oracle** (`docs/v0p2_pillarEF_*`):
rational KdV `u = 2/(x−x_1)² + 2/(x−x_2)²`, poles obey
`ẍ_i = 4/(x_i−x_j)³`; with `x(0) = ±1, ẋ(0) = 0` the exact
trajectories are `x_{1,2}(t) = ±½√(1 + 4t²)`.

## Acceptance (PRD §v0.2 north star)

- [ ] Shared-Q `robust_pade` — mutation-proof vs 3 independent oracles (V1e).
- [ ] Vector `Coefficients`/`PadeStepper`/`Problems` GREEN on an
      `A_2⁽¹⁾` smoke test; v0.1 PIV self-validation bit-identical (V3, V5b).
- [ ] An `A_n⁽¹⁾` pole-field figure, `n ≥ 4`, reproducible from `figures/` (V8a).
- [ ] A `P_I⁽²⁾` tritronquée pole-field figure vs Kapaev–Klein–Grava 2015 (V8b).
- [ ] All v0.1 tests still GREEN (continuous).
- [ ] ADRs 0019 (shared-Q), 0020 (vector-y storage), 0021 (vector
      step-control), 0022 (new vector problem types).

## Out of scope

Garnier `G_n` (multi-`z` — a PDE / multi-time problem, a different
project); performance work and registry registration (v0.3); tau-function
evaluation and rigorous Arb-ball enclosures (PRD §out-of-scope).
