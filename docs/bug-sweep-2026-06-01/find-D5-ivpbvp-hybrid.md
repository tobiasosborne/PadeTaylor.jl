# Bug sweep D5 — IVPBVPHybrid (FFW 2017 §3 PFS↔BVP coupling + asymptotic-series IC)

## Area

`src/IVPBVPHybrid.jl` — the Tier-3 hybrid driver composing
`PathNetwork.path_network_solve` (PFS / Padé-Taylor IVP) with
`BVP.bvp_solve` (Chebyshev spectral) on the pole-free sector of the
tronquée P̃_III problem (FFW 2017 Figure 5). Special focus per the
assignment: (1) the PFS↔BVP matching condition at the sector boundary,
and (2) the asymptotic-series IC helper `pIII_asymptotic_ic` and its
hard-coded a_2 coefficient.

## References checked

- `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:31`
  — canonical P_III: `u'' = (u')²/u − u'/z + (αu²+β)/z + γu³ + δ/u`.
  Note `γu³` is **not** divided by z.
- `…FFW2017…md:43` — modified P̃_III RHS:
  `w'' = (1/w)(w')² + (1/4)(αw² + γw³ + βeᶻ + δe^{2ζ}/w)`.
- `…FFW2017…md:222` — sector geometry `−3π/4 < arg z < 9π/4` ↔
  `−3π/2 < Im ζ < 9π/2`; existence-theorem family `γ=0, α=1=−δ, β`
  arbitrary; `u ~ z^{1/3}`.
- `…FFW2017…md:230` — asymptotic ansatz
  `u ~ z^{1/3}[1 + Σ_{n≥1} a_n (z^{1/3})^{-2n}]`.
- `…FFW2017…md:232` — "generating the coefficients a_n recursively"
  (FFW use optimal truncation with many terms).
- `…FFW2017…md:243` (eq. 6) — published Fig 5 boundary ICs
  `z₂ = 30 e^{(-3π/4 + π/12)i}`, `u(z₂) = 2.384379236170 − 1.993845650158i`.
- `…FFW2017…md:245-247` — BVP coupling: PFS values on the two curved
  boundaries are the Dirichlet BCs; derivative-match convergence.
- `…FFW2017…md:252-264` — κ_r condition number derivation.
- `src/CoordTransforms.jl:120-156` — `pIII_z_to_ζ`, `pIII_ζ_to_z`,
  `pIII_transformed_rhs`.
- `src/BVP.jl:185-342` — `BVPSolution` struct + `bvp_solve`; node
  ordering `u_nodes[1]=u_b` (t=+1↔z_b), `u_nodes[N+1]=u_a` (t=−1↔z_a).
- `docs/adr/0014-ivp-bvp-hybrid.md`, `docs/worklog/039-ivp-bvp-hybrid.md`.
- `test/ivp_bvp_hybrid_test.jl:25-51` (IB.1.1 reference value).

All symbolic series algebra was verified independently with sympy
(see "Code behavior" cites below).

## Findings

### [HIGH] `a_2` asymptotic coefficient is wrong: ships ≈ −0.22208, correct value for the Fig 5 family is exactly 0

**Location**: `src/IVPBVPHybrid.jl:332-342` (the `a[2] = ((4/9) + δ·a[1]²)/(2δ)`
formula), with the flawed derivation in the comment block lines 267-331.

**Ground truth (cited)**: FFW md:31 (canonical P_III) + md:222 (family
`γ=0, α=1, δ=−1`) + md:230 (ansatz). Substituting
`u = s + a_1 s^{-1} + a_2 s^{-3} + …` (`s = z^{1/3}`) into the canonical
equation and collecting powers of `1/s` gives the order-by-order
balance (sympy, exact):

```
s^{-1}: −δ − 1                              ⇒ requires δ = −1
s^{-3}: a_1·δ − 2a_1 − β                    ⇒ a_1 = β/(δ−2)  = −β/3 at δ=−1
s^{-5}: −a_1²δ − a_1² + a_2·δ − 2a_2        ⇒ a_2 = a_1²(δ+1)/(δ−2) = 0 at δ=−1
s^{-7}: a_1³δ − 2a_1a_2δ − 2a_1a_2 + 4a_1/9 + a_3δ − 2a_3   ⇒ determines a_3
```

For the Fig 5 family (δ=−1): **a_2 = 0** and a_3 = 533/216000 ≈ 0.0024676.

**Code behavior**: `_pIII_asymptotic_coeffs` matches a_2 at the **wrong
order**. The comment (lines 267-269) states a_2 is "the a_2 coefficient
match at order z^{-7/3}" (i.e. s^{-7}) — but a_2 is in fact fixed at the
s^{-5} order, which the code never writes down. The code jumps from a_1
(s^{-3}) straight to s^{-7}, skipping s^{-5} entirely. Worse, the s^{-7}
balance it does write (line 322: `4a_1/9 = −δa_1³ + 2δa_1a_2`) is itself
missing the `−2a_1a_2` term: the comment at lines 285 and 318 asserts
"u²/z contribution at z^{-7/3}: 0", but `u²/z` has an s^{-7} term equal
to `2a_1a_2` (sympy: `u²/z = s^{-1}+2a_1 s^{-3}+(a_1²+2a_2)s^{-5}+2a_1a_2 s^{-7}+…`).
Net: the shipped formula `a_2 = ((4/9)+δa_1²)/(2δ)` evaluates to
−0.2220833… for δ=−1, a_1=1/60.

Cross-validation against FFW's *published* IC (md:243, eq. 6) at the
principal-branch-reachable point `z₂` (arg = −120°):
- correct series (a_2 = 0): `u(z₂) = 2.3843871 − 1.9938428i`, matching
  FFW's `2.384379236170 − 1.993845650158i` to **8.4e-6** (the residual
  is exactly the a_3·30^{−5/3} ≈ 8.5e-6 truncation, confirming a_2=0);
- shipped series (a_2 = −0.22208): `u(z₂) = 2.3880885 − 2.0002537i`,
  off from FFW by **7.4e-3** — three orders of magnitude worse.

The spurious shift in u from the wrong a_2 is `|a_2|·z^{-1} ≈ 0.2221/30
≈ 7.4e-3` (and ≈0.22 in the ζ-frame `w = z·u`, since `w ≈ z·u`).

**Mechanism (intermittent discontinuity)**: the wrong a_2 corrupts both
(a) the two boundary ICs that seed the PFS walks
(`solve_pole_free_hybrid` lines 491-496) and (b) the per-collocation-node
BVP initial guess (`_bvp_solve_on_slice` line 662, `u,_ =
asymptotic_ic_fn(z)` → `w = z·u`). Worklog 039 Friction 2 (lines 250-265)
documents that the BVP discretisation has **two Newton basins** — one
matching the true tronquée solution, one a different branch — and that
which basin Newton lands in depends on the initial guess. A per-node
initial guess that is wrong by ≈0.22 in w can tip a slice whose true
solution sits near the basin boundary into the *wrong* basin while
adjacent slices stay in the correct basin. Because `sol(ζ)` linearly
interpolates `w` across bracketing `Re ζ` slices (lines 762-774), one
slice landing in the wrong basin produces a jump between neighbouring
slices → an intermittent discontinuity in the computed sheet. The
data/order dependence (which slices flip) is precisely the "intermittent"
signature the maintainer reports.

**Intermittent?**: Yes — the always-on part is a fixed ~7e-3 IC bias; the
discontinuity is conditional on a near-basin-boundary slice being tipped,
which is order/N/ζ dependent (Friction 2's N-sweep already showed
`w(centre)` jumping between basins at different N).

**Confidence**: 0.9 (formula mismatch proven symbolically + cross-validated
against FFW's published IC eq. 6; the intermittency mechanism is
plausible-but-not-directly-demonstrated, hence not 0.95).

### [MEDIUM] Test IB.1.1 enshrines the wrong a_2 — it cannot catch the bug

**Location**: `test/ivp_bvp_hybrid_test.jl:38-42` (and worklog 039:174).

**Ground truth (cited)**: FFW md:243 publishes the real IC; CLAUDE.md
Rule 5 requires a test assert against a *known-correct* value, not a
value recomputed from the same code path.

**Code behavior**: the test pins `real(u(1000; n_terms=2)) =
10.001444583333333` and the comment (lines 38-41) *derives that number
from the same wrong formula* `a_2 = (4/9 − 1/3600)/(−2) = −0.22208`. The
correct value (a_2 = 0) is `10.001666666…`. So the reference is the bug,
not an independent oracle. The IB.1.2 ratio test and the M2 mutation
(perturb `u_sum` by +1e-3, worklog line 206) also cannot detect a wrong
a_2 because none of them compares against FFW eq. 6 or an independent
recurrence.

**Mechanism**: not itself a runtime discontinuity, but it is *why the
a_2 bug survived* — the GREEN suite gives false confidence. Listed
because the assignment asks to judge the truncation against FFW 2017.

**Intermittent?**: No (test-time, deterministic).

**Confidence**: 0.92 (the pinned constant is arithmetically identical to
the wrong-formula output; verified by computation).

### [LOW] `a_1 = −β/3` hard-coded for all δ, but is only valid at δ = −1

**Location**: `src/IVPBVPHybrid.jl:258` (`a[1] = −Float64(β)/3`), with the
helper advertising "δ and β are free" (docstring line 174).

**Ground truth (cited)**: FFW md:222 fixes the family at `α = 1 = −δ`;
the general s^{-3} balance gives `a_1 = β/(δ−2)`, which equals −β/3 only
at δ=−1. (The s^{-1} balance `−δ−1=0` already forces δ=−1 for `u ~ z^{1/3}`
to be a valid leading term at all.)

**Code behavior**: the helper validates `(α,γ)=(1,0)` but leaves δ
unchecked while computing `a_1 = −β/3`. A caller passing δ ≠ −1 silently
gets a wrong a_1 (and the leading ansatz is itself inconsistent). FFW's
family is δ=−1 only, so this is a latent loose-contract issue, not a
live Fig 5 bug.

**Mechanism**: would corrupt ICs (same propagation path as the a_2
finding) only if a caller supplies δ ≠ −1; for the shipped Fig 5 usage
(δ=−1) a_1 is correct.

**Intermittent?**: No (deterministic, and not exercised on Fig 5).

**Confidence**: 0.85 (general formula proven symbolically; severity LOW
because the live family pins δ=−1).

### [LOW] `sol(ζ)` does not perform the advertised PFS↔BVP glue continuity check

**Location**: `src/IVPBVPHybrid.jl:709-775` (docstring lines 720-723 claim
"asserts their values agree to within `glue_tol`"; the body never
evaluates the PFS side).

**Ground truth (cited)**: FFW md:245-247 — the coupling requires the
BVP boundary derivative to match the PFS boundary derivative (the
convergence criterion). ADR-0014 §"Step 4 — Glue" (md doc lines 61-64)
says "Continuity across the boundary is asserted to `glue_tol`; violation
fails loud."

**Code behavior**: `sol(ζ)` uses `glue_tol` only as slop in the
inside/outside membership test (lines 742-743) and returns the BVP
interpolant; it never compares against the PFS boundary value/derivative,
so a genuine matching violation passes silently. Worklog 039 Friction 4
(lines 302-310) acknowledges glue_tol is reduced to a membership tolerance,
but the docstring/ADR were not updated in lockstep (CLAUDE.md Law 2).
This is a docs↔code mismatch and a missing safety check — not a
computation-discontinuity source by itself, but it removes the guard that
would have *caught* the a_2-induced boundary mismatch at runtime.

**Mechanism**: a silent boundary mismatch (e.g. the ~7e-3/0.22 from the
a_2 bug, or a wrong-basin BVP slice) is never reported, so a
discontinuity at the PFS↔BVP seam can go undetected.

**Intermittent?**: No (the missing check is always missing).

**Confidence**: 0.8 (the absence of any PFS evaluation in the callable is
plain from the code; classifying impact as low because it is a guard, not
a producer).

## Areas verified correct

- **P̃_III RHS in `_bvp_solve_on_slice`** (`src/IVPBVPHybrid.jl:644-652`)
  matches FFW md:43 exactly: `f = wp²/w + (αw² + γw³ + βeᶻ + δe^{2ζ}/w)/4`
  with `eζ = exp(ζ)`, `e2ζ = exp(2ζ)`. The analytic Jacobians
  `∂f_w = −wp²/w² + (2αw + 3γw² − δe^{2ζ}/w²)/4` and `∂f_wp = 2wp/w` are
  the correct partials (verified by hand differentiation; the `−δe^{2ζ}/w²`
  sign in ∂f_w is right).
- **`pIII_transformed_rhs`** (`src/CoordTransforms.jl:150-156`) is byte-
  identical to the slice RHS and to FFW md:43.
- **Coordinate transforms** `pIII_z_to_ζ` / `pIII_ζ_to_z`
  (`src/CoordTransforms.jl:120-137`): `ζ=2log z`, `w=z·u`,
  `w'=(z·u+z²·u')/2`, and the inverse `u=w/z`, `u'=(2w'−w)/z²` are exact
  inverses and match the worklog derivation (md doc lines 85-90).
- **`a_1 = −β/3`** is the correct s^{-3} coefficient *for the shipped
  δ = −1 family* (= 1/60 at β = −1/20). Confirmed symbolically.
- **Termwise derivative of the asymptotic series** (lines 220-229): the
  `up_sum -= (2n−1)·a[n]·pow` with the `/(3 s²)` prefactor reproduces
  `d/dz[s(1 + Σ a_n s^{-2n})]` exactly (sympy: difference = 0). The minus
  sign and the (2n−1) factor are correct.
- **`s = z^{1/3}` exponent** (line 215) is correct (FFW md:230); the M4
  mutation guards it.
- **BVP BC↔endpoint wiring** in `_bvp_solve_on_slice` (za=im_lo↔w_b,
  zb=im_hi↔w_t) is internally consistent with `bvp_solve`'s node ordering
  (`src/BVP.jl:290-291`).
- **κ_r condition-number prose** (module docstring lines 18-25) matches
  FFW md:254 (eq. 7) and md:262/264 (`27/16·30^{4/3} ≈ 157`).
- **`_harvest_at_re` bracket/interp** (lines 587-619): the
  `searchsortedfirst(...) − 1` + `clamp(.,1,n−1)` + endpoint guards are a
  correct piecewise-linear lookup with no off-by-one; the conjugate-
  transpose class is N/A here (no `'` adjoints on complex data in this
  module).
