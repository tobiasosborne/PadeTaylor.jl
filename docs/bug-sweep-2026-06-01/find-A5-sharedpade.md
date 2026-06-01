# Bug sweep A5 — SharedPade (shared-denominator type-II Hermite–Padé lift)

## Area

`src/SharedPade.jl` — the shared-denominator (type-II Hermite–Padé /
simultaneous Padé) lift of robust Padé to `d`-component vectors. Audited
against the scalar GGT path (`src/RobustPade.jl`), `external/chebfun/padeapprox.m`,
GGT 2013 (markdown), Mano–Tsuda 2017 (eq. 2.2 / 2.6, extracted from PDF),
and the pillar-A findings spec. Downstream consumers: `VectorStepper.jl`,
`VectorProblems.jl` (the vector-ODE path/segment solver — exactly where
intermittent discontinuities would surface).

## References checked

- `src/SharedPade.jl:114-141` — `_toeplitz_block` / `_upper_block`.
- `src/SharedPade.jl:185-268` — main solve: degree-reduction loop, SVD
  null vector, QR-reweighting, normalisation, trimming.
- `src/RobustPade.jl:196-202` — `_lower_tri_toeplitz` (`Z[i,j]=c[i-j+1]`).
- `src/RobustPade.jl:413-453` — scalar `:svd` path: `C = Z[m+2:m+n+1,:]`,
  null vector, QR-reweighting `qr(adjoint(C*D))`, `a = Z[1:m+1,1:n+1]*b`.
- `src/RobustPade.jl:459-498` — `_trim_and_normalise` (leading-zⁿ cancel +
  trailing trim + `b[1]=1`).
- `src/LinAlg.jl:42,75-109` — `pade_svd` returns `F.Vt` (= V*, the adjoint);
  null vector is `Vt[end,:]`.
- `external/chebfun/padeapprox.m:89-120` — `Z=toeplitz(col(1:m+n+1),row(1:n+1))`,
  `C=Z(m+2:m+n+1,:)`, `b=V(:,n+1)`, `qr((C*D)')` (line 113 documents the
  historical `.'`-vs-`'` adjoint bug), `a=Z(1:m+1,1:n+1)*b`.
- `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/...md:57-58`
  (eq. 2.10: `C̃` top-left entry is `c_{m+1}`), `:73-75` (b₀=0 / common
  factor zᵏ subcase), `:87-93` (Algorithm 1 steps).
- `references/hermite_pade/ManoTsuda2017_..._MathZ285.pdf` p.10 (eq. 2.2:
  `A^i_j(k,l)[m,n]=a^i_{j+m-n}`, top-left `a^i_j`), p.12 (eq. 2.6 with block
  `A^i_m(n,m+1)`; condition `[f_i P^{(0)}_0]_m^{m+n-1}=0`, numerator degree
  `≤ m-1`).
- `docs/v0p2_pillarA_hermite_pade_findings.md:69-114` (block construction),
  `:402-461` (algorithm steps 1-8 + failure modes).
- `docs/adr/0019-shared-denominator-pade.md:69-89` (bit-identical d=1 claim,
  QR-reweighting port).
- `test/shared_pade_test.jl:72-92,138-185,263-269,290-308` (SP.1.1/1.2/1.6/2.1
  all feed EXACT rationals via `_ratio_jet`).

## Findings

### [CRITICAL] Toeplitz denominator block is shifted one coefficient — d=1 does NOT reduce to the scalar GGT `(m,m)` path; it solves a Mano–Tsuda `(m-1,m)` system

**Location**: `src/SharedPade.jl:117` (`idx = m + rr - cc + 1`), used at
`:192-193`; interacts with `_upper_block` at `:262`.

**Ground truth (cited)**:
- GGT eq. 2.10 (`...GGT2013_...md:57-58`): the scalar denominator matrix
  `C̃` has **top-left entry `c_{m+1}`** and `C̃[r,c]=c_{m+r-c+1}`. The
  scalar path realises exactly this: `C=Z[m+2:m+n+1,:]` with
  `Z[i,j]=cv[i-j+1]` gives `C[r,c]=cv[(m+1+r)-c+1]=cv[m+r-c+2]=c_{m+r-c+1}`
  (`src/RobustPade.jl:413-414,199`). Top-left `C[1,1]=cv[m+2]=c_{m+1}`. ✓
- The matching condition `P=fQ+O(z^{m+n+1})` makes the *below-the-line*
  rows enforce vanishing of `[z^{m+1}],…,[z^{m+n}]` while `a₀…a_m` stay
  free (GGT eq. 2.6/2.8 split at row `m+1`).
- Mano–Tsuda eq. 2.2 (`...MathZ285.pdf` p.10): `A^i_j(k,l)[m,n]=a^i_{j+m-n}`,
  top-left `a^i_j`. The simultaneous-Padé system eq. 2.6 (p.12) uses
  `A^i_m(n,m+1)` → entry `a^i_{m+r-c}`, **top-left `a^i_m=c_m`**, and the
  condition `[f_i P^{(0)}_0]_m^{m+n-1}=0` enforces vanishing of
  `[w^m],…,[w^{m+n-1}]` with numerator degree `≤ m-1`.

**Code behavior**: `_toeplitz_block[rr,cc] = c[m+rr-cc+1]`, so top-left
`blk[1,1]=c[m+1]=c_{m-1+1}=c_m` (0-based `c_m`), i.e. it builds the
**Mano–Tsuda** block (matching window `[z^m]…[z^{m+n-1}]`), NOT the GGT
block (window `[z^{m+1}]…[z^{m+n}]`). The two differ by exactly one
coefficient row. To equal the scalar `C̃` the index must be
`m + rr - cc + 2`.

Worse, this is **internally inconsistent** with the numerator step: the
`_upper_block` (`:132-141`, `up[r,c]=c_{r-c}`) recovers `a₀…a_m` exactly
as the GGT "above-the-line" rows (it correctly mirrors
`Z[1:m+1,1:n+1]`). But the denominator block's FIRST row already forces
`[z^m]=Σ_k c_{m-k} b_k = 0`, which is the SAME quantity `_upper_block`
returns as `a_m`. So row `z^m` is double-counted: the recovered top
numerator coefficient `a_m` is driven to ≈0 by the shared null vector,
collapsing the result to numerator degree `m-1`. The split point between
"numerator rows" and "constraint rows" is off by one relative to a clean
`(m,m)` system.

Hand check (m=1, jet=[c₀,c₁,c₂]):
- scalar GGT: `C̃=[c₂ c₁]`, enforces `[z²]=c₂b₀+c₁b₁=0` → genuine type (1,1).
- SharedPade: `blk=[c₁ c₀]`, enforces `[z¹]=c₁b₀+c₀b₁=0` → type (0,1);
  `_upper_block` then yields `a₁=c₁b₀+c₀b₁≈0`. Different approximant.

Consumes one fewer Taylor coefficient than the scalar path: SharedPade
reads `c₀…c_{2m-1}` (jet-length requirement is `≥ m+1` at `:175`), the
scalar `robust_pade(jet,m,m)` reads `c₀…c_{2m}` (`padeapprox.m:62-63`,
`src/RobustPade.jl:377-381`). The highest available coefficient `c_{2m}`
is silently discarded.

**Mechanism (intermittent discontinuity)**: For an EXACT rational input
of type `(≤m, ≤m)` (every test uses `_ratio_jet`, an exact rational —
`test/shared_pade_test.jl:72-92,142,173-174,266,294`), the true `Q`
annihilates the matching equations at *every* power, so both windows
recover the same exact `Q`; the off-by-one is invisible. For a
**truncated transcendental jet** — the real workload in `VectorStepper`
stepping vector ODEs (`src/VectorStepper.jl:242`,
`src/VectorProblems.jl:290`) — the GGT `(m,m)` approximant and the
realised `(m-1,m)` approximant genuinely differ. The shared `Q` (whose
roots are the tracked poles and feed step control) is then a different,
one-order-lower-matched polynomial. Because `m` is chosen adaptively per
step and the jets differ per step, the size and even direction of the
error in the pole estimate varies step-to-step → the recovered pole
positions / solution values jump between steps that happen to be
near-rational (off-by-one harmless) and steps that are strongly
transcendental (off-by-one bites). That is the textbook signature of an
**intermittent** discontinuity in the computed trajectory, and it lives
at the heart of the core arithmetic — consistent with the maintainer's
suspicion.

**Intermittent?**: Yes — data-dependent. Invisible on exact-rational
inputs (hence all current tests pass); bites on generic/transcendental
jets, varying per step.

**Confidence**: 0.9. The block index is demonstrably one off from GGT
eq. 2.10 / the scalar `C̃` (cited line-for-line), and the module docstring
+ ADR-0019 both assert bit-identical d=1 reduction, which the code does
not deliver on non-rational inputs. The one residual uncertainty (→ not
1.0): the code DOES faithfully implement Mano–Tsuda eq. 2.6, so this could
be argued a deliberate convention choice rather than a transcription slip
— but then the upper/lower split is internally inconsistent (double-counts
`[z^m]`) and the stated scalar-equivalence contract is violated, so it is
a genuine defect either way.

### [MEDIUM] Degree-reduction loop accepts a full-column-rank stack (no genuine null space) for d≥2 and returns a spurious denominator

**Location**: `src/SharedPade.jl:196,208-209` (`ρ = count(s->s>τ,S)`;
`if ρ ≥ m_cur break`) and the guard at `:224-225` (`n_near = count(s->s≤τ,S)`;
throws only if `n_near > 1`).

**Ground truth (cited)**: scalar GGT breaks only on `ρ == n`
(`src/RobustPade.jl:420`, `padeapprox.m:95`, GGT Algorithm step 5
`...GGT2013_...md:90`) — i.e. exactly one singular value at/below τ, a
genuine isolated 1-D null space. Uniqueness needs rank exactly `m`
(Mano–Tsuda p.12, "unique iff rank = m"; `SharedPade.jl` docstring
`:75-78`).

**Code behavior**: `A_full` is `(d·m_cur)×(m_cur+1)`, so `S` has
`m_cur+1` singular values (for `d≥2`). The break test `ρ ≥ m_cur` admits
`ρ = m_cur+1` (ALL σ above τ, full column rank, NO null vector below
noise). The loop still falls through; `Vt[end,:]` is then the smallest
*signal* singular vector, not a null vector, and the QR step produces a
`b` that is not a real annihilator. The failure-mode-4 guard counts
`n_near = #{σ ≤ τ}`; when full-rank `n_near = 0`, so `n_near > 1` is
false and nothing throws. (For d=1 the stack is `m×(m+1)`, `S` has only
`m` values, `ρ ≥ m_cur ⇔ ρ = m_cur`, so this pathology cannot occur —
consistent with the scalar `ρ==n`.)

**Mechanism (intermittent discontinuity)**: For d≥2 vector steps where
no shared `Q` of degree `m` actually exists at tolerance (a common
regime far from singularities, where the components are not yet rational
to order `m`), the routine should reduce `m` or throw; instead it returns
a denominator extracted from a non-null singular direction. Its roots are
near-arbitrary, sensitive to rounding and SVD sign/ordering — so pole
estimates and step control jump intermittently exactly on the d≥2 steps
that hit this branch. Should at minimum require `ρ == m_cur` (one σ at/below
τ) before accepting, and throw / reduce otherwise.

**Intermittent?**: Yes — only on d≥2 steps where the stack is full column
rank at tolerance; order/rounding sensitive.

**Confidence**: 0.6. The branch is real and reachable; the precise harm
depends on whether real VectorStepper inputs ever present a full-rank
stack at the working tol (I could not run code to confirm frequency).
The divergence from the scalar `ρ==n` contract is certain.

### [LOW] Leading-zᵏ common-factor cancellation is omitted (diverges from scalar `_trim_and_normalise`; throws where GGT would succeed)

**Location**: `src/SharedPade.jl:240-245` (Q(0)≈0 → throw) and the
explicit note `:250-252` ("No leading-z^λ cancellation is needed here").

**Ground truth (cited)**: GGT step 6/7 (`...GGT2013_...md:73-75,92`),
`padeapprox.m:122-127`, and the scalar `_trim_and_normalise`
(`src/RobustPade.jl:463-476`) cancel a common factor `zᵏ` when
`b₀=…=b_{λ-1}=0` (the `C` singular / b₀=0 subcase) and then normalise.
This is a legitimate success path, not an error.

**Code behavior**: SharedPade throws "pole at expansion centre" (failure
mode 3) whenever `b[1]` is below threshold, instead of cancelling the
leading factor. For inputs where the GGT scalar path would cancel `zᵏ`
and return a valid lower-degree approximant, the d=1 SharedPade call
throws — a behavioural divergence from the asserted bit-identical
reduction (docstring `:46-53`, ADR-0019 `:69-77`).

**Mechanism (intermittent discontinuity)**: Less a discontinuity than an
intermittent hard failure — a step whose shared `Q` legitimately has a
factor `z` (b₀=0) aborts instead of returning, so the path solver hits an
exception on some steps and not others. Lower severity because it fails
loud (Rule 1) rather than returning a silent wrong value.

**Intermittent?**: Yes — only fires when recovered `b₀≈0`.

**Confidence**: 0.55. The omission vs the scalar path is certain
(line-cited); whether it ever fires on real inputs is unverified (no
execution allowed).

### [LOW] Failure-mode-3 guard uses an OR of two thresholds (looser than either alone)

**Location**: `src/SharedPade.jl:240` —
`abs(b[1]) > τ || abs(b[1]) > tol_t || throw(...)`.

**Ground truth (cited)**: scalar leading-zero detection compares against
a single threshold `tol` (`padeapprox.m:123`, `src/RobustPade.jl:463`).

**Code behavior**: throws only if `abs(b[1])` is below BOTH `τ=tol·‖c‖`
and `tol_t`. When `‖c‖>1`, `τ>tol_t` so the effective gate is `τ`; when
`‖c‖<1`, the gate is `tol_t`. The guard thus uses the *larger* of the two
thresholds (most permissive), letting a genuinely tiny `b[1]` through in
the regime where the other threshold would have caught it. Then `b ./ b[1]`
(`:253,262`) divides by a near-zero, amplifying noise into both Q and the
numerators.

**Mechanism (intermittent discontinuity)**: On steps where `b[1]` lands
in the gap between `tol_t` and `τ`, a near-singular normalisation injects
large, rounding-sensitive coefficients into `Q` and `P_i`, producing
spurious near-origin poles that jump per step. Narrow regime, hence low.

**Intermittent?**: Yes — only when `b[1]` falls in the inter-threshold gap.

**Confidence**: 0.4. The OR is clearly looser than a single threshold;
the practical impact depends on `‖c‖` and is unverified.

## Areas verified correct

- **Conjugate-transpose / adjoint hazard (the A1 class) — CORRECT.** The
  QR-reweighting uses `qr(adjoint(A_full * D))` (`src/SharedPade.jl:234`),
  matching the scalar `qr(adjoint(C*D))` (`src/RobustPade.jl:443`) and
  `padeapprox.m:113` `(C*D)'` — the post-July-2018 *adjoint* form, not the
  erroneous `.'` plain transpose. `D` is a real diagonal
  (`abs(bk)+eps_T`, `:233`), so no hidden transpose-vs-adjoint ambiguity.
  `b_init = Vt[end,:]` (`:219`) is used only to size the real diagonal `D`
  (magnitudes are conjugation-invariant), exactly as the scalar path; the
  V1a mutation record (`test/...:242-249`) independently confirms `b_init`
  is conditioning-only. No `transpose()` / `permutedims` on complex data
  anywhere in the module.
- **QR null-column selection — CORRECT.** `F.Q[:, m_cur+1]` (`:235`) is
  the right generalisation of the scalar `F.Q[:, n+1]`
  (`src/RobustPade.jl:446`): `adjoint(A_full*D)` is
  `(m_cur+1)×(d·m_cur)`, its Householder `Q` is `(m_cur+1)×(m_cur+1)`, and
  column `m_cur+1` is the null direction of `A_full*D`. The V1a mutation A′
  (`test/...:251-254`) proves selecting `m_cur` instead is caught.
- **`_upper_block` numerator recovery — CORRECT (matches scalar).**
  `up[r,c]=c_{r-c}` for `r≥c` (`:132-141`) reproduces `Z[1:m+1,1:n+1]` at
  `n=m` (`src/RobustPade.jl:450,199`); `a = upper·b` is the GGT/Mano–Tsuda
  "above-the-line" product. (It is the *denominator* block that is shifted,
  not this one — see the CRITICAL finding.)
- **Trailing-zero trim thresholds — CONSISTENT with scalar.** Denominator
  trimmed at `tol_t` (`:254`) ↔ scalar `tol_typed` (`src/RobustPade.jl:479`);
  numerator trimmed at `τ` (`:263`) ↔ scalar `ts_typed`
  (`src/RobustPade.jl:484`). Matches `padeapprox.m:130,134`.
- **`b[1]=1` normalisation sign-robustness — CORRECT.** `b ./= norm(b)`
  then divide a and b by `b[1]` (`:236,253,262`); global sign cancels in
  P/Q, confirmed by V1a mutation B (`test/...:256-259`).
- **τ definition — CORRECT.** `τ = tol·‖vcat(jets)‖₂` (`:182-183`) is the
  vector analogue of GGT `τ = tol·‖c‖₂` (`...GGT2013_...md` Algorithm 2
  step 2; `padeapprox.m:66`; scalar `src/RobustPade.jl:384-385`).
- **`pade_svd` null-vector convention — CONSISTENT.** `Vt[end,:]`
  (`:219`) matches the scalar usage (`src/RobustPade.jl:435`) and the
  documented `F.Vt = V*` return (`src/LinAlg.jl:42`).
