# Bug sweep E2 — off-by-one / index-shift errors across `src/`

Auditor focus: Taylor-coefficient index (0-based math vs 1-based Julia),
degree-vs-length confusion, MATLAB 1-indexed slices ported verbatim,
companion-matrix subscripts, cyclic-mod indices, segment/breakpoint
indexing. Read-only; no Julia executed.

## Area

WHOLE `src/` — every `.jl` file read for index/off-by-one transcription
errors, with primary attention on the MATLAB→Julia ports:

  - `src/RobustPade.jl` (port of `external/chebfun/padeapprox.m`, GGT 2013)
  - `src/SharedPade.jl` (vector type-II Hermite–Padé, mirrors RobustPade)
  - `src/LinAlg.jl` (SVD dispatcher — null-vector column indexing)
  - `src/Coefficients.jl`, `src/VectorCoefficients.jl` (Taylor recurrence)
  - `src/PadeStepper.jl`, `src/VectorStepper.jl` (rescale, Horner, FFW T(h))
  - `src/StepControl.jl`, `src/VectorStepControl.jl` (Jorba–Zou, ported from
    `external/TaylorIntegration.jl/src/integrator/stepsize.jl`)
  - `src/BVP.jl`, `src/Laplace2D.jl` (Chebyshev `cheb.m` / barycentric port)
  - `src/PathNetwork.jl`, `src/VectorPathNetwork.jl`,
    `src/VectorPathNetworkStage2.jl`, `src/VectorWedgeStep.jl`
    (wedge directions, nearest-visited, dense-output `t=(z-z_v)/h`)
  - `src/NoumiYamada.jl` (cyclic `mod` indices), `src/PoleField.jl`
  - `src/Problems.jl` (segment-breakpoint dense output)

## References checked

  - `external/chebfun/padeapprox.m:62` — `c = [f(:); zeros(m+n+1-len,1)]; c=c(1:m+n+1)` (pad/truncate).
  - `external/chebfun/padeapprox.m:76-77` — `row=[c(1) zeros(1,n)]; col=c`.
  - `external/chebfun/padeapprox.m:89` — `Z = toeplitz(col(1:m+n+1), row(1:n+1))`.
  - `external/chebfun/padeapprox.m:92` — `C = Z(m+2:m+n+1,:)` (the C̃ slice — KEY).
  - `external/chebfun/padeapprox.m:106-120` — `[U,S,V]=svd(C,0); b=V(:,n+1)`; QR-reweight `Q(:,n+1)`; `a=Z(1:m+1,1:n+1)*b`.
  - `external/chebfun/padeapprox.m:113` — documented historical `.'`-vs-`'` adjoint bug.
  - `external/chebfun/padeapprox.m:122-138` — leading/trailing trim + `b(1)` normalise.
  - `references/markdown/TrefethenSMIM_2000_book/TrefethenSMIM_2000_book.md:445-451` — `cheb.m` listing.
  - `references/markdown/TrefethenSMIM_2000_book/TrefethenSMIM_2000_book.md:458` — eq. (6.6) negative-row-sum diagonal.
  - `external/TaylorIntegration.jl/src/integrator/stepsize.jl:17-90` — Jorba–Zou `stepsize`/`_second_stepsize`.
  - `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:74-91` — T(h) error coefficient `ε_k = c_k + Σ_{r=1..ν} b_r c_{k-r}` and rescale eq. (2).
  - `references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex:85-88` — `\eqref{A2n}` RHS.
  - `docs/v0p2_pillarA_hermite_pade_findings.md:18-21` — GGT C̃ entry `c_{m+i-j}`; block null-space spec.
  - `docs/v0p2_pillarA_hermite_pade_findings.md:407-412` — pillar-A step 2: block (r,c) entry `jets[i][m+r-c+1]`.

## Findings

### [CRITICAL] `SharedPade._toeplitz_block` builds the C̃ block one Taylor-coefficient too low (off-by-one vs `padeapprox.m`)

**Location:** `src/SharedPade.jl:114-123` (the `idx = m + rr - cc + 1` line, 117).

**Ground truth (cited):**
- The denominator matching equations of a type-(m,m) Padé are the
  coefficients of `z^{m+1}, …, z^{2m}` of `f(z)Q(z)`: for `s = 1..m`,
  `Σ_{j=0}^{m} b_j c_{m+s-j} = 0`. Row `s`, column `j` (0-based j) holds
  `c_{m+s-j}`; in 1-based `(r,c)` that is the coefficient of order
  `m + r - c + 1`.
- The canonical implementation confirms this exactly:
  `external/chebfun/padeapprox.m:89,92`. `Z = toeplitz(col, row)` with
  `col = c` gives `Z(i,j) = c_{i-j}` (0-based), and the slice
  `C = Z(m+2:m+n+1,:)` (`padeapprox.m:92`) makes row `r` originate from
  MATLAB row `m+1+r`, so `C(r,c) = c_{(m+1+r)-c} = c_{m+r-c+1}`. The
  faithful scalar port reproduces this: `src/RobustPade.jl:413-414`
  (`Z = _lower_tri_toeplitz(cv, m+n+1, n+1)` with
  `Z[i,j]=c[i-j+1]=c_{i-j}` per `src/RobustPade.jl:199`, sliced
  `C = Z[(m+2):(m+n+1), :]`), giving the same `C[r,c] = c_{m+r-c+1}`.

**Code behavior:** `_toeplitz_block` sets `idx = m + rr - cc + 1` and
`blk[rr,cc] = c[idx]`. Since `c[1] = c_0`, the stored entry is the
Taylor coefficient of order `idx - 1 = m + rr - cc` — i.e. `c_{m+r-c}`,
exactly **one order below** the canonical `c_{m+r-c+1}`. Concrete check
at `m = 1`: the scalar `C̃ = [c_2, c_1]` (`padeapprox.m:92`/RobustPade),
but `_toeplitz_block(c,1) = [c_1, c_0]`. The block therefore encodes the
matching equations for `z^{m}, …, z^{2m-1}` instead of `z^{m+1}, …,
z^{2m}`; its right null vector is a denominator that matches the series
at the wrong orders.

The module's own docstring (`src/SharedPade.jl:108-112`) and the pillar-A
spec it cites (`docs/v0p2_pillarA_hermite_pade_findings.md:410`,
"(r,c) entry of block i is `jets[i][m + r - c + 1]`") both carry the same
off-by-one — the spec transcription of GGT eq. (2.10) / ManoTsuda eq.
(2.2) dropped the `+1` when converting the 0-based `c_{m+i-j}` formula
into a 1-based array index, and the code faithfully implemented the buggy
spec. (The companion numerator builder `_upper_block`,
`src/SharedPade.jl:132-141`, IS correct: `idx = rr-cc+1` ⇒ `c_{r-c}`,
matching `a_k = Σ_j b_j c_{k-j}` and RobustPade's `Z[1:m+1,1:n+1]`. Only
the denominator block is shifted, so the recovered `Q` — and hence the
pole set and every component value — is wrong, while the numerator
formula is correct but consumes the wrong `b`.)

**Why the test suite did not catch it:** `test/shared_pade_test.jl`'s
d=1 reduction oracle (SP.1.1, `test/shared_pade_test.jl:138-158`) and the
d≥2 cases (SP.1.2) feed **exact rationals** `f = P/Q` whose true `Q` has
degree ≤ the requested `m`. For an exact rational the true denominator
annihilates *all* the higher matching equations, so it lies in the null
space of BOTH the correct C̃ (orders `m+1..2m`) and the off-by-one block
(orders `m..2m-1`). The off-by-one is invisible to any exact-rational
oracle; it only bites on a genuinely non-rational jet (a real ODE Taylor
jet), which never reproduces a finite Padé exactly.

**Mechanism (intermittent discontinuity):** `_toeplitz_block` is the only
denominator-block builder for the entire vector pipeline:
`SharedPade.shared_denominator_pade` → `VectorStepper`
(`src/VectorStepper.jl:242`) → `VectorWedgeStep._candidate_pole_disc`
/ `_select_wedge` (`src/VectorWedgeStep.jl:345-355,406`) → the per-node
canonical store `_canonical_pade` (`src/VectorPathNetwork.jl:825,953-963`)
→ Stage-2 dense output `t=(z_f-z_v)/h_v`
(`src/VectorPathNetworkStage2.jl:436-438`) and `VectorPoleField`. Because
the recovered shared `Q` matches the series one order too low, its roots
(the system's "poles") are placed slightly wrong and its value `Q(1)`
carries an order-`h` extra error. The size of that error scales with the
truncation tail `c_{2m}·h^{2m}` of the jet, so it is **negligible where
the step is comfortably inside the radius of convergence and grows
sharply as a node approaches a pole** — precisely the near-pole region
the path-network is designed to probe. Two neighbouring grid cells served
by different visited nodes (different `(z_v, h_v)`, hence different
mis-placed `Q`) can therefore disagree by far more than the nominal
tolerance, producing the reported intermittent jumps in the computed
vector solution. The dependence on which node is "nearest" (a function of
the shuffle/visit order through `_select_wedge` and `_nearest_visited`)
makes the discontinuity data/order-dependent rather than uniform.

**Intermittent?** Yes — error is jet-tail-magnitude-dependent (vanishes
for exact rationals and far-from-pole steps, large near poles) and
nearest-node/visit-order dependent.

**Confidence:** 0.9. The matrix mismatch versus `padeapprox.m:92` and
`RobustPade.jl:413-414` is demonstrated arithmetically and by the m=1
worked example; the matching-equation derivation pins which side is
correct. The 0.1 reservation is that I could not execute Julia to
numerically confirm a non-rational jet produces a divergent `Q`, and a
subtle compensating convention elsewhere in the vector stack (none found
on read) cannot be 100% excluded without a run. Recommended confirmation
(no tolerance relaxation): feed a non-exact-rational jet (e.g. the
Taylor jet of `exp` or a real `A_2^(1)` ODE jet) to
`shared_denominator_pade([jet], m)` and to `robust_pade(jet, m, m;
method=:svd)` and compare `Q` — they should match and currently will not.

### [LOW] `PoleField._extract_poles_core` clusters after a non-stable sort keyed only on `|t*|`

**Location:** `src/PoleField.jl:156-167`.

**Ground truth (cited):** The module docstring (`src/PoleField.jl:64-68`)
specifies "Clustering is therefore *greedy in increasing* `|t*|`: the
first (best-placed) candidate to land in a cluster becomes its
representative."

**Code behavior:** `sort!(candidates; by = c -> c[2])` sorts by `|t*|`
only; Julia's default `sort!` is not stable. When two candidates share an
identical `|t*|` (e.g. a symmetric double-pole doublet straddling the true
location at equal radius), their relative order — and thus which becomes
the cluster representative — is unspecified.

**Mechanism:** could make the reported pole *representative* (and the
`min_support` accounting at the cluster boundary) depend on sort
internals rather than data. In practice ties at identical float `|t*|`
are rare and both tied candidates fall within `cluster_atol` of the same
physical pole, so the effect on the *reported pole location* is bounded
by `cluster_atol`, not a true solution discontinuity.

**Intermittent?** Yes (tie-dependent) but very rare and bounded.

**Confidence:** 0.3. This is a reproducibility nit, not a demonstrated
discontinuity in a computed trajectory; flagged for completeness.

## Areas verified correct

- **`src/RobustPade.jl` C̃ slice, null-vector, QR-reweight, numerator,
  trim/normalise** — `_lower_tri_toeplitz` (`:199`, `Z[i,j]=c[i-j+1]`),
  `C = Z[(m+2):(m+n+1), :]` (`:414`), `b = Vt[end,:]` (`:435`),
  `F.Q[:, n+1]` (`:446`), `a = Z[1:(m+1),1:(n+1)]*b` (`:450`), and the
  `_trim_and_normalise` `findfirst`/`findlast`/`b[1]` logic
  (`:463-495`) all reproduce `padeapprox.m:89,92,109,116,120,122-138`
  faithfully. The `adjoint(C*D)` at `:443` correctly uses the conjugate
  transpose, matching `padeapprox.m:113`'s post-July-2018 `'` (the
  documented historical `.'` bug is NOT present here).
- **`src/LinAlg.jl` null-vector column** — `Vt` from a full SVD of the
  `n×(n+1)` C̃ is `(n+1)×(n+1)`; `Vt[end,:]` is the right-singular vector
  of the smallest singular value (S descending), equal to MATLAB
  `V(:,n+1)` (`padeapprox.m:109`). `src/LinAlg.jl:38-44` documents this
  correctly.
- **`src/SharedPade.jl` numerator builder `_upper_block`** (`:132-141`):
  `c_{r-c}`, correct (`a_k = Σ_j b_j c_{k-j}`), matches RobustPade's
  `Z[1:m+1,1:n+1]`. Only the denominator block is shifted (finding 1).
- **`src/Coefficients.jl`** — first-order `y[k]=f_t[k-1]/k` (`:129`) and
  second-order `u[j]=f_t[j-2]/(j*(j-1))`, `up[j-1]=j*u[j]` (`:198-199`)
  reproduce FW 2011 §2.1.2 method (b)
  (`references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:96-107`)
  with the correct 0-based(TaylorSeries getindex)/1-based(loop) mapping.
- **`src/VectorCoefficients.jl`** — `y[i][k]=f_t[i][k-1]/k` (`:166`),
  the component-wise lift; identical recurrence indexing.
- **`src/PadeStepper.jl`** — `_rescale_by_powers` (`out[k]=c[k]·h^{k-1}`
  ⇒ `c̃_k=h^k c_k`, `:342-350`); Horner `_evaluate_pade` (`:367-382`)
  and derivative coefficients `(k-1)*P.a[k]`/`(k-1)*P.b[k]`
  (`:402-403`); FFW `ε̃_{n+1}` reads `coefs_u_scaled[n_int+2]` and
  `[n_int+2-r]` with `P_u.b[r+1]` (`:472-477`) — exactly FFW eq. at
  `FFW2017…md:76`. All indices correct.
- **`src/VectorStepper.jl`** — `_rescale_by_powers` (`:286-294`),
  `_eval_poly` Horner (`:305-311`), `Q(1)` evaluation and per-component
  `num(1)/Q(1)` (`:256-264`). Correct.
- **`src/StepControl.jl` `step_jorba_zou`** — loop `for k in (p-1,p)`
  reading `coefs[k+1]` (order-k) and `(ε/|·|)^(1/k)` (`:180-183`); the
  `_second_stepsize` fallback `for j in 1:(p-2)` reading `coefs[j+1]`
  (`:194-197`). Bit-for-bit equivalent to
  `external/TaylorIntegration.jl/src/integrator/stepsize.jl:23-33,77-89`
  given `Taylor1` getindex is 0-based-by-order. `step_pade_root` forward
  projection (`:253-261`) correct.
- **`src/VectorStepControl.jl` `vector_step_jorba_zou`** — same
  `(p-1,p)` loop and `_coef_vector(jets,k,d)=jets[i][k+1]`
  (`:190-193,225-226`); reduces to the scalar selector at d=1.
- **`src/BVP.jl` `_chebyshev_D1`** (`:525-547`) — off-diagonal
  `(c_i/c_j)·(-1)^{i+j}/(t_i-t_j)` with `c=2` at the two endpoints, plus
  the negative-row-sum diagonal; matches `cheb.m`
  (`TrefethenSMIM…md:448-451`) and eq. (6.6) (`…md:458`) exactly
  (`(-1)^{i+j}=(-1)^{i-j}`). Nodes `cos(jπ/N)` j=0..N (`:272`), affine
  map and `D2`-scale `((z_b-z_a)/2)^2` chain-rule factor (`:294,310`)
  correct. `_barycentric_eval` weights `(-1)^{j-1}` with endpoint-halving
  (`:577-579`) correct.
- **`src/Laplace2D.jl`** — `_cheb_D1` copied verbatim from BVP (same
  formula); not re-derived here beyond confirming it is the identical
  body.
- **`src/Problems.jl`** — dense-output segment scan
  `while k<length(h) && z_T>z[k+1]` (`:227-229`), `t=(z-z[k])/h[k]`,
  `up=P'(t)/h[k]` (`:231-235`); breakpoint count `n+1` vs `n` segments
  handled correctly; `taylor_eval` Horner (`:252-258`). Correct.
- **`src/PathNetwork.jl`** — wedge `θ=goal_dir+wedge_angles[k]`,
  `h_step=h(cosθ+i sinθ)` (`:943-945`); `@threads` writes `evals[k]`
  per-k with read-only shared inputs (no aliasing); `_select_candidate`
  `argmin` deterministic (`:991-998`); Stage-2 `t=(z_f-z_v)/h_v`,
  `up=P'(t)/h_v` (`:689-691`); canonical real-h Padé stored
  (`:605-618`); `_nearest_visited` lexicographic tiebreak (`:882-894`).
  Correct. (The known shuffle-dependent Schwarz asymmetry at
  `:751-769` is documented with a cure, not a transcription bug.)
- **`src/VectorPathNetwork.jl` / `VectorWedgeStep.jl`** — canonical store
  rebuilt with real `h_cur` at `z_new` (`:822-826`), matching the
  real-`h_v` Stage-2 eval; `_adaptive_h` maps root `t*`→z-distance
  `h_prev·|t*|` (`:550`); `_select_wedge` sequential, deterministic
  lexicographic max (`:401-451`). Indexing correct (the shared-Q itself
  is wrong via finding 1, but the *plumbing* indices are right).
- **`src/VectorPathNetworkStage2.jl`** — `t=(z_f-z_v)/h_v`,
  `num(t)/q_t` (`:436-438`); `_validity_radius` rescale + root clamp
  `0.5·h_v·min|t*|` (`:315-332`); `_rescale_by_powers`/`_eval_poly`
  copies match. Correct.
- **`src/NoumiYamada.jl`** — RHS `f_j' = f_j(Σ_{r=1..n} f_{j+2r-1} −
  Σ f_{j+2r}) + α_j` with `slot(j)=mod(j,d)+1`, `d=2n+1` (`:142-148`)
  is the verbatim `\eqref{A2n}`
  (`references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex:85-88`);
  cyclic mod indexing correct.
