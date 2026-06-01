# Bug sweep C1 — BVP + VectorBVP: Chebyshev D matrix + Newton collocation

## Area

`src/BVP.jl` (scalar 2nd-order Chebyshev-Newton BVP, 2-arg and 3-arg
RHS overloads, barycentric dense output) and `src/VectorBVP.jl`
(first-order vector-system Chebyshev collocation + τ-method + Newton +
`Taylor1` autodiff Jacobian). Focus per assignment: the Chebyshev
first-derivative matrix `D1`, `D2 = D1·D1`, boundary-row handling, the
Newton Jacobian sign/assembly, the affine-scale factor, the
barycentric weights, and any transpose/adjoint hazard on complex data.

Verdict: **no faithfully-transcribed-equation bug found.** Every
load-bearing formula in both modules matches the cited ground truth
(Trefethen SMIM Theorem 7 / `cheb.m`; Weideman–Reddy Eq. 14 / 15 +
`chebdif.m` / `chebint.m`; FW 2011 Eq. 3.2; the `bvp_recipe.md`
derivations). The D-matrix is additionally cross-validated to Octave
`chebdif` output in `test/_oracle_bvp.jl`, and the high-suspicion
τ-method pairing / scale factor / autodiff index are all mutation-proven
(`test/vector_bvp_test.jl` VB.7.1 A–D). One low-severity API-contract
laxity (complex-segment domain check uses only `real(t_star)`) is the
only deviation, and it cannot perturb the computed nodal solution.

## References checked

- `references/markdown/TrefethenSMIM_2000_book/TrefethenSMIM_2000_book.md:360`
  — Chebyshev points `x_j = cos(jπ/N), j=0…N` (descending +1…−1).
- `.../TrefethenSMIM_2000_book.md:407-429` — Theorem 7: off-diagonal
  `(D_N)_{ij} = (c_i/c_j)(−1)^{i+j}/(x_i−x_j)`, `c_i=2` at `i=0,N` else 1;
  corner `(D_N)_{00}=(2N²+1)/6`, `(D_N)_{NN}=−(2N²+1)/6`.
- `.../TrefethenSMIM_2000_book.md:444-461` — `cheb.m` source +
  negative-row-sum diagonal identity (6.6) `(D_N)_{ii}=−Σ_{j≠i}(D_N)_{ij}`.
- `references/markdown/WeidemanReddy2000_DMSUITE_ACMTOMS26/...md:419-425`
  — Eq. 14, `D^(1)_{k,j}=c_k(−1)^{j+k}/(c_j(x_k−x_j))`, `c_1=c_N=2`.
- `.../WeidemanReddy2000...md:414-417` — Eq. 15 barycentric, weight
  `(−1)^j/c_j` (endpoints halved).
- `references/markdown/FW2011_painleve_methodology_JCP230/...md:176-194`
  — affine map (3.1 preamble), Eq. 3.2 `D_2 u = ¼(z_b−z_a)²(6u²+…)`,
  BCs `u_0=u_b, u_N=u_a`, barycentric dense output, `D_1 u` endpoint
  derivative diagnostic.
- `references/bvp_recipe.md:36-93` (nodes, affine map, D1 off-diag +
  diagonal), `:140-207` (Newton residual/Jacobian + scale factor),
  `:221-263` (barycentric weights + back-map).
- `external/DMSUITE/chebdif.m:38,45-58` — `x=sin(...)` descending nodes,
  `C=toeplitz((-1).^k)` with endpoint doubling, `Z=1/(x_k−x_j)`,
  off-diag recursion, `D(L)=-sum(D')` negative-row-sum.
- `external/DMSUITE/chebint.m:26-36` — `w=(-1).^[0:N-1]`, `w(1)/2`,
  `w(N)/2`, coincident-node guard `D+eps*(D==0)`.
- `docs/adr/0014-ivp-bvp-hybrid.md:159-162` — 3-arg-RHS Jacobian
  `J = D₂_ii − scale·diag(∂f/∂u) − (scale/half_diff)·diag(∂f/∂u')·D₁_ii`.
- `docs/adr/0023-vector-bvp-solver.md:22-71` — discretisation, τ-method
  row replacement at node `t_0`, autodiff-Jacobian convention.
- `test/_oracle_bvp.jl` — pinned Octave `chebdif`/`chebint` outputs
  (D1/D2 cross-validation; corner `(D_N)_{00}=5.5` at N=4 = (2·16+1)/6).
- `test/vector_bvp_test.jl:398-455` — mutation-proof footer (A: scale,
  B: τ-pairing, C/D: autodiff coefficient index).

## Findings

### [LOW] Complex-segment dense-output domain check ignores Im(t*), allowing silent off-segment extrapolation

- **Location**: `src/BVP.jl:491` and `src/VectorBVP.jl:318`.
- **Ground truth**: the module docstrings promise a `DomainError` for
  any `z` "outside the segment" (`src/BVP.jl:133-135`, `:482-484`;
  `src/VectorBVP.jl:309-311`). The `bvp_recipe.md:259-263` back-map is
  `t* = (2z*−z_a−z_b)/(z_b−z_a)`; on a genuinely complex segment a
  query `z` off the segment line yields a `t*` with nonzero imaginary
  part.
- **Code behavior**: the guard tests only `real(t_star) ≤ 1+100eps &&
  real(t_star) ≥ −1−100eps`. For a complex segment, an off-segment `z`
  whose preimage has `|Im t*| ≫ 0` but `Re t* ∈ [−1,1]` passes the
  guard and is barycentric-extrapolated rather than rejected.
- **Mechanism**: this is a deterministic API-contract laxity in *dense
  output only*. It does not touch the Newton-solved `u_nodes` /
  `Y_nodes`, so it cannot introduce an intermittent discontinuity into
  the computed solution at the collocation nodes; at most it returns a
  silently-extrapolated value where the contract said it would throw.
  Listed for completeness, not as the suspected bug.
- **Intermittent?**: No (deterministic; query-geometry dependent only).
- **Confidence**: 0.6 that this is a real contract deviation; ~0.05
  that it is the maintainer's "intermittent discontinuity".

## Areas verified correct

- **`_chebyshev_D1` off-diagonal sign + endpoint weights**
  (`src/BVP.jl:528-537`, `src/VectorBVP.jl:375-383`). `ci/cj` with
  `c=2` at `i==1||i==np1` (0-indexed `0,N`), `sign=(-1)^{i+j}` via
  `iseven(i+j)`, divided by `t[i]-t[j]`. Exactly Trefethen Theorem 7
  Eq. 6.5 (`TrefethenSMIM...md:417,422`) and Weideman–Reddy Eq. 14
  (`...md:421`). `(-1)^{i+j}=(-1)^{j+k}` (parity), so commuted indices
  are harmless.
- **Diagonal negative-row-sum** (`src/BVP.jl:538-545`,
  `src/VectorBVP.jl:385-392`). `D[i,i]=-Σ_{j≠i}D[i,j]` is Trefethen
  Eq. 6.6 (`...md:458`) / `chebdif.m:56` `D(L)=-sum(D')`. Reproduces
  the corner `(2N²+1)/6` (oracle `test_d2_N4_D1[1,1]=5.5`).
- **Node ordering descending** `t_j=cos(jπ/N), j=0…N` →
  `t_0=+1,…,t_N=−1` (`src/BVP.jl:272`, `src/VectorBVP.jl:214`). Matches
  `TrefethenSMIM...md:360`, `bvp_recipe.md:19-24`, `chebdif.m:38`
  (which is descending), and oracle `test_d2_N4_nodes`.
- **Affine map orientation** `z(t)=half_sum+half_diff·t`,
  `half_diff=(z_b−z_a)/2` ⇒ `t=+1↦z_b`, `t=−1↦z_a`
  (`src/BVP.jl:275-277`, `src/VectorBVP.jl:215-217`). Matches FW 2011
  Eq. 3.1 (`FW2011...md:178`) and `bvp_recipe.md:41-44`.
- **Dirichlet BC enforcement** `u_nodes[1]=u_b` (node 0 ↔ z_b),
  `u_nodes[N+1]=u_a` (node N ↔ z_a) (`src/BVP.jl:290-291`); `bc_col`
  uses columns `1` (z_b,u_b) and `N+1` (z_a,u_a) (`:297`). Matches FW
  2011 `u_0=u_b, u_N=u_a` (`FW2011...md:190`) given the descending map.
- **Affine scale factor** `scale=(z_b−z_a)²/4` and residual `R =
  D2_ii·u + bc_col − scale·F` (`src/BVP.jl:294,310`). Exactly FW 2011
  Eq. 3.2 (`FW2011...md:186`) / `bvp_recipe.md:182`.
- **2-arg Newton Jacobian** `J = D2_ii − scale·Diagonal(∂f/∂u)`
  (`src/BVP.jl:315`). Exactly `bvp_recipe.md:189`. Sign correct.
- **3-arg Newton Jacobian** `J = D2_ii − scale·diag(∂f/∂u) −
  (scale·inv_hd)·diag(∂f/∂u')·D1_ii`, `inv_hd=1/half_diff`
  (`src/BVP.jl:426,439-440`). Matches ADR-0014 (`:159-162`); derivation
  checked: `scale·inv_hd = (z_b−z_a)/2 = half_diff`; chain-rule
  `∂u'_int/∂u_int = inv_hd·D1_ii` correct; `bc_col_D1 = D1[int,[1,N+1]]·
  [u_b,u_a]` pairs z_b-column with u_b (`:415,428`). Sign correct.
- **Barycentric weights + coincident-node guard**
  (`src/BVP.jl:560-585`, `src/VectorBVP.jl:484-503`). `sign_w =
  iseven(j-1)?+1:-1` = `(-1)^{j-1}` = 0-indexed `(-1)^j`; endpoints
  halved; second-form `num/den`. Matches Berrut–Trefethen / WR Eq. 15
  (`...md:416`) and `chebint.m:30-36`. Guard returns the node value at
  the continuous limit, so no discontinuity is introduced near nodes.
- **VectorBVP first-order residual** `R = (D1⊗I)·Y − s·f`, `s=(z_b−z_a)/2`
  (`src/VectorBVP.jl:216,419-424`). Correct chain rule `dY/dt=s·f`
  (ADR-0023 `:24-31`). Mutation A proves the `s` factor is load-bearing.
- **VectorBVP τ-method endpoint↔node pairing**
  (`src/VectorBVP.jl:235-239,427`): `B_a` → node-N block (z_a), `B_b` →
  node-0 block (z_b), BC rows replace node-0's `d` rows. Matches
  ADR-0023 (`:38-48`); Mutation B confirms the *swapped* pairing
  produces the reflected solution, so the shipped pairing is correct.
- **`Dop = kron(D1, I_d)`** node-major stack operator
  (`src/VectorBVP.jl:221`) — correct ordering for `[Y_0;…;Y_N]` layout.
- **`_autodiff_jacobian`** (`src/VectorBVP.jl:441-453`): fresh `yt` per
  column (no cross-column aliasing), seed component `i` with unit linear
  series, read coefficient `[1]` ⇒ `Jf[r,i]=∂f_r/∂y_i`. Matches
  ADR-0023 (`:56-65`). Mutation C/D confirm `[1]` (not `[0]`) is correct.
- **No transpose/adjoint hazard.** Neither file applies a postfix `'`
  (adjoint) or `transpose()` to `D1`/`D2`/`Dop` or any complex array;
  the only `D'` token is `D(L)=-sum(D')` *inside the MATLAB reference*
  `chebdif.m`, not the Julia port. The Julia diagonal is the explicit
  loop `D[i,i]=-Σ D[i,j]`, sidestepping the conjugate-transpose class
  entirely (the historical hazard noted at `padeapprox.m:113` does not
  recur here).
- **No in-place buffer aliasing / order dependence.** `u_int =
  u_nodes[int]` with `int=2:N` copies (range indexing), so `u_int.-=Δu`
  does not mutate `u_nodes`; `J = copy(Dop)` before in-place row edits
  (`src/VectorBVP.jl:261`); `J_bc` rebuilt from constants. No RNG, no
  iteration-order-dependent accumulation that could produce intermittency.
