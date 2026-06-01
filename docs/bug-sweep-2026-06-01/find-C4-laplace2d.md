# Bug sweep C4 — Laplace2D: tensor-product Chebyshev Dirichlet Laplacian

Auditor: read-only sweep, 2026-06-01. Hunt: subtle equation/algorithm
mis-transcription producing intermittent discontinuities.

## Area

`src/Laplace2D.jl` — in-house 2D-Chebyshev spectral Dirichlet Laplace
solver on a rectangle `[a,b]×[c,d]` (voter (2) of the ADR-0024
triple-method harmonic-extension fill). Specifically:

- the affine-scaled second-derivative matrices `D2x`, `D2y`;
- the Kronecker-sum interior Laplacian `L`;
- the interior/boundary index split and Dirichlet-into-RHS absorption;
- the column-major `vec`/`reshape` consistency;
- the tensor-barycentric dense-output callable.

## References checked

- `references/markdown/TrefethenSMIM_2000_book/TrefethenSMIM_2000_book.md:445-461`
  — `cheb.m`: `x = cos(pi*(0:N)/N)'` (descending, node 1 = +1), weight
  vector `c = [2; ones(N-1,1); 2].*(-1).^(0:N)'`, off-diagonal
  `D=(c*(1./c)')./dX`, diagonal by negative-row-sum identity (eq. 6.6).
- `references/markdown/TrefethenSMIM_2000_book/TrefethenSMIM_2000_book.md:747`
  — discrete Laplacian (7.5) of dimension `(N-1)²×(N-1)²` as a Kronecker
  sum on the interior grid in lexicographic ordering (p. 68).
- `references/markdown/TrefethenSMIM_2000_book/TrefethenSMIM_2000_book.md:754-760`
  — Program 23 (compare p16.m): `D2 = D^2; D2 = D2(2:N,2:N); I = eye(N-1);
  L = -kron(I,D2) - kron(D2,I);`. NOTE: ch. 7 Program 16 itself (p. 68-70,
  eq. 7.5, the lexicographic-ordering paragraph) is **not present** in this
  markdown extraction; lines 754-760 (Program 23, the eigen-problem variant)
  are the only Kronecker-sum Laplacian in the file. The sign-vs-eigenvalue
  framing differs (Program 23 negates for `-Δ`), so it had to be matched
  carefully against the actual elliptic solve.
- `references/markdown/BerrutTrefethen2004_barycentric_SIAMReview/...md:130`
  (second-form formula eq. 4.2) and `:171` (Chebyshev-2 weights eq. 5.3:
  `w_j=(-1)^j δ_j`, `δ_j=1/2` at `j=0,n`), `:143,:175` (interval invariance
  of the weights).
- `src/BVP.jl:525-547` — `_chebyshev_D1` (the verbatim source of `_cheb_D1`).
- `src/BVP.jl:279-297` — `D2 = D1*D1`, affine scale, and the `bc_col`
  boundary-contribution idiom.
- `src/BVP.jl:560-585` — `_barycentric_eval` (the verbatim source of `_bary`).
- `docs/adr/0024-laplace-harmonic-extension.md:133-142` — the ADR statement
  of `L = kron(D²x[int,int], I_y) + kron(I_x, D²y[int,int])` and the
  BC-into-RHS split.
- `test/laplace2d_test.jl` — closed-form harmonic oracles O1/O2/O3, the
  non-square-rectangle test L2D.3, and the L2D.8 mutation-proof block
  toggling `swap_kron` and `drop_scale`.

## Findings

### [LOW] Top-docstring states the Kronecker factor order swapped vs the (correct) code

**Location**: `src/Laplace2D.jl:43` (top docstring) vs `src/Laplace2D.jl:215-216`
(code) and `src/Laplace2D.jl:212` (inline comment).

**Ground truth (cited)**: The discrete Laplacian must satisfy
`L·vec(Φ) = vec(D2x·Φ + Φ·D2yᵀ)`. For column-major `vec` the Kronecker
identity `vec(A·M·B) = (Bᵀ⊗A)·vec(M)` gives
`vec(D2x·Φ) = kron(I_niy, D2x_ii)·vec(Φ)` and
`vec(Φ·D2y_iiᵀ) = kron(D2y_ii, I_nix)·vec(Φ)`
(Trefethen Kronecker-sum Laplacian,
`references/markdown/TrefethenSMIM_2000_book/TrefethenSMIM_2000_book.md:747,754-760`;
ADR-0024:133-142).

**Code behavior**: The code at lines 215-216 builds
`L = kron(I_niy, D2x_ii) + kron(D2y_ii, I_nix)` — this is **correct**.
The inline comment at line 212 (`kron(D²y_ii, I_x) + kron(I_y, D²x_ii)`)
matches the code. But the **top-of-file docstring at line 43** writes
`L = kron(D²x_ii, I_y) + kron(I_x, D²y_ii)`, which is the *swapped* order
(the `D2x` factor placed first). The L2D.8 mutation-proof
(`test/laplace2d_test.jl:172-174,197-201`) confirms that the swapped order
is the WRONG operator and breaks the oracle. So the docstring describes a
known-wrong operator; the code is right.

**Mechanism**: Documentation-only. No runtime effect — the executed code
is correct. The risk is purely to a future maintainer who edits the
assembly to "match the docstring" and silently introduces the
grid-transpose bug. Not a source of the reported intermittent
discontinuity.

**Intermittent?**: No (doc only; no runtime path).

**Confidence**: 0.97 that the docstring/code disagree as stated; 0.97 that
the *code* is the correct order (cross-checked by the Kronecker identity,
the inline comment, `PadeTaylor.jl:199`, and the L2D.8 mutation-proof).

### [LOW] Coincident-node guard in `_bary` uses an absolute `eps(T)` tolerance that does not scale with the (mapped) interval

**Location**: `src/Laplace2D.jl:328-331` (`atol = eps(T)`; the
`abs(x_star - x_nodes[j]) ≤ atol` short-circuit).

**Ground truth (cited)**: Berrut-Trefethen 2004
(`references/markdown/BerrutTrefethen2004_barycentric_SIAMReview/...md:179`,
section-7 remark) handles the `x = x_j` case by a guard; the formula
itself (eq. 4.2, `:130`) is numerically stable away from nodes, and the
weights are interval-invariant (`:143,:175`). The reference does not
mandate a relative tolerance; the standard `chebint`/Chebfun guard is a
near-exact-equality test.

**Code behavior**: On a wide affine-mapped grid (e.g. `Re w ∈ [log r₀,
log r₁]` with large magnitudes), `eps(T)` is an absolute threshold that
does not scale with node spacing. The only consequences are: (a) a query
point extremely close to (but not equal to) a node is interpolated by the
stable barycentric ratio instead of being snapped — which is correct and
accurate; (b) an exactly-equal node is caught. There is **no** branch
where the guard produces a wrong value: the barycentric ratio is the
correct limit. Matches `BVP._barycentric_eval:565` exactly.

**Mechanism**: No discontinuity mechanism. The barycentric ratio is
continuous through the near-node region; the guard only avoids a literal
0/0 at the node. Listed for completeness, not as a defect.

**Intermittent?**: No.

**Confidence**: 0.9 that this is not a bug; flagged only because an
absolute `eps` tolerance on mapped data is a recurring smell.

## Areas verified correct

- **Chebyshev grid + descending convention** (`Laplace2D.jl:179-183`):
  `tx[j] = cos(j·π/Nx)` for `j∈0:Nx` ⇒ `tx[1]=+1`, `tx[Nx+1]=-1`,
  matching `cheb.m` `x = cos(pi*(0:N)/N)'`
  (`TrefethenSMIM_2000_book.md:447`). Affine map `xs = (a+b)/2 +
  (b-a)/2·t` ⇒ `xs[1]=b` (hi), `xs[Nx+1]=a` (lo). Self-consistent and
  matched by the line-204 comment.

- **`_cheb_D1` first-derivative matrix** (`Laplace2D.jl:294-316`): byte-for-byte
  identical to `BVP._chebyshev_D1` (`BVP.jl:525-547`). Endpoint factors
  `c=2` at `i=1,np1`; off-diagonal `(c_i/c_j)·(-1)^{i+j}/(t_i-t_j)`;
  diagonal by negative-row-sum (eq. 6.6,
  `TrefethenSMIM_2000_book.md:454-459`). The `iseven(i+j)` sign matches
  `cheb.m`'s `c = [...].*(-1).^(0:N)'` after the 0-index→1-index shift
  (the 0-indexed `(-1)^(i+j)` and the 1-indexed `(-1)^(i+j)` agree because
  the exponent differs by the even constant 2). Verified, no sign error,
  no off-by-one.

- **Affine D² scaling** (`Laplace2D.jl:188-189`): `D2 = (D1·D1)·(2/(b-a))²`.
  `dt/dx = 2/(b-a)`, second derivative ⇒ squared chain-rule factor; matches
  `BVP.jl:294` (`scale = (z_b-z_a)²/4`, the reciprocal placement). The
  square is always positive, so no sign dependence on `b<a`/`d<c`. Pinned
  by the non-square test L2D.3 (`[-2,3]×[-1,4]`) and by L2D.8 MUTATION 1
  (`drop_scale`), which bites only because the widths 5 vs 2 are unequal.

- **Kronecker-sum interior Laplacian** (`Laplace2D.jl:213-216`):
  `L = kron(I_niy, D2x_ii) + kron(D2y_ii, I_nix)` with `D2x_ii =
  D2x[2:Nx,2:Nx]`, `D2y_ii = D2y[2:Ny,2:Ny]`. Re-derived via
  `vec(A·M·B)=(Bᵀ⊗A)vec(M)` for column-major `vec` (x fast); both terms
  reproduce `φ_xx` and `φ_yy` exactly. Cross-checked against the Program 23
  Kronecker sum (`TrefethenSMIM_2000_book.md:760`, modulo the `-Δ`
  eigen-sign that is irrelevant for `Δφ=0`) and the L2D.8 `swap_kron`
  mutation-proof. Asymmetric oracles O1 (`x²−y²`) and O3 (`x³−3xy²`) on a
  non-square grid would expose any x↔y transpose; tests are GREEN.

- **Interior/boundary index split** (`Laplace2D.jl:192-193`): `ix=iy=2:N`
  interior, `bx=by=[1,N+1]` boundary, matching MATLAB `2:N` interior
  (`TrefethenSMIM_2000_book.md:758-759`) under the 1-index shift. `nix=Nx-1`,
  `niy=Ny-1`. No off-by-one.

- **Boundary trace assignment** (`Laplace2D.jl:200-208`): `bc_x_lo`(x=a)
  → `phi[Nx+1,:]`; `bc_x_hi`(x=b) → `phi[1,:]`; `bc_y_lo`(y=c) →
  `phi[:,Ny+1]`; `bc_y_hi`(y=d) → `phi[:,1]`. Each lands on the edge whose
  node value equals the stated coordinate (`xs[1]=b`, `xs[Nx+1]=a`, etc.).
  Corner nodes are overwritten by the y-edges after the x-edges but are
  never read into `rhs` (the RHS loops use interior-only columns/rows), so
  the overwrite is harmless. Verified.

- **Dirichlet-into-RHS absorption + sign** (`Laplace2D.jl:222-233`):
  `rhs[:,col] += D2x[ix,bx]·phi[bx,j]` (φ_xx boundary part, for each
  interior y-column `j`) and `rhs[row,:] += D2y[iy,by]·phi[i,by]` (φ_yy
  boundary part, for each interior x-row `i`). The `bx`/`by` column order
  of `D2x[ix,bx]` is paired consistently with `phi[bx,j]`/`phi[i,by]` —
  no index mismatch. The equation `L·vec(Φ) + vec(rhs) = 0` is solved as
  `φ_int = L \ (-vec(rhs))` (line 233), so the minus sign is correct. This
  is the 2D analogue of `BVP.jl:297` `bc_col`. Closed-form oracles confirm.

- **`vec`/`reshape` round-trip** (`Laplace2D.jl:233-234`): `rhs` is
  `nix×niy`, `vec` column-major (x fast) matches `L`'s block structure;
  `reshape(φ_int, nix, niy)` is the exact inverse, written back to
  `phi[ix,iy]`. Consistent.

- **`_sample` expected lengths + fail-fast** (`Laplace2D.jl:200-203,246-256`):
  x-edges expect `niy+2 = Ny+1` (vary over y); y-edges expect `nix+2 = Nx+1`
  (vary over x). Length-mismatch and `Nx/Ny<3` throw `ArgumentError` with a
  suggestion (Rule 1). Pinned by L2D.7.

- **`_bary` tensor-barycentric callable** (`Laplace2D.jl:270-342`):
  second-form barycentric (eq. 4.2,
  `BerrutTrefethen2004...md:130`) with Chebyshev-2 weights
  `w_j=(-1)^{j-1}·(1/2 at endpoints)` (eq. 5.3, `...md:171`); identical to
  `BVP._barycentric_eval:560-585`. Interval-invariance of the weights
  (`...md:143,:175`) justifies using `[-1,1]` weights on the mapped nodes.
  Tensor order (interpolate each `phi[:,j]` in x, then the resulting column
  in y, lines 274-275) is consistent with `phi[i,j]=φ(xs[i],ys[j])`. Pinned
  by L2D.4 (off-node) including the exact-node short-circuit assertion.

- **No intermittency surface in this module**: `rhs` is freshly zeroed
  (line 222); the two accumulation loops write disjoint contributions of a
  read-only `phi`; there is no RNG, no in-place buffer aliasing, no
  order-dependent mutation, and the single `L \ (-vec(rhs))` solve is
  deterministic. Any intermittent discontinuity in the assembled figure is
  therefore unlikely to originate inside `Laplace2D.jl`; the conformal
  `w=log x` map and tritronquée boundary-data supply (F3, outside this
  module) and the per-grid-point majority vote are the more plausible
  homes.
