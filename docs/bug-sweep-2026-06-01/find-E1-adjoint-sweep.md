# Bug sweep E1 — conjugate-transpose vs plain-transpose on complex data

Date: 2026-06-01
Auditor: subagent "find-E1-adjoint-sweep" (read-only)
Assignment: WHOLE `src/` — adjoint vs transpose on complex data (the
generalized prime suspect for intermittent discontinuities).

## Area

Every `src/*.jl` file, scanned for the postfix-quote adjoint operator,
`adjoint(`, `transpose(`, `permutedims(`, dot-quote `.'`, `dot(` / `⋅`,
backslash linear solves, and SVD/QR/eigen calls — with each occurrence
on data that *can be complex* (Padé coefficients, Toeplitz/Hermite-Padé
blocks, Chebyshev matrices applied to complex BVP solutions, the
Laplace tensor system, the branch-cut Cramer solve) adjudicated for
whether the math needs a PLAIN transpose or a CONJUGATE transpose, and
whether the code matches.

## References checked

- `external/chebfun/padeapprox.m:69` — special-case r=0 test uses the
  inf-norm: `norm(c(1:m+1),inf) <= tol*norm(c,inf)`.
- `external/chebfun/padeapprox.m:89-93` — Toeplitz `Z`, `C = Z(m+2:m+n+1,:)`,
  `rho = sum(svd(C) > ts)`.
- `external/chebfun/padeapprox.m:106-120` — `[U,S,V]=svd(C,0)`,
  `b=V(:,n+1)`, reweighting `D=diag(abs(b)+sqrt(eps))`,
  **`[Q,R]=qr((C*D)')` with comment "until July 2018 there was an
  erroneous `.'` here"**, `b=D*Q(:,n+1)`, `b=b/norm(b)`,
  `a = Z(1:m+1,1:n+1)*b`.
- `external/chebfun/padeapprox.m:122-138` — leading/trailing-zero trim
  with `tol` (for b) and `ts` (for a), normalise to `b(1)=1`.
- `external/chebfun/padeapprox.m:150` — `poles = roots(b(end:-1:1))`.
- `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/GGT2013_robust_pade_via_SVD_SIREV55.md:102`
  — "b must be a (right) null vector of the n×(n+1) matrix".
- GGT2013 md:117-122 — definition of `C̃` (Toeplitz) and `C` (square,
  first column deleted).
- GGT2013 md:128 — SVD `C̃ = UΣV*`.
- GGT2013 md:132 — "the final column of V provides … the null vector b".
- GGT2013 md:226-232 — Algorithm 2: `τ = tol·‖c‖₂`; rank = #σ **> τ**;
  null right singular vector b; trim with `tol`/`τ`.
- GGT2013 md:236 — three "reweighting" lines go beyond Algorithm 2:
  QR of a column-reweighted matrix.
- **GGT2013 md:279 — the PAPER'S Figure-1 listing shows the OLD line
  `[Q,R] = qr((C*D).');` with the DOT-QUOTE (plain transpose)** — i.e.
  the erroneous pre-July-2018 form that `padeapprox.m:113` later fixed.
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:330-336`
  — classical Toeplitz system eqs. (5.4)/(5.5).

## Findings

### [LOW] Comment/docstring glyph says "ᵀ" (transpose) where the code correctly uses `adjoint` — `SharedPade.jl`

- **Location:** `src/SharedPade.jl:234` (code), and the doc references
  `src/SharedPade.jl:59` ("`QR((A_full·D)ᵀ)`") and the inline comment on
  line 234 ("`# (A_full·D)ᵀ in MATLAB`").
- **Ground truth (cited):** `external/chebfun/padeapprox.m:113`
  `[Q,R] = qr((C*D)')` — the postfix-quote **conjugate transpose**, with
  the comment "until July 2018 there was an erroneous `.'` here". The
  scalar sibling `src/RobustPade.jl:443` writes the same operation and
  documents it correctly as `# (C*D)' in MATLAB`.
- **Code behavior:** The CODE is correct — `F = qr(adjoint(A_full * D))`
  uses `adjoint` (conjugate transpose), matching the corrected chebfun
  source. Only the *prose* (`ᵀ` glyph, twice) is imprecise: `ᵀ` denotes
  the plain transpose, contradicting the actual `adjoint`.
- **Mechanism (why this matters even though code is correct):** This is
  the exact documented historical-bug class. A future maintainer reading
  the `ᵀ` glyph and "matching the comment to the code" could change
  `adjoint` → `transpose`, silently introducing the July-2018 bug on
  complex Hermite-Padé blocks (vector solver). That WOULD produce
  intermittent discontinuities (only complex jets, only certain matrices).
  No bug exists today; the risk is latent in the documentation.
- **Intermittent?** Not currently (no runtime defect). The *latent* bug
  it could seed would be intermittent (complex-data only).
- **Confidence:** 0.9 that the glyph is a doc inaccuracy and the code is
  mathematically correct (verified against padeapprox.m:113 and the
  null-space derivation below).

### [LOW] LinAlg docstring says null vector is `transpose(Vt[end,:])` (no conjugate) — harmless in current callers, latent trap

- **Location:** `src/LinAlg.jl:42` (docstring prose only — not code).
- **Ground truth (cited):** GGT2013 md:128 `C̃ = UΣV*`, md:132 "the final
  column of V provides … the null vector b". With Julia's `svd` returning
  `Vt = V*` (conjugate transpose of V), the last right singular vector is
  `V[:,end] = conj(Vt[end,:])` = `adjoint`/`Vt[end,:]'`, **not**
  `transpose(Vt[end,:])` (which omits the conjugate on complex data).
- **Code behavior:** Both callers — `RobustPade.jl:435` and
  `SharedPade.jl:219` — read `b_init = Vector{T}(Vt[end, :])` (the raw
  row, no conjugate). Crucially, `b_init` is then used ONLY to form the
  magnitude-diagonal `D = Diagonal([abs(bk) + √ε for bk in b_init])`
  (`RobustPade.jl:442`, `SharedPade.jl:233`). `abs(·)` is
  conjugation-invariant, so the missing conjugate has NO effect on the
  result: the final `b` is recomputed from the QR null vector
  (`b = D * F.Q[:, n+1]`), never from `b_init` directly.
- **Mechanism:** No current defect. Latent trap: if a future change ever
  used `b_init` as the denominator directly (phase-sensitive), the missing
  conjugate would corrupt the complex denominator → intermittent
  discontinuities on complex data.
- **Intermittent?** Not currently. Latent form would be intermittent.
- **Confidence:** 0.88 (verified b_init feeds only `abs`; doc claim of
  `transpose` is the imprecise part).

## Areas verified correct

- **`src/RobustPade.jl:443` — `F = qr(adjoint(C * D))`. CORRECT.** This
  is the load-bearing QR-reweighting null-vector computation. It uses
  `adjoint` (conjugate transpose), matching the CORRECTED chebfun source
  `external/chebfun/padeapprox.m:113` (`qr((C*D)')`) and NOT the erroneous
  plain-transpose form preserved in the GGT2013 paper listing
  (`GGT2013…md:279` `qr((C*D).')`). Derivation: `Q(:,n+1)` from
  `qr((C*D)')` is orthogonal to `range((C*D)')` = `null(C*D)`, so it lies
  in the (reweighted) denominator null space; using a plain transpose on
  complex `C` would instead target `conj(range((C*D)'))`, corrupting the
  complex denominator. The codebase got the exact historical-bug case
  right. The inline comment `# (C*D)' in MATLAB` is accurate.
- **`src/SharedPade.jl:234` — `F = qr(adjoint(A_full * D))`. CORRECT
  (code).** Same `adjoint` operation as the scalar path; the `d=1` block
  reduces bit-for-bit to the scalar `C̃` (module docstring lines 46-53),
  so the same null-space argument applies. Only the surrounding prose
  glyph `ᵀ` is imprecise (see LOW finding above).
- **`b_init` magnitude-only usage** (`RobustPade.jl:435,442`;
  `SharedPade.jl:219,233`) — fed exclusively into `abs(·)` for the
  reweighting diagonal `D`, so the SVD `Vt[end,:]` conjugation convention
  is immaterial to the result. Matches padeapprox.m:112 `diag(abs(b)+…)`.
- **`src/LinAlg.jl:75-108` — `pade_svd` backends** return `(F.U, F.S, F.Vt)`
  unchanged for LAPACK and GenericLinearAlgebra; no transpose/adjoint
  manipulation of the factors. Consistent across real and complex element
  types.
- **`src/RobustPade.jl:413-424` — Toeplitz `Z`, `C = Z[m+2:m+n+1,:]`,
  rank via `count(s -> s > ts_typed, S)`.** Matches padeapprox.m:89-93
  and GGT2013 md:228 (strict `> τ`). Strict `>` is correct vs the
  reference; no off-by-one in the row slice (`m+2 : m+n+1`).
- **`src/RobustPade.jl:259-318` — `classical_pade_diagonal`** Toeplitz
  build `T_mat[i,j] = cv[m+i-j+1]` = `c_{m+i-j}` and RHS `-c_{m+i}` match
  FW2011 eq. (5.4) (FW2011…md:330-336); numerator
  `a_k = Σ_{j=0..min(k,m)} b[j]·c_{k-j}` matches eq. (5.5). Solve is a
  plain `F \ rhs` (no transpose). Correct.
- **`src/RobustPade.jl:459-498` — `_trim_and_normalise`** uses `tol` for
  b leading/trailing and `ts` for a trailing, normalise to `b[1]=1`,
  matching padeapprox.m:122-138 and GGT2013 md:231-232.
- **`src/BVP.jl:525-547` / `src/VectorBVP.jl:372` — `_chebyshev_D1`.**
  Off-diagonal `(c_i/c_j)·(-1)^{i+j}/(t_i-t_j)` with `iseven(i+j)`
  (1-indexed) = `(-1)^{(i-1)+(j-1)}` (0-indexed); diagonal by
  negative-row-sum. Built in natural orientation, applied as `D1 * u`
  (no transpose). Real matrix on complex data → plain matrix-multiply,
  which is what the code does. Matches Trefethen `cheb.m` convention.
- **`src/BVP.jl:280-325` — Newton solve** `J = D2_ii − scale·Diag(∂f/∂u)`,
  `R = D2_ii*u_int + bc_col − scale·F`, `J \ R`. Plain solve; no
  transpose; `D2 = D1*D1`. No conjugation hazard.
- **`src/VectorBVP.jl:221-271` — `Dop = kron(D1, I_d)`** applied as
  `Dop*Y`; Jacobian built in the same orientation; solved via `_solve`.
  No transpose/adjoint. (The `kron` left/right ordering is an
  ordering-consistency question, not a conjugation question, and is
  internally consistent between residual and Jacobian.)
- **`src/Laplace2D.jl:215-234` — tensor Laplacian** `L = kron(I,D2x) +
  kron(D2y,I)`, solved `L \ (-vec(rhs))`. Real-valued harmonic problem;
  plain solve; `vec`/`reshape` column-major consistent. The `ᵀ` in the
  line-221 comment is descriptive of a length-2 vector slice, not a code
  transpose.
- **`src/BranchTracker.jl:112-121` — `segment_crosses_cut` Cramer's
  rule.** `det = imag(d·conj(u))`, `t = imag(δ·conj(u))/det`,
  `s = imag(δ·conj(d))/det`. Verified: solving `t·d − s·u = δ` with the
  2D cross product `a×b = imag(a·conj(b))` gives exactly these formulas
  (`t = (δ×u)/(d×u)`, `s = (δ×d)/(d×u)`). The `conj` is used consistently
  as the cross-product/inner-product convention, not as a stray adjoint.
  Correct.
- **No naked postfix-quote (`'`) adjoint operator exists in code
  anywhere in `src/`** — every `'` matched by grep is a prose apostrophe
  (`u'`, `f_j'`, possessive). The codebase uniformly uses the explicit
  `adjoint(...)` function (2 sites only). **No `.'` (dot-quote) anywhere.**
  **No `transpose(`/`permutedims(`/`dot(`/`⋅`/`eigen`/`pinv`/`nullspace`
  in code.** This makes the complex-data linear-algebra surface small and
  fully auditable: it is exactly the two `adjoint(qr(...))` sites, both
  verified correct.
- **`src/PoleField.jl:139-143` / `src/StepControl.jl:249` — root finding**
  via `Polynomials.roots(Polynomial(P.b))`. `Polynomial` takes coeffs
  low-to-high (b[1]=b_0), the correct Polynomials.jl convention; no
  MATLAB-style `b(end:-1:1)` reversal needed (the reversal in
  padeapprox.m:150 is a MATLAB `polyval`/`roots` artifact). No transpose.
