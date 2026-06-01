# Bug sweep A1 — RobustPade + LinAlg (GGT 2013 Algorithm 2 / padeapprox.m port)

Date: 2026-06-01
Auditor: parallel agent A1 (read-only)
Scope: `src/RobustPade.jl`, `src/LinAlg.jl`

## Area

The GGT 2013 Algorithm 2 SVD-based Padé path (`robust_pade(...; method=:svd)`)
and the FW 2011 §5.1.4 classical Toeplitz-`\` path
(`classical_pade_diagonal`), plus the SVD dispatcher `LinAlg.pade_svd`.
The maintainer flagged this as the prime suspect for a transcription bug
producing intermittent discontinuities in computed solutions.

Verdict up front: **I found no transcription bug in the core RobustPade/LinAlg
algorithm.** Every flagged suspicion item (the historical adjoint-vs-transpose
class, the `tol`/`ts` threshold split, the 1-based `C` extraction shift, the
null-vector index, the pole/coefficient order, the normalisation) was diffed
line-by-line against `padeapprox.m` and the GGT 2013 markdown and matches the
ground truth. The two genuinely-low-confidence observations below are
behavioural notes about *where* the SVD path runs (rarely, on fallback), not
demonstrated mismatches.

## References checked

- `external/chebfun/padeapprox.m:62,66,69,76-77,89,92-93,100-101,106,109,112-113,116-117,120,123,126-127,130,134,137-138,146,150` — canonical reference impl.
- `references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/GGT2013_robust_pade_via_SVD_SIREV55.md:204` (eq. 2.5 `||b||=1` normalisation), `:230-235` (eq. 2.10/2.11 `C̃`/`C`, null vector = last col of V), `:248-258` (Algorithm 2 steps 2,4,5,7,8 with the `tol` vs `τ=tol·||c||₂` split), `:260-285` (Figure 1 embedded code).
- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:330-336` (eqs. 5.4, 5.5), `:346` (singular fallback), `:350` ("method of choice ... Toeplitz approach + backslash").
- `docs/adr/0002-bigfloat-svd-via-genericlinalg.md:11-19,33-39,98-104` — SVD dispatch + relative-accuracy guarantee.
- `docs/adr/0005-classical-pade-default-at-float64.md:11-28,56-85` — element-type dispatch.
- `~/.julia/packages/GenericLinearAlgebra/WvuVM/src/svd.jl:327-380` (bidiagonalize wide case), `:638-668` (`svd!` `full=true` builds `(n+1)×(n+1)` Vᴴ), `:666,672-676` (default `atol=rtol=0` zeroes nothing) — confirms `pade_svd(C; full=true)` returns a full `(n+1)×(n+1)` Vt so `Vt[end,:]` is the genuine null vector at `BigFloat`.

## Findings

### [LOW] Hypothesis #1 (adjoint vs transpose in the QR-reweighting) is VERIFIED CORRECT — not a bug

- Location: `src/RobustPade.jl:443` (and the mirror `src/SharedPade.jl:234`).
- Ground truth: `external/chebfun/padeapprox.m:113`
  `[Q, R] = qr((C*D)');           % until July 2018 there was an erroneous .' here`.
  In MATLAB the postfix quote `'` is the **conjugate** transpose; the comment
  records that the plain transpose `.'` was the historical bug for complex `c`.
- Code behavior: the Julia port writes `F = qr(adjoint(C * D))`. Julia's
  `adjoint` is the conjugate transpose — the **correct** operator, matching the
  fixed MATLAB. `D = Diagonal([abs(bk) + eps_T ...])` is real-diagonal, so it
  commutes with conjugation cleanly; `b = D * F.Q[:, n + 1]` then satisfies
  `C·b = (C·D)·Q[:,n+1] = 0` because `Q[:,n+1]` spans the left null space of
  `(C·D)* = adjoint(C·D)`. Algebra checks out.
- Mechanism: N/A — no defect. (If it *had* been `transpose`, it would have bitten
  only complex-coefficient near-rank-deficient blocks, i.e. intermittently. It
  does not.)
- Intermittent: N/A.
- Confidence: 0.97 that this is correct (demonstrated equality with the cited
  fixed reference line).

### [LOW] `b_init = Vt[end, :]` is `conj(V[:,end])` for complex data, but only `abs()` of it is used — harmless

- Location: `src/RobustPade.jl:435,442`.
- Ground truth: `padeapprox.m:109` `b = V(:,n+1)`; GGT 2013 md:231-233 "the final
  column of V provides ... the null vector b".
- Code behavior: `pade_svd` returns `F.Vt = V*`, so `Vt[end,:] = conj(V[:,end])`
  for complex input — a conjugated null vector relative to MATLAB. **However**
  `b_init` is consumed only as `D = Diagonal([abs(bk)+eps_T for bk in b_init])`,
  and `abs(conj(x)) == abs(x)`. The actual returned `b` comes from the QR step
  (line 446), not from `b_init`. So the conjugation cannot affect the result.
- Mechanism: would only matter if `b_init` were used un-conjugated downstream; it
  is not.
- Intermittent: no.
- Confidence: 0.9 this is harmless as written.

### [LOW] The SVD path is exercised only on rare fallbacks for the dominant `Complex{Float64}` walker — narrows where any latent SVD-path defect could hide, but is not itself a defect

- Location: dispatch at `src/RobustPade.jl:368-375` + defaults `:176-180`; hot-path
  caller `src/PadeStepper.jl:314,465` (`robust_pade(c, m, n)` with no `method`).
- Ground truth: `docs/adr/0005-classical-pade-default-at-float64.md:15-18` —
  `Float64`/`Complex{Float64}` default to `:classical`.
- Code behavior: the path-network walker uses `Complex{Float64}` and diagonal
  `m=n`, so it runs `classical_pade_diagonal` (LU Toeplitz), **not** the SVD/QR
  path. The SVD path runs only when (a) `m≠n`, (b) the Toeplitz is exactly
  singular and the try/catch at `:369-374` falls through, (c) the user forces
  `:svd`, or (d) the element type is `BigFloat`/`Complex{BigFloat}`. I verified
  the classical path matches FW eqs. (5.4)/(5.5) exactly (see "verified correct"
  below), and the SVD path matches padeapprox.m exactly. The observation only
  matters as guidance: if an intermittent discontinuity is real and lives in
  Padé arithmetic, the *classical* path is the one running on essentially every
  step, and it is a single LU solve — deterministic, no rank-counting, no
  tie-break. An intermittent symptom is therefore more likely to originate
  outside RobustPade (step/direction selection, nearest-node lookup, branch
  tracking) than inside it.
- Mechanism: not a defect; routing note only. The singular-fallback at
  `:369-374` catches `SingularException` only — if a near-singular (not exactly
  singular) Toeplitz produces an `Inf`/`NaN` from `\` without `lu` reporting
  `!issuccess`, it would propagate; but that is FW's documented accept
  (`FW2011 md:346`) and matches `padeapprox.m`'s reliance on backslash too.
- Intermittent: the *routing* is data-dependent (singular fallback fires only on
  some inputs), but no wrong value is produced.
- Confidence: 0.6 as a routing observation; 0.05 that a bug lives here.

## Areas verified correct

Each item below was diffed against the cited reference line and matches.

1. **Toeplitz `Z` construction.** `src/RobustPade.jl:196-202`
   `Z[i,j]=c[i-j+1]` (i≥j) vs `padeapprox.m:89`
   `toeplitz(col(1:m+n+1),row(1:n+1))` with `col=c`, `row=[c(1) zeros(1,n)]`,
   which yields exactly the lower-triangular Toeplitz `Z[i,j]=c_{i-j}`. Size
   `(m+n+1)×(n+1)`. Matches.

2. **`C̃` extraction / 1-based index shift.** `src/RobustPade.jl:414`
   `C = Z[(m+2):(m+n+1), :]` vs `padeapprox.m:92` `C = Z(m+2:m+n+1,:)` — byte-
   identical 1-based slice. `C[1,1]=Z[m+2,1]=c_{m+1}` matches GGT eq. (2.10)
   (md:230). `C` is `n×(n+1)`.

3. **Special case r≡0.** `src/RobustPade.jl:396-399`
   `maximum(abs, cv[1:m+1]) ≤ tol·maximum(abs, cv)` vs `padeapprox.m:69`
   `norm(c(1:m+1),inf) <= tol*norm(c,inf)` (and GGT Algorithm 2 step 2,
   md:249). `‖·‖∞ = maximum(abs,·)`; uses `tol` (not `ts`); uses `≤`. Matches,
   including the `typemin(Int)` = `μ=-∞` convention (GGT §2, md:204 "r=0 ⇒
   μ=−∞, ν=0").

4. **Rank-count threshold.** `src/RobustPade.jl:418`
   `ρ = count(s -> s > ts_typed, S)` with `ts = tol·‖c‖₂` vs `padeapprox.m:93`
   `rho = sum(svd(C) > ts)` and GGT step 4 (md:251, "greater than τ=tol·‖c‖₂").
   Strictly-greater on both sides; `ts` (absolute), not `tol`. Matches.

5. **Diagonal hopping update.** `src/RobustPade.jl:423-424`
   `m -= n-ρ; n = ρ` vs `padeapprox.m:100-101` `m = m-(n-rho); n = rho` and
   GGT step 5 (md:252). Matches.

6. **Null-vector / QR-reweighting block.** `src/RobustPade.jl:434-447` vs
   `padeapprox.m:106-117`: full SVD → last right singular vector;
   `D=diag(|b|+√eps)`; `qr(adjoint(C·D))`; `b=D·Q[:,n+1]`; `b/=‖b‖`. Matches,
   including the `adjoint` (conjugate-transpose) fix (see finding above).

7. **Numerator recovery.** `src/RobustPade.jl:450`
   `a = Z[1:m+1, 1:n+1]·b` vs `padeapprox.m:120` `a = Z(1:m+1,1:n+1)*b`
   (GGT step 6, md:253). `Z` is the rebuilt-on-final-(m,n) matrix in both.
   Matches.

8. **Leading `z^λ` cancellation.** `src/RobustPade.jl:463-476`
   `lam = findfirst(|b|>tol)-1`, drop `b[1:lam]` and `a[1:lam]` vs
   `padeapprox.m:123,126-127` `lam = find(abs(b)>tol,1)-1; b=b(lam+1:end);
   a=a(lam+1:end)` (GGT step 7, md:254). Uses `tol`. Matches. The all-below-tol
   branch throws loudly (Rule 1) instead of MATLAB's undefined `find([])`.

9. **Trailing-`b` trim uses `tol`; trailing-`a` trim uses `ts`.**
   `src/RobustPade.jl:479` `findlast(|b|>tol)` and `:484` `findlast(|a|>ts)`
   vs `padeapprox.m:130` (`>tol`) and `:134` (`>ts`), and GGT step 8 (md:255,
   "`|b...|≤tol` remove last of q; `|a...|≤τ` remove last of p"). The two
   thresholds are correctly distinct — this was a flagged easy-to-mis-port
   item and is ported correctly.

10. **Normalisation `b₀=1`.** `src/RobustPade.jl:493-495` `a/=b1; b/=b1` vs
    `padeapprox.m:137-138` `a=a/b(1); b=b/b(1)` (GGT step 9, md:256). The
    intermediate `‖b‖₂=1` (line 447 / padeapprox:117, GGT eq. 2.5 md:204) is
    correctly washed out by the final `/b(1)`. Evaluation in
    `PadeStepper._evaluate_pade` (`src/PadeStepper.jl:367-382`) is consistent
    with `b[1]=1`. Matches.

11. **`n==0` branch.** `src/RobustPade.jl:407-411` returns `a=cv[1:m+1]`,
    `b=[1]`, then trims trailing `a` with `ts` and normalises — equivalent to
    `padeapprox.m:82-86` `break` skipping the `n>0` block, then line 134/137-138.
    Verified equivalent (the leading-zero and trailing-`b` trims in
    `_trim_and_normalise` are no-ops on `b=[1]`). Matches.

12. **Coefficient padding.** `src/RobustPade.jl:378-381` pad/truncate `cv` to
    `m+n+1` vs `padeapprox.m:62-63` `c=[f(:);zeros(...)]; c=c(1:m+n+1)`. `ts`
    computed once on the original-(m,n) `cv` (`:384-385`) as MATLAB does.
    Matches.

13. **Classical Padé FW eqs. (5.4)/(5.5).** `src/RobustPade.jl:277-314`.
    `T_mat[i,j]=cv[m+i-j+1]=c_{m+i-j}` matches FW (5.4) (`FW md:330`);
    `rhs[i]=-cv[m+i+1]=-c_{m+i}` matches the `-[c_{v+1}..c_{2v}]` RHS;
    numerator `a_k=Σ_{j=0..k} b_j c_{k-j}` (with `b_0=1`) is algebraically
    identical to FW (5.5) `a_0=c_0`, `a_k=Σ_{j=1..k} c_{k-j}b_j + c_k`
    (`FW md:336`). `lu(...; check=false)`+`issuccess` singular detection maps to
    FW's "outright singular" case (`FW md:346`). Matches. (Mutation P1
    sign-flip on RHS is covered by CP.1 tests per ADR-0005:147-150.)

14. **`LinAlg.pade_svd` dispatch + `full` semantics.** `src/LinAlg.jl:75-109`.
    LAPACK path for `Float32/64`+complex (matches `padeapprox.m:93,106` DK SVD);
    `GenericLinearAlgebra.svd` one-sided-Jacobi for `BigFloat`/generic +
    complex (ADR-0002:17,33-39). Confirmed by reading
    `GenericLinearAlgebra/src/svd.jl:638-668`: for the wide `n×(n+1)` input with
    `full=true`, `Vᴴ` is built `(n+1)×(n+1)`, so `Vt[end,:]` is the true null
    vector — matching the LinAlg docstring (`src/LinAlg.jl:38-44`) and the
    MATLAB `V(:,n+1)`. Default `atol=rtol=0` zeroes no singular values
    (`svd.jl:666,672-676`), so it does not perturb the `s>ts` rank count.
    The relative-accuracy guarantee that makes the `σᵢ<tol·‖c‖₂` threshold
    load-bearing at arb-prec (ADR-0002:33-39) is honoured by routing `BigFloat`
    to Jacobi. Matches.

15. **`default_tol`.** `src/RobustPade.jl:154-158`: `1e-14` at Float64 (matches
    `padeapprox.m:36` and GGT md:248 `tol=10⁻¹⁴`); `2^(-precision+10)` at higher
    precision (matches ADR-0002 rationale, ADR notes 10 bits slack). Matches.

16. **Poles ordering (informational).** `padeapprox.m:146,150` reverse `b`/`a`
    (`polyval(a(end:-1:1),z)`, `roots(b(end:-1:1))`) because MATLAB `polyval`
    is high-to-low while these coeff vectors are low-to-high. The Julia port
    does **not** reproduce the pole/residue branch in RobustPade (it returns a
    `PadeApproximant` and evaluates via low-to-high Horner in PadeStepper:367),
    so there is no order-reversal to mis-port here. Evaluation order is internally
    consistent (low-to-high in both `a` and `b`). Correct by construction.
