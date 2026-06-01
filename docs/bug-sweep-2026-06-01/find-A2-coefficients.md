# Bug sweep 2026-06-01 — A2: Coefficients + VectorCoefficients

Auditor focus: the Taylor / automatic-differentiation recurrence that turns
`y' = f(z, y)` (and `u'' = f(z, u, u')`) into Taylor coefficients.

## Area

Files audited (every relevant line read):

- `src/Coefficients.jl` — scalar `taylor_coefficients_1st` (1st order) and
  `taylor_coefficients_2nd` (2nd order).
- `src/VectorCoefficients.jl` — `vector_taylor_coefficients` (1st-order vector
  system, component-wise lift).

Special-focus items demanded by the assignment:
1. coefficient-vs-derivative (factorial) confusion;
2. the 2nd-order reduction seeding both `y(0)` and `y'(0)` and recurring both;
3. the AD bootstrap integration step `c_{k+1} = .../(k+1)`;
4. the `d = 1` vector reduction being bit-identical to the scalar routine;
5. in-place buffer aliasing / order-dependence (the intermittency vector).

## References checked (path:line)

- **FW 2011 §2.1.2 method (b)** —
  `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:96-107`
  ("at first, use only one term of (2.3), substitute into (2.4) and integrate,
  giving two correct terms. For each similar step, one gains one correct
  coefficient"). Eq. (2.3) `y(t+h)=c_0+c_1 h+...` and eq. (2.4)
  `dy(t+h)/dh = f(t+h, y(t+h))` are the integration relation the loop discretises.
- **FW 2011 §2.1.1** (line 88-94): `c_n = y^{(n)}(t)/n!` — the explicit statement
  that the stored object is the *normalised* derivative (the Taylor coefficient),
  not the raw derivative. Used to rule out factorial confusion.
- **Jorba & Zou 2005 eq. (5)** —
  `references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/JorbaZou2005_taylor_IVP_package_ExpMath14.md:276`:
  `a^{[j]}(t) = (1/j!) a^{(j)}(t)` — the "normalised j-th derivative" is exactly
  the Taylor coefficient.
- **Jorba & Zou 2005 recurrence** —
  `.../JorbaZou2005_taylor_IVP_package_ExpMath14.md:405-416`:
  `x^{[n+1]} = (1/(n+1)) u_2^{[n]}` with `u_2 = x'` the RHS, and the explicit
  remark "The factor 1/(n+1) comes from the definition given in equation (5)."
  Also line 437: `x_n(t) = Σ_j x_n^{[j]}(t_n) h^j` (coefficient = `x^{[j]}`).
- **TaylorIntegration.jl canonical `solcoeff!`** —
  `external/TaylorIntegration.jl/src/TaylorIntegration.jl:26-32`:
  `a[k] = b[k-1] / k`. The canonical Julia integration step; the loop body is a
  byte-for-byte match. Driven by `jetcoeffs!`
  (`external/TaylorIntegration.jl/src/integrator/jetcoeffs.jl:23-43` scalar,
  `:66-93` vector), which advances `ord = 0:order-1`, `ordnext = ord+1`.
- **TaylorSeries.jl Taylor1 indexing** —
  `external/TaylorSeries.jl/src/auxiliary.jl:107` `getindex(a::Taylor1, n) =
  a.coeffs[n+1]` and `:113` `setindex!(a::Taylor1{T}, x, n) = a.coeffs[n+1] = x`.
  Confirms `Taylor1` indexing is 0-based at the math level (maps to 1-based
  internal storage). Version pinned at `[0.21]` in `Project.toml:38`; external
  clone is `0.21.7` (`external/TaylorSeries.jl/Project.toml:3`).
- **TaylorSeries.jl Taylor1 constructor (zero-padding)** —
  `external/TaylorSeries.jl/src/constructors.jl:40-50`: `Taylor1(coeffs, order)`
  pads with `zero` up to `order+1`. Confirms the seeds `T[y0]`, `T[y0,y1]`, `T[y1]`
  are correctly zero-extended.
- **ADR-0001** — `docs/adr/0001-four-layer-architecture.md:13` (Coefficients is the
  jet layer; "A bug in the coefficient layer surfaces as wrong Taylor coefficients
  ... independently checkable against analytic derivatives", line 37-39).
- **Tests** — `test/coefficients_test.jl` (3.1.1 exp, 3.1.2 Bessel-ratio, 3.1.3
  Weierstrass-℘ 2nd-order oracle to rtol 1e-12, 3.1.4 BigFloat-256, mutation
  procedure lines 100-125); `test/vector_coefficients_test.jl` (VC.1.1 d=1
  bit-identical oracle, VC.1.2/1.3 decoupled+coupled closed forms, VC.1.4
  BigFloat-256, mutation procedure lines 135-155).

## Findings

No bugs found in this area. Every recurrence index, divisor, sign, and seed maps
exactly onto the cited ground truth. Detail below under "Areas verified correct".

(One sub-200-confidence observation about `@inbounds` is recorded so the parent
agent has it, but it is verified *safe*, not a bug — see VC-OK-5.)

## Areas verified correct

### CO-OK-1 — 1st-order integration step `y[k] = f_t[k-1] / k` is correct (confidence 0.97)

`src/Coefficients.jl:127-130`:

```julia
@inbounds for k in 1:order
    f_t = f(z, y)
    y[k] = f_t[k-1] / k
end
```

Ground truth: Jorba-Zou eq. (5) + recurrence (`JorbaZou2005...md:276,405-416`)
gives `x^{[n+1]} = (1/(n+1)) u_2^{[n]}` where `u_2 = x' = f`. With the loop index
`k = n+1`, that is `c_k = f_{k-1} / k` — exactly the code. Identical to the
canonical `solcoeff!` `a[k] = b[k-1]/k`
(`external/TaylorIntegration.jl/src/TaylorIntegration.jl:31`). Index map verified:
`f_t[k-1]` reads math-order `k-1` via `getindex(...,k-1)=coeffs[k]`
(`auxiliary.jl:107`); `y[k]` writes math-order `k` via `setindex!(...,k)=coeffs[k+1]`
(`auxiliary.jl:113`). **No factorial confusion**: the stored object is the
normalised derivative (= coefficient), so the only division needed is `/k`, which
is present and correct; there is no spurious `*k!` or `/k!`. Mutation A
(`test/coefficients_test.jl:106-110`, `/k → /(k+1)`) was confirmed to RED test
3.1.1. The exp / Bessel-ratio oracles (3.1.1, 3.1.2) pin the values.

### CO-OK-2 — independent-variable seed `z = z0 + Taylor1(T, order)` is correct (confidence 0.97)

`src/Coefficients.jl:120` (and `:177`, and `VectorCoefficients.jl:148`).
`Taylor1(T, order)` sets `coeffs[2] = one(T)` (`constructors.jl:89`), i.e. the
math-order-1 coefficient is 1, so `z = z0 + h`. Matches eq. (2.4)'s expansion
variable `h` (`FW2011...md:101`). `order ≥ 1` is guaranteed by the `order < 1`
guards, so the `order == 0` degenerate branch of `Taylor1(::Type,order)`
(`constructors.jl:86`) is never hit.

### CO-OK-3 — 2nd-order reduction seeds and recurs BOTH `u` and `u'` correctly (confidence 0.95)

`src/Coefficients.jl:183-200`:

```julia
u  = Taylor1(T[y0, y1], order)   # u[0]=y0, u[1]=y1
up = Taylor1(T[y1],     order)   # up[0]=y1
@inbounds for j in 2:order
    f_t = f(z, u, up)
    u[j] = f_t[j-2] / (j * (j-1))
    up[j-1] = j * u[j]
end
```

This is the highest-suspicion item in the assignment; it checks out. Derivation
against ground truth (Jorba-Zou eq. (5), `...md:276`): write the companion system
`u' = up`, `up' = f`. Jorba-Zou gives `u^{[j]} = up^{[j-1]}/j` and
`up^{[j]} = f^{[j-1]}/j`. Eliminating `up`: `u^{[j]} = (f^{[j-2]}/(j-1))/j =
f^{[j-2]}/(j(j-1))` — exactly line 198. The resync `up[j-1] = j*u[j]` is the
formal-derivative identity `up^{[j-1]} = ((j-1)+1) u^{[(j-1)+1]} = j u^{[j]}`
(equivalently `up^{[j-1]} = f^{[j-2]}/(j-1)`, consistent with the Jorba-Zou form).
- **Seed of `u'(0)` is correct**: `up[0] = y1 = u'(z0)`, and `u[1] = y1` so the
  formal-derivative invariant `up[0] = 1·u[1]` already holds at entry; that is why
  the `j=2` pass needs no top-of-loop `up` resync (docstring lines 164-165 / inline
  189-195 match the code). No off-by-one in the derivative seed.
- **No `up` lag bug**: `f_t[j-2]` depends on inputs up to order `j-2`; pass `j-1`
  set `up[j-2]` via `up[(j-1)-1]`, so `up[0..j-2]` are all available on entry to
  pass `j`. The last coefficient `up[order]` is never set, but `up` is internal
  scaffolding — only `u.coeffs` is returned (line 202) and the highest needed `up`
  index across the whole loop is `order-2`. Harmless.
- Mutation C (`test/coefficients_test.jl:112-119`, `/(j(j-1)) → /(j(j+1))`) was
  confirmed to RED test 3.1.3. The Weierstrass-℘ closed-form oracle (3.1.3) pins
  all 31 coefficients to rtol 1e-12, double-validated by AsymptoticDSolveValue and
  Series[WeierstrassP] agreeing to 3.78e-15.

### CO-OK-4 — `order` guards / fail-loud and type stability are correct (confidence 0.95)

`order < 1` (1st, `Coefficients.jl:115`) and `order < 2` (2nd, `:172`) throw
`ArgumentError` with `suggestion`/`detail` text (Rule 1). Seeds use explicit
`T[...]` so the coefficient `Vector` is `Vector{T}` not `Vector{Float64}`
(`:124`, `:183-184`); confirmed correct against the silent-widening hazard the
docstring calls out. Tests 3.1.4 (BigFloat-256) asserts `eltype == BigFloat`.

### VC-OK-5 — vector loop `y[i][k] = f_t[i][k-1] / k` is correct and aliasing-free (confidence 0.95)

`src/VectorCoefficients.jl:159-168`:

```julia
@inbounds for k in 1:order
    f_t = f(z, y)
    length(f_t) == d || throw(ArgumentError(...))
    for i in 1:d
        y[i][k] = f_t[i][k-1] / k
    end
end
```

Component-wise lift of CO-OK-1; the per-component arithmetic is identical
(`y[i][k] = f_t[i][k-1]/k`), matching `solcoeff!` per component. **Aliasing /
intermittency check — the prime suspect for "intermittent discontinuity" — comes
back clean**:
- The seed array is built by comprehension
  `Taylor1{T}[Taylor1(T[y0[i]], order) for i in 1:d]` (`:153`), so each component
  is a *distinct* `Taylor1` object — no shared backing buffer between components.
- `f_t = f(z, y)` is evaluated **once per pass** and returns a fresh array; all
  `f_t[i][k-1]` reads use a single consistent snapshot of `y` (whose first `k`
  coefficients are correct) *before* any `y[i][k]` write. So coupling order
  `y[1]` then `y[2]` cannot let a half-updated `y[1]` leak into `f_t[2]`. No
  order/RNG dependence, no in-place aliasing.
- The `@inbounds` is safe: `f_t[i][k-1]` has max index `order-1` (valid for a
  series of order `order`), `y[i][k]` max index `order`; both within `0..order`.
  Cannot mask an out-of-range read.
- d=1 reduction is bit-identical to the scalar routine (same single-component
  arithmetic, same constructor, same copy-out). Test VC.1.1
  (`test/vector_coefficients_test.jl:32-48`) asserts `vec_jet[1] == scalar_jet`
  with `==` (not `isapprox`); the FW Bessel-ratio and exp ODEs are the oracles.
  Mutation A (`/k → /(k+1)`) and Mutation B (`for i in 1:1` instead of `1:d`) were
  both confirmed to RED the suite (lines 135-155), proving the recurrence and the
  component-wise loop are genuinely exercised.

### VC-OK-6 — vector fail-loud guards correct (confidence 0.93)

`order < 1`, empty `y0` (`d == 0`), and RHS dimension mismatch
(`length(f_t) != d`) all throw `ArgumentError` (`VectorCoefficients.jl:138-145,
161-164`). The dimension-mismatch check is *inside* the loop after the first
`f(z,y)` call, so a mis-shaped RHS fails loudly on pass 1 rather than silently
truncating. Test VC.1.5 covers all three.

---

## Summary

The coefficient layer (`src/Coefficients.jl`, `src/VectorCoefficients.jl`) is a
faithful, correctly-indexed port of the FW 2011 §2.1.2 method-(b) /
Jorba-Zou eq.(5) bootstrap, identical to the canonical TaylorIntegration.jl
`solcoeff!` step. No sign error, no off-by-one, no factorial/derivative
confusion, no branch-cut issue (none present — pure rational arithmetic), no
in-place aliasing or order dependence. The intermittent-discontinuity bug the
maintainer is hunting is **very unlikely to live here**; this layer is among the
better-pinned in the repo (Weierstrass-℘ and Bessel-ratio closed-form oracles
plus a d=1 bit-identical cross-check, all mutation-proven). Recommend the search
move to the Padé/SVD conversion (data-dependent block structure → classic source
of intermittency) and the path-network / step-control layers.
