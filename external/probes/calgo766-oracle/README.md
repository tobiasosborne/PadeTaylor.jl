# Calgo 766 shared-Q oracle (`calgo766-oracle`)

An **independent** reference implementation of the type-II Hermite–Padé
(simultaneous-Padé / shared-denominator) approximant, for mutation-proving
`PadeTaylor.SharedPade.shared_denominator_pade`. This is the v0.2 analogue
of `padeapprox-oracle/` (which validates `RobustPade` against Chebfun's
`padeapprox.m` in v0.1).

The oracle is ACM TOMS **Algorithm 766**, `VECTOR_PADE`, from Cabay, Jones
& Labahn, *"Algorithm 766: Experiments with a Weakly Stable Algorithm for
Computing Padé–Hermite and Simultaneous Padé Approximants"*, ACM TOMS 23(1),
91–110 (March 1997), **DOI 10.1145/244768.244790**. Its role here is the
role `docs/v0p2_pillarA_hermite_pade_findings.md:330–340` (§6) assigns it:
a published, weakly stable FORTRAN routine for the exact problem
`SharedPade` solves, written by a different team with a different algorithm
(Beckermann–Labahn striped-Sylvester recurrence, not a block-Toeplitz SVD).

## Files

| File           | Role |
|----------------|------|
| `fetch.sh`     | Downloads + extracts the Calgo 766 FORTRAN from the netlib TOMS mirror into `src-0766/` (gitignored). |
| `build.sh`     | Compiles `src-0766/` + `wrapper.f90` with `gfortran` into `libcalgo766.so`. |
| `wrapper.f90`  | The `working_area_VECTOR_PADE` module Calgo needs, plus a `bind(C)` shim `vpade_type2` with a flat, `ccall`-friendly ABI. **Committed.** |
| `capture.jl`   | `ccall`s `vpade_type2` on the SharedPade test inputs, pins outputs, sanity-checks vs `shared_denominator_pade`. **Committed.** |
| `oracles.txt`  | Pinned Calgo-766 type-II outputs (shared `Q`, numerators `P_i`). **Committed.** |
| `src-0766/`, `libcalgo766.so`, `*.o`, `*.mod` | Gitignored — reproducible via `fetch.sh` + `build.sh`. |

## Reproduce

```sh
./fetch.sh        # netlib -> src-0766/
./build.sh        # src-0766/ + wrapper.f90 -> libcalgo766.so
julia --project=../../.. capture.jl    # ccall, pin oracles.txt, sanity-check
```

## Why a hand-written `bind(C)` wrapper

Calgo 766's `VECTOR_PADE` has a pure-Fortran-90 interface: every array
argument is *assumed-shape* (passed as a compiler-private array descriptor,
not a bare pointer) and it `USE`s a host module of allocatable scratch
arrays the caller must supply. Neither is `ccall`-able. `wrapper.f90`
therefore (1) defines that `working_area_VECTOR_PADE` module and (2)
exposes one `bind(C)` entry point `vpade_type2` whose entire ABI is plain
`real*8` / `int32` arrays passed by reference in column-major order. The
FORTRAN core (`vector_pade.f90`, `linpack.f`, `blas1.f`) is used verbatim
from netlib — only the *calling shell* is ours, so the oracle remains
genuinely independent of `SharedPade`.

`ccall` conventions handled in `capture.jl` / `wrapper.f90`:
`bind(C, name="vpade_type2")` fixes the symbol name (no gfortran `_`
mangling); all scalars are `Ref{Cint}` / `Ref{Float64}`, all arrays
`Ptr{...}` passed column-major; `wrapper.f90` marshals the caller's
`Float64` data into Calgo's default `real` and back.

## Convention mapping: Calgo 766 type-II ↔ `SharedPade`

This is the load-bearing section for bead V1e. Read it before wiring the
oracle into `test/`.

### Calgo 766 type-II output

`VECTOR_PADE(k, n, A, tau, S, gamma, S_star, gamma_star, kappa, ...)`
takes `k+1` power series `A(0..k)` and a degree vector `n(0..k)`, and
returns the *scaled* simultaneous-Padé system `S_star`. **Row 0** of
`S_star` is the type-II approximant; from `build_delta_T_star` (driver2,
`alpha = 0`):

```
S_star(:,0,beta) * A(0)  -  S_star(:,0,0) * A(beta)  =  O(z**(||n||+1))
```

so with `A(0)` set to the constant series `1` and `A(beta) := f_beta`:

```
Q       :=  S_star(:,0,0)        degree ||n|| - n(0)
P_beta  :=  S_star(:,0,beta)     degree ||n|| - n(beta)
f_beta  ~=  P_beta / Q           (to order z**(||n||+1))
```

The wrapper writes `Q` to `Qout` and `P_1..P_k` to `Pout`.

### Mapping table

| Aspect | Calgo 766 type-II | `SharedPade` | Bridge |
|---|---|---|---|
| Components | `d = k` components → `k+1` series `A(0..k)` | `d` jets | set `A(0)=1`, `A(beta)=jets[beta]`, `k=d` |
| Denominator | `Q = S_star(:,0,0)` | `denominator` | both **low-to-high** |
| Numerators | `P_beta = S_star(:,0,beta)` | `numerators[beta]` | both **low-to-high** |
| Coefficient order | `z**0` first | `z**0` first (`b[1]=Q(0)`) | identical |
| Normalisation | `S_star` is **scaled** — divide row `i` by `gamma_star(i)`; Calgo does **not** force `Q(0)=1` | `Q(0)=b[1]=1` exactly | `wrapper.f90` divides row 0 by `gamma_star(0)`; `capture.jl` then divides `Q,P` by `Q[1]` |
| Degree control | per-series `n(j)`; `Q` deg `= ||n||-n(0)`, `P_beta` deg `= ||n||-n(beta)` | single shared `m` | choose `n` so type-II degrees ≥ true degrees (see below) |
| Matching order | `O(z**(||n||+1))` | `O(z**(m+n+1)) = O(z**(2m+1))` at `n=m` | differ — compare *function values*, not raw coeffs |

### The degree subtlety V1e must know

`SharedPade` uses a single shared degree `m` for `Q` and every `P_i`.
Calgo's type-II degrees are `||n||-n(j)`, set by the per-series `n`. They
**cannot in general be made equal to `m`** for arbitrary `d`. The oracle
sidesteps this: a rational function that *is* exactly `P/Q` of low type
is recovered exactly by **any** type-II Padé of sufficiently high order —
the extra degrees of freedom are absorbed by a common spurious factor
that divides `Q` and every `P_i` alike, so the *reduced rational function*
`P_i/Q` is still the true one.

Concretely, `oracles.txt` case `SP_1_2_d2_shared_Q` has `n = [2,3,2]`,
giving a degree-5 `Q` — the true degree-2 `Q = [1,-0.4,0.2]` times a
degree-3 spurious factor that also divides `P_1, P_2`. The pinned `Q`/`P`
coefficients therefore look nothing like the true low-degree coefficients,
**by design**.

**Consequence for V1e:** compare the oracle and `SharedPade` by
**rational-function values at sample points** and by **pole locations**
(roots of the *reduced* denominator), never by raw coefficient vectors.
This is the same value-invariant `test/shared_pade_test.jl` already uses
(`_ratval` at several `z`), and the same lesson the SP.1.1 docstring states
("compare rational-function VALUES … not just raw coefficients"). `capture.jl`'s
sanity check follows this rule.

## Accuracy floor

The Calgo `Sp/` build uses default `real` = **single precision**. The
oracle is therefore good to ~`1e-6`–`1e-7`, not to `1e-12`. The sanity
check in `capture.jl` uses a `1e-5` tolerance; V1e assertions against this
oracle must use a comparably loose tolerance. (A `Dp/` double-precision
build is not in the netlib archive; if V1e needs tighter agreement, the
block-Toeplitz determinant oracle V1c or the AAA oracle V1d are the
higher-precision routes — see findings §6.)

## Sanity-check result

`capture.jl` run 2026-05-20, gfortran 13.3.0, Julia 1.10:

```
SP_1_1_d1_scalar     flag = 0   max |Calgo - SharedPade| & |Calgo - truth| = 2.73e-08   AGREE
SP_1_2_d2_shared_Q   flag = 1   max |Calgo - SharedPade| & |Calgo - truth| = 2.69e-08   AGREE
```

Calgo 766 and `SharedPade.shared_denominator_pade` agree to **~3e-8** on
both the SP.1.1 scalar rational jet and the SP.1.2 `d=2` shared-`Q` jets —
at the single-precision floor of the Calgo build. The oracle is sound.

`flag = 1` on SP.1.2 is **not** a failure: it reports that an *intermediate*
striped-Sylvester matrix along the diagonal is ill-conditioned (expected —
embedding an exact low-degree rational in a high-degree type makes those
matrices near-singular). The Calgo header (`vector_pade.f90` lines ~152–171)
guarantees the **row-0 type-II approximant is still valid** in that case;
the value-agreement above confirms it. `flag = 2` (singular) or `flag = 3`
(bad input) *would* be failures — `capture.jl` errors on `flag = 3`.
