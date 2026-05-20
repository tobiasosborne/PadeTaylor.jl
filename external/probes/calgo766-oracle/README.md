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
| `build.sh`     | Compiles `src-0766/` + `wrapper.f90` with `gfortran` into `libcalgo766_dp.so` (double precision) and `libcalgo766.so` (single precision). |
| `wrapper.f90`  | The `working_area_VECTOR_PADE` module Calgo needs, plus a `bind(C)` shim `vpade_type2` with a flat, `ccall`-friendly ABI. **Committed.** |
| `capture.jl`   | `ccall`s `vpade_type2` (from `libcalgo766_dp.so`) on the SharedPade test inputs, pins outputs, sanity-checks vs `shared_denominator_pade`. **Committed.** |
| `oracles.txt`  | Pinned Calgo-766 type-II outputs (shared `Q`, numerators `P_i`). **Committed.** |
| `src-0766/`, `libcalgo766*.so`, `*.o`, `*.mod` | Gitignored — reproducible via `fetch.sh` + `build.sh`. |

## Reproduce

```sh
./fetch.sh        # netlib -> src-0766/
./build.sh        # src-0766/ + wrapper.f90 -> libcalgo766_dp.so + libcalgo766.so
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
`Float64` data into Calgo's `real` arrays and back. In the
double-precision build those marshalling conversions are double→double
no-ops, and the `ccall` ABI is unchanged because the `bind(C)` interface
is written with `iso_c_binding`'s `c_double` / `c_int` named constants —
which always denote the C `double` / `int` kinds and are *not* affected
by `gfortran`'s `-fdefault-real-8` (see "Accuracy floor" below).

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

The netlib TOMS mirror (`https://www.netlib.org/toms/766`) ships **only**
a single-precision `Sp/` source tree — there is no `Dp/` double-precision
variant. But the Calgo 766 FORTRAN is **fully precision-agnostic**: every
floating declaration in `vector_pade.f90`, `linpack.f`, `blas1.f` and
`wrapper.f90` is a bare `real` (no `kind` parameter, no
`selected_real_kind`, no `real*8`), and the only numeric literal of any
note — `-1.0e0` in `linpack.f` — is a default-`real`-kind literal.
Compiling with

```
gfortran -fdefault-real-8 -fdefault-double-8 ...
```

therefore promotes the entire algorithm to **IEEE double precision** with
no source edits. `build.sh` produces both builds:

| Library             | Precision | Oracle accuracy floor |
|---------------------|-----------|-----------------------|
| `libcalgo766_dp.so` | double (`-fdefault-real-8`) | ~`1e-15` (machine epsilon) |
| `libcalgo766.so`    | single (default `real`)     | ~`1e-6`–`1e-7` |

`capture.jl` loads `libcalgo766_dp.so`. The single-precision library is
kept only for comparison. The `ccall` ABI is identical for both: the
`bind(C)` interface in `wrapper.f90` uses `iso_c_binding`'s `c_double` /
`c_int`, which are unaffected by `-fdefault-real-8`, so `capture.jl`'s
`ccall` signature stays `Float64` / `Cint`.

V1e may therefore assert agreement against this oracle at **`1e-12`** —
on par with the block-Toeplitz determinant oracle V1c and the AAA oracle
V1d, rather than the loose `1e-5` a single-precision build would force.

## Sanity-check result

`capture.jl` run 2026-05-20, gfortran 13.3.0, Julia 1.10, double-precision
build (`libcalgo766_dp.so`):

```
SP_1_1_d1_scalar     flag = 0   max |Calgo - SharedPade| & |Calgo - truth| = 2.48e-16   AGREE
SP_1_2_d2_shared_Q   flag = 2   max |Calgo - SharedPade| & |Calgo - truth| = 8.94e-16   AGREE
```

Calgo 766 and `SharedPade.shared_denominator_pade` agree to **~1e-15** —
machine epsilon — on both the SP.1.1 scalar rational jet and the SP.1.2
`d=2` shared-`Q` jets. The sanity check asserts `< 1e-12`. The oracle is
sound.

`flag = 2` on SP.1.2 is **not** a failure. The Calgo header
(`vector_pade.f90` lines 158–172) defines `flag = 2` as "the Sylvester
matrix at the point `n` is numerically singular", but explicitly adds:
"the first row of `S_star` still yields a simultaneous Padé approximant
of type `n`; the remaining rows and columns are meaningless." The wrapper
extracts **exactly** row 0 of `S_star` — the type-II result — so `flag = 2`
is benign for this oracle. The singularity is genuine: embedding an exact
degree-2 rational at the high degree vector `n = [2,3,2]` makes the
Sylvester matrix rank-deficient by construction. At single precision the
old `Sp/` build masked this as mere ill-conditioning (`flag = 1`); double
precision detects the true singularity. The row-0 approximant is
unaffected — the ~`1e-15` value-agreement above confirms it. Only
`flag = 3` (bad input) is a genuine error; `capture.jl` errors on it.
