#!/usr/bin/env python3
# =============================================================================
# capture.py — high-precision (mpmath) oracle generator for
# test/corpus_robust_pade_test.jl (bead padetaylor-7pav, epic
# padetaylor-p1v0).
#
# Two oracles are produced, both independent of the Julia implementation:
#
#   (A) adr0002-rank-boundary-svd
#       f(z) = 1/(1-z) + 1e-13/(1-z/2)  =>  c_k = 1 + 1e-13*(1/2)^k.
#       The (m=10,n=10) Padé builds the n×(n+1) = 10×11 lower-triangular
#       Toeplitz block  C[i,j] = c[m+1+i-j]  (GGT 2013 eq.2.6, the bottom
#       n rows of the (m+n+1)×(n+1) Toeplitz Z, matching RobustPade.jl's
#       `Z[(m+2):(m+n+1), :]`).  We compute its singular values at 50 dps.
#       The genuine rank is 2 (two poles, z=1 and z=2); sigma_2 sits a
#       hair (~3.9%) above the Float64 GGT threshold 1e-14*||c||_2.  This
#       is the ADR-0002 stress case: relative-accuracy Jacobi SVD must
#       resolve sigma_2; absolute-accuracy Demmel-Kahan (LAPACK) might
#       not.  Catalogue: docs/test_corpus/01_corpus_catalogue.md:1527-1542.
#
#   (B) geom-noisy-1on1z (deterministic, noise-free reduction variant)
#       The catalogue's frozen Octave oracle (test/_oracles.jl 2.1.6)
#       uses MATLAB randn(state,42), unreproducible here.  But the
#       *structural* claim — tol above noise collapses (10,10) to (0,1)
#       with denominator root ≈ 1 — is regenerable with mpmath-built
#       deterministic noise.  We emit the exact noise-free reduction
#       (c_k = 1 exactly => exact type (0,1), b=[1,-1]) as the
#       byte-stable oracle, plus a deterministic tiny-noise vector so the
#       Julia test can assert the tol>noise collapse without depending on
#       a particular RNG.  Catalogue lines 1493-1508.
#
# Reproduce:  python3 external/probes/corpus-oracles/robust-pade/capture.py
# =============================================================================
import mpmath

mpmath.mp.dps = 50

# ---- (A) ADR-0002 rank-boundary Toeplitz singular values --------------------
delta = mpmath.mpf('1e-13')
c = [1 + delta * mpmath.mpf('0.5') ** k for k in range(21)]   # c_0..c_20
m, n = 10, 10
# Bottom-n-rows Toeplitz block C (n x (n+1)), GGT 2013 eq.2.6 /
# RobustPade.jl Z[(m+2):(m+n+1), :].  C[i,j] = c[m+1+i-j], 0-indexed i,j.
C = mpmath.matrix(n, n + 1)
for i in range(n):
    for j in range(n + 1):
        idx = m + 1 + i - j
        if 0 <= idx < len(c):
            C[i, j] = c[idx]
U, S, V = mpmath.svd(C)                 # S is descending min(n,n+1)=10 vector
print("# (A) adr0002-rank-boundary-svd")
print("# 10x11 lower-tri Toeplitz of c_k = 1 + 1e-13*(0.5)^k, mpmath dps=50")
cnorm = mpmath.sqrt(sum(x * x for x in c[m:m + n + 1]))  # ||c[m..m+n]||_2 used in GGT
# GGT thresholds on the FULL c[1..m+n+1] 2-norm (RobustPade uses norm(cv), cv length m+n+1):
cnorm_full = mpmath.sqrt(sum(x * x for x in c[:m + n + 1]))
for k in range(min(4, len(S))):
    print(f"sigma_{k+1} = {mpmath.nstr(S[k], 25)}")
print(f"||c[1:m+n+1]||_2  = {mpmath.nstr(cnorm_full, 25)}")
thr_f64 = mpmath.mpf('1e-14') * cnorm_full
print(f"Float64 GGT threshold 1e-14*||c||_2 = {mpmath.nstr(thr_f64, 12)}")
print(f"sigma_2 / threshold = {mpmath.nstr(S[1] / thr_f64, 8)}")
# Demmel-Kahan absolute error bound on sigma_2 at Float64:
dk = 10 * mpmath.mpf(2) ** (-52) * S[0]
print(f"DK abs-error bound (10*2^-52*sigma_1) = {mpmath.nstr(dk, 8)} "
      f"= {mpmath.nstr(100 * dk / S[1], 6)}% of sigma_2")

# ---- (B) geom-noisy / noise-free exact reduction ---------------------------
print()
print("# (B) geom-noisy-1on1z structural reduction")
# Exact noise-free: c_k = 1 => 1/(1-z) is exactly type (0,1): a=[1], b=[1,-1].
print("# noise-free c_k=1 exact Pade type (0,1): a=[1], b=[1,-1], root z=1")
# Deterministic tiny noise (mpmath, reproducible): c_k = 1 + eps*sin(k+1),
# eps = 1e-9 (well below tol=1e-5).  Emit as Float64-printable literals so
# the Julia test pins an identical vector (no RNG dependence).
eps = mpmath.mpf('1e-9')
cnoisy = [1 + eps * mpmath.sin(mpmath.mpf(k + 1)) for k in range(21)]
print("noisy_c = [" + ", ".join(mpmath.nstr(x, 18) for x in cnoisy) + "]")
