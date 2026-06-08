#!/usr/bin/env python3
"""capture.py — mpmath/sympy oracle for corpus_out_of_class_test.jl.

Family G (docs/test_corpus/02_corpus_extension_plan.md:292-306): OUT-OF-CLASS /
FAIL-LOUD.  The package is meromorphic-only (Rule 1: fail loud, never silently
lie).  The probe ODE has an ESSENTIAL singularity — NOT a pole — at z = 0:

    u'' = u·(1 + 2z) / z^4       exact solution   u(z) = e^{1/z}

z = 0 is an essential singularity of e^{1/z} (the Laurent series 1 + 1/z +
1/(2z^2) + ... has infinitely many negative-power terms).  The local jet's
Taylor coefficients grow super-geometrically (|c_k|^{1/k} drifts upward without
bound), which is the signature of non-meromorphy that the meromorphic-only
contract should detect and refuse.  This script is the residual==0 GATE plus
the e^{1/z} pin source for the EMPIRICAL behaviour test.

This is the candidate Rule-1 bug `padetaylor-v1ub`.  Whether the Julia solver
SILENTLY BRIDGES this (returns finite plausible values, no throw/NaN — bug
CONFIRMED) or FAILS LOUD (step-control collapse / DomainError / NaN — bug
REFUTED) is determined EMPIRICALLY by driving the REAL solver; this script only
supplies the independent e^{1/z} ground truth to measure the relerr curve
against.

Run:  python3 external/probes/corpus-oracles/out-of-class/capture.py
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 50

# ---------------------------------------------------------------------------
# (1) RESIDUAL==0 GATE — sympy, symbolic.  u = e^{1/z} must satisfy
#     u'' = u·(1 + 2z)/z^4   exactly.
# ---------------------------------------------------------------------------
z = sp.symbols('z')
u = sp.exp(1 / z)
res = sp.simplify(sp.diff(u, z, 2) - u * (1 + 2 * z) / z**4)
assert res == 0, f"residual != 0: {res}"

# Confirm z=0 is an ESSENTIAL singularity (not a pole): the principal part of
# the Laurent series of e^{1/z} about z=0 has infinitely many terms.  We check
# that the first several negative-power coefficients (1/z^k) are all nonzero —
# a pole would terminate.  e^{1/z} = sum_{n>=0} 1/(n! z^n).
for n in range(1, 8):
    coeff = sp.Rational(1, sp.factorial(n))         # coefficient of z^{-n}
    assert coeff != 0, f"Laurent coeff of z^-{n} is zero (would be a pole)"
print("# residual==0 GATE PASSED; z=0 confirmed ESSENTIAL (not a pole).")
print()


def u_cf(zz):
    """Exact closed form u(z) = e^{1/z}, mpmath."""
    return mp.e ** (1 / zz)


def up_cf(zz):
    """u'(z) = -e^{1/z}/z^2."""
    return -mp.e ** (1 / zz) / zz**2


def emit(label, val):
    if isinstance(val, mp.mpc):
        print(f"const {label} = complex({mp.nstr(val.real, 45)}, "
              f"{mp.nstr(val.imag, 45)})")
    else:
        print(f"const {label} = {mp.nstr(val, 45)}")


# ---------------------------------------------------------------------------
# (2) INITIAL CONDITIONS for the negative-axis approach run (real, z<0).
#     We integrate in the +z (real-increasing) direction toward 0 from the
#     left, starting at z0 = -1:  u(-1) = e^{-1},  u'(-1) = -e^{-1}/1 = -e^{-1}.
#     On z<0, u = e^{1/z} is REAL, finite, and ->0 as z->0^-.
# ---------------------------------------------------------------------------
print("# Negative-axis approach: IC at z0 = -1 (real).")
emit("CFAIL_U0_AT_NEG1",  u_cf(mp.mpf('-1.0')))      # e^{-1}
emit("CFAIL_UP0_AT_NEG1", up_cf(mp.mpf('-1.0')))     # -e^{-1}
print()

# Reference u(z) = e^{1/z} on the FAR side of the start (z in (-1, 0)),
# approaching the essential singularity from the left.  The relerr-vs-distance
# curve is measured against these.  z<0 => e^{1/z} in (0, e^{-1}), well-behaved.
print("# u(z)=e^{1/z} on the negative axis, approaching z=0^- from z0=-1.")
emit("CFAIL_U_AT_NEG0p5",  u_cf(mp.mpf('-0.5')))     # e^{-2}
emit("CFAIL_U_AT_NEG0p2",  u_cf(mp.mpf('-0.2')))     # e^{-5}
emit("CFAIL_U_AT_NEG0p1",  u_cf(mp.mpf('-0.1')))     # e^{-10}
emit("CFAIL_U_AT_NEG0p05", u_cf(mp.mpf('-0.05')))    # e^{-20}
print()

# ---------------------------------------------------------------------------
# (3) ACROSS-ZERO window IC + far target.  Start at z0 = -0.2 (u = e^{-5}),
#     target z = +0.2 (exact u = e^{+5}).  Crossing z=0 means crossing the
#     essential singularity; e^{1/z} blows up to +inf as z->0^+.
# ---------------------------------------------------------------------------
print("# Across-zero window: IC at z0 = -0.2, far target z = +0.2.")
emit("CFAIL_U0_AT_NEG0p2",  u_cf(mp.mpf('-0.2')))    # e^{-5} ~= 0.006738
emit("CFAIL_UP0_AT_NEG0p2", up_cf(mp.mpf('-0.2')))   # -e^{-5}/0.04
emit("CFAIL_U_AT_POS0p2",   u_cf(mp.mpf('0.2')))     # e^{+5} ~= 148.413
print()

# A milder positive-axis target just inside z>0 (e^{1/z} for small z>0 is huge):
emit("CFAIL_U_AT_POS0p5",   u_cf(mp.mpf('0.5')))     # e^{2}
emit("CFAIL_U_AT_POS1p0",   u_cf(mp.mpf('1.0')))     # e^{1}
print()
print("# end of pins.")
