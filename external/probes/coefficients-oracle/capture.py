"""capture.py -- Phase 3 Coefficients oracle (Python cross-check).

Independent oracle in Python (sympy + mpmath) for cases 3.1.1, 3.1.2, 3.1.4.
Case 3.1.3 (Weierstrass P) skipped: mpmath has no native WeierstrassP and
sympy's series expansion of WeierstrassP is unreliable; wolframscript is
the primary for that case.

Usage:  python3 capture.py
Outputs: oracles_python.txt (Julia-readable literals).

Run AFTER capture.wl, then run verify.jl to assert agreement to 1e-14
(Float64 cases) and string-match (BigFloat-256 case).
"""

import sympy as sp
import mpmath
from sympy import Symbol, Rational, exp, besselj, factorial, Integer, series

OUT = open('oracles_python.txt', 'w')
OUT.write("# Phase 3 Coefficients oracle (Python sympy/mpmath cross-check).\n")
OUT.write("# Captured by: python3 external/probes/coefficients-oracle/capture.py\n")
OUT.write("# Compared against oracles.txt (wolframscript) by verify.jl.\n\n")


def write_float64_vec(name, values):
    OUT.write(f"{name} = [")
    OUT.write(", ".join(repr(float(v)) for v in values))
    OUT.write("]\n\n")


def write_bigfloat_str_vec(name, values_str):
    OUT.write(f"{name} = [\n  ")
    OUT.write(",\n  ".join(f'"{s}"' for s in values_str))
    OUT.write("\n]\n\n")


# --- 3.1.1: exp(z) Taylor coeffs, order 10 ----------------------------------
OUT.write("# Case 3.1.1 -- exp(z) Taylor coeffs, order 10. Source: c_k = 1/k! exact rationals.\n")
exp_coefs = [Rational(1, sp.factorial(k)) for k in range(11)]
write_float64_vec("c_3_1_1_exp_o10_py", exp_coefs)

# --- 3.1.2: t * J_{3/4}(t^2/2) / J_{-1/4}(t^2/2}, order 14 ------------------
# sympy's series() can struggle on Bessel ratios with fractional order at 0.
# Approach: use the explicit power-series expansion of J_nu(z),
#   J_nu(z) = sum_{k>=0} (-1)^k / (k! Gamma(nu + k + 1)) * (z/2)^{nu+2k}
# For our function f(t) = t * J_{3/4}(t^2/2) / J_{-1/4}(t^2/2):
#   J_{3/4}(t^2/2)  = sum_k (-1)^k / (k! Gamma(k+7/4)) * (t/sqrt(2))^{3/2+4k} / 2^{...}
#   J_{-1/4}(t^2/2) = sum_k (-1)^k / (k! Gamma(k+3/4)) * (t/sqrt(2))^{-1/2+4k} / 2^{...}
# The leading t^{-1/2} of J_{-1/4} cancels the t^{1/2} prefactor in
# t * J_{3/4} / J_{-1/4}, leaving a Taylor series in t^4 starting at t^3.
#
# Concretely, factor out t^{3/2} from J_{3/4}(t^2/2) and t^{-1/2} from
# J_{-1/4}(t^2/2): both become Taylor series in t^4. The ratio is a Taylor
# series in t^4, multiplied by t^3 (= t * t^{3/2} / t^{-1/2}), so f(t) is a
# Taylor series in t starting at t^3 with only every-fourth term nonzero.
#
# We compute symbolically and trust sympy to handle this.
t = Symbol('t', positive=True, real=True)


def J_series(nu, arg, n_terms):
    """Power-series expansion of BesselJ[nu, arg] up to n_terms. arg is symbolic."""
    nu = Rational(nu).limit_denominator(1000)
    s = Integer(0)
    for k in range(n_terms):
        coef = Rational((-1) ** k) / (sp.factorial(k) * sp.gamma(nu + k + 1))
        s += coef * (arg / 2) ** (nu + 2 * k)
    return s


# Expand both Bessels with enough terms so that the ratio is correct to t^14.
# J_{3/4}(t^2/2) needs terms up through t^{3/2 + 4k}; for ratio accuracy to
# t^14 we need k up to ~4. Pad generously.
n_keep = 8
num = J_series(Rational(3, 4), t**2 / 2, n_keep)
den = J_series(Rational(-1, 4), t**2 / 2, n_keep)
ratio = sp.simplify(t * num / den)

# Series expand around t=0 to order 15 (so we capture coefficients of t^0..t^14).
expanded = sp.series(ratio, t, 0, 15).removeO()
expanded = sp.expand(expanded)

bessel_coefs = []
for k in range(15):
    coef = expanded.coeff(t, k)
    bessel_coefs.append(coef)

OUT.write("# Case 3.1.2 -- Bessel ratio Taylor coeffs, order 14.\n")
OUT.write("# Source: explicit J_nu power series, sympy series of the ratio.\n")
write_float64_vec("c_3_1_2_bessel_o14_py", bessel_coefs)

# --- 3.1.4: 1/k! at 80 decimal digits, k=0..60 -------------------------------
mpmath.mp.dps = 80
OUT.write("# Case 3.1.4 -- 1/k! at 80 decimal digits (mpmath), k=0..60.\n")
str_coefs = [mpmath.nstr(mpmath.mpf(1) / mpmath.factorial(k), 78,
                         strip_zeros=False) for k in range(61)]
write_bigfloat_str_vec("c_3_1_4_exp_o60_py_str", str_coefs)

OUT.close()
print("Wrote oracles_python.txt")
