#!/usr/bin/env python3
# capture.py -- INDEPENDENT mpmath oracle for the Taylor-jet corpus
# (bead padetaylor-jwhf, epic padetaylor-p1v0, bucket B2).
#
# WHAT THIS SCRIPT IS
# -------------------
# A *second, independent* implementation of the scalar 2nd-order Taylor
# recurrence that src/Coefficients.jl::taylor_coefficients_2nd computes
# via TaylorSeries.jl.  The whole point of the B2 corpus gap is that NO
# BigFloat coefficient-level oracle existed for taylor_coefficients_2nd
# (the routine under the entire scalar Painleve/Weierstrass pipeline).
# This script *re-derives the recurrence from the ODE in Python* and
# pins the coefficient vector; the Julia test then asserts the Julia
# routine reproduces it.  Asserting the routine against itself would be
# worthless (CLAUDE.md Rule 5) -- so the recurrence below is written out
# by hand from the convolution, NOT by calling any Julia code.
#
# THE RECURRENCE (re-derived, with citation)
# ------------------------------------------
# For u'' = f(z, u, u') we expand u(z0+h) = sum_n c[n] h^n.  Then
#   u''(z0+h) = sum_{m>=0} (m+2)(m+1) c[m+2] h^m.
#
#   * Case u'' = 6 u^2  (wp-taylor-recurrence, catalogue md:1565-1580):
#       u^2 = sum_m ( sum_{k=0}^{m} c[k] c[m-k] ) h^m  (Cauchy product).
#       Matching h^{n-2} on both sides (m = n-2):
#         n(n-1) c[n] = 6 * sum_{k=0}^{n-2} c[k] c[n-2-k]
#       => c[n] = 6/(n(n-1)) * sum_{k=0}^{n-2} c[k] c[n-2-k],  n >= 2.
#       The convolution index range is k = 0 .. n-2 (n-1 terms);
#       derived above from u''=6u^2 and CITED to catalogue closed_form.
#
#   * Case u'' = 6 u^2 + z  (pi1-taylor-recurrence, catalogue md:1582-1597):
#       the inhomogeneity z = z0 + h contributes its Taylor coeffs:
#       [z]_m = z0 for m=0, 1 for m=1, 0 otherwise.  With z0=0 that adds
#       +0 at n=2 (so c[2] is UNCHANGED = 3 u0^2) and +1 at n=3:
#         c[3] = (6 * (2 c0 c1) + 1) / (3*2) = (12 c0 c1 + 1)/6 = 2 c0 c1 + 1/6.
#       For n>=4 the z-term has vanished, recurrence is identical to WP.
#
# INITIAL CONDITIONS -- shared with the Julia test (load-bearing)
# --------------------------------------------------------------
# The oracle is only meaningful if the SAME ICs are fed to both the
# Python recurrence and the Julia routine.  We therefore print the ICs
# as 80-dps strings; the Julia test parses the IDENTICAL strings into
# BigFloat (and round-trips to Float64).  The WP ICs are obtained from
# the closed form u = WP(z-1; 0, 2) so they are a genuine Weierstrass
# solution (catalogue wp-taylor-coeff-recurrence-b2 ic_bc); the PI ICs
# are the FW2011 eq.(4.1) tritronquee values (catalogue md:1586).
#
# Run:  python3 capture.py        (mpmath 1.3, dps = 80)

import mpmath as mp

mp.mp.dps = 80

# 80-dps strings consumed verbatim by both the Python recurrence and the
# Julia test (parse(BigFloat, ...) and parse(Float64, ...)).  Anything
# the Julia side does NOT parse from these exact strings would break the
# shared-IC contract and make the oracle a lie.
WP_U0_STR = None   # filled below from the closed form
WP_UP0_STR = None

ORDER = 30          # B2 catalogue order (30+); we pin 0..30.


# --- Weierstrass closed-form ICs (genuine WP(-1; 0, 2) solution) ---------
# Same equianharmonic construction as the scalar-pole-bridge capture.py:
# WP(w; 0, 2) = e3 + (e1 - e3)/sn^2(sqrt(e1-e3) w, k^2), w = z - 1.
pi = mp.pi
e1 = (mp.mpf(1) / 2) ** (mp.mpf(1) / 3)
e2 = e1 * mp.exp(mp.mpc(0, 2 * pi / 3))
e3 = e1 * mp.exp(mp.mpc(0, 4 * pi / 3))
k_sq = (e2 - e3) / (e1 - e3)


def wp(z, c1=-1):
    w = z + c1
    return (e3 + (e1 - e3) / mp.ellipfun('sn', mp.sqrt(e1 - e3) * w, k_sq) ** 2).real


def wp_deriv(z, c1=-1):
    return mp.diff(lambda zz: wp(zz, c1), z).real


WP_U0 = wp(mp.mpf(0))
WP_UP0 = wp_deriv(mp.mpf(0))
WP_U0_STR = mp.nstr(WP_U0, 78)
WP_UP0_STR = mp.nstr(WP_UP0, 78)

# --- PI tritronquee ICs (FW2011 eq.(4.1), catalogue md:1586) -------------
PI_U0_STR = "-0.18755430834049490"
PI_UP0_STR = "0.30490556026122890"


def recurrence_2nd(u0, u1, order, inhomog):
    """INDEPENDENT 2nd-order Taylor recurrence for u'' = 6 u^2 + g(z),
    z0 = 0.  `inhomog` maps Taylor index m of g(z0+h) to its coefficient;
    here g=0 (WP) or g=z so [z]_1 = 1 (PI).  Re-derived above, not ported
    from Julia."""
    c = [mp.mpf(u0), mp.mpf(u1)]
    for n in range(2, order + 1):
        conv = sum(c[k] * c[n - 2 - k] for k in range(n - 1))  # k=0..n-2
        rhs = 6 * conv + inhomog(n - 2)                        # g coeff at h^{n-2}
        c.append(rhs / (n * (n - 1)))
    return c


def wp_inhomog(m):
    return mp.mpf(0)


def pi_inhomog(m):
    # g(z) = z = z0 + h, z0 = 0 -> coeff 1 at h^1, else 0.
    return mp.mpf(1) if m == 1 else mp.mpf(0)


def emit(name, coeffs):
    print(f"# {name}: c[0..{len(coeffs) - 1}] (80-dps strings)")
    print(f"const {name} = [")
    for c in coeffs:
        print(f'    "{mp.nstr(c, 78)}",')
    print("]")
    print()


print("# ===== Taylor-jet oracle (mpmath dps=80, INDEPENDENT recurrence) =====")
print(f'const CTJ_WP_U0_STR  = "{WP_U0_STR}"')
print(f'const CTJ_WP_UP0_STR = "{WP_UP0_STR}"')
print(f'const CTJ_PI_U0_STR  = "{PI_U0_STR}"')
print(f'const CTJ_PI_UP0_STR = "{PI_UP0_STR}"')
print()

wp_c = recurrence_2nd(WP_U0, WP_UP0, ORDER, wp_inhomog)
pi_c = recurrence_2nd(mp.mpf(PI_U0_STR), mp.mpf(PI_UP0_STR), ORDER, pi_inhomog)

# Spot-check invariants the catalogue demands (md:1580, md:1597).
assert abs(wp_c[2] - 3 * WP_U0 ** 2) < mp.mpf(10) ** -70, "WP c[2] != 3 u0^2"
u0 = mp.mpf(PI_U0_STR)
u1 = mp.mpf(PI_UP0_STR)
assert abs(pi_c[2] - 3 * u0 ** 2) < mp.mpf(10) ** -70, "PI c[2] != 3 u0^2 (z0=0)"
assert abs(pi_c[3] - (2 * u0 * u1 + mp.mpf(1) / 6)) < mp.mpf(10) ** -70, \
    "PI c[3] != 2 u0 u1 + 1/6"

emit("CTJ_WP_C", wp_c)
emit("CTJ_PI_C", pi_c)

# Cross-check the WP recurrence against the FW2011 Table 5.1 keystone by
# summing the series at a point WELL inside the convergence radius (z=0.3,
# radius ~1 to the pole at z=1) and comparing to the closed form.
zc = mp.mpf('0.3')
series_val = sum(wp_c[n] * zc ** n for n in range(ORDER + 1))
closed_val = wp(zc)
print(f"# WP series(0.3) = {mp.nstr(series_val, 30)}")
print(f"# WP closed(0.3) = {mp.nstr(closed_val, 30)}")
print(f"# |series - closed| at z=0.3 = {mp.nstr(abs(series_val - closed_val), 6)}")
print()

# =====================================================================
# B6 -- orphaned step-control branches (jorba_zou fallback + relative).
# Closed-form Jorba-Zou formulae re-derived from TI.jl stepsize.jl, NOT
# read off the existing oracle.  catalogue md:1650-1682.
# =====================================================================
print("# ===== B6 Jorba-Zou step-control oracle (dps=50) =====")
mp.mp.dps = 50

# jorba-zou-fallback-trigger-b6 (md:1650): even-only cos jet through c[28],
# c[29]=c[30]=0 -> primary formula yields Inf (both candidates zero) ->
# _second_stepsize scan: for each Taylor index j it solves the "term
# magnitude = 1" criterion |c[j]| h^j = 1  =>  h = (1/|c[j]|)^(1/j), the
# SAME j in coefficient AND exponent (src/StepControl.jl:194-197, the
# verbatim TI.jl _second_stepsize port; confirmed bit-for-bit faithful in
# docs/bug-sweep-2026-06-01/find-A3-stepcontrol.md:134-145).  The largest
# nonzero coefficient is c[28], so the max is at j=28:
#     h = (1/|c[28]|)^(1/28) = (28!)^(1/28) = 11.2980899840442.
#
# *** CATALOGUE ERROR (md:1656) ***  The catalogue pins (28!)^(1/27) =
# 12.3596..., mismatching coefficient index 28 with exponent 1/27.  That
# violates the |c[j]| h^j = 1 criterion (coefficient and exponent indices
# must agree) and disagrees with the verbatim TI.jl port.  The corpus
# exists to root out bugs; here the defect is in the ORACLE RECIPE, not the
# code.  We pin the code-faithful, TI.jl-faithful value (28!)^(1/28); see
# the test header for the full triage.
h_fallback = (mp.factorial(28)) ** (mp.mpf(1) / 28)
h_fallback_cat = (mp.factorial(28)) ** (mp.mpf(1) / 27)  # the WRONG catalogue value
print(f"const CTJ_JZ_FALLBACK_H = {mp.nstr(h_fallback, 45)}")
print(f"# catalogue-claimed (WRONG) (28!)^(1/27) = {mp.nstr(h_fallback_cat, 30)}")

# jorba-zou-relative-mode-b6 (md:1667): full cos jet (c[0]=1, c[30]!=0),
# eps_abs=1e-15 < eps_rel*|c0|=1e-12 -> eps_eff = 1e-12.  Only c[30] nonzero
# among {c[29],c[30]} -> h = (eps_eff/|c[30]|)^(1/30) = (1e-12 * 30!)^(1/30).
eps_eff = mp.mpf('1e-12')
h_rel = (eps_eff / (mp.mpf(1) / mp.factorial(30))) ** (mp.mpf(1) / 30)
print(f"const CTJ_JZ_RELATIVE_H = {mp.nstr(h_rel, 45)}")

# Counter-value for the relative-mode mutation check: if eps_abs=1e-15 were
# (wrongly) used, h would be (1e-15 * 30!)^(1/30) -- a DIFFERENT number.
eps_abs = mp.mpf('1e-15')
h_rel_wrong = (eps_abs / (mp.mpf(1) / mp.factorial(30))) ** (mp.mpf(1) / 30)
print(f"# (mutation counter-value, eps_abs misused) = {mp.nstr(h_rel_wrong, 20)}")

# sin-jet odd-parity value (== existing h_4_1_1_TI): (1e-12 * 29!)^(1/29).
h_sin = (mp.mpf('1e-12') * mp.factorial(29)) ** (mp.mpf(1) / 29)
print(f"# sin-jet h (== h_4_1_1_TI, NOT re-pinned) = {mp.nstr(h_sin, 20)}")
