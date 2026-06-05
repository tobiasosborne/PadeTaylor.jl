#!/usr/bin/env python3
# capture.py -- mpmath regeneration oracle for the scalar pole-bridge corpus.
#
# Regenerates the transcendental / Bessel / elliptic constants pinned in
#   test/_oracle_corpus_scalar_pole_bridge.jl
# at 50 decimal places.  Each printed line is a labelled Julia-readable
# literal: copy the value into the oracle file with this script + dps cited.
#
# Ground truth recipes are verbatim from
#   docs/test_corpus/01_corpus_catalogue.md
# "Full candidate records" for:
#   wp-equianharmonic-g0g2, riccati-tan-1storder / tan-2nd-order,
#   simple-pole-riccati, coth-simple-pole, tanh-complex-pole,
#   fw-riccati-bessel, jacobi-sn-2ndorder.
#
# Run:  python3 capture.py     (mpmath 1.3 installed)
#
# Tolerances in the test are justified by the METHOD's accuracy, not fitted;
# these constants are the *exact* oracle the method is measured against.

import mpmath as mp

mp.mp.dps = 50


def show(label, val):
    # Print enough digits that Float64 and BigFloat-256 parses are exact.
    print(f"{label} = {mp.nstr(val, 45)}")


print("# --- wp-equianharmonic (u'' = 6u^2, u = WP(z-1; 0, 2)) ---")
pi = mp.pi
e1 = (mp.mpf(1) / 2) ** (mp.mpf(1) / 3)
e2 = e1 * mp.exp(mp.mpc(0, 2 * pi / 3))
e3 = e1 * mp.exp(mp.mpc(0, 4 * pi / 3))
k_sq = (e2 - e3) / (e1 - e3)


def wp(z, c1=-1):
    # WP(w; 0, 2) = e3 + (e1-e3) / sn^2(sqrt(e1-e3)*w, k_sq), w = z + c1.
    # The sn is SQUARED (catalogue closed_form for wp-equianharmonic-g0g2).
    w = z + c1
    return (e3 + (e1 - e3) / mp.ellipfun('sn', mp.sqrt(e1 - e3) * w, k_sq) ** 2).real


def wp_deriv(z, c1=-1):
    return mp.diff(lambda zz: wp(zz, c1), z).real


# IC at z=0 (catalogue ic_bc), and pole-bridge eval points around pole z=1.
show("WP_U0", wp(mp.mpf(0)))
show("WP_UP0", wp_deriv(mp.mpf(0)))
show("WP_U_AT_0p5", wp(mp.mpf('0.5')))
show("WP_UP_AT_0p5", wp_deriv(mp.mpf('0.5')))
show("WP_U_AT_1p05", wp(mp.mpf('1.05')))
show("WP_UP_AT_1p05", wp_deriv(mp.mpf('1.05')))
show("WP_U_AT_1p4", wp(mp.mpf('1.4')))
show("WP_UP_AT_1p4", wp_deriv(mp.mpf('1.4')))
# Cross-check against FW2011 Table 5.1 published values (eq 5.3):
#   u(30) = 1.095098255959744, u(1e4) = 21.02530339471055.
show("WP_CHECK_U30", wp(mp.mpf('30')))
show("WP_CHECK_U1e4", wp(mp.mpf('10000')))

print()
print("# --- tan (u' = 1+u^2 -> tan(z); 2nd-order u''=2u+2u^3) ---")
# IC u(0)=0, u'(0)=1; pole at pi/2; bridge to z=pi where tan(pi)=0.
for zz in ['0.5', '1.0', '2.0', '2.5', '3.0']:
    show(f"TAN_U_AT_{zz.replace('.', 'p')}", mp.tan(mp.mpf(zz)))
# derivative u' = 1 + tan^2 = sec^2
for zz in ['2.0', '2.5', '3.0']:
    t = mp.tan(mp.mpf(zz))
    show(f"TAN_UP_AT_{zz.replace('.', 'p')}", 1 + t * t)

print()
print("# --- coth (u' = 1-u^2 -> coth(z); 2nd-order u''=-2u+2u^3) ---")
# IC at z=-1: u(-1)=coth(-1), u'(-1)=1-coth^2(-1); pole at z=0; bridge to z=+1.
show("COTH_U0_M1", mp.coth(mp.mpf('-1')))
show("COTH_UP0_M1", 1 - mp.coth(mp.mpf('-1')) ** 2)
for zz in ['0.5', '1.0']:
    c = mp.coth(mp.mpf(zz))
    show(f"COTH_U_AT_{zz.replace('.', 'p')}", c)
    show(f"COTH_UP_AT_{zz.replace('.', 'p')}", 1 - c * c)

print()
print("# --- fw-riccati-bessel (u'=t^2+u^2; 2nd-order u''=2t+2t^2 u+2u^3) ---")
nu14 = mp.mpf(1) / 4


def Jm14(x):
    return mp.cos(pi / 4) * mp.besselj(nu14, x) - mp.sin(pi / 4) * mp.bessely(nu14, x)


def y_FW(t):
    t = mp.mpf(t)
    return t * mp.besselj(3 * nu14, t ** 2 / 2) / Jm14(t ** 2 / 2)


def yprime_FW(t):
    t = mp.mpf(t)
    return t ** 2 + y_FW(t) ** 2  # exact via the ODE itself


z1 = mp.findroot(Jm14, 2.0)
show("FW_FIRST_POLE_T", mp.sqrt(2 * z1))
for zz in ['0.5', '1.0', '1.5']:
    show(f"FW_U_AT_{zz.replace('.', 'p')}", y_FW(zz))
    show(f"FW_UP_AT_{zz.replace('.', 'p')}", yprime_FW(zz))
# Value PAST the first pole (t1 = 2.0031): the closed form continues.
show("FW_U_AT_2p1", y_FW('2.1'))
show("FW_UP_AT_2p1", yprime_FW('2.1'))

print()
print("# --- tanh imaginary-pole (u'=1-u^2 -> tanh(z); u''=-2u+2u^3) ---")
# IC u(0)=0, u'(0)=1; poles at i*(pi/2+k pi); imaginary bridge to z=i*pi/4: u=i.
for zz in ['0.25', '0.5', '0.75', '1.0']:
    z = mp.mpc(0, float(mp.pi) * float(zz))  # z = i*(zz*pi)
    t = mp.tanh(z)
    show(f"TANH_U_AT_I_{zz.replace('.', 'p')}PI", t)
    show(f"TANH_UP_AT_I_{zz.replace('.', 'p')}PI", 1 - t * t)

print()
print("# --- jacobi-sn complex bridge (u''=-(5/4)u+(1/2)u^3, u=sn(z,1/4)) ---")
m = mp.mpf(1) / 4
K = mp.ellipk(m)
Kp = mp.ellipk(1 - m)
show("SN_K", K)
show("SN_KP", Kp)
# Exact algebraic ICs at z = i*Kp/2: u = i*sqrt(2), u' = 3/sqrt(2).
zic = mp.mpc(0, float(Kp) / 2)
show("SN_U_IC", mp.ellipfun('sn', zic, m))
# derivative sn' = cn*dn
cn = mp.ellipfun('cn', zic, m)
dn = mp.ellipfun('dn', zic, m)
show("SN_UP_IC", cn * dn)
# real-axis (no-pole) sanity values
for zz in ['0.5', '1.0', '1.5']:
    s = mp.ellipfun('sn', mp.mpf(zz), m)
    cn = mp.ellipfun('cn', mp.mpf(zz), m)
    dn = mp.ellipfun('dn', mp.mpf(zz), m)
    show(f"SN_U_AT_{zz.replace('.', 'p')}", s)
    show(f"SN_UP_AT_{zz.replace('.', 'p')}", cn * dn)
# bridge target: z = i*3*Kp/2 -> u = -i*sqrt(2) (exact)
zt = mp.mpc(0, 3 * float(Kp) / 2)
show("SN_U_AT_I_3KP_2", mp.ellipfun('sn', zt, m))
