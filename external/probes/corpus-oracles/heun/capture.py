#!/usr/bin/env python3
# capture.py -- INDEPENDENT mpmath oracle for the Heun corpus
# (bead padetaylor-h5ho, epic padetaylor-p1v0).
#
# WHAT THIS IS
# ------------
# A from-scratch re-implementation of the HeunG / HeunC three-term Frobenius
# recurrences (DLMF 31.3.2-31.3.4 for HeunG; Motygin 2015 eqs (3)-(4) for
# HeunC) in mpmath at 50 decimal places.  It does NOT call the Julia heun_g /
# heun_c -- that would be self-referential.  It is the EXTERNAL ground truth
# against which test/corpus_heun_test.jl pins.
#
# It also re-derives / cross-checks the closed-form reductions:
#   * heung-2f1-reduction : HeunG(2,2,1,2,3/2,1;z) =?= 2F1(1/2,1;3/2;2z-z^2)
#   * heung-pole-at-z1    : HeunG(2,4,1,4,2,2;z)   =?= 1/(1-z)   (all c_n=1)
#   * heung-polynomial    : HeunG(2,q*,-1,2,3,4;z) =?= 1 + (q*/6) z, q*=-6+-2sqrt6
#   * heung-trivial       : HeunG(a,0,0,b,g,d;z)   =?= 1
# against the mpmath hypergeometric native + the recurrence, so each closed
# form is verified BEFORE being used as an oracle (Law 1: never trust the
# catalogue's numeric value without an independent check).
#
# For the OFF-monoculture numeric Heun values (complex-a, |a|<1, |a|>1,
# eps=-3/+2, fractional gamma) the recurrence IS the oracle; we additionally
# print a wolframscript cross-check command at the bottom (run separately).
#
# Run:  python3 capture.py    (mpmath 1.3 installed)

import mpmath as mp

mp.mp.dps = 50


def show(label, val):
    if isinstance(val, mp.mpc):
        print(f"{label} = {mp.nstr(val.real, 40)} {'+' if val.imag >= 0 else '-'} "
              f"{mp.nstr(abs(val.imag), 40)}j")
    else:
        print(f"{label} = {mp.nstr(val, 40)}")


# ---------------------------------------------------------------------------
# HeunG Frobenius recurrence (DLMF 31.3.2-31.3.4), INDEPENDENT re-derivation.
#   c0 = 1, c1 = q/(a*gamma)
#   R_j c_{j+1} = (Q_j + q) c_j - P_j c_{j-1},   eps = a+b+1-g-d
#   P_j = (j-1+alpha)(j-1+beta)
#   Q_j = j*((j-1+gamma)(1+a) + a*delta + eps)
#   R_j = a (j+1)(j+gamma)
# ---------------------------------------------------------------------------
def heung_series(a, q, al, be, ga, de, z, N=300):
    a, q, al, be, ga, de, z = map(mp.mpmathify, (a, q, al, be, ga, de, z))
    ep = al + be + 1 - ga - de
    c = [mp.mpf(0)] * (N + 1)
    c[0] = mp.mpf(1)
    c[1] = q / (a * ga)
    for j in range(1, N):
        P = (j - 1 + al) * (j - 1 + be)
        Q = j * ((j - 1 + ga) * (1 + a) + a * de + ep)
        R = a * (j + 1) * (j + ga)
        c[j + 1] = ((Q + q) * c[j] - P * c[j - 1]) / R
    return sum(c[k] * z**k for k in range(N + 1)), c


# ---------------------------------------------------------------------------
# HeunC Frobenius recurrence (Motygin 2015 eqs (3)-(4) == src/Heun.jl form).
#   b0 = 1, b1 = -q/gamma
#   P_n b_n = Q_n b_{n-1} + R_n b_{n-2}
#   P_n = n(n+gamma-1)
#   Q_n = -q + (n-1)(gamma+delta-eps+n-2)
#   R_n = alpha + (n-2) eps
# ---------------------------------------------------------------------------
def heunc_series(q, al, ga, de, ep, z, N=400):
    q, al, ga, de, ep, z = map(mp.mpmathify, (q, al, ga, de, ep, z))
    b = [mp.mpf(0)] * (N + 1)
    b[0] = mp.mpf(1)
    b[1] = -q / ga
    for n in range(2, N + 1):
        P = n * (n + ga - 1)
        Q = -q + (n - 1) * (ga + de - ep + n - 2)
        R = al + (n - 2) * ep
        b[n] = (Q * b[n - 1] + R * b[n - 2]) / P
    return sum(b[k] * z**k for k in range(N + 1)), b


print("# === CHN.1  heung-2f1-reduction: HeunG(2,2,1,2,3/2,1;z) = 2F1(1/2,1;3/2;2z-z^2) ===")
half = mp.mpf(1) / 2
for zt in ('0.1', '0.3', '0.5', '0.7'):
    z = mp.mpf(zt)
    hv, _ = heung_series(2, 2, 1, 2, mp.mpf(3) / 2, 1, z)
    fv = mp.hyp2f1(half, 1, mp.mpf(3) / 2, 2 * z - z * z)
    show(f"HEUNG_2F1_z{zt}", hv)
    print(f"#   2F1 oracle  z={zt}: {mp.nstr(fv, 40)}   |diff|={mp.nstr(abs(hv-fv),5)}")

print("# === CHN.2  heung-pole-at-z1: HeunG(2,4,1,4,2,2;z) = 1/(1-z), all c_n=1 ===")
_, cpole = heung_series(2, 4, 1, 4, 2, 2, mp.mpf('0.3'))
print(f"#   first 10 coeffs (should all be 1): {[mp.nstr(cpole[k],4) for k in range(10)]}")
for zt in ('0.3', '0.5', '0.9'):
    z = mp.mpf(zt)
    hv, _ = heung_series(2, 4, 1, 4, 2, 2, z)
    show(f"HEUNG_POLE_z{zt}", hv)
    print(f"#   1/(1-z)     z={zt}: {mp.nstr(1/(1-z), 40)}   |diff|={mp.nstr(abs(hv-1/(1-z)),5)}")

print("# === CHN.3a heung-polynomial-degree1: HeunG(2,q*,-1,2,3,4;z)=1+(q*/6)z, q*=-6+-2sqrt6 ===")
for sgn, tag in ((+1, 'p'), (-1, 'm')):
    qstar = -6 + sgn * 2 * mp.sqrt(6)
    show(f"HEUNG_POLY_q{tag}", qstar)
    for zt in ('0.5', '0.9'):
        z = mp.mpf(zt)
        hv, cc = heung_series(2, qstar, -1, 2, 3, 4, z)
        cf = 1 + (qstar / 6) * z
        show(f"HEUNG_POLY_{tag}_z{zt}", hv)
        print(f"#   linear      z={zt}: {mp.nstr(cf,40)}  |diff|={mp.nstr(abs(hv-cf),5)} "
              f" c2={mp.nstr(cc[2],5)} (should be 0)")

print("# === CHN.3b heung-trivial-alpha0-q0: HeunG(a,0,0,b,g,d;z)=1 ===")
for zt in ('0.3', '0.7'):
    hv, _ = heung_series(2, 0, 0, 2, 3, 4, mp.mpf(zt))
    show(f"HEUNG_TRIVIAL_z{zt}", hv)

print("# === CHN.4a heung-complex-a: HeunG(1/2+i/2,1,1,2,3,4;0.3) [EXACT-RATIONAL a] ===")
hv, _ = heung_series(mp.mpc(mp.mpf(1) / 2, mp.mpf(1) / 2), 1, 1, 2, 3, 4, mp.mpf('0.3'))
show("HEUNG_COMPLEX_A_z0p3", hv)

print("# === CHN.4b heung-small-a: HeunG(1/2,1,1,2,3,4;z) (|a|<1, R=1/2) ===")
for zt in ('0.1', '0.3', '0.45'):
    hv, _ = heung_series(mp.mpf(1) / 2, 1, 1, 2, 3, 4, mp.mpf(zt))
    show(f"HEUNG_SMALL_A_z{zt}", hv)

print("# === CHN.4c heung-large-q-a3: HeunG(3,10,1,2,2,3;z) (|a|>1, large q) ===")
for zt in ('0.25', '0.5'):
    hv, _ = heung_series(3, 10, 1, 2, 2, 3, mp.mpf(zt))
    show(f"HEUNG_LARGE_Q_z{zt}", hv)

print("# === CHN.5  heunc-epsilon-neg3 / -pos2 / -fractional-gamma ===")
for tag, (q, al, ga, de, ep), zs in (
        ('NEG3', (1, 1, 2, 3, -3), ('0.1', '0.5', '0.9')),
        ('POS2', (1, 1, 2, 3, 2), ('0.1', '0.5')),
        ('FRACG', (1, 1, mp.mpf(1) / 2, 3, -1), ('0.1', '0.5'))):
    for zt in zs:
        hv, _ = heunc_series(q, al, ga, de, ep, mp.mpf(zt))
        show(f"HEUNC_{tag}_z{zt}", hv)

print("# === CSF special-fn IVP oracle constants (mpmath natives) ===")
show("AIRY_AI0", mp.airyai(0))
show("AIRY_AIP0", mp.airyai(0, derivative=1))
for zt in ('1', '2', '3', '-1', '-2'):
    show(f"AIRY_z{zt}", mp.airyai(mp.mpf(zt)))
show("BESSEL_J0_1", mp.besselj(0, 1))
show("BESSEL_J1_1", mp.besselj(1, 1))
for zt in ('2', '3', '5'):
    show(f"BESSEL_J0_z{zt}", mp.besselj(0, mp.mpf(zt)))
a2, b2, c2 = mp.mpf(1) / 2, mp.mpf(1) / 3, mp.mpf(5) / 4
for zt in ('0.4', '0.7', '0.9'):
    show(f"GAUSS2F1_z{zt}", mp.hyp2f1(a2, b2, c2, mp.mpf(zt)))
show("GAUSS2F1_cplx", mp.hyp2f1(a2, b2, c2, mp.mpc('0.5', '0.5')))
# 2f1-elementary-reduction: F(1/2,3;3;z) = (1-z)^(-1/2)
for zt in ('0.4', '0.9'):
    z = mp.mpf(zt)
    fv = mp.hyp2f1(half, 3, 3, z)
    print(f"#   2F1-elem z={zt}: hyp2f1={mp.nstr(fv,40)} (1-z)^-1/2={mp.nstr((1-z)**(-half),40)}"
          f" |diff|={mp.nstr(abs(fv-(1-z)**(-half)),5)}")

print()
print("# --- wolframscript cross-checks (run separately to confirm trust) ---")
print('# wolframscript -code "N[HeunG[0.5+0.5*I,1,1,2,3,4,0.3],40]"')
print('# wolframscript -code "N[HeunG[1/2,1,1,2,3,4,0.3],40]"')
print('# wolframscript -code "N[HeunG[3,10,1,2,2,3,0.25],40]"')
print('# wolframscript -code "N[HeunC[1,1,2,3,-3,0.1],40]"')
print('# wolframscript -code "N[HeunC[1,1,2,3,2,0.1],40]"')
print('# wolframscript -code "N[HeunC[1,1,1/2,3,-1,0.1],40]"')
print('# wolframscript -code "N[HeunG[2,2,1,2,3/2,1,0.3],40]"')
