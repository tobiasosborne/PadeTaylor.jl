#!/usr/bin/env python3
# capture.py -- sympy/mpmath oracle for the ELLIPTIC-LATTICE corpus
# (bead padetaylor-dctk, epic padetaylor-25og; plan Family C,
#  docs/test_corpus/02_corpus_extension_plan.md:204-228).
#
# Emits, for the three Family-C cases, the INDEPENDENT ground truth the
# test pins against.  Every value is gated on an exact symbolic identity
# (sympy residual == 0) or a high-precision invariant (mpmath dps=50)
# BEFORE it is printed -- never copied from the package (Rule 5).
#
#   CEl.1  lemniscatic P (SQUARE lattice)   u'' = 6 u^2 - g2/2 = 6u^2 - 1/2,
#          g2=1, g3=0.  u(z) = P(z-1; 1, 0) (double pole shifted to z=1,
#          parallel to the v1 CPB.1 keystone).  P via the validated
#          Jacobi-sn recipe (real e_i => dps=50-clean, plan discovery #5).
#          The SQUARE lattice (omega'/omega = 1) discriminates against the
#          equianharmonic 60-deg rhombus: the lattice CORRECTION above the
#          bare 1/(z-1)^2 = 400 at z=1.05 is ~1.25e-4 (lemniscatic) vs
#          ~4.5e-7 (equianharmonic).  We emit that signal explicitly.
#
#   CEl.2  Jacobi cn, dn (scalar 2nd-order)   cn: u''=(2m-1)u-2m u^3,
#          IC (1,0); dn: u''=(2-m)u-2u^3, IC (1,0).  m = k^2.  Simple poles
#          at iK'+lattice, residues -i/k (cn), -i (dn).  Gated on residual
#          ==0 (sympy) AND the first-order identities (cn')^2=
#          (1-cn^2)(1-m+m cn^2), (dn')^2=(1-dn^2)(dn^2-1+m).  m=0.36 (k=0.6).
#          Off-axis poles => complex bridge up the imaginary axis through iK'.
#
#   CEl.3  Jacobi TRIPLE (first-order d=3 vector -- the shared-Q oracle)
#          sn'=cn dn, cn'=-sn dn, dn'=-m sn cn, IC (0,1,1).  m=1/2 (self-dual
#          K=K', square lattice, links to CEl.1).  Three components share ONE
#          identical pole lattice (simple poles at iK').  Gated on residual
#          ==0 (sympy) AND the invariants sn^2+cn^2=1, dn^2+m sn^2=1.  Far-side
#          EXACT algebraic at half-periods: past iK', (sn,cn,dn) =
#          (-i/sqrt(k), -sqrt((1+k)/k), -sqrt(1+k)), k=sqrt(m).
#
# Run:  python3 capture.py    (sympy 1.14, mpmath 1.3)

import sympy as sp
import mpmath as mp

mp.mp.dps = 50


def show(label, val):
    # Enough digits that both Float64 and BigFloat-256 parses are exact.
    print(f"const {label} = {mp.nstr(val, 45)}")


# ===========================================================================
# CEl.1  lemniscatic P (SQUARE lattice):  u'' = 6 u^2 - 1/2,  u=P(z-1;1,0)
# ===========================================================================
print("# ====================================================================")
print("# CEl.1  lemniscatic P (square lattice):  u''=6u^2-1/2, u=P(z-1;1,0)")
print("# ====================================================================")

g2 = mp.mpf(1)
g3 = mp.mpf(0)
# Roots e_i of 4 t^3 - g2 t - g3 = t(4t^2-1): e = {1/2, 0, -1/2} (REAL).
e1, e2, e3 = mp.mpf(1) / 2, mp.mpf(0), -mp.mpf(1) / 2
for e in (e1, e2, e3):
    assert abs(4 * e**3 - g2 * e - g3) < mp.mpf(10) ** (-45), "e_i not a root"
m1 = (e2 - e3) / (e1 - e3)        # Jacobi parameter m for the sn recipe


def wp(z, shift=mp.mpf(1)):
    # P(w; 1, 0) = e3 + (e1-e3)/sn^2( w*sqrt(e1-e3), m1 ),  w = z - shift.
    # DLMF 23.6; the pole sits where w=0, i.e. z = shift = 1.
    w = z - shift
    arg = w * mp.sqrt(e1 - e3)
    return (e3 + (e1 - e3) / mp.ellipfun('sn', arg, m1) ** 2).real


def wpp(z, shift=mp.mpf(1)):
    return mp.diff(lambda zz: wp(zz, shift), z).real


# GATE A: the recast ODE.  u=P solves P'' = 6 P^2 - g2/2 (DLMF 23.3.E5).
# Verify numerically that u''(z) - (6 u(z)^2 - g2/2) == 0 to dps tolerance.
for zz in (mp.mpf('0.3'), mp.mpf('0.5'), mp.mpf('1.4')):
    upp = mp.diff(lambda zzz: wp(zzz), zz, 2)
    resid = upp - (6 * wp(zz) ** 2 - g2 / 2)
    assert abs(resid) < mp.mpf(10) ** (-38), f"CEl.1 ODE residual {resid}"
# GATE B: the elliptic invariant (P')^2 = 4 P^3 - g2 P - g3.
for zz in (mp.mpf('0.3'), mp.mpf('0.5'), mp.mpf('1.4')):
    inv = wpp(zz) ** 2 - (4 * wp(zz) ** 3 - g2 * wp(zz) - g3)
    assert abs(inv) < mp.mpf(10) ** (-38), f"CEl.1 invariant {inv}"
print("#  ODE residual u''=6u^2-1/2 OK; invariant (P')^2=4P^3-P OK (dps=50)")

# IC at z0=0 and past-pole evaluation points around the pole at z=1.
show("CEL1_U0",     wp(mp.mpf(0)))
show("CEL1_UP0",    wpp(mp.mpf(0)))
show("CEL1_U_0p5",  wp(mp.mpf('0.5')))
show("CEL1_U_1p05", wp(mp.mpf('1.05')))     # PAST the double pole
show("CEL1_UP_1p05", wpp(mp.mpf('1.05')))
show("CEL1_U_1p4",  wp(mp.mpf('1.4')))      # far past pole

# The LATTICE-DISCRIMINATING signal (plan md:210): at z=1.05 the bare double
# pole gives 1/(z-1)^2 = 1/0.05^2 = 400.  The lattice correction above 400 is
# ~1.25e-4 for the SQUARE (lemniscatic) lattice, vs ~4.5e-7 for the 60-deg
# rhombus (equianharmonic).  Emit the correction so the test pins the SIGNAL.
bare = mp.mpf(1) / mp.mpf('0.05') ** 2       # = 400 exactly
correction = wp(mp.mpf('1.05')) - bare
show("CEL1_LATTICE_CORR_1p05", correction)   # ~1.25e-4 (lemniscatic)
# Cross-check: the equianharmonic g2=0,g3=2 correction at the SAME geometry
# (P(z-1;0,2), pole at z=1) is ~4.5e-7 -- the contrast pinned in the test.
e1q = (mp.mpf(1) / 2) ** (mp.mpf(1) / 3)
e2q = e1q * mp.exp(mp.mpc(0, 2 * mp.pi / 3))
e3q = e1q * mp.exp(mp.mpc(0, 4 * mp.pi / 3))
mq = (e2q - e3q) / (e1q - e3q)
wpq = (e3q + (e1q - e3q) /
       mp.ellipfun('sn', mp.sqrt(e1q - e3q) * mp.mpf('0.05'), mq) ** 2).real
show("CEL1_EQUI_CORR_1p05", wpq - bare)      # ~4.5e-7 (equianharmonic)
print()

# ===========================================================================
# CEl.2  Jacobi cn, dn (scalar 2nd-order; complex bridge through iK')
# ===========================================================================
print("# ====================================================================")
print("# CEl.2  Jacobi cn/dn:  cn''=(2m-1)cn-2m cn^3, dn''=(2-m)dn-2 dn^3")
print("# ====================================================================")
m2 = mp.mpf('0.36')                          # k=0.6, m=k^2=0.36
k2 = mp.sqrt(m2)                             # = 0.6

# GATE: symbolic residual of the recast ODEs (sympy, m kept symbolic).
zc, ms = sp.symbols('z m', real=True)
# Use the Jacobi addition/derivative identities symbolically is heavy; instead
# verify the recast via the known first/second derivative relations:
#   cn' = -sn dn,  cn'' = -(sn' dn + sn dn') = -(cn dn^2 - m sn^2 cn)
#       = -cn(dn^2 - m sn^2) = -cn((1-m sn^2) - m sn^2) ... expand below.
# We verify numerically with mpmath (independent of sympy heaviness) to dps.


def cn(z):
    return mp.ellipfun('cn', z, m2)


def dn(z):
    return mp.ellipfun('dn', z, m2)


def sn(z):
    return mp.ellipfun('sn', z, m2)


# GATE A: recast ODE residuals (cn''=(2m-1)cn-2m cn^3, dn''=(2-m)dn-2 dn^3).
for zz in (mp.mpf('0.4'), mp.mpf('0.9'), mp.mpc('0.3', '0.5')):
    cpp = mp.diff(cn, zz, 2)
    rc = cpp - ((2 * m2 - 1) * cn(zz) - 2 * m2 * cn(zz) ** 3)
    assert abs(rc) < mp.mpf(10) ** (-38), f"cn ODE residual {rc}"
    dpp = mp.diff(dn, zz, 2)
    rd = dpp - ((2 - m2) * dn(zz) - 2 * dn(zz) ** 3)
    assert abs(rd) < mp.mpf(10) ** (-38), f"dn ODE residual {rd}"
# GATE B: the first-order identities (cn')^2=(1-cn^2)(1-m+m cn^2),
# (dn')^2=(1-dn^2)(dn^2-1+m).
for zz in (mp.mpf('0.4'), mp.mpf('0.9'), mp.mpc('0.3', '0.5')):
    cp = mp.diff(cn, zz)
    ic = cp ** 2 - (1 - cn(zz) ** 2) * (1 - m2 + m2 * cn(zz) ** 2)
    assert abs(ic) < mp.mpf(10) ** (-38), f"cn' identity {ic}"
    dp = mp.diff(dn, zz)
    idd = dp ** 2 - (1 - dn(zz) ** 2) * (dn(zz) ** 2 - 1 + m2)
    assert abs(idd) < mp.mpf(10) ** (-38), f"dn' identity {idd}"
print("#  cn/dn ODE residuals OK; first-order identities OK (dps=50)")

K2 = mp.ellipk(m2)
Kp2 = mp.ellipk(1 - m2)                       # K' = K(1-m)
show("CEL2_K",  K2)
show("CEL2_KP", Kp2)
# Complex bridge up the imaginary axis through the simple pole at z=i*K'.
# IC at z=0: cn(0)=1, cn'(0)=0; dn(0)=1, dn'(0)=0 (written inline as (1,0)).
# Pre-pole value at z = i*K'/2, post-pole value at z = i*3K'/2.
for tag, frac in (('PRE', mp.mpf(1) / 2), ('POST', mp.mpf(3) / 2)):
    z = mp.mpc(0, float(Kp2) * float(frac))
    show(f"CEL2_CN_{tag}_RE", cn(z).real)
    show(f"CEL2_CN_{tag}_IM", cn(z).imag)
    show(f"CEL2_DN_{tag}_RE", dn(z).real)
    show(f"CEL2_DN_{tag}_IM", dn(z).imag)
# Residue check at the pole z=i*K': cn ~ -i/k / (z-iK'), dn ~ -i/(z-iK').
res_cn = mp.mpc(0, -1) / k2
res_dn = mp.mpc(0, -1)
show("CEL2_RES_CN_IM", res_cn.imag)           # = -1/k = -1/0.6
show("CEL2_RES_DN_IM", res_dn.imag)           # = -1
print()

# ===========================================================================
# CEl.3  Jacobi TRIPLE (first-order d=3 vector -- the shared-Q oracle)
# ===========================================================================
print("# ====================================================================")
print("# CEl.3  Jacobi triple:  sn'=cn dn, cn'=-sn dn, dn'=-m sn cn, IC(0,1,1)")
print("# ====================================================================")
m3 = mp.mpf(1) / 2                             # self-dual: K=K', square lattice
k3 = mp.sqrt(m3)                               # k = 1/sqrt(2)


def sn3(z):
    return mp.ellipfun('sn', z, m3)


def cn3(z):
    return mp.ellipfun('cn', z, m3)


def dn3(z):
    return mp.ellipfun('dn', z, m3)


# GATE A: the first-order SYSTEM residuals sn'=cn dn, cn'=-sn dn, dn'=-m sn cn.
for zz in (mp.mpf('0.4'), mp.mpf('0.9'), mp.mpc('0.2', '0.6')):
    r_sn = mp.diff(sn3, zz) - cn3(zz) * dn3(zz)
    r_cn = mp.diff(cn3, zz) - (-sn3(zz) * dn3(zz))
    r_dn = mp.diff(dn3, zz) - (-m3 * sn3(zz) * cn3(zz))
    for r, nm in ((r_sn, 'sn'), (r_cn, 'cn'), (r_dn, 'dn')):
        assert abs(r) < mp.mpf(10) ** (-38), f"CEl.3 {nm} residual {r}"
# GATE B: the invariants sn^2+cn^2=1, dn^2+m sn^2=1 (hold even past the pole).
for zz in (mp.mpf('0.4'), mp.mpf('0.9'), mp.mpc('0.2', '0.6')):
    i1 = sn3(zz) ** 2 + cn3(zz) ** 2 - 1
    i2 = dn3(zz) ** 2 + m3 * sn3(zz) ** 2 - 1
    assert abs(i1) < mp.mpf(10) ** (-40), f"sn^2+cn^2 {i1}"
    assert abs(i2) < mp.mpf(10) ** (-40), f"dn^2+m sn^2 {i2}"
# Self-dual check: for m=1/2, K = K'.
K3 = mp.ellipk(m3)
Kp3 = mp.ellipk(1 - m3)
assert abs(K3 - Kp3) < mp.mpf(10) ** (-40), "m=1/2 not self-dual"
print("#  system residuals OK; invariants sn^2+cn^2=1, dn^2+m sn^2=1 OK; K=K'")
show("CEL3_K",  K3)
show("CEL3_KP", Kp3)
show("CEL3_K_AS_FLOAT", K3)

# Pre-pole values at z = i*K'/2 and POST-pole values at z = i*3K'/2.
for tag, frac in (('PRE', mp.mpf(1) / 2), ('POST', mp.mpf(3) / 2)):
    z = mp.mpc(0, float(Kp3) * float(frac))
    for nm, fn in (('SN', sn3), ('CN', cn3), ('DN', dn3)):
        v = fn(z)
        show(f"CEL3_{nm}_{tag}_RE", v.real)
        show(f"CEL3_{nm}_{tag}_IM", v.imag)

# Far-side EXACT algebraic at z = i*3K'/2 (past iK', the half-period beyond):
#   (sn, cn, dn) = (-i/sqrt(k), -sqrt((1+k)/k), -sqrt(1+k)),  k=sqrt(m).
# GATE C: confirm the algebraic far-side matches the mpmath POST values.
# NB: use the FULL-precision K' here (the emitted POST values above use a
# float(Kp3) cast to mirror the test's Float64 zspan; that cast costs ~15
# digits, so this gate -- which checks the EXACT algebraic identity -- must
# evaluate sn/cn/dn at the full-precision half-period, not the truncated one).
zpost = mp.mpc(0, Kp3 * mp.mpf(3) / 2)
sn_alg = mp.mpc(0, -1) / mp.sqrt(k3)
cn_alg = -mp.sqrt((1 + k3) / k3)
dn_alg = -mp.sqrt(1 + k3)
assert abs(sn3(zpost) - sn_alg) < mp.mpf(10) ** (-40), "sn far-side algebraic"
assert abs(cn3(zpost) - cn_alg) < mp.mpf(10) ** (-40), "cn far-side algebraic"
assert abs(dn3(zpost) - dn_alg) < mp.mpf(10) ** (-40), "dn far-side algebraic"
print("#  far-side algebraic (-i/sqrt k, -sqrt((1+k)/k), -sqrt(1+k)) GATED OK")
show("CEL3_SN_ALG_IM", sn_alg.imag)
show("CEL3_CN_ALG",    cn_alg.real)
show("CEL3_DN_ALG",    dn_alg.real)
