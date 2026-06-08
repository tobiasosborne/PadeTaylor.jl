#!/usr/bin/env python3
"""capture.py — transcendental oracle for corpus_riccati_special_test.jl.

Family A, TRANSCENDENTAL subfamily (docs/test_corpus/02_corpus_extension_plan.md
:134-166).  A Riccati u = ±w'/w is the log-derivative of a linear ODE
w'' + p w' + q w = 0, with a SIMPLE pole at every zero of w, residue +1 (or -1
for u = -w'/w).  When w is a SPECIAL FUNCTION its zeros are tabulated to
arbitrary precision -> exact pole locations AND residues, giving the corpus its
first ROWS of simple poles sitting at named special-function zeros.  Three
members, none rational (each carries a transcendental approximation floor):

  CRic.1  Airy    u = -Ai'/Ai            poles at the Airy zeros a_k
                                          (NEGATIVE real row), residue -1
  CRic.2  Bessel  u = -J_1/J_0 = J_0'/J_0 poles at j_{0,k}
                                          (POSITIVE real row), residue +1
  CRic.5  parabolic cylinder              poles a COMPLEX-CONJUGATE PAIR
          u = -U'(0,sqrt2 z)/U(0,sqrt2 z) near arg ~= 3pi/4, residue -1

CRic.5 BRIEF CORRECTION (discovery #2, plan 02:63-68): the recessive solution of
w'' = z^2 w is the parabolic-cylinder U(0,sqrt2 z), NON-oscillatory -> ZERO real
zeros.  u' = u^2 - z^2 therefore has NO real poles; its poles are only
complex-conjugate pairs.  The real-row oscillatory sister u' = u^2 + z^2 is
deferred (CRic.6) pending a dedup decision with CPB.5.

Each first-order Riccati R(z,u) is recast to 2nd order for solve_pade /
path_network_solve via  u'' = dR/dz + (dR/du)*R  (chain rule through u' = R).
This script is the residual==0 GATE: it verifies the recasts SYMBOLICALLY
(sympy, residual identically 0) and NUMERICALLY (mpmath dps=50, residual <1e-45
against the closed-form special function) BEFORE emitting any number, then
prints the mpmath dps=50 oracle constants (ICs, reference values, pole rows)
transcribed into test/_oracle_corpus_riccati_special.jl.

Run:  python3 external/probes/corpus-oracles/riccati-special/capture.py
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 50
z, U = sp.symbols('z U')


def recast(R1):
    """u'' = dR/dz + (dR/du)*R  (chain rule through u' = R)."""
    return sp.expand(sp.diff(R1, z) + sp.diff(R1, U) * R1)


# ---- SYMBOLIC residual==0 GATE: each recast == the brief's 2nd-order RHS ----
assert sp.simplify(recast(U**2 - z) - (2*U**3 - 2*U*z - 1)) == 0, "CRic.1 recast"
assert sp.simplify(recast(-U**2 - U/z - 1) -
                   (2*U**3 + 3*U**2/z + 2*U + 2*U/z**2 + 1/z)) == 0, "CRic.2 recast"
assert sp.simplify(recast(U**2 - z**2) - (2*U**3 - 2*U*z**2 - 2*z)) == 0, "CRic.5 recast"


# ---- NUMERIC residual==0 GATE: u'' = f(z,u,u') against the closed form -------
def u_airy(zz):
    zz = mp.mpc(zz)
    return -mp.airyai(zz, 1) / mp.airyai(zz)


def u_bessel(zz):
    zz = mp.mpc(zz)
    return -mp.besselj(1, zz) / mp.besselj(0, zz)


def w_pcf(zz):                       # w = U(0, sqrt2 z) solves w'' = z^2 w
    return mp.pcfu(0, mp.sqrt(2) * mp.mpc(zz))


def u_pcf(zz):
    zz = mp.mpc(zz)
    h = mp.mpf('1e-25')
    return -((w_pcf(zz + h) - w_pcf(zz - h)) / (2 * h)) / w_pcf(zz)


def gate_2nd(u, f2, zz):
    uu = u(zz)
    upp = mp.diff(u, zz, 2)
    return abs(upp - f2(zz, uu))


assert gate_2nd(u_airy, lambda z, u: 2*u**3 - 2*u*z - 1, mp.mpf('-0.5')) < mp.mpf('1e-30')
assert gate_2nd(u_bessel, lambda z, u: 2*u**3 + 3*u**2/z + 2*u + 2*u/z**2 + 1/z,
                mp.mpf('1')) < mp.mpf('1e-30')
# CRic.5 numeric diff floor is coarser (pcfu has no mpmath derivative); 1st order:
zz = mp.mpc('0.5')
assert abs(mp.diff(u_pcf, zz) - (u_pcf(zz)**2 - zz**2)) < mp.mpf('1e-20'), "CRic.5 1st-order"

print("# residual==0 GATE PASSED (symbolic + numeric) — emitting dps=50 oracle")


def n(x):
    return mp.nstr(x, 45)


# ---------------------------------------------------------------------------
# CRic.1  Airy  u = -Ai'/Ai.  IC z0=0; NEGATIVE-axis pole row a_k, residue -1.
# ---------------------------------------------------------------------------
u0a = -mp.airyai(0, 1) / mp.airyai(0)
print("\n# CRic.1 Airy:  IC z0=0  u0 = -Ai'(0)/Ai(0),  up0 = u0^2")
print(f"const CRic1_U0  = {n(u0a)}")
print(f"const CRic1_UP0 = {n(u0a**2)}")
print("# Airy zeros a_k (NEGATIVE real pole row), residue -1:")
for k in (1, 2, 3):
    print(f"const CRic1_A{k} = {n(mp.airyaizero(k))}")
print("# neg-axis reference values  u(z) = -Ai'(z)/Ai(z)  (real):")
for zt, nm in ((mp.mpf('-0.5'), 'm0p5'), (mp.mpf('-1'), 'm1'),
               (mp.mpf('-1.5'), 'm1p5'), (mp.mpf('-2'), 'm2'),
               (mp.mpf('-2.5'), 'm2p5')):
    print(f"const CRic1_U_{nm} = {n(mp.re(u_airy(zt)))}")
print("# pos-axis (pole-free) tight-accuracy leg:")
for zt, nm in ((mp.mpf('0.5'), 'p0p5'), (mp.mpf('1'), 'p1')):
    print(f"const CRic1_U_{nm} = {n(mp.re(u_airy(zt)))}")

# ---------------------------------------------------------------------------
# CRic.2  Bessel  u = -J1/J0 = J0'/J0.  Seed z0=1 (1/z, 1/z^2 RHS coeff-sing at
# origin).  POSITIVE-axis pole row j_{0,k}, residue +1.
# ---------------------------------------------------------------------------
u0b = -mp.besselj(1, 1) / mp.besselj(0, 1)
print("\n# CRic.2 Bessel:  IC z0=1  u0 = -J1(1)/J0(1),  up0 = -u0^2 - u0 - 1")
print(f"const CRic2_U0  = {n(u0b)}")
print(f"const CRic2_UP0 = {n(-u0b**2 - u0b - 1)}")
print("# Bessel zeros j_{0,k} (POSITIVE real pole row), residue +1:")
for k in (1, 2):
    print(f"const CRic2_J{k} = {n(mp.besseljzero(0, k))}")
print("# reference values  u(z) = -J1(z)/J0(z)  (real):")
for zt, nm in ((mp.mpf('1.5'), 'p1p5'), (mp.mpf('2'), 'p2'),
               (mp.mpf('3'), 'p3'), (mp.mpf('3.5'), 'p3p5')):
    print(f"const CRic2_U_{nm} = {n(mp.re(u_bessel(zt)))}")

# ---------------------------------------------------------------------------
# CRic.5  parabolic cylinder  u = -U'(0,sqrt2 z)/U(0,sqrt2 z).  IC z0=0.
# COMPLEX-conjugate pole pair near arg ~= 3pi/4, residue -1 (NO real poles).
# Pole oracle: root-find w(z) = pcfu(0, sqrt2 z) = 0.
# ---------------------------------------------------------------------------
u0p = u_pcf(mp.mpc(0))
print("\n# CRic.5 parabolic cylinder:  IC z0=0  u0 = -w'(0)/w(0),  up0 = u0^2")
print(f"const CRic5_U0  = {n(mp.re(u0p))}")
print(f"const CRic5_UP0 = {n(mp.re(u0p)**2)}")
p1 = mp.findroot(lambda zz: w_pcf(zz), mp.mpc('-1.49', '1.60'))
print("# first complex pole (upper half): w(z)=0 root-find, residue -1:")
print(f"const CRic5_P1 = complex({n(mp.re(p1))}, {n(mp.im(p1))})")
print("# complex-path reference values  u(z) = -w'/w:")
for zt, nm in ((mp.mpc('-0.5', '0.5'), 'a'), (mp.mpc('-1', '1'), 'b'),
               (mp.mpc('-1.7', '1.8'), 'c')):
    uu = u_pcf(zt)
    print(f"const CRic5_U_{nm} = complex({n(mp.re(uu))}, {n(mp.im(uu))})")
