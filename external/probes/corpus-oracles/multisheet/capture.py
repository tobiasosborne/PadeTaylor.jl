#!/usr/bin/env python3
# capture.py -- mpmath regeneration oracle for the multi-sheet corpus
# (bead padetaylor-5qqr, epic padetaylor-p1v0).
#
# WHY THIS SCRIPT EXISTS (break the self-reference)
# -------------------------------------------------
# The catalogue's P~III / P~V "sheet oracle" values were produced by an
# RK4-in-mpmath recipe; re-deriving them with mpmath here is NOT circular
# against the *package*, whose transformed-RHS closures are an entirely
# independent Float64 implementation of the FFW formulae.  This script is the
# non-self-referential oracle.  It produces:
#
#   1. The sheet-0 P~III / P~V solution values via an mpmath RK4 integration
#      of the transformed ODE along the catalogue winding path, AND the
#      Schwarz-reflection (conjugate-symmetry) residual along the mirror path
#      -- reproducing the catalogue's published sheet-0 numbers to ~15 digits
#      (cross-validation, not self-reference: independent integrator here vs
#      Padé-Taylor stepper in the package).
#   2. Direct mpmath evaluations of the four transformed RHS closures
#      (P~III, P~V, PVI-zeta eq.3, PVI-eta eq.5) at clean interior points.
#      The Julia test asserts the PACKAGE's RHS evaluator equals THESE
#      numbers -- mpmath formula vs Julia formula, two independent typings of
#      the same FFW equation.  This is the airtight B10 RHS check.
#   3. The PVI-eta double-exp cancellation reference (E-1 = e^{e^eta}-1 near
#      the k=1 fixed singularity), used to prove BigFloat-256 recovers what
#      Float64 loses (gap #3).
#
# Ground-truth recipes verbatim from docs/test_corpus/01_corpus_catalogue.md
# "Full candidate records": pIII-tilde-sheet-oracle, pV-tilde-sheet-oracle,
# pVI-eta-cancellation, pVI-eta-subexpr-rhs, pIII-tilde-conjugate-symmetry.
# Transformed equations cross-checked against FFW2017
# references/.../FFW2017_painleve_riemann_surfaces_preprint.md:43-49 (PIII/PV
# eqs), :144 (PVI zeta eq.3), :154 (PVI eta eq.5).
#
# Run:  python3 capture.py     (mpmath 1.3; ~30 s total)

import cmath
import mpmath as mp

pi = mp.pi
j = mp.mpc(0, 1)


def show(label, val):
    print(f"{label} = {mp.nstr(val, 40)}")


def rk4_path(f, y0, legs, n):
    """y' = f(zeta, y) along polyline `legs` (complex endpoints) with n RK4
    steps per leg.  y = [w, w']; f returns [w', w'']."""
    y = [mp.mpc(c) for c in y0]
    for a, b in zip(legs[:-1], legs[1:]):
        a, b = mp.mpc(a), mp.mpc(b)
        h = (b - a) / n
        z = a
        for _ in range(n):
            k1 = f(z, y)
            k2 = f(z + h / 2, [y[i] + h / 2 * k1[i] for i in range(2)])
            k3 = f(z + h / 2, [y[i] + h / 2 * k2[i] for i in range(2)])
            k4 = f(z + h, [y[i] + h * k3[i] for i in range(2)])
            y = [y[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i])
                 for i in range(2)]
            z += h
    return y


# ===========================================================================
# (1) P~III / P~V sheet-0 + conjugate-symmetry  (dps=40, n=4000)
# ===========================================================================
mp.mp.dps = 40

print("# --- P~III sheet-0 (FFW Fig 1: a=b=-1/2, g=1, d=-1; IC w0=1/4,wp0=9/16) ---")
a3, b3, g3, d3 = mp.mpf(-1) / 2, mp.mpf(-1) / 2, mp.mpf(1), mp.mpf(-1)


def f_pIII(z, y):
    w, wp = y
    return [wp, (1 / w) * wp ** 2 + (mp.mpf(1) / 4) * (
        a3 * w ** 2 + g3 * w ** 3 + b3 * mp.e ** z + d3 * mp.e ** (2 * z) / w)]


y0_3 = [mp.mpc(1) / 4, mp.mpc(9) / 16]
s0 = rk4_path(f_pIII, y0_3, [0, 0.5 * j, -1 + 0.5 * j], 4000)        # sheet 0
show("PIII_W_SHEET0", s0[0])
sm = rk4_path(f_pIII, y0_3, [0, -0.5 * j, -1 - 0.5 * j], 4000)       # mirror
show("PIII_CONJ_SYM_ERR", abs(sm[0] - mp.conj(s0[0])))

print()
print("# --- P~V sheet-0 (FFW Fig 6: a=1,b=-1,g=1,d=-1/2; IC w0=2,wp0=-1) ---")
a5, b5, g5, d5 = mp.mpf(1), mp.mpf(-1), mp.mpf(1), mp.mpf(-1) / 2


def f_pV(z, y):
    w, wp = y
    return [wp, (1 / (2 * w) + 1 / (w - 1)) * wp ** 2
            + (w - 1) ** 2 * (a5 * w + b5 / w)
            + g5 * mp.e ** z * w
            + d5 * mp.e ** (2 * z) * w * (w + 1) / (w - 1)]


y0_5 = [mp.mpc(2), mp.mpc(-1)]
v0 = rk4_path(f_pV, y0_5, [0, 0.5 * j, -1 + 0.5 * j], 4000)
show("PV_W_SHEET0", v0[0])
vm = rk4_path(f_pV, y0_5, [0, -0.5 * j, -1 - 0.5 * j], 4000)
show("PV_CONJ_SYM_ERR", abs(vm[0] - mp.conj(v0[0])))

# ===========================================================================
# (2) Transformed-RHS interior-point oracles  (dps=50, no integration)
# ===========================================================================
mp.mp.dps = 50
print()
print("# --- transformed-RHS interior values (independent FFW-formula evals) ---")

# P~III RHS at zeta=0.3+0.7i, w=0.5+0.2i, w'=0.1-0.05i (FFW Fig 1 params)
ez = mp.e ** mp.mpc('0.3', '0.7')
w, wp = mp.mpc('0.5', '0.2'), mp.mpc('0.1', '-0.05')
rhs3 = (1 / w) * wp ** 2 + (mp.mpf(1) / 4) * (
    a3 * w ** 2 + g3 * w ** 3 + b3 * ez + d3 * ez ** 2 / w)
show("PIII_RHS_INTERIOR", rhs3)

# P~V RHS at zeta=0.4+0.6i, w=1.5+0.3i, w'=-0.2+0.1i (FFW Fig 6 params)
ez = mp.e ** mp.mpc('0.4', '0.6')
w, wp = mp.mpc('1.5', '0.3'), mp.mpc('-0.2', '0.1')
rhs5 = (1 / (2 * w) + 1 / (w - 1)) * wp ** 2 + (w - 1) ** 2 * (a5 * w + b5 / w) \
    + g5 * ez * w + d5 * ez ** 2 * w * (w + 1) / (w - 1)
show("PV_RHS_INTERIOR", rhs5)

# PVI-eta RHS at eta=log(log 10)+0.5i (FFW Fig 2 params (4,-4,8,-8))
A, B, G, D = mp.mpf(4), mp.mpf(-4), mp.mpf(8), mp.mpf(-8)
eta = mp.mpc(mp.log(mp.log(mp.mpf(10))), mp.mpf('0.5'))
v, vdot = mp.mpc('0.429534600325223'), mp.mpc('-0.037235820650106481')
E, e_eta = mp.e ** (mp.e ** eta), mp.e ** eta
Em1, vE = E - 1, v - E
tA = (mp.mpf(1) / 2) * (1 / v + 1 / (v - 1) + 1 / vE) * vdot ** 2
tB = -(e_eta * E / Em1 + e_eta * E / vE - 1) * vdot
Pe = A + B * E / v ** 2 + G * Em1 / (v - 1) ** 2 + D * E * Em1 / vE ** 2
tC = v * (v - 1) * vE * e_eta ** 2 / Em1 ** 2 * Pe
show("PVI_ETA_RHS_INTERIOR", tA + tB + tC)

# PVI-zeta RHS at zeta=0.7+0.3i, w=0.42+0.1i, w'=-0.05-0.02i (same params)
zeta = mp.mpc('0.7', '0.3')
w, wp = mp.mpc('0.42', '0.1'), mp.mpc('-0.05', '-0.02')
ez = mp.e ** zeta
ezm1, wm1, wmez = ez - 1, w - 1, w - ez
t1 = (mp.mpf(1) / 2) * (1 / w + 1 / wm1 + 1 / wmez) * wp ** 2
t2 = -(ez / ezm1 + ez / wmez) * wp
Pz = A + B * ez / w ** 2 + G * ezm1 / wm1 ** 2 + D * ez * ezm1 / wmez ** 2
t3 = w * wm1 * wmez / ezm1 ** 2 * Pz
show("PVI_ZETA_RHS_INTERIOR", t1 + t2 + t3)

# ===========================================================================
# (3) PVI-eta double-exp cancellation reference  (dps=60)  [gap #3]
# ===========================================================================
mp.mp.dps = 60
print()
print("# --- PVI-eta double-exp cancellation E-1 near k=1 singularity ---")
eta_sing = mp.log(2 * pi) + j * pi / 2
eta_test = eta_sing + mp.mpf('1e-10')
Em1 = mp.e ** (mp.e ** eta_test) - 1
show("ETA_EM1_REAL", Em1.real)
show("ETA_EM1_IMAG", Em1.imag)
ef = complex(float(eta_test.real), float(eta_test.imag))
Em1_f64 = cmath.exp(cmath.exp(ef)) - 1
print(f"# Float64 cmath E-1 = {Em1_f64!r}")
print(f"# Float64 real-part rel error vs mpmath = "
      f"{abs((Em1_f64.real - float(Em1.real)) / float(Em1.real)):.4e}")
