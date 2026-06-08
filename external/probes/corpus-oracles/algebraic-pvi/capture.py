#!/usr/bin/env python3
"""capture.py — sympy/mpmath oracle for corpus_algebraic_pvi_test.jl.

Family F (docs/test_corpus/02_corpus_extension_plan.md:274-289): the FIRST
NON-self-referential PVI oracle (v1 PVI figure values are all self-pinned,
capability-map gap #9).  Two EXTERNAL algebraic solutions of the standard
sixth-Painleve (PVI) ODE, each gated by a sympy symbolic residual==0 check
BEFORE any value is emitted, then pinned numerically by mpmath at dps=50.

THE EXACT PVI ODE USED (auditable wiring; standard form, e.g. Iwasaki-Kimura-
Shimomura-Yoshida, or FFW2017 eq. preceding the z=e^zeta transform):

  u'' = (1/2)(1/u + 1/(u-1) + 1/(u-z)) (u')^2
        - (1/z + 1/(z-1) + 1/(u-z)) u'
        + u(u-1)(u-z)/(z^2 (z-1)^2) *
          ( alpha + beta*z/u^2 + gamma*(z-1)/(u-1)^2 + delta*z(z-1)/(u-z)^2 )

  CPvi.1  u(z) = z^(-1/2)        solves PVI with (a,b,g,d)=(0, 0, 1/8, -5/8).
  CPvi.2  u(z) = (sqrt(z)+1)/2   solves PVI with (a,b,g,d)=(1/2, 0, 1/8, -5/8).

Both have a branch point ONLY at z=0 (which coincides with PVI's OWN fixed
singularity) -> two sheets +-sqrt; NO generic (movable, off-{0,1,inf}) branch.

THE PACKAGE FRAME (what corpus_algebraic_pvi_test.jl actually pins against).
`pVI_transformed_rhs(a,b,g,d)` is the zeta-PLANE closure `(zeta,w,wp) -> w''`
(src/SheetTracker.jl:155, FFW2017 eq. 3, md:144) under z=e^zeta, w(zeta)=u(z):

      zeta = log z,   w = u,   wp = dw/dzeta = z*u'(z),
      w'' = d^2 w / d zeta^2 = z*u'(z) + z^2 * u''(z).

So the test feeds the package the algebraic state (zeta=log z, w=u(z),
wp=z*u'(z)) and asserts the package's w'' equals the ANALYTIC
w'' = z*u' + z^2*u'' computed here independently from the closed form.  That
pins the (a,b,g,d) wiring + the e^zeta Jacobian against an external algebraic
solution -- exactly the gap-#9 non-self-reference check.

Probe points: z=4 (real) and z=2+3i (complex), BOTH sheets (+- branch of sqrt).
At z=4: u=z^(-1/2) -> +1/2 (other sheet -1/2); the complex point likewise has
its negated partner for CPvi.1 (u is odd in sqrt(z)).

Run:  python3 external/probes/corpus-oracles/algebraic-pvi/capture.py
"""
import sympy as sp
import mpmath as mp

mp.mp.dps = 50


def emit(label, val):
    if isinstance(val, mp.mpc):
        print(f"const {label} = complex({mp.nstr(val.real, 45)}, "
              f"{mp.nstr(val.imag, 45)})")
    else:
        print(f"const {label} = {mp.nstr(val, 45)}")


# ===========================================================================
# Symbolic PVI residual gate.  Build the standard PVI ODE operator and confirm
# each algebraic solution annihilates it IDENTICALLY (sympy simplify -> 0).
# ===========================================================================
z = sp.symbols('z')
a, b, g, d = sp.symbols('alpha beta gamma delta')


def pvi_residual(u_expr, av, bv, gv, dv):
    """u'' - RHS(u,u',z; a,b,g,d) for the standard PVI ODE, fully substituted."""
    u = u_expr
    up = sp.diff(u, z)
    upp = sp.diff(u, z, 2)
    kinetic = sp.Rational(1, 2) * (1 / u + 1 / (u - 1) + 1 / (u - z)) * up**2
    friction = -(1 / z + 1 / (z - 1) + 1 / (u - z)) * up
    param = (av + bv * z / u**2 + gv * (z - 1) / (u - 1)**2 +
             dv * z * (z - 1) / (u - z)**2)
    nonlinear = u * (u - 1) * (u - z) / (z**2 * (z - 1)**2) * param
    rhs = kinetic + friction + nonlinear
    return sp.simplify(upp - rhs)


# CPvi.1 : u = z^(-1/2), (a,b,g,d) = (0, 0, 1/8, -5/8).
u1 = z**sp.Rational(-1, 2)
res1 = pvi_residual(u1, sp.Integer(0), sp.Integer(0),
                    sp.Rational(1, 8), sp.Rational(-5, 8))
assert res1 == 0, f"CPvi.1 PVI residual != 0: {res1}"

# CPvi.2 : u = (sqrt(z)+1)/2, (a,b,g,d) = (1/2, 0, 1/8, -5/8).
u2 = (sp.sqrt(z) + 1) / 2
res2 = pvi_residual(u2, sp.Rational(1, 2), sp.Integer(0),
                    sp.Rational(1, 8), sp.Rational(-5, 8))
assert res2 == 0, f"CPvi.2 PVI residual != 0: {res2}"

print("# PVI residual==0 GATE PASSED for BOTH algebraic solutions:")
print("#   CPvi.1  u=z^(-1/2)      (a,b,g,d)=(0,0,1/8,-5/8)   residual == 0")
print("#   CPvi.2  u=(sqrt(z)+1)/2 (a,b,g,d)=(1/2,0,1/8,-5/8) residual == 0")
print("# PVI ODE used: u'' = (1/2)(1/u+1/(u-1)+1/(u-z))u'^2")
print("#               - (1/z+1/(z-1)+1/(u-z))u'")
print("#               + u(u-1)(u-z)/(z^2(z-1)^2)*")
print("#                 (a + b z/u^2 + g(z-1)/(u-1)^2 + d z(z-1)/(u-z)^2)")
print()


# ===========================================================================
# mpmath dps=50 numeric pins.  For each case + each sheet + each probe point we
# emit u, u', and the ZETA-FRAME second derivative wpp = z*u' + z^2*u'' that the
# package's pVI_transformed_rhs must reproduce, PLUS the zeta-frame state
# (wp = z*u') the test feeds in.  sheet = +1 (principal +sqrt) or -1 (-sqrt).
# ===========================================================================
def state_cpvi1(zv, sheet):
    """u=z^(-1/2) on the given sqrt sheet.  Returns (u, up, upp, wp, wpp)."""
    s = mp.mpf(1) if sheet > 0 else mp.mpf(-1)
    rt = s * mp.sqrt(zv)                 # chosen sqrt branch
    u = 1 / rt                           # z^(-1/2) = 1/sqrt(z)
    up = -sp.Rational(1, 2)              # placeholder; compute explicitly below
    up = mp.mpf(-1) / 2 * rt / zv**2     # d/dz z^(-1/2) = -1/2 z^(-3/2) = -1/(2 z^(3/2))
    upp = mp.mpf(3) / 4 / (rt * zv**2)   # 3/4 z^(-5/2)
    wp = zv * up
    wpp = zv * up + zv**2 * upp
    return u, up, upp, wp, wpp


def state_cpvi2(zv, sheet):
    """u=(sqrt(z)+1)/2 on the given sqrt sheet."""
    s = mp.mpf(1) if sheet > 0 else mp.mpf(-1)
    rt = s * mp.sqrt(zv)
    u = (rt + 1) / 2
    up = mp.mpf(1) / 4 / rt              # 1/2 * 1/(2 sqrt z) = 1/(4 sqrt z)
    upp = mp.mpf(-1) / 8 / (rt * zv)     # -1/8 z^(-3/2)
    wp = zv * up
    wpp = zv * up + zv**2 * upp
    return u, up, upp, wp, wpp


# Self-check the closed-form u',u'' against sympy-derived values at dps high.
for (uexpr, fn) in ((u1, state_cpvi1), (u2, state_cpvi2)):
    up_sym = sp.lambdify(z, sp.diff(uexpr, z), 'mpmath')
    upp_sym = sp.lambdify(z, sp.diff(uexpr, z, 2), 'mpmath')
    zt = mp.mpf('1.7')                   # +sheet test point
    _, up_n, upp_n, _, _ = fn(zt, +1)
    assert abs(up_n - up_sym(zt)) < mp.mpf('1e-45'), "u' self-check"
    assert abs(upp_n - upp_sym(zt)) < mp.mpf('1e-45'), "u'' self-check"
print("# closed-form u', u'' cross-checked vs sympy.diff (dps=50, <1e-45).")
print()

PROBES = [("Z4", mp.mpf('4.0')), ("Z2P3I", mp.mpc('2.0', '3.0'))]

for case, fn in (("CPVI1", state_cpvi1), ("CPVI2", state_cpvi2)):
    print(f"# --- {case} pins (u, zeta-frame wp, zeta-frame wpp) ---")
    for plab, zv in PROBES:
        for sheet, slab in ((+1, "PLUS"), (-1, "MINUS")):
            u, up, upp, wp, wpp = fn(zv, sheet)
            emit(f"{case}_{plab}_{slab}_U", mp.mpc(u))
            emit(f"{case}_{plab}_{slab}_WP", mp.mpc(wp))
            emit(f"{case}_{plab}_{slab}_WPP", mp.mpc(wpp))
    print()

# Sanity: brief-stated headline values for CPvi.1.
u_z4 = state_cpvi1(mp.mpf('4.0'), +1)[0]
assert abs(u_z4 - mp.mpf('0.5')) < mp.mpf('1e-45'), f"CPvi.1 u(4)!=1/2: {u_z4}"
u_c = state_cpvi1(mp.mpc('2.0', '3.0'), +1)[0]
# brief: 0.4643254526508149... - 0.2484994409113033...i
assert abs(u_c.real - mp.mpf('0.4643254526508149')) < mp.mpf('1e-15')
assert abs(u_c.imag - mp.mpf('-0.2484994409113033')) < mp.mpf('1e-15')
print("# Headline CPvi.1 spot values reproduced: u(4)=+1/2 ; "
      "u(2+3i)=0.46432545...-0.24849944...i (brief-verified).")
print("# All residual==0 GATES PASSED; values above are mpmath dps=50 pins.")
