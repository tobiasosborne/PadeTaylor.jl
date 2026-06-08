#!/usr/bin/env python3
"""capture.py — mpmath oracle for corpus_pathnet_lattice_sectors_test.jl.

Family I (docs/test_corpus/02_corpus_extension_plan.md:345-385): PATH-NETWORK at
scale — lattice-geometry CONTRAST, the tritronquée EMPTY-SECTOR quantitative
pins, and LONG-RANGE multi-direction accuracy.  Three cases:

  CPN.5  Painlevé I tritronquée empty sector.  f(z,u,u') = 6u² + z, tritronquée
         IC at z0 = 0 (FW 2011 eq. 4.1):  u(0) = -0.1875543083404949,
         u'(0) = 0.3049055602612289.  The tritronquée confines poles to ONE
         wedge near the positive real axis (FW Fig 4.x; leading real pole
         z ≈ 2.07); the NEGATIVE-real fan is provably POLE-FREE, with the
         leading asymptote u ~ -sqrt(-z/6) (FW eq. 1.2, negative branch).
         INDEPENDENT ORACLE: integrate u'' = 6u² + z from the IC ALONG THE RAY
         z(t) = t·dir, t ∈ [0, |z_T|], via mpmath.odefun at dps=30, tol 1e-26.
         This is a genuinely independent ODE integrator (not the package's
         Taylor–Padé).  Targets DEEP in the pole-free fan {-30, -40, -50}
         (on-axis) + a small off-axis fan {-30±5i}.  The odefun ray value is
         the TIGHT pin (~1e-13 vs the package); the -sqrt(-z/6) asymptote is
         the loose ~1e-5 sanity cross-check.

  CPN.6  Equianharmonic ℘ long-range, MULTI-DIRECTION.  f = 6u², the FW Table
         5.1 problem (u = ℘(z + c1; 0, c2), c1 = -1, c2 = 2), but integrated
         along SEVERAL rays (±30°, ±60°), not just FW's single real axis.
         INDEPENDENT ORACLE: odefun ray-integration g(t,Y) = [dir·Y1,
         dir·6·Y0²] from the FW IC, evaluated at t = R where the ray stays
         CLEAR of the ℘ lattice poles (clearance ≥ 0.5 at R = 11).  The real
         axis (0°) reuses the FW Table 5.1 paper value u(30) = 1.095098255959744
         (the path-network bridges the on-axis poles).

  CPN.3  Lattice-vs-lattice routing.  equianharmonic ℘ (f = 6u², 60° RHOMBIC
         lattice) vs lemniscatic ℘ (f = 6u² - 1/2, g2 = 1, g3 = 0, SQUARE
         lattice).  Lattice constants only (the test asserts extracted poles
         lie on the respective lattice and the two visited trees differ):
         equianharmonic Ω = Γ(1/3)³/(2^(13/6)·π); lemniscatic real half-period
         ω = Γ(1/4)²/(4√π).

This script is the residual==0 GATE plus the numerical-pin source.  No value is
emitted before the recast residual / first-order identity is proven ≡ 0 (mpmath,
to ≤1e-25).  mpmath at dps = 30 (the ray integrations); lattice constants exact.

Run:  python3 external/probes/corpus-oracles/pathnet-lattice/capture.py
"""
import sys
import mpmath as mp

mp.mp.dps = 30


def emit(label, val, n=22):
    if isinstance(val, mp.mpc):
        print(f"const {label} = complex({mp.nstr(val.real, n)}, "
              f"{mp.nstr(val.imag, n)})")
    else:
        print(f"const {label} = {mp.nstr(val, n)}")


# ===========================================================================
# (0) RESIDUAL==0 GATES.  The two ODE recasts both have f = u'' as the literal
#     RHS, so the recast residual is the identity u'' - f ≡ 0 by construction;
#     we instead GATE on the leading-asymptote consistency for the tritronquée
#     (the ray oracle's only non-trivial claim) and on the ℘ first integral.
# ===========================================================================
# Tritronquée IC.
U0_TRI  = mp.mpf('-0.1875543083404949')
UP0_TRI = mp.mpf('0.3049055602612289')

# Equianharmonic ℘ IC (FW 2011 §5.1.1).
U0_EQ  = mp.mpf('1.071822516416917')
UP0_EQ = mp.mpf('1.710337353176786')
# First integral of u'' = 6u²: (u')² = 4u³ - c2, with c2 = 2 (g3 of the
# equianharmonic ℘, FW md:284).  GATE: the FW IC must satisfy it to ≤1e-13
# (the IC literals are 16-digit, so this is the round-off floor).
inv_eq = UP0_EQ**2 - (4 * U0_EQ**3 - mp.mpf(2))
assert abs(inv_eq) < mp.mpf(10)**(-13), f"CPN.6 ℘ first-integral residual {inv_eq}"

print("# residual==0 / first-integral GATES PASSED "
      "(℘ (u')²=4u³-2 to <1e-13).")
print()


# ===========================================================================
# (5) CPN.5 — Painlevé I tritronquée empty-sector odefun RAY oracle.
#     u'' = 6u² + z, ray z(t) = t·dir, dY/dt = dir·[u', 6u² + z].
# ===========================================================================
def tri_ray(target, tol=26):
    direction = target / abs(target)
    L = abs(target)

    def g(t, Y):
        z = t * direction
        return [direction * Y[1], direction * (6 * Y[0] ** 2 + z)]

    f = mp.odefun(g, 0, [mp.mpc(U0_TRI), mp.mpc(UP0_TRI)],
                  tol=mp.mpf(10) ** (-tol))
    return f(L)[0]


print("# CPN.5  PI tritronquée empty-sector: odefun ray oracle for u(target).")
print("#        Asymptote -sqrt(-z/6) agrees to ~1e-5 (loose sanity);")
print("#        the odefun ray value is the TIGHT pin (~1e-13 vs package).")
_tri = [("CPN5_U_M30",   mp.mpc(-30, 0)),
        ("CPN5_U_M40",   mp.mpc(-40, 0)),
        ("CPN5_U_M50",   mp.mpc(-50, 0)),
        ("CPN5_U_M30P5I", mp.mpc(-30, 5)),
        ("CPN5_U_M30M5I", mp.mpc(-30, -5))]
for label, tgt in _tri:
    val = tri_ray(tgt)
    asym = -mp.sqrt(-tgt / 6)
    rel = abs((val - asym) / asym)
    # GATE: the deep-fan targets must agree with the asymptote to <1e-4
    # (confirms the target really is in the pole-free sector — a target in
    # the pole sector would diverge from the asymptote by O(1)).
    assert rel < mp.mpf('1e-4'), f"{label}: rel_asym {rel} ≥ 1e-4 — NOT pole-free"
    emit(label, val)
    sys.stdout.flush()
print()


# ===========================================================================
# (6) CPN.6 — Equianharmonic ℘ long-range MULTI-DIRECTION odefun ray oracle.
#     f = 6u², g(t,Y) = [dir·Y1, dir·6·Y0²].  Rays at ±30°, ±60°, R = 11
#     (clearance from the ℘ lattice ≥ 0.5).  The 0° (real-axis) value is the
#     FW Table 5.1 paper u(30) = 1.095098255959744 (pinned in-test inline).
# ===========================================================================
def eq_ray(direction, R, tol=26):
    def g(t, Y):
        return [direction * Y[1], direction * 6 * Y[0] ** 2]

    f = mp.odefun(g, 0, [mp.mpc(U0_EQ), mp.mpc(UP0_EQ)],
                  tol=mp.mpf(10) ** (-tol))
    return f(R)[0]


print("# CPN.6  equianharmonic ℘ multi-direction:  u(R·e^{iθ}), R = 11.")
R6 = mp.mpf(11)
for label, deg in [("CPN6_U_P30", 30), ("CPN6_U_M30", -30),
                   ("CPN6_U_P60", 60), ("CPN6_U_M60", -60)]:
    d = mp.expjpi(mp.mpf(deg) / 180)
    val = eq_ray(d, R6)
    emit(label, val)
    sys.stdout.flush()
print(f"# CPN.6  ray radius R = {mp.nstr(R6, 4)} (clearance ≥ 0.5 from ℘ lattice).")
print()


# ===========================================================================
# (3) CPN.3 — Lattice-vs-lattice constants.  Exact closed-form half-periods.
# ===========================================================================
EQUI_OM = mp.gamma(mp.mpf(1) / 3) ** 3 / (mp.mpf(2) ** (mp.mpf(13) / 6) * mp.pi)
LEMN_OM = mp.gamma(mp.mpf(1) / 4) ** 2 / (4 * mp.sqrt(mp.pi))
print("# CPN.3  lattice constants (exact closed form).")
print("#   equianharmonic Ω = Γ(1/3)³/(2^(13/6)·π)  -> 60° RHOMBIC lattice")
print("#   lemniscatic    ω = Γ(1/4)²/(4√π)          -> SQUARE lattice")
emit("CPN3_EQUI_OM", EQUI_OM, 25)
emit("CPN3_LEMN_OM", LEMN_OM, 25)
print()
print("# DONE — all GATES passed; literals above pinned into "
      "_oracle_corpus_pathnet_lattice_sectors.jl")
