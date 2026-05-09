"""capture.py -- Phase 5 PadeStepper oracle (Python cross-check).

Independent integration of u'' = 6u^2 (equianharmonic ℘ companion) and
u'' = 6u^2 + z (PI) from FW 2011 ICs via mpmath.odefun (Taylor-method
arbitrary-precision ODE solver).

Usage:  python3 capture.py
Outputs: oracles_python.txt (Julia-readable literals).

mpmath has no native Weierstrass-P, so we do not produce a closed-form
℘ cross-check; instead we cross-check the NDSolve numbers from
wolframscript by integrating the same ODEs in Python at the same
working precision.  Three-source triangulation for the test cases:

  ℘ at z=0.5:   wolframscript closed-form ≡ wolframscript NDSolve ≡ mpmath odefun
  ℘ at z=0.9:   same three sources
  ℘ at z=0.95:  same three sources
  PI at z=0.5:  wolframscript NDSolve ≡ mpmath odefun
"""

import mpmath as mp

mp.mp.dps = 40  # 40 decimal digits >> Float64 precision

u0  = mp.mpf('1.071822516416917')
up0 = mp.mpf('1.710337353176786')

# u'' = 6 u^2  → as a first-order system [u, up]:
#   d/dt [u, up] = [up, 6 u^2]
def F_weier(t, y):
    u_, up_ = y[0], y[1]
    return [up_, 6 * u_**2]

# u'' = 6 u^2 + z  → similar:
def F_PI(t, y):
    u_, up_ = y[0], y[1]
    return [up_, 6 * u_**2 + t]

print("Integrating u'' = 6u^2 with mpmath.odefun ...")
solW = mp.odefun(F_weier, mp.mpf('0'), [u0, up0], tol=mp.mpf('1e-30'))

print("Integrating u'' = 6u^2 + z (PI) with mpmath.odefun ...")
solP = mp.odefun(F_PI, mp.mpf('0'), [u0, up0], tol=mp.mpf('1e-30'))

with open('oracles_python.txt', 'w') as f:
    f.write("# Phase 5 PadeStepper oracle (Python mpmath cross-check).\n")
    f.write("# Cross-validates capture.wl (wolframscript NDSolve) at 40 dps.\n\n")

    # ---- 5.1.1: u(0.5), u'(0.5) on u'' = 6 u^2 ----
    u_at, up_at = solW(mp.mpf('0.5'))
    f.write("# Case 5.1.1 -- mpmath.odefun on u'' = 6u^2 to z=0.5.\n")
    f.write(f"u_5_1_1_at_05_py  = {float(u_at)!r}\n")
    f.write(f"up_5_1_1_at_05_py = {float(up_at)!r}\n\n")

    # ---- 5.1.2: u(0.9), u'(0.9) on u'' = 6 u^2 ----
    u_at, up_at = solW(mp.mpf('0.9'))
    f.write("# Case 5.1.2 -- mpmath.odefun on u'' = 6u^2 to z=0.9.\n")
    f.write(f"u_5_1_2_at_09_py  = {float(u_at)!r}\n")
    f.write(f"up_5_1_2_at_09_py = {float(up_at)!r}\n\n")

    # ---- 5.1.3: PI to z=0.5 ----
    u_at, up_at = solP(mp.mpf('0.5'))
    f.write("# Case 5.1.3 -- mpmath.odefun on u'' = 6u^2 + z (PI) to z=0.5.\n")
    f.write(f"u_5_1_3_PI_at_05_py  = {float(u_at)!r}\n")
    f.write(f"up_5_1_3_PI_at_05_py = {float(up_at)!r}\n\n")

    # ---- 5.1.4: u(0.95), u'(0.95) on u'' = 6 u^2 (near pole at z=1) ----
    u_at, up_at = solW(mp.mpf('0.95'))
    f.write("# Case 5.1.4 -- mpmath.odefun on u'' = 6u^2 to z=0.95 (near pole at z=1).\n")
    f.write(f"u_5_1_4_at_095_py  = {float(u_at)!r}\n")
    f.write(f"up_5_1_4_at_095_py = {float(up_at)!r}\n")

print("Wrote oracles_python.txt")
