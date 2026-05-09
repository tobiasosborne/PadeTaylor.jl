"""capture.py -- Phase 6 Problems oracle (Python cross-check).

Independent integration of u'' = 6u^2 (equianharmonic ℘ companion)
from FW 2011 ICs via mpmath.odefun (Taylor-method arbitrary-precision
ODE solver), used as the THIRD source in the pole-bridge demo
verification (worklog 004).

Usage:  python3 capture.py
Outputs: oracles_python.txt (Julia-readable literals).

Scope notes (revised for the Phase-6 v1 pivot):

  - We integrate from FW ICs at z=0 ONCE at 40 dps and query at
    z = 0.5 and z = 0.95.  Do NOT integrate past z=0.95 — the pole at
    z=1 cannot be reliably crossed by mpmath.odefun (no Padé
    extension), and the previous capture.py's range to z=30 across 12
    poles was the source of long hangs / non-convergence.
  - For z > 1 (z = 1.05, z = 1.4) the closed-form WeierstrassP from
    capture.wl is the sole primary source; mpmath cannot help past
    the pole.  Padé bridges it; that's the whole point of the test.
  - The Phase-5 padestepper-oracle established 3-source consensus
    (closed-form ≡ NDSolve ≡ mpmath.odefun) on this same ODE at
    z=0.5 to <1e-15.  Replicating it here gives independent
    confirmation that nothing has drifted between probes.
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

print("Integrating u'' = 6u^2 from z=0 to z=0.95 with mpmath.odefun ...")
sol = mp.odefun(F_weier, mp.mpf('0'), [u0, up0], tol=mp.mpf('1e-30'))

with open('oracles_python.txt', 'w') as f:
    f.write("# Phase 6 Problems oracle (Python mpmath cross-check).\n")
    f.write("# Range integrated: z in [0, 0.95] at 40 dps.  Past z=1 (the\n")
    f.write("# pole), mpmath.odefun cannot continue; closed-form is sole source.\n\n")

    # ---- z = 0.5 ----
    print("  evaluating at z=0.5 ...")
    u_at, up_at = sol(mp.mpf('0.5'))
    f.write("# u(0.5) and u'(0.5) via mpmath.odefun (BEFORE pole at z=1).\n")
    f.write(f"u_at_0_5_py  = {float(u_at)!r}\n")
    f.write(f"up_at_0_5_py = {float(up_at)!r}\n\n")

    # ---- z = 0.95 ----
    print("  evaluating at z=0.95 ...")
    u_at, up_at = sol(mp.mpf('0.95'))
    f.write("# u(0.95) and u'(0.95) via mpmath.odefun (NEAR pole, BEFORE z=1).\n")
    f.write(f"u_at_0_95_py  = {float(u_at)!r}\n")
    f.write(f"up_at_0_95_py = {float(up_at)!r}\n")

print("Wrote oracles_python.txt")
