"""capture.py -- Phase 6 Problems oracle (Python cross-check).

Cross-validates the wolframscript closed-form WeierstrassP values
against an independent ODE integration (mpmath.odefun, Taylor-method
arbitrary precision) for the Phase-6 test cases that lie within
practical integration range.

Usage:  python3 capture.py
Outputs: oracles_python.txt

Scope notes:
  - We integrate u'' = 6 u^2 from FW 2011 ICs over z ∈ [0, 30] once
    at 40 decimal digits, then query at z = 7.123, 28.261, 30.0.
    This covers tests 6.1.1, 6.1.3, 6.1.5.
  - z = 10^4 (test 6.1.2) is not feasible with mpmath.odefun at any
    useful precision; we rely on the FW Table 5.1 paper-pinned
    reference + wolframscript closed-form for that case.
  - BigFloat-256 (test 6.1.4) at 80 dps via mpmath would take many
    minutes; the wolframscript 78-digit pin is sole source (same
    closed-form, two languages: an honest one-source pin for that
    very-high-precision corner).
"""

import mpmath as mp

mp.mp.dps = 40

u0  = mp.mpf('1.071822516416917')
up0 = mp.mpf('1.710337353176786')

def F_weier(t, y):
    u_, up_ = y[0], y[1]
    return [up_, 6 * u_**2]

print("Integrating u'' = 6u^2 from z=0 to z=30 with mpmath.odefun ...")
sol = mp.odefun(F_weier, mp.mpf('0'), [u0, up0], tol=mp.mpf('1e-25'))

with open('oracles_python.txt', 'w') as f:
    f.write("# Phase 6 Problems oracle (Python mpmath cross-check).\n")
    f.write("# Range integrated: z in [0, 30] at 40 dps; z=10^4 and BF-256 not covered.\n\n")

    for z_query, label in [
        (mp.mpf('7.123'),  '7_123'),
        (mp.mpf('28.261'), '28_261'),
        (mp.mpf('30'),     '30'),
    ]:
        print(f"  evaluating at z={z_query} ...")
        u_at, up_at = sol(z_query)
        f.write(f"# u({z_query}) and u'({z_query}) via mpmath.odefun.\n")
        f.write(f"u_at_{label}_py  = {float(u_at)!r}\n")
        f.write(f"up_at_{label}_py = {float(up_at)!r}\n\n")

print("Wrote oracles_python.txt")
