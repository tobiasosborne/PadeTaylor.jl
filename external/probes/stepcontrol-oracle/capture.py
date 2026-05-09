"""capture.py -- Phase 4 StepControl oracle, Python cross-check.

Independent computation of (eps/|x[k]|)^(1/k) at high precision via
mpmath, plus an independent root + projection computation via sympy
for the Pade-pole tests.

Usage:  python3 capture.py
Outputs: oracles_python.txt

Cross-validates with capture.jl (TI.jl primary) and capture.wl
(wolframscript) via verify.jl.
"""

import mpmath
import sympy as sp

OUT = open('oracles_python.txt', 'w')

OUT.write("# Phase 4 StepControl oracle (Python mpmath/sympy cross-check).\n")
OUT.write("# Cross-validates capture.jl (TI.jl) and capture.wl (wolframscript).\n\n")

# ---- 4.1.1: step_jorba_zou expected at high precision -----------------
# Formula consensus (paper §3.3.1 + TI.jl stepsize):
#   h = min over k in {p-1, p} of (eps / |x[k]|)^(1/k)
# For x[k] = 1/k!: |x[k]| = 1/k!, so (eps/|x[k]|)^(1/k) = (eps · k!)^(1/k).
mpmath.mp.dps = 50
eps = mpmath.mpf('1e-12')
p = 30

def candidate(k):
    """(eps * k!)^(1/k) at 50 dps."""
    fk = mpmath.factorial(k)
    return (eps * fk) ** (mpmath.mpf(1) / k)

h_p_minus_1 = candidate(p - 1)
h_p         = candidate(p)
h_min       = min(h_p_minus_1, h_p)

OUT.write("# Case 4.1.1 -- (eps · k!)^(1/k) for k in {29, 30}, eps = 1e-12.\n")
OUT.write("# Source: mpmath at 50 decimal digits.\n")
OUT.write(f"h_4_1_1_at_29 = {float(h_p_minus_1)!r}\n")
OUT.write(f"h_4_1_1_at_30 = {float(h_p)!r}\n")
OUT.write(f"h_4_1_1_min   = {float(h_min)!r}\n")
OUT.write(f'h_4_1_1_min_50dps = "{mpmath.nstr(h_min, 48)}"\n')
OUT.write("\n")

# ---- 4.1.3: roots of 1 - z/2 ------------------------------------------
OUT.write("# Case 4.1.3 -- roots of Q(z) = 1 - z/2.  Pole at z = 2.\n")
z = sp.Symbol('z')
roots_4_1_3 = sp.solve(1 - z/2, z)
OUT.write(f"roots_4_1_3 = ComplexF64[{', '.join(f'{complex(r)!r}' for r in roots_4_1_3)}]\n")
OUT.write("step_4_1_3_expected = 2.0\n\n")

# ---- 4.1.4: roots of 1 - 0.6 z + 0.1 z^2 (after b[1]=1 normalization) -
OUT.write("# Case 4.1.4 -- roots of Q(z) = 1 - 0.6 z + 0.1 z^2.  Poles at z = 3 ± i.\n")
roots_4_1_4 = sp.solve(sp.Rational(1,1) - sp.Rational(6,10)*z + sp.Rational(1,10)*z**2, z)
OUT.write(f"roots_4_1_4 = ComplexF64[{', '.join(f'{complex(r)!r}' for r in roots_4_1_4)}]\n")
OUT.write("step_4_1_4_expected = 3.0\n")

OUT.close()
print("Wrote oracles_python.txt")
