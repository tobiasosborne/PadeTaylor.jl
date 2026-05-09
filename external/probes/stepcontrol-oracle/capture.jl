# capture.jl -- Phase 4 StepControl oracle, Julia/TI.jl primary.
#
# Pin step-size oracle values for tests 4.1.1, 4.1.2 by running
# TaylorIntegration.jl's `stepsize` on the test inputs and writing
# Julia-readable literals to oracles_julia.txt.
#
# Usage (from project root):
#   julia --project=external/probes/stepcontrol-oracle external/probes/stepcontrol-oracle/capture.jl
#
# Cross-validated by capture.py (Python mpmath, independent computation
# of (eps/|x[k]|)^(1/k)) and by capture.wl (wolframscript, same).

using TaylorIntegration: stepsize
using TaylorSeries: Taylor1
using Polynomials: Polynomial, roots

const OUT = joinpath(@__DIR__, "oracles_julia.txt")

open(OUT, "w") do io
    println(io, "# Phase 4 StepControl oracle (TI.jl primary).")
    println(io, "# Captured by: julia --project=. capture.jl")
    println(io, "# Format: Julia-readable literals (top-level assignments).")
    println(io)

    # ---- Case 4.1.1 / 4.1.2: exp Taylor coefs, order 30, eps = 1e-12 ----
    #
    # Per the consensus of (a) Jorba-Zou 2005 §3.3.1 eq. 11 and (b)
    # TI.jl src/integrator/stepsize.jl:17-35, the step-size formula at
    # fixed order p is
    #
    #     h = min over k in {p-1, p} of (eps / |x[k]|)^(1/k).
    #
    # No e^2 divisor and no 0.7-safety factor: those appear in DESIGN.md
    # but have no source in either the paper or TI.jl (grep "0.7" in
    # both: zero matches).  The paper uses /e^2 only for the
    # asymptotically-optimal-(p,h) case where p is also a free variable;
    # PadeTaylor fixes p (= 30 per FW 2011), so the eps-driven form is
    # the appropriate reduction of the paper formula.
    #
    # We capture TI.jl's output as the load-bearing oracle for 4.1.1 and
    # 4.1.2; capture.py reproduces this independently in mpmath.
    println(io, "# Case 4.1.1 / 4.1.2 -- exp Taylor coefs c[k] = 1/k!, k=0..30; eps = 1e-12.")
    println(io, "# Source: TaylorIntegration.stepsize(Taylor1(c, 30), 1e-12).")

    coefs_4_1_1 = [1.0 / factorial(big(k)) for k in 0:30]
    # Convert to Float64 explicitly to match the test environment.
    coefs_4_1_1_f64 = Float64.(coefs_4_1_1)
    taylor_4_1_1 = Taylor1(coefs_4_1_1_f64, 30)
    eps_4_1_1 = 1e-12

    h_TI = stepsize(taylor_4_1_1, eps_4_1_1)
    println(io, "# Direct TI.jl output, 1e-15 round-trip precision:")
    println(io, "h_4_1_1_TI = ", repr(h_TI))
    # Also pin the exact coefs (so we can reconstruct in test without
    # depending on factorial-precision in Julia 1.10 vs 1.12).
    print(io, "coefs_4_1_1 = [")
    for (i, c) in enumerate(coefs_4_1_1_f64)
        i > 1 && print(io, ", ")
        print(io, repr(c))
    end
    println(io, "]")
    println(io, "eps_4_1_1 = ", repr(eps_4_1_1))
    println(io)

    # ---- Case 4.1.3: P = 1/(1 - z/2), z_current=0, target=5 ----
    #
    # Denominator polynomial Q(z) = 1 - z/2.  Polynomials.jl convention
    # is lowest-order coefficient first, so b = [1.0, -0.5].  Roots:
    # solve 1 - z/2 = 0 -> z = 2.
    # Projection of (z_pole - z_current) onto direction (target -
    # z_current)/|target - z_current|: Re((2 - 0)·conj(1))/|1| = 2.
    # Capped at |target - z_current| = 5; min(2, 5) = 2.
    println(io, "# Case 4.1.3 -- P = 1/(1 - z/2); pole at z = 2.")
    println(io, "# Source: Polynomials.roots([1.0, -0.5]).")
    Q_4_1_3 = Polynomial([1.0, -0.5])
    roots_4_1_3 = roots(Q_4_1_3)
    print(io, "roots_4_1_3 = ComplexF64[")
    for (i, r) in enumerate(roots_4_1_3)
        i > 1 && print(io, ", ")
        print(io, repr(ComplexF64(r)))
    end
    println(io, "]")
    println(io, "step_4_1_3_expected = 2.0")
    println(io)

    # ---- Case 4.1.4: P = 1/(1 + (z - 3)^2) = 1/(z^2 - 6z + 10) ----
    #
    # Normalize so b[1] = 1 (GGT 2013 step 7): divide numerator and
    # denominator by 10.  b = [1.0, -0.6, 0.1].  Roots of
    # 0.1 z^2 - 0.6 z + 1: z = (0.6 ± sqrt(0.36 - 0.4))/0.2 = 3 ± i.
    # Path direction = +1 (real); projection of 3+i onto +1 is Re(3+i)
    # = 3.  Same for 3-i.  min over forward poles = 3 < 5 = cap.
    println(io, "# Case 4.1.4 -- P = 1/(1 + (z - 3)^2); poles at z = 3 ± i.")
    println(io, "# Source: Polynomials.roots([1.0, -0.6, 0.1]) (b normalized so b[1] = 1).")
    Q_4_1_4 = Polynomial([1.0, -0.6, 0.1])
    roots_4_1_4 = roots(Q_4_1_4)
    print(io, "roots_4_1_4 = ComplexF64[")
    for (i, r) in enumerate(roots_4_1_4)
        i > 1 && print(io, ", ")
        print(io, repr(ComplexF64(r)))
    end
    println(io, "]")
    println(io, "step_4_1_4_expected = 3.0")
end

println("Wrote ", OUT)
