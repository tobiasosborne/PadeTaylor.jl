using TaylorSeries
using Arblib

println("="^70)
println("PadeTaylor.jl Stage-0 Probe: TaylorSeries.jl + Arblib.jl")
println("Julia version: ", VERSION)
println("="^70)

results = Dict{String, String}()

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

"""Return true if two Arb values overlap (balls intersect)."""
function arb_overlaps(a::Arb, b::Arb)
    return Arblib.overlaps(a, b)
end

"""Return true if Arb a contains Arb b (ball a contains ball b)."""
function arb_contains(a::Arb, b::Arb)
    return Arblib.contains(a, b)
end

# ─────────────────────────────────────────────────────────────────────────────
# PROBE 1: Construction
# ─────────────────────────────────────────────────────────────────────────────
println("\n--- PROBE 1: Construction ---")
probe1_pass = false
try
    t1 = Taylor1(Arb, 60)
    println("  Taylor1(Arb, 60) => ", typeof(t1), " order=", get_order(t1))

    t2 = Taylor1([Arb(1), Arb(2), Arb(0.5)], 60)
    println("  Taylor1([Arb(1), Arb(2), Arb(0.5)], 60) => ", typeof(t2), " order=", get_order(t2))

    @assert t2[0] == Arb(1) "coeff[0] mismatch"
    @assert t2[1] == Arb(2) "coeff[1] mismatch"
    @assert t2[2] == Arb(0.5) "coeff[2] mismatch"

    println("  PASS: Both constructors work; coefficient indexing verified.")
    probe1_pass = true
    results["probe1"] = "PASS"
catch e
    println("  FAIL: ", typeof(e), ": ", e)
    results["probe1"] = "FAIL: $(typeof(e)): $e"
end

# ─────────────────────────────────────────────────────────────────────────────
# PROBE 2: Arithmetic
# ─────────────────────────────────────────────────────────────────────────────
println("\n--- PROBE 2: Arithmetic ---")
probe2_pass = false
try
    ORDER = 60
    # x is the Taylor variable around 0 at order ORDER
    x = Taylor1([Arb(0), Arb(1)], ORDER)

    # (1 + x)^2 two ways
    one_plus_x = 1 + x
    sq1 = one_plus_x^2
    sq2 = (1 + x) * (1 + x)

    println("  (1+x)^2 via ^2:   ", sq1)
    println("  (1+x)^2 via mult: ", sq2)

    # Assert coefficients are equal for all stored terms
    for k in 0:ORDER
        c1 = sq1[k]
        c2 = sq2[k]
        @assert arb_overlaps(c1, c2) "Coeff[$k] mismatch: $c1 vs $c2"
    end
    println("  All $(ORDER+1) coefficients of (1+x)^2 match between the two methods.")

    # Basic arithmetic
    a = Taylor1([Arb(1), Arb(2), Arb(3)], ORDER)
    b = Taylor1([Arb(4), Arb(5), Arb(6)], ORDER)
    _ = a + b
    _ = a - b
    _ = a * b
    # Division: divide a by something with non-zero constant term
    c = Taylor1([Arb(2), Arb(1)], ORDER)
    d_result = a / c
    println("  Addition, subtraction, multiplication, division: all executed.")

    probe2_pass = true
    results["probe2"] = "PASS"
catch e
    println("  FAIL: ", typeof(e), ": ", e)
    showerror(stdout, e, catch_backtrace())
    results["probe2"] = "FAIL: $(typeof(e)): $e"
end

# ─────────────────────────────────────────────────────────────────────────────
# PROBE 3: Transcendentals
# ─────────────────────────────────────────────────────────────────────────────
println("\n--- PROBE 3: Transcendentals (order 60) ---")
probe3_pass = false
try
    ORDER = 60
    x = Taylor1([Arb(0), Arb(1)], ORDER)

    # exp(x)
    println("  Computing exp(x)...")
    ex = exp(x)
    println("  exp(x) computed, type: ", typeof(ex))

    exp_failures = Int[]
    for k in 0:ORDER
        expected = Arb(1) / Arb(factorial(big(k)))
        got = ex[k]
        if !arb_overlaps(got, expected)
            push!(exp_failures, k)
            println("    MISMATCH at k=$k: got=$got expected=$expected")
        end
    end
    if isempty(exp_failures)
        println("  exp(x): all $(ORDER+1) coefficients correct (overlap with 1/k!).")
    else
        println("  exp(x): $(length(exp_failures)) coefficient mismatches at orders: $exp_failures")
    end

    # sin(x)
    println("  Computing sin(x)...")
    sx = sin(x)
    sin_failures = Int[]
    for k in 0:ORDER
        # sin coeff at k: (-1)^((k-1)/2) / k! for odd k, 0 for even k
        if iseven(k)
            expected = Arb(0)
        else
            m = (k - 1) ÷ 2
            sign = iseven(m) ? Arb(1) : Arb(-1)
            expected = sign / Arb(factorial(big(k)))
        end
        got = sx[k]
        if !arb_overlaps(got, expected)
            push!(sin_failures, k)
            println("    MISMATCH at k=$k: got=$got expected=$expected")
        end
    end
    if isempty(sin_failures)
        println("  sin(x): all $(ORDER+1) coefficients correct.")
    else
        println("  sin(x): $(length(sin_failures)) mismatches at: $sin_failures")
    end

    # cos(x)
    println("  Computing cos(x)...")
    cx = cos(x)
    cos_failures = Int[]
    for k in 0:ORDER
        if isodd(k)
            expected = Arb(0)
        else
            m = k ÷ 2
            sign = iseven(m) ? Arb(1) : Arb(-1)
            expected = sign / Arb(factorial(big(k)))
        end
        got = cx[k]
        if !arb_overlaps(got, expected)
            push!(cos_failures, k)
            println("    MISMATCH at k=$k: got=$got expected=$expected")
        end
    end
    if isempty(cos_failures)
        println("  cos(x): all $(ORDER+1) coefficients correct.")
    else
        println("  cos(x): $(length(cos_failures)) mismatches at: $cos_failures")
    end

    # log(1+x) — coefficients are (-1)^(k+1)/k for k>=1, 0 for k=0
    println("  Computing log(1+x)...")
    lx = log(1 + x)
    log_failures = Int[]
    for k in 0:ORDER
        if k == 0
            expected = Arb(0)
        else
            sign = iseven(k+1) ? Arb(1) : Arb(-1)   # (-1)^(k+1): +1 for even k+1 (odd k), -1 for odd k+1 (even k)
            expected = sign / Arb(k)
        end
        got = lx[k]
        if !arb_overlaps(got, expected)
            push!(log_failures, k)
            println("    MISMATCH at k=$k: got=$got expected=$expected")
        end
    end
    if isempty(log_failures)
        println("  log(1+x): all $(ORDER+1) coefficients correct.")
    else
        println("  log(1+x): $(length(log_failures)) mismatches at: $log_failures")
    end

    if isempty(exp_failures) && isempty(sin_failures) && isempty(cos_failures) && isempty(log_failures)
        probe3_pass = true
        results["probe3"] = "PASS"
    else
        results["probe3"] = "FAIL: coefficient mismatches in exp/sin/cos/log"
    end
catch e
    println("  FAIL: ", typeof(e), ": ", e)
    showerror(stdout, e, catch_backtrace())
    results["probe3"] = "FAIL: $(typeof(e)): $e"
end

# ─────────────────────────────────────────────────────────────────────────────
# PROBE 4: High-order stress (order 80)
# ─────────────────────────────────────────────────────────────────────────────
println("\n--- PROBE 4: High-order stress (order 80) ---")
probe4_pass = false
try
    ORDER = 80
    x = Taylor1([Arb(0), Arb(1)], ORDER)
    two = Taylor1([Arb(2)], ORDER)  # constant 2

    println("  Computing exp(x)^2 at order $ORDER...")
    ex2 = exp(x)^2

    println("  Computing exp(2x) at order $ORDER...")
    # 2x as a Taylor1
    twox = Taylor1([Arb(0), Arb(2)], ORDER)
    e2x = exp(twox)

    # Coefficients of exp(2x): 2^k / k!
    failures = Int[]
    for k in 0:ORDER
        c1 = ex2[k]
        c2 = e2x[k]
        expected = Arb(2)^k / Arb(factorial(big(k)))
        if !arb_overlaps(c1, c2)
            push!(failures, k)
            println("    exp(x)^2 vs exp(2x) MISMATCH at k=$k: exp(x)^2=$(c1) exp(2x)=$(c2)")
        end
        if !arb_overlaps(c1, expected)
            println("    exp(x)^2 vs 2^k/k! MISMATCH at k=$k: got=$(c1) expected=$(expected)")
        end
    end

    if isempty(failures)
        println("  exp(x)^2 and exp(2x) agree on all $(ORDER+1) coefficients at order $ORDER.")
        probe4_pass = true
        results["probe4"] = "PASS"
    else
        println("  $(length(failures)) mismatches between exp(x)^2 and exp(2x).")
        results["probe4"] = "FAIL: $(length(failures)) coefficient mismatches"
    end
catch e
    println("  FAIL: ", typeof(e), ": ", e)
    showerror(stdout, e, catch_backtrace())
    results["probe4"] = "FAIL: $(typeof(e)): $e"
end

# ─────────────────────────────────────────────────────────────────────────────
# PROBE 5: Arb precision propagation
# ─────────────────────────────────────────────────────────────────────────────
println("\n--- PROBE 5: Arb precision propagation (256-bit) ---")
probe5_pass = false
try
    setprecision(Arb, 256)
    ORDER = 60

    x = Taylor1([Arb(0), Arb(1)], ORDER)
    ex = exp(x)

    println("  Arb precision set to 256 bits.")
    println("  exp(x) computed at 256-bit precision.")

    # Check radii of coefficients — should be << 2^-200 (a slack of ~56 bits)
    threshold = Float64(2.0^(-200))
    max_radius = 0.0
    for k in 0:ORDER
        c = ex[k]
        r = Float64(Arblib.radius(c))
        max_radius = max(max_radius, r)
        if r > threshold
            println("  WARNING: coeff[$k] radius=$r exceeds threshold=$threshold")
        end
    end
    println("  Max radius across all exp(x) coefficients: ", max_radius)
    println("  Threshold (2^-200): ", threshold)

    # Also check a specific known coefficient: exp coeff[60] = 1/60!
    c60 = ex[60]
    expected_c60 = Arb(1) / Arb(factorial(big(60)))
    println("  exp coeff[60] = ", c60)
    println("  1/60! (Arb)   = ", expected_c60)
    println("  Overlap: ", Arblib.overlaps(c60, expected_c60))

    if max_radius < threshold
        println("  PASS: All radii < 2^-200; precision is genuinely 256-bit.")
        probe5_pass = true
        results["probe5"] = "PASS: max_radius=$(max_radius) < 2^-200=$(threshold)"
    else
        println("  DEGRADED: Some radii exceed 2^-200 threshold.")
        results["probe5"] = "DEGRADED: max_radius=$(max_radius) >= 2^-200=$(threshold)"
    end
catch e
    println("  FAIL: ", typeof(e), ": ", e)
    showerror(stdout, e, catch_backtrace())
    results["probe5"] = "FAIL: $(typeof(e)): $e"
end

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
println("\n", "="^70)
println("SUMMARY")
println("="^70)
for k in ["probe1", "probe2", "probe3", "probe4", "probe5"]
    status = get(results, k, "NOT RUN")
    println("  $k: $status")
end

all_pass = all(startswith(get(results, k, ""), "PASS") for k in ["probe1","probe2","probe3","probe4","probe5"])
println("\nOverall: ", all_pass ? "ALL PASS" : "SOME FAILURES")
println("="^70)
