# battery.jl — the audition battery: entire + meromorphic + high-d systems,
# each with a jet builder (via the package's own vector_taylor_coefficients,
# rescaled by h^k exactly as VectorStepper.jl:237) and a ground-truth oracle.
#
# Constants/oracles cited from the test suite (see the scout digest):
#   harmonic  test/vector_stepper_test.jl:92,108  (VS.1.2 h=0.7, ~5e-9 @+2)
#   exp pair  test/vector_stepper_test.jl:149-162  (~4e-8 @+2)
#   ℘        test/_oracle_problems.jl:32-33,83-84  (FW Table 5.1 pinned)
#   CM       test/calogero_moser_test.jl:150-163
#   NY A₄    figures/_noumi_yamada_a4_helpers.jl:91-109  (conservation Σf=t)

using PadeTaylor: vector_taylor_coefficients, noumi_yamada_rhs

struct Sys
    name   :: String
    kind   :: Symbol          # :entire | :meromorphic
    d      :: Int
    m      :: Int             # denominator degree = production order ÷ 2
    h      :: Float64
    rhs                        # f(z, y) -> ẏ
    y0     :: Vector{ComplexF64}
    z0     :: ComplexF64
    oracle                     # () -> Vector{ComplexF64}, or :conservation
end

# Rescaled jets c̃_k = h^k c_k, built to order 2m+pad so the held-out
# coefficient at order K+1 exists for diagnostics.  Cell degree stays m.
function rescaled_jets(s::Sys; pad = 4)
    order = 2 * s.m + pad
    raw = vector_taylor_coefficients(s.rhs, s.z0, s.y0, order)
    return [[s.h^(k - 1) * raw[i][k] for k = 1:length(raw[i])] for i = 1:s.d]
end

eval_at_one(num, Q) = sum(num) / sum(Q)        # P(1)/Q(1)

function step_error(s::Sys, cell)
    yv = [eval_at_one(cell.numerators[i], cell.denominator) for i = 1:s.d]
    if s.oracle === :conservation
        return abs(sum(yv) - (s.z0 + s.h))     # Σf(t) = t drift
    else
        tru = s.oracle()
        return maximum(abs(yv[i] - tru[i]) for i = 1:s.d)
    end
end

function build_battery()
    sys = Sys[]

    # 1. harmonic [cos, -sin] — entire, d=2, the VS.1.2 setup
    push!(sys, Sys("harmonic", :entire, 2, 15, 0.7,
        (z, y) -> [y[2], -y[1]], ComplexF64[1, 0], 0.0 + 0im,
        () -> ComplexF64[cos(0.7), -sin(0.7)]))

    # 2. exp pair [exp, z·exp] — entire, d=2
    push!(sys, Sys("exp-pair", :entire, 2, 15, 0.2,
        (z, y) -> [y[1], y[1] + y[2]], ComplexF64[1, 0], 0.0 + 0im,
        () -> ComplexF64[exp(0.2), 0.2 * exp(0.2)]))

    # 3. Weierstrass-℘ companion [u, u'] — meromorphic (regular step h=0.5,
    #    pole at z=1 not crossed), d=2; FW Table 5.1 pinned oracle at z=0.5
    push!(sys, Sys("weierstrass-℘", :meromorphic, 2, 15, 0.5,
        (z, y) -> [y[2], 6 * y[1]^2],
        ComplexF64[1.071822516416917, 1.710337353176786], 0.0 + 0im,
        () -> ComplexF64[4.0044646690030875, 15.964278048239492]))

    # 4. Calogero–Moser N=2 — meromorphic (off-axis poles), d=4; closed form
    cm = (z, y) -> begin
        dd = y[1] - y[2]; a = 4 / dd^3
        [y[3], y[4], a, -a]
    end
    cmoracle = () -> begin
        t = 0.1; x1 = sqrt(1 + t^2 / 2); v1 = t / (2 * x1)
        ComplexF64[x1, -x1, v1, -v1]
    end
    push!(sys, Sys("calogero-moser", :meromorphic, 4, 15, 0.1, cm,
        ComplexF64[1, -1, 0, 0], 0.0 + 0im, cmoracle))

    # 5. Noumi–Yamada A₄ — meromorphic high-d, d=5; conservation oracle Σf=t.
    #    THE headline high-d target (~6.4e-6 @+2 vs ~2e-9 @+1).
    nyα = ComplexF64[0.30, 0.10, 0.25, 0.15, 0.20]
    nyrhs = noumi_yamada_rhs(2; α = nyα)
    push!(sys, Sys("noumi-yamada-A₄", :meromorphic, 5, 12, 0.3, nyrhs,
        ComplexF64[0.7, -0.3, 0.5, -0.55, -0.35], 0.0 + 0im, :conservation))

    return sys
end
