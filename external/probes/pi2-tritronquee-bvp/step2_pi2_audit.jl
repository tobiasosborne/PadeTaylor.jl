# Step 2 — reproduce the P_I^(2) tritronquee BVP attempt and audit the BCs.
#
# Step 1 proved vector_bvp_solve is sound at d=4 / 2+2-split (solErr ~1e-15,
# BCres = 0). So if the P_I^(2) attempt converges to u POSITIVE while the BC
# pins u NEGATIVE, the question is: does the converged solution actually
# satisfy its own BC?  Audit that here.

using PadeTaylor.PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP: vector_bvp_solve
using LinearAlgebra

t   = 0.0
f   = painleve_hierarchy(:I, 2; t = t)
Jf  = painleve_hierarchy_jacobian(:I, 2; t = t)

# KKG asymptotic on the negative real axis: u ~ -cbrt(6)*|x|^(1/3)
u_kkg(x) = -cbrt(6.0) * abs(x)^(1//3)

println("="^70)
println("STEP 2 — P_I^(2) tritronquee BVP audit  (t = $t)")
println("="^70)

function run_case(za, zb, N; n_terms = 2, guess_seed = true, label = "")
    ya = pI2_tritronquee_ic(za; t = t, n_terms = n_terms)   # at z_a (more neg)
    yb = pI2_tritronquee_ic(zb; t = t, n_terms = n_terms)   # at z_b
    # 2+2 split: pin (u,u') at each end
    Ba = zeros(4,4); Bb = zeros(4,4)
    Ba[1,1] = 1; Ba[2,2] = 1
    Bb[3,1] = 1; Bb[4,2] = 1
    g = [ya[1], ya[2], yb[1], yb[2]]
    guess = guess_seed ?
        (z -> pI2_tritronquee_ic(z; t = t, n_terms = n_terms)) : nothing
    sol = vector_bvp_solve(f, za, zb, Ba, Bb, g;
                           N = N, maxiter = 40, jacobian = Jf,
                           initial_guess = guess)
    ya_num = sol(za); yb_num = sol(zb)
    bcres  = norm(Ba*ya_num + Bb*yb_num - g, Inf)
    # interior u sample
    xs = range(za, zb; length = 9)
    uvals = [real(sol(x)[1]) for x in xs]
    kkgvals = [u_kkg(x) for x in xs]
    println("\n[$label]  segment [$za, $zb], N=$N, n_terms=$n_terms, ",
            guess_seed ? "seeded guess" : "ZERO guess")
    println("  iters=$(sol.iterations)  ODEres=$(sol.residual_inf)")
    println("  BC g            = $g")
    println("  y(z_a) numeric  = $ya_num")
    println("  y(z_b) numeric  = $yb_num")
    println("  BC residual ‖B_a y_a + B_b y_b - g‖_inf = $bcres")
    println("  interior u (numeric): ", round.(uvals, digits=3))
    println("  interior u (KKG)    : ", round.(kkgvals, digits=3))
    println("  seed u(z_a)=$(ya[1])  seed u(z_b)=$(yb[1])  (both should be <0)")
    return sol
end

try
    run_case(-20.0, -2.0, 60;  label = "A: [-20,-2] N=60 seeded")
catch e
    println("\n[A] threw: ", sprint(showerror, e))
end

try
    run_case(-20.0, -12.0, 80; label = "B: [-20,-12] N=80 seeded")
catch e
    println("\n[B] threw: ", sprint(showerror, e))
end

try
    run_case(-20.0, -18.0, 60; label = "C: [-20,-18] N=60 seeded (short)")
catch e
    println("\n[C] threw: ", sprint(showerror, e))
end

try
    run_case(-20.0, -2.0, 60;  guess_seed = false, label = "D: [-20,-2] N=60 ZERO guess")
catch e
    println("\n[D] threw: ", sprint(showerror, e))
end
