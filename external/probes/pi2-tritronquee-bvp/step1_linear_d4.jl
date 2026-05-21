# Step 1 — d=4 LINEAR sanity test of vector_bvp_solve.
#
# ODE  u'''' = u  has the closed-form basis {cosh, sinh, cos, sin}.
# Pick a specific solution and check vector_bvp_solve recovers it AND
# satisfies its own boundary condition B_a y(z_a) + B_b y(z_b) = g.
#
# Companion system  y = (u, u', u'', u'''),  y' = (y2, y3, y4, y1).

using PadeTaylor.VectorBVP: vector_bvp_solve
using LinearAlgebra

# closed-form solution: u = cosh(x) + 0.5*sin(x)  (a generic combination)
uexact(x)    =  cosh(x) + 0.5*sin(x)
u1(x)        =  sinh(x) + 0.5*cos(x)
u2(x)        =  cosh(x) - 0.5*sin(x)
u3(x)        =  sinh(x) - 0.5*cos(x)
yexact(x)    = [uexact(x), u1(x), u2(x), u3(x)]

f_lin(x, y)  = [y[2], y[3], y[4], y[1]]      # u'''' = u

za, zb = -1.0, 2.0

println("="^70)
println("STEP 1 — d=4 linear ODE  u'''' = u,  segment [$za, $zb]")
println("="^70)

function bc_residual(sol, Ba, Bb, g)
    ya = sol(sol.z_a); yb = sol(sol.z_b)
    norm(Ba*ya + Bb*yb - g, Inf)
end

function maxerr(sol)
    xs = range(za, zb; length = 41)
    maximum(norm(sol(x) - yexact(x), Inf) for x in xs)
end

# --- variant A: 2+2 split, pin (y1,y2) at each end -----------------------
println("\n--- Variant A: 2+2 split, pin (u,u') at BOTH ends ---")
let
    Ba = zeros(4,4); Bb = zeros(4,4)
    # rows 1,2: pin y1,y2 at z_a   ; rows 3,4: pin y1,y2 at z_b
    Ba[1,1] = 1; Ba[2,2] = 1
    Bb[3,1] = 1; Bb[4,2] = 1
    g = [yexact(za)[1], yexact(za)[2], yexact(zb)[1], yexact(zb)[2]]
    for N in (24, 40, 64)
        sol = vector_bvp_solve(f_lin, za, zb, Ba, Bb, g; N=N, maxiter=20)
        println("  N=$N  iters=$(sol.iterations)  ODEres=$(sol.residual_inf)  ",
                "BCres=$(bc_residual(sol,Ba,Bb,g))  solErr=$(maxerr(sol))")
    end
end

# --- variant B: 1+3 split ------------------------------------------------
println("\n--- Variant B: 1+3 split, pin u at z_a; u,u',u'' at z_b ---")
let
    Ba = zeros(4,4); Bb = zeros(4,4)
    Ba[1,1] = 1
    Bb[2,1] = 1; Bb[3,2] = 1; Bb[4,3] = 1
    g = [yexact(za)[1], yexact(zb)[1], yexact(zb)[2], yexact(zb)[3]]
    for N in (24, 40, 64)
        sol = vector_bvp_solve(f_lin, za, zb, Ba, Bb, g; N=N, maxiter=20)
        println("  N=$N  iters=$(sol.iterations)  ODEres=$(sol.residual_inf)  ",
                "BCres=$(bc_residual(sol,Ba,Bb,g))  solErr=$(maxerr(sol))")
    end
end

# --- variant C: 3+1 split ------------------------------------------------
println("\n--- Variant C: 3+1 split, pin u,u',u'' at z_a; u at z_b ---")
let
    Ba = zeros(4,4); Bb = zeros(4,4)
    Ba[1,1] = 1; Ba[2,2] = 1; Ba[3,3] = 1
    Bb[4,1] = 1
    g = [yexact(za)[1], yexact(za)[2], yexact(za)[3], yexact(zb)[1]]
    for N in (24, 40, 64)
        sol = vector_bvp_solve(f_lin, za, zb, Ba, Bb, g; N=N, maxiter=20)
        println("  N=$N  iters=$(sol.iterations)  ODEres=$(sol.residual_inf)  ",
                "BCres=$(bc_residual(sol,Ba,Bb,g))  solErr=$(maxerr(sol))")
    end
end

# --- variant D: 2+2 split, pin (y1,y3) at z_a, (y2,y4) at z_b ------------
println("\n--- Variant D: 2+2 split, pin (u,u'') at z_a; (u',u''') at z_b ---")
let
    Ba = zeros(4,4); Bb = zeros(4,4)
    Ba[1,1] = 1; Ba[2,3] = 1
    Bb[3,2] = 1; Bb[4,4] = 1
    g = [yexact(za)[1], yexact(za)[3], yexact(zb)[2], yexact(zb)[4]]
    for N in (24, 40, 64)
        sol = vector_bvp_solve(f_lin, za, zb, Ba, Bb, g; N=N, maxiter=20)
        println("  N=$N  iters=$(sol.iterations)  ODEres=$(sol.residual_inf)  ",
                "BCres=$(bc_residual(sol,Ba,Bb,g))  solErr=$(maxerr(sol))")
    end
end

# --- endpoint check at the node level (bypass barycentric eval) ----------
println("\n--- Variant A re-run: print raw node-0 / node-N states ---")
let
    Ba = zeros(4,4); Bb = zeros(4,4)
    Ba[1,1] = 1; Ba[2,2] = 1
    Bb[3,1] = 1; Bb[4,2] = 1
    g = [yexact(za)[1], yexact(za)[2], yexact(zb)[1], yexact(zb)[2]]
    sol = vector_bvp_solve(f_lin, za, zb, Ba, Bb, g; N=40, maxiter=20)
    Y0 = sol.Y_nodes[1]; YN = sol.Y_nodes[end]
    println("  node 0 (t=+1, z=$(sol.nodes_z[1]))     Y = $Y0")
    println("  node N (t=-1, z=$(sol.nodes_z[end]))   Y = $YN")
    println("  yexact(z_a=$za) = $(yexact(za))")
    println("  yexact(z_b=$zb) = $(yexact(zb))")
    println("  => node 0 should match z_b, node N should match z_a")
    println("  BC rows: B_a*y(z_a)+B_b*y(z_b)-g = ",
            Ba*sol(za) + Bb*sol(zb) - g)
end
