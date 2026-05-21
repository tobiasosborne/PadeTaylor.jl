# Step 4 — the BC must SELECT the tritronquee uniquely.
#
# Step 3b proved: the (u,u')-only 2+2 BC admits >=2 P_I^(2) solutions with
# identical (u,u') boundary data; continuation cannot escape the wrong one
# because the wrong branch satisfies the SAME BVP. The fix is a BC that
# the wrong branch does NOT satisfy. Options:
#   S1: 2+2 split pinning DIFFERENT derivative orders at the two ends so
#       the 4 pinned components together over-specify enough to exclude
#       the wrong branch:   (u,u'') at z_a  +  (u,u'') at z_b   etc.
#   S2: pin (u,u',u'',u''') split 2+2 as (u,u') at z_a + (u'',u''') at z_b
#   S3: pin all 4 at z_a + all 4 at z_b  -> 8 eqs, d=4 -> LEAST SQUARES
#       (would need solver change; here we test via a 4+4 that the SOLVER
#       cannot do — so instead probe whether ANY 2+2 selection works).
#
# Also: verify the n_terms=2 asymptotic seed accuracy for u'' and u'''.

using PadeTaylor.PainleveHierarchy: painleve_hierarchy, painleve_hierarchy_jacobian,
                                    pI2_tritronquee_ic
using PadeTaylor.VectorBVP: vector_bvp_solve
using LinearAlgebra

t   = 0.0
f   = painleve_hierarchy(:I, 2; t = t)
Jf  = painleve_hierarchy_jacobian(:I, 2; t = t)
u_kkg(x)   = -cbrt(6.0) * abs(x)^(1//3)
seed(z; n) = pI2_tritronquee_ic(z; t = t, n_terms = n)

function kkg_err(sol, za, zb)
    xs = range(za, zb; length = 25)
    maximum(abs(real(sol(x)[1]) - u_kkg(x)) for x in xs)
end

println("="^70)
println("STEP 4 — BC selection: which 2+2 split picks the tritronquee?")
println("="^70)

# all 2+2 splits: choose 2 component indices pinned at z_a, the other 2 at z_b
# (a "complementary" 2+2) AND non-complementary (same 2 at both ends).
function trial(name, a_idx, b_idx, za, zb, N, nt; tol=1e-9)
    ya = seed(za;n=nt); yb = seed(zb;n=nt)
    Ba = zeros(4,4); Bb = zeros(4,4)
    g  = zeros(4)
    row = 1
    for i in a_idx; Ba[row,i]=1; g[row]=ya[i]; row+=1; end
    for i in b_idx; Bb[row,i]=1; g[row]=yb[i]; row+=1; end
    try
        sol = vector_bvp_solve(f, za, zb, Ba, Bb, g;
                               N=N, tol=tol, maxiter=120, jacobian=Jf,
                               initial_guess = z->seed(z;n=nt))
        e = kkg_err(sol, za, zb)
        println("  [$name] a=$a_idx b=$b_idx N=$N -> iters=$(sol.iterations) ",
                "ODEres=", round(sol.residual_inf,sigdigits=3),
                " KKGerr=", round(e,sigdigits=4), "  ",
                e < 5e-3 ? "*** TRITRONQUEE ***" : "wrong")
        return (sol, e)
    catch ex
        println("  [$name] a=$a_idx b=$b_idx N=$N threw: ",
                first(sprint(showerror,ex),110))
        return (nothing, Inf)
    end
end

# short interval first (best seed accuracy)
za, zb, nt = -20.0, -16.0, 2
println("\n### Segment [$za,$zb], N=80, all distinct 2+2 splits")
splits = [
  ("uu'|uu'      ", [1,2],[1,2]),
  ("uu'|u'u''    ", [1,2],[2,3]),
  ("uu'|u''u'''  ", [1,2],[3,4]),
  ("uu''|uu''    ", [1,3],[1,3]),
  ("uu''|u'u'''  ", [1,3],[2,4]),
  ("uu'''|uu'''  ", [1,4],[1,4]),
  ("u'u''|u'u''  ", [2,3],[2,3]),
  ("uu'u''|u (3+1)",[1,2,3],[1]),
  ("u|uu'u'''(1+3)",[1],[1,2,4]),
  ("uu'u'''|u'' ", [1,2,4],[3]),
  ("uu''u'''|u  ", [1,3,4],[1]),
]
for (nm, ai, bi) in splits
    trial(nm, ai, bi, za, zb, 80, nt)
end

# also probe seed accuracy: integrate the ODE with the seed as IC over a
# tiny step using a high-order check is overkill; instead just print the
# seed vs what the leading KKG asymptotic predicts for u'' , u'''.
println("\n### Seed accuracy (n_terms=2) vs leading-order asymptotic")
for x in (-20.0, -16.0, -8.0, -4.0)
    s1 = seed(x; n=1); s2 = seed(x; n=2)
    println("  x=$x  n1=", round.(s1,sigdigits=5), "  n2=", round.(s2,sigdigits=5))
end
