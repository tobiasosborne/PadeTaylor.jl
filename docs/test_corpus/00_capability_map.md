# PadeTaylor.jl Capability and Coverage Map

This is the capability and coverage map produced in Phase 1 of epic **padetaylor-p1v0** (comprehensive ground-truth test corpus), to drive a corpus that exercises every PadeTaylor capability and roots out correctness bugs.
It was produced by 8 read-only slice readers plus a synthesis pass.

## Capability buckets (B1–B22)

| ID | Bucket | Priority | Coverage | Ground-truth strategy |
|----|--------|----------|----------|-----------------------|
| B1 | B1 — Robust Padé core (scalar): SVD, GGT-2013 Algorithm 2 (:svd), FW-2011 classical Toeplitz (:classical), tol/method dispatch | P1 | strong | closed-form (Padé of exp/geometric) + cross-solver (a 256-bit mpmath/Mathematica re-run of padeapprox-equivalent for the same 6 cases to pin the BigFloat :svd path, currently pinned only against closed-form exp); Mathematica PadeApproximant for a coefficient-level (m,m) classical oracle at m>=10 on a non-exp function |
| B2 | B2 — Scalar Taylor-jet generation: 1st-order (method b bootstrap) and 2nd-order (u,u' joint evolution + up-resync) | P0 | moderate | published-table / Mathematica: Series[WeierstrassP] and Series[PainleveI] at WorkingPrecision>=50 to order 30+ for u''=6u^2 AND u''=6u^2+z at BigFloat-256, mirroring _oracle_coefficients.jl; plus a ComplexF64 1st-order jet vs Mathematica Series; direct assertion on returned up (differentiate the u-series) |
| B3 | B3 — Scalar Padé-Taylor stepper (fixed-h): one FW 5-step recipe step + analytic dense-output evaluation | P1 | strong | closed-form (known rational 1/(1-z/2) for standalone _evaluate_pade/_deriv at several z + a z hitting a denominator root for DomainError) + cross-solver (mpmath.odefun / NDSolve WP=50 for a direct BigFloat-256 pade_step! on Weierstrass, and an odd-order step) |
| B4 | B4 — Scalar adaptive (FFW-2017) controller: truncation-error estimate, q-rescale, adaptive step loop | P1 | moderate | closed-form (independently-derived FFW eps_{n+1}*h^{n+1}/a(h) on a non-exp function via mpmath to pin ffw_truncation_error directly) + crafted pathological RHS/h forcing max_rescales exhaustion (fail-loud) + BigFloat adaptive step vs mpmath.odefun |
| B5 | B5 — Scalar IVP driver: PadeTaylorProblem builder + solve_pade + PadeTaylorSolution dense callable | P0 | moderate | published-table (FW Table 5.1, already pinned) + closed-form Weierstrass for a MULTI-SEGMENT solve (h_max=0.5 over [0,1.5]) exercising the segment loop + k-scan at several z; cheap throw-tests (h_max<=0, scalar-y0 rejection, max_steps=1 exhaustion, out-of-range dense eval); a within-radius taylor_eval value pin (not just divergence) |
| B6 | B6 — Orphaned scalar step-control: Jorba-Zou first-step + FW pole-distance heuristic (DEAD in scalar pipeline) | P2 | weak | cross-solver (TaylorIntegration.stepsize on an odd-parity sin-like jet to trigger _second_stepsize; eps_rel*\|c0\|>eps_abs relative-mode case) + closed-form (hand-computed projection for a complex-direction step with poles either side; all-poles-behind full-distance return; constant-denominator nu=0) |
| B7 | B7 — Path-network whole-plane IVP solver (Stage-1 wedge walk + Stage-2 dense fill) and per-point accessors | P1 | strong | closed-form (℘ off-node eval_at value pin; hand-built visited_z fixture with exact-distance ties for the sheet-blind tiebreak; constructed evals fixture where a forbidden Inf-u candidate is the S^1-closest wedge to confirm :steepest_descent + refuse-mode interaction) + published-table (high-precision Painlevé pole/value table from Mathematica NDSolve / mpmath for PI/PII/PIV/PVI to extend beyond the single ℘ closed-form — the slice's biggest ground-truth hole) |
| B8 | B8 — Pole-field detection: Laplacian residual, mask classifier, scalar/network pole extraction (Froissart filter) | P1 | strong | closed-form (hand-built PadeApproximant with a deliberate near-coincident pole-zero doublet to pin the min_residue drop; a known double pole asserting exactly 1 cluster; \|∇²(1/z)\|=2/\|z\|³ magnitude pin) + published-table (Mathematica/mpmath Painlevé pole tables to extend beyond the ℘ lattice) |
| B9 | B9 — Edge-gated region-growing solver + morphological operators (dilate/erode/open/flood-fill) | P1 | weak | closed-form (hand-built small BitMatrix fixtures with textbook erosion/dilation/opening/connected-component results — trivial, currently absent) + published-table (high-precision PI tritronquée pole positions via Mathematica NDSolve / Boutroux-tritronquée pole table to pin extracted pole accuracy, not just sector counts) |
| B10 | B10 — Multi-sheet machinery: Painlevé coordinate transforms (PIII/PV/PVI zeta + PVI eta double-exp) and transformed RHS closures | P0 | moderate | closed-form (algebraic transform at multiple complex z incl non-principal sheet and complex-z eta round-trip; recompute at setprecision(256) to bound the E-1 cancellation precision loss) + cross-solver (independent mpmath/Mathematica PIII/PV/PVI value on the transformed plane to break the loose-vs-tight self-reference) |
| B11 | B11 — Riemann-sheet winding bookkeeping: winding_delta, accumulate_winding/sheet_index, segment_crosses_cut, multi-branch step_sheet_update | P0 | moderate | closed-form (oblique-angle segment_crosses_cut at alpha=π/4 with branch at 1+1i hand-computed; a step with \|dθ\|>=π asserting either fail-loud per Rule 1 or an explicit pin of the documented-wrong behaviour; a hand-built step crossing two cuts confirming both counters bump; a realistic many-step CCW circumambulation with known winding number) |
| B12 | B12 — Painlevé problem/solution wrappers: PainleveProblem builder (6 eqs + 2 PVI frames), solver forwarding, PainleveSolution access surface (callable/poles/grid_values) | P1 | moderate | closed-form (hand-built PathNetworkSolution wrapped as :transformed PIII and as :transformed_eta PVI asserting poles map via exp(zeta/2)/exp(exp(eta)) and grid_values remap up — isolates wrapper frame-mapping from solver; near-branch-point IC z0=branch+1e-12 to pin the exact-equality guard blind spot) + cross-solver (end-to-end :transformed_eta path_network_solve whose grid_values round-trip the z-frame grid and match a direct PVI solve) |
| B13 | B13 — Vector Taylor jets, vector step-control, vector Padé step (shared-Q cell A/B) | P0 | moderate | closed-form (complex vector ODE y'=i*y coeffs; d>=2 shared-pole complex-h step to a known z) + cross-solver (determinant/AAA at d=2,3,4,m=4,BigFloat,Rational for cell-B with a TRANSCENDENTAL d=1 discriminator + 4-mutation proof; an mpmath/Mathematica complex-jet coefficient table for a coupled d>=2 meromorphic system) |
| B14 | B14 — Vector IVP/path-network driver, wedge selector, adaptive-h, Stage-2 fill, shared-Q pole extraction | P0 | moderate | closed-form (hand-built 5-candidate wedge fixture with closed-form discs + planted doublet not shrinking disc; hand-built sol with one genuine + one doublet asserting Froissart-filtered; pole-on-chord geometry to deterministically fire :all_candidates_failed/:step_collapse) + cross-solver (a Calogero-Moser or coupled meromorphic system whose trajectory passes THROUGH a shared complex pole, mpmath/NDSolve-pinned — the v0.2 shared-Q pole-bridging keystone is currently never validated against an external vector transcendent value) |
| B15 | B15 — Scalar spectral BVP: Chebyshev-Newton collocation (2-arg and 3-arg-u' RHS), barycentric callable, D1/D2 matrices | P0 | moderate | closed-form (3-arg: linear variable-coefficient u''=a(z)u'+b(z)u or nonlinear u''=(u')^2/u with u=A*e^{kz}, asserting node values to spectral floor — isolates the D1-coupled Jacobian term; complex-oblique u''=u with z_b=1+1i and cosh BCs to test the real(t*) guard; a rank-deficient/singular-Jacobian fixture) |
| B16 | B16 — Vector spectral BVP (general linear two-point BC) and 2D Laplace harmonic fill | P1 | moderate | closed-form (vector BVP with a GENUINELY general endpoint-coupling/periodic BC and known solution — e.g. periodic BC on y'=(y2,-y1) giving sin/cos, or a mixing BC — to exercise B_a*Y_N+B_b*Y_0 assembly non-trivially; rank-deficient-B fixture for the singular-Jacobian throw; inconsistent-corner Laplace2D data to confirm fail-fast or document silent behaviour) |
| B17 | B17 — IVP<->BVP chain/lattice dispatch and pole-free-sector IVP+BVP hybrid | P0 | weak | closed-form (nonlinear chain weierstrass-P/PI through IVP->BVP->IVP; a BVP-first-segment chain; cosh value-check on the edge-gated lattice AUTO path; a discontinuous slice stack confirming the promised glue assertion fires — currently it does NOT, so this test would expose the unimplemented guard) + cross-solver (independent BF-256+many-term-series mpmath/Mathematica PIII tronquée sector reference to replace the circular self-check) |
| B18 | B18 — Special-function evaluators: Heun (general + confluent), Frobenius recurrences, RHS closures | P1 | moderate | published-table (consume the full 42-record oracles.txt spanning a in {complex, \|a\|<1, \|a\|>1}, complex/large params; Mathematica SeriesCoefficient[HeunG/HeunC] for direct element-by-element coefficient-vector oracles; a BigFloat-256 Mathematica run for the ADR-0018 #4 high-precision path; connection-formula reference for HeunC z>=1 to tighten from 1e-3) |
| B19 | B19 — Painlevé closed-form families (PII rational/Airy, PIV entire) — self-validating oracles + pole-crossing acceptance | P1 | strong | closed-form (pii_airy(0; θ=π/2) IC + solver cross-check where both Ai/Bi contribute; a direct _u_pii_airy(1) derivative assertion vs finite-difference-of-closed-form; integrate u_3 from (0,0) ACROSS z≈1.508 asserting match on the far side — the Padé-bridge stress case) + cross-solver (Mathematica symbolic confirmation that u=-(2/3)z forces (α,β)=(0,-2/9) uniquely) |
| B20 | B20 — Named transcendents (tritronquée PI, Hastings-McLeod PII) + PI-hierarchy (P_I^(2)) builders, Jacobian, asymptotic IC | P1 | moderate | published-table (FW2011 published tritronquée pole locations / first-pole-on-negative-real-axis, OR a high-precision asymptotic u(z)~-sqrt(-z/6) match integrated FROM the constructor IC — to confirm the IC truly selects the tritronquée; an independent-of-FW2014 literature value for the Hastings-McLeod u(0) connection constant) + cross-solver (mpmath Taylor-method P_I^(2) integration from a verified asymptotic IC to pin one downstream u(x); a θ-boundary test for the (2π,4π] half-open sheet window) |
| B21 | B21 — Noumi-Yamada A_{2n}^(1) systems: even-parity RHS builder, rational solutions, affine-Weyl Bäcklund symmetry | P2 | weak | cross-solver (mpmath/Mathematica integration of A_2^(1) or A_4^(1) from a fixed IC to a downstream t — an independent numerical NON-rational NY transcendent value to validate the vector walk against something external) + closed-form (n=3 A_6 group-relation spot-checks; a longer Weyl-word orbit to a known rational solution) |
| B22 | B22 — Cross-cutting infra: oracle-regeneration harnesses (verify.jl), shared-Padé live oracles, duplicated _chebyshev_D1 (3 copies) | P1 | moderate | cross-solver (extend the live-independent-impl pattern to vector stepper and BVP — run a small in-suite mpmath-equivalent / RK reference rather than only frozen Float64 literals; regenerate the empty pathnetwork-long-range BF-256 oracle) + closed-form (a trivial cross-module test asserting BVP._chebyshev_D1 == VectorBVP._chebyshev_D1 == Laplace2D._cheb_D1 elementwise to catch copy drift) |

## Highest-value gaps (most dangerous first)

Gap #1 (B17) has been CONFIRMED and filed as bead padetaylor-q0yq.

1. CONFIRMED LIVE DOC-VS-CODE DISCREPANCY (most dangerous, B17): solve_pole_free_hybrid / IVPBVPSolution docstrings (src/IVPBVPHybrid.jl:64, :726-729, :736) promise the callable 'asserts BVP and PFS values agree to within glue_tol; Failure throws', but the callable body (src/IVPBVPHybrid.jl:738-781) performs ONLY inside/outside dispatch + linear interpolation — there is NO dual-side evaluation and NO glue-continuity check anywhere. A discontinuous PFS-vs-BVP glue would pass silently. The fail-loud Rule-1 guard the user believes exists does not. A test feeding a discontinuous slice stack expecting a throw would expose this immediately.

2. BigFloat 2nd-order Taylor recurrence has NO coefficient-level oracle (B2): taylor_coefficients_2nd (src/Coefficients.jl:171) — the recurrence under the ENTIRE scalar Painlevé/Weierstrass pipeline — is oracle-pinned ONLY at Float64 (test 3.1.3); BigFloat is checked solely end-to-end via solve_pade value (6.1.6). A precision-driven sign/off-by-one specific to BigFloat, or a corrupted up-resync (up[j-1]=j*u[j]) that self-corrects in u, is invisible.

3. The (e^zeta-1) and double-exp (E-1) catastrophic-cancellation hotspots near the PVI branch lattice have ZERO BigFloat coverage (B10): pVI_transformed_rhs (src/SheetTracker.jl:155) and pVI_eta_transformed_rhs (src/SheetTracker.jl:206) cancel near zeta=2πik / eta near log(2π) where the entire multi-sheet machinery operates. A cancellation bug invisible at Float64 (the only precision tested across the whole multivalued slice) would corrupt every sheet-1/+-2 value — and those values are themselves only loosely self-pinned (see next).

4. winding_delta |dtheta|>=pi precondition is NOT enforced and NOT tested (B11): src/SheetTracker.jl:283 silently loses a full revolution if a single step encloses a branch (CONFIRMED: no guard in source). accumulate_winding/sheet_index/step_sheet_update all inherit this. A walker producing too-large steps near a branch corrupts every downstream sheet index with no test catching it — the single most likely place a real sheet-tracking correctness bug hides.

5. Vector shared-Q pole-BRIDGING (the v0.2 keystone) is never validated against an external transcendent through a complex pole (B14): vector_path_network_solve / vector_pade_step_with_pade! (src/VectorStepper.jl:225, src/VectorPathNetwork.jl:642) are pinned only against closed-form poles (Riccati d=3, Weierstrass which is a scalar-oracle re-encoding) and the smooth-only Calogero-Moser. No oracle exercises a vector step landing a shared Padé STRADDLING a complex pole for a coupled system with a known pole location. extract_poles_shared_q's Froissart filter (src/VectorPoleField.jl:294) is not directly tested.

6. Scalar BVP 3-arg (u'-dependent) overload has NO closed-form oracle (B15): bvp_solve 3-arg (src/BVP.jl:370) — its D1-coupled NON-diagonal Jacobian term (scale/half_diff)*diag(df/dup)*D1 could be mis-scaled by an O(1) factor and survive, since it is exercised ONLY through PIII hybrid self-consistency at 5e-2 thresholds. A linear variable-coefficient or u=A*e^{kz} closed-form test would isolate it.

7. The real(t*)-only domain guard silently extrapolates off oblique-complex segments (B15/B16): CONFIRMED identical at src/BVP.jl:491 and src/VectorBVP.jl:318 — an off-axis complex z with |real(t*)|<=1 but large imag(t*) PASSES the guard and the barycentric interpolant silently extrapolates off-segment. No scalar oracle uses an oblique complex segment; the figures (FW Fig 4.1) actually do.

8. solve_pade multi-segment loop + callable k-scan are NEVER exercised by solve_pade itself (B5): every problems_test case is single-segment (h_max=1.5>=window). The multi-segment while-loop, multi-push of h_vec/pade_vec, and the callable's while-k-scan (src/Problems.jl:217) are untested for the headline scalar driver. A k-scan off-by-one or wrong-segment selection in dense interpolation would pass the entire scalar suite.

9. FFW multi-sheet figure values are self-referentially pinned (B10/B22): the ffw_fig_1..7 PRIMARY oracle is a loose-vs-tight in-house cross-check; published-number comparison is deliberately one order looser than FFW; sheet +-1/+-2 thresholds reach 8e-1. A systematic multi-sheet bug present in BOTH the loose and tight runs passes undetected. No external high-precision multi-sheet oracle exists anywhere in the suite.

10. VectorBVP general/coupling/periodic linear BC is untested despite being ADR-0023's entire rationale (B16): vector_bvp_solve (src/VectorBVP.jl:180) is tested ONLY with Dirichlet SELECTOR B matrices ([1 0;0 0]); a bug in the B_a*Y_N + B_b*Y_0 tau-method assembly that manifests only for an endpoint-coupling/mixing/periodic B (the stated purpose) would survive. Singular-Jacobian throw paths (scalar AND vector) also untested.

11. Edge-gated morphology (dilate/erode/open/flood-fill) has ZERO direct unit tests (B9): src/EdgeGatedSolve.jl:140-213 covered only transitively through ONE tritronquée fixture. An 8-vs-4-connectivity error or an erosion grid-edge boundary bug (lines 172-174) would be invisible, AND LatticeDispatcher._dilate_one (src/LatticeDispatcher.jl:447) assumes byte-matching connectivity for the up_grid derivative scatter with no cross-module test — a 4-vs-8 mismatch silently misplaces derivative values at diagonal frontier cells.

12. :steepest_descent + refuse-mode Inf-sentinel interaction is an untested latent bug (B7): _select_candidate (src/PathNetwork.jl:987) :steepest_descent selects purely on angle θ_sd ignoring u, so it can pick a forbidden candidate (whose u was nullified to Inf by the refuse-mode filter) that only fails at the downstream isfinite guard — silently stepping through a branch cut. Never tested in combination.

13. SharedPadeCellB +1-window m_eff index has no transcendental discriminator (B13): build_square_cell (src/SharedPadeCellB.jl:109) — the sole cell-B test (SPD.1) uses an EXACT RATIONAL where an off-by-one in the m_eff=ceil(m/d)*d+1 window index survives; d==m boundary, cnorm<=tol, post-QR, and collapse->nothing branches all untested.

14. Three verbatim _chebyshev_D1 copies can silently diverge (B22): CONFIRMED at src/BVP.jl:525, src/VectorBVP.jl, src/Laplace2D.jl:186 (_cheb_D1) — only BVP's copy is DMSUITE-oracled directly; a fix or bug in one copy will not propagate and no test asserts byte-identity.

15. PainleveSolution PIII and eta frame-mappings are never exercised through the access surface (B12): all :transformed callable/poles/grid_values tests use the PV single-exp map; PIII (z=exp(zeta/2), u=w/z) and eta double-exp (z=exp(exp(eta)), up=vp/(z*zeta)) wrapper remaps (src/PainleveSolution.jl:159/194) are untested, AND poles() uses _coord projection (dropping up) while grid_values uses the full from_frame triple — a bug in the eta/PIII up-component is untested. :transformed_eta is NEVER forwarded through any solver (src/Painleve.jl:340/365).

## Corpus search directives (Phase 2 seed)

1. Find high-precision pole-location and on-axis value tables for the Painlevé-I tritronquée and the P_I^(2) second-member tritronquée (KKG 2015 / Boutroux): the first poles on the negative real axis and a few interior u(x) values, computable via Mathematica NDSolve WorkingPrecision>=40 or mpmath from a verified asymptotic IC. Targets B7, B9, B20 (currently the ONLY pole oracle is the closed-form Weierstrass lattice; no Painlevé pole table exists anywhere).

2. Find / derive closed-form solutions of 3-arg (u'-dependent) scalar 2nd-order BVPs with Dirichlet data: a linear variable-coefficient u''=a(z)u'+b(z)u with elementary solution, and a nonlinear u''=(u')^2/u with u=A*e^{kz}, on real AND oblique-complex segments (z_b-z_a with both Re and Im nonzero). Targets B15 — isolates the untested D1-coupled Jacobian term and the real(t*)-only domain guard.

3. Find a coupled (d>=2) meromorphic vector ODE with a KNOWN closed-form shared pole location whose trajectory passes THROUGH the pole in the complex plane — e.g. an integrable Calogero-Moser / Ruijsenaars system with colliding poles, or a 2x2 Riccati with a shared denominator — to validate vector shared-Q pole-bridging against an external transcendent. Targets B13, B14 (the v0.2 keystone, currently never externally pinned through a pole).

4. Find closed-form or high-precision-tabulated solutions of the transformed-plane PIII/PV/PVI equations (FFW 2017 eqs in zeta and eta) at specific complex (zeta,w,wp) / (eta,v,vp) points on sheets 0, +-1, +-2, computable via mpmath / Mathematica integration along the SAME winding path FFW uses. Targets B10 — replaces the loose-vs-tight self-referential multi-sheet oracle and exposes the (e^zeta-1)/(E-1) cancellation at BigFloat.

5. Derive Mathematica Series[WeierstrassP] and Series[PainleveI] coefficient tables at WorkingPrecision>=50 to order 30+ for u''=6u^2 and u''=6u^2+z, AND an mpmath/Mathematica Taylor-coefficient table for a coupled d>=2 meromorphic vector system, both at BigFloat-256. Targets B2, B13 — the 2nd-order scalar recurrence and the vector jet recurrence have NO BigFloat coefficient-level oracle.

6. Find the full multi-parameter Heun oracle space (consume the existing 42-record external/probes/heun-oracle/oracles.txt plus extend): HeunG/HeunC at a in {complex, |a|<1, |a|>1}, complex and large-magnitude (q,α,β,γ,δ,ε), Mathematica SeriesCoefficient[HeunG/HeunC,{z,0,k}] for direct coefficient-vector oracles, and a connection-formula reference for HeunC on the z>=1 second sheet (to tighten from 1e-3). Targets B18 — breaks the single-tuple ε=-3 monoculture.

7. Find an independent numerical value (mpmath / Mathematica WP>=50) for a GENUINE Painlevé transcendent integrated to a downstream point: tritronquée PI at u(z) on the pole-free sector, Hastings-McLeod PII at u(x>0) with the published u(0) connection constant from a source independent of FW2014, and the FFW-published PVI value u(10)=0.429534600325223 asserted via an END-TO-END solve from a different IC (currently used only as an IC seed). Targets B20, B12, B19 — converts self-referential IC tautologies into solver oracles.

8. Find / construct a closed-form vector first-order BVP with a GENUINELY general (non-Dirichlet) linear two-point BC and known solution — periodic B_a=I,B_b=-I on y'=(y2,-y1) giving sin/cos, or an endpoint-mixing/coupling B — plus a rank-deficient-B fixture. Targets B16 — exercises the B_a*Y_N+B_b*Y_0 tau-method assembly and the singular-Jacobian throw, the untested core of ADR-0023.

9. Construct textbook morphological-operator fixtures: small hand-built BitMatrices with known 8-connected dilation, erosion (including grid-edge windows hanging off the lattice), opening, and connected-component flood-fill outputs; plus a fixture confirming LatticeDispatcher._dilate_one connectivity matches EdgeGatedSolve._solve_targets. Targets B9, B17 — the morphology modules have zero direct tests.

10. Construct hand-built oracle fixtures for the orphaned and weakly-covered control paths: a TaylorIntegration.stepsize value on an odd-parity sin-like jet triggering _second_stepsize fallback, a complex-direction step_pade_root projection hand-computed, a 5-candidate vector wedge fixture with closed-form pole-free discs and a planted Froissart doublet (residue<1e-4) that must NOT shrink the disc, and a hand-built PadeApproximant with a near-coincident pole-zero to pin the scalar/vector min_residue Froissart drop. Targets B6, B8, B14.

11. Construct sheet-tracking edge-case fixtures: a single walker step with |dtheta|>=pi enclosing a branch (assert fail-loud or pin documented-wrong behaviour), an oblique-cut-angle segment_crosses_cut at alpha=pi/4 with branch at 1+1i hand-computed, a step crossing TWO branch cuts simultaneously, and a realistic fine-step CCW circumambulation with known winding number. Targets B11 — the highest-risk untested sheet-tracking branches.

12. Find a complex-segment (off-real-axis, e.g. vertical imaginary-axis) 2nd-order BVP with an mpmath.odefun reference shot from a known IC, and a BigFloat-256 BVP reference value, to cover the off-real-axis and high-precision BVP regimes the figures actually exercise but no oracle pins. Targets B15, B22 — _oracle_bvp.jl is Float64 and real-segment only (its '_complex_N24' case is actually real [0.1,0.5]).

## Per-slice capability detail

### core-pade

**pade_svd** (internal-algorithm) — `src/LinAlg.jl:75 (LAPACK), src/LinAlg.jl:94 (generic real), src/LinAlg.jl:106 (generic complex)`

SVD of the n x (n+1) lower-triangular Toeplitz that GGT 2013 Algorithm 2 (eq 2.6) needs for rank counting + null-vector extraction. ADR-0002: relative-accuracy c*2^-p*sigma_i per SV via Demmel-Veselic Jacobi is load-bearing for the sigma_i < tol*||c||_2 threshold at arb-prec.

Regimes: Float64/Float32 (LAPACK absolute-accuracy DK); Complex{Float64/32}; BigFloat-256 (Jacobi relative-accuracy); Complex{BigFloat}; full=false thin SVD; full=true (needed for the (n+1)-th null column); rank-deficient wide n x (n+1) matrices; Arb path lives only in ext/PadeTaylorArblibExt.jl (out of this slice)

Current coverage: test/linalg_test.jl 1.1.1-1.1.4 + 1.1.3b. STRONG (closed-form/published): 1.1.1 vs svdvals atol 1e-14; 1.1.2 Hilbert-10 F64 recon; 1.1.3 Hilbert-10 BigFloat-256 recon <1e-65 + eltype===BigFloat; 1.1.3b full=true (n+1)x(n+1) Vt and A*v_null~0; 1.1.4 rank-deficient [1 2;2 4] BigFloat S[2]/S[1]<2^-100. Mutation 1.1.5 documented.

Gaps: No Complex{BigFloat} SVD test (only real BigFloat). No near-block-boundary Toeplitz where a small genuine SV sits just above tol*||c|| (the actual GGT motivation - 1.1.4 uses a CLEAN zero, not a near-threshold SV). full=true at BigFloat untested. No degenerate-shape (1x2) SVD test.

Ground truth needed: A hand-built Toeplitz with a known SV at ~10*tol*||c|| where F64 DK misclassifies but BigFloat Jacobi keeps it - proves the ADR-0002 relative-accuracy claim. Mathematica/mpmath SVD of a Complex{BigFloat} matrix for the complex generic path.

**robust_pade** (public-api) — `src/RobustPade.jl:356`

Type-(m,n) Pade approximant via GGT 2013 Algorithm 2 (:svd, port of external/chebfun/padeapprox.m lines 60-143 incl QR-reweighting 111-117 beyond the paper) OR FW 2011 5.1.4 Toeplitz-backslash (:classical). Diagonal-hopping rank reduction, Froissart-doublet removal, common-z^lambda cancellation, b[1]=1 normalisation.

Regimes: :classical (m==n,n>0 -> classical_pade_diagonal; SingularException -> :svd) vs :svd; off-diagonal m!=n with :classical routes to :svd; r==0 special case (max|c[1:m+1]|<=tol*||c||_inf -> a=[0],b=[1],mu=typemin(Int)); diagonal-stripe reduction; Froissart removal (mu!=nu); defect>0 ill-posed boundary; Float64 / BigFloat-256 / Complex{Float64}; tol override; zero-pad when length(c)<m+n+1; n==0 trivial denominator

Current coverage: test/robustpade_test.jl 2.1.1-2.1.8 (:svd) + classical_pade_test.jl CP.1.1-CP.1.9 (:classical & dispatch). STRONG: 2.1.1 exp(2,2) closed-form + Octave oracle; 2.1.2/2.1.3 published GGT reductions (20->7, 20->10) Octave mu/nu oracle + functional eval; 2.1.4 tan(z^4) Froissart (20,16) Octave-FFT oracle; 2.1.5 1+z^2 defect-1 collapse; 2.1.6 noisy tol-thresholded (0,1); 2.1.8 BigFloat-256 no-reduction + eltype + 1e-40 functional. Mutations P1-P4 + QR-reweighting documented.

Gaps: Negative m / negative n ArgumentErrors UNTESTED (only :bogus). _trim_and_normalise error branches (all-b<tol; post-trim b empty; all-a<ts->a=[0]) UNTESTED (unreachable-by-design, unproven). n==0 early-return UNTESTED. No Complex{BigFloat} test. 2.1.3 coefficient match deliberately NOT asserted (functional only) - coefficient correctness near defect>0 rests on functional eval alone.

Ground truth needed: Add Complex{BigFloat} case vs Mathematica; explicit n==0 (pure-polynomial) case; crafted inputs hitting each _trim_and_normalise branch to confirm fail-loud vs silent-wrong.

**classical_pade_diagonal** (internal-algorithm) — `src/RobustPade.jl:259`

FW 2011 5.1.4 eqs (5.4)/(5.5): m x m Toeplitz T[i,j]=c_{m+i-j}, LU solve (\), numerator a_k=sum b_j c_{k-j}. Faster+more accurate than SVD on smooth trajectory inputs (worklog 020: 580-1000x per-step). Throws SingularException on zero pivot (FW line 346).

Regimes: m=0 trivial r==c_0; well-conditioned diagonal; exactly-singular Toeplitz (rank-1 geometric series, defect-1 1+z^2) -> SingularException; Float64 / Complex{Float64} / BigFloat (via :classical override); zero-pad/truncate to 2m+1; m<0 ArgumentError

Current coverage: classical_pade_test.jl CP.1.1 (exp(2,2) closed form vs SVD), CP.1.2 (exp(15,15) functional 1e-14), CP.1.3 (log(1.2-z)(10,10) functional), CP.1.4+CP.1.7 (singular->SingularException), CP.1.8 (m<0 ArgumentError), CP.1.9 (Complex{Float64} eltype + functional). STRONG closed-form/functional. Mutation P1 propagates ~60 RED across whole IVP stack.

Gaps: No BigFloat coefficient-level oracle (only mu/nu shape via override in CP.1.6). m=0 trivial path UNTESTED. eq-(5.5) numerator recurrence for m>2 only validated functionally (CP.1.2/1.3), not at coefficient level.

Ground truth needed: Coefficient-level (m,m) oracle at m>=10 from Mathematica PadeApproximant for a non-exp function to pin the eq-(5.5) recurrence beyond functional eval. BigFloat classical vs high-precision reference.

**default_tol / _default_pade_method** (internal-algorithm) — `src/RobustPade.jl:154 (default_tol), src/RobustPade.jl:176 (_default_pade_method)`

Singular-value threshold + method selection per element type (ADR-0002/ADR-0005). Float64 matches Chebfun 1e-14; arb-prec 10-bit slack above the working-precision floor.

Regimes: Float64; Float32; BigFloat (precision-dependent); Complex element types delegate to real

Current coverage: NONE direct. Implicit: robustpade 2.1.8 relies on BigFloat default tol giving no reduction; classical CP.1.6 relies on _default_pade_method (mutation P3 proves dispatch load-bearing).

Gaps: No unit test asserts default_tol(BigFloat) value at a given setprecision, nor Float32 (the 1f-6 path is entirely unexercised - Float32 appears nowhere in the corpus). _default_pade_method(Float32)/(Complex{Float32}) untested.

Ground truth needed: Cheap closed-form pins currently missing: default_tol(Float64)==1e-14, default_tol BigFloat-256 ==2.0^(-246), the Float32 path.

**taylor_coefficients_1st** (public-api) — `src/Coefficients.jl:114`

FW 2011 2.1.2 method (b) bootstrap for y'=f(z,y): seed y=Taylor1(T[y0]); for k=1..order y[k]=f(z,y)[k-1]/k. Delegates Taylor arithmetic to TaylorSeries.jl Taylor1{T} (RESEARCH 3.3). No silent widening.

Regimes: Float64; BigFloat-256 (order 60); order<1 ArgumentError; type-stability eltype(out)===T; arbitrary analytic RHS

Current coverage: coefficients_test.jl 3.1.1 (exp order 10, wolframscript+sympy oracle rtol 1e-14), 3.1.2 (z^2+y^2 Bessel-ratio order 14 oracle), 3.1.4 (exp order 60 BigFloat-256 vs 80-dps Mathematica/mpmath tol 2^-200). STRONG closed-form/published. Mutation A (k+1 off-by-one) REDs 3.1.1.

Gaps: order<1 ArgumentError UNTESTED. No ComplexF64 1st-order test here (complex jets only in vector slice). Type-stability eltype only checked at F64 and BigFloat - Complex untested.

Ground truth needed: order<1 throw test; a ComplexF64 1st-order jet vs Mathematica Series.

**taylor_coefficients_2nd** (public-api) — `src/Coefficients.jl:171`

FW 2011 2.1.2 method (b) for u''=f(z,u,u'): evolve u and formal derivative up jointly; pass j: u[j]=f(z,u,up)[j-2]/(j(j-1)), resync up[j-1]=j*u[j]. The recurrence underpinning the entire scalar Painleve/Weierstrass IVP pipeline.

Regimes: Float64; BigFloat (transitively via solve_pade 6.1.6, NOT direct); order<2 ArgumentError; Painleve-I RHS (6u^2+z) and Weierstrass RHS (6u^2); near-pole large-coefficient growth (c_30~3e33)

Current coverage: coefficients_test.jl 3.1.3 (u''=6u^2 FW ICs order 30 Float64, wolframscript AsymptoticDSolveValue + Series[WeierstrassP] cross-check rtol 1e-12). STRONG at Float64. Mutation C (j(j+1) off-by-one) REDs 3.1.3.

Gaps: BIGFLOAT 2ND-ORDER COEFFICIENT-LEVEL ORACLE MISSING - only direct at Float64; BigFloat only end-to-end via solve_pade (6.1.6 value, not coefficients). order<2 ArgumentError UNTESTED. +z (Painleve-I) jet only at stepper level (5.1.3), not bare-jet. up-resync invariant only indirectly checked via u output - a bug that corrupts up but self-corrects in u is invisible.

Ground truth needed: BigFloat-256 coefficient-level oracle for u''=6u^2 (and 6u^2+z) at order 30+ from Mathematica Series[WeierstrassP/Painleve] WP>=50, mirroring _oracle_coefficients. A direct assertion on the returned up (differentiate the u-series).

**step_jorba_zou** (internal-algorithm) — `src/StepControl.jl:162`

Jorba-Zou 2005 3.3.1 eq.11 first-step-size in the fixed-order TI.jl form: h=min over k in {p-1,p} of (eps/|c[k+1]|)^(1/k), eps=eps_abs or eps_rel*|c0|; _second_stepsize fallback scan j=1..p-2 when both candidates zero. Verbatim port of TaylorIntegration.jl stepsize.jl.

Regimes: normal last-two-coefficient path; eps_abs vs eps_rel*|c0| selection; zero candidate skip; both c[p-1],c[p] zero -> _second_stepsize fallback; all-zero -> ArgumentError throw; length<2 ArgumentError; Float64 only in tests

Current coverage: stepcontrol_test.jl 4.1.1/4.1.2 (exp coefs order 30 eps 1e-12 -> 4.501206370338986, 3-source consensus TI.jl==mpmath==Mathematica 47-digit, round-trip vs TI.stepsize). STRONG published-number on PRIMARY formula only. Mutation A (1/(k+1) exponent) REDs it.

Gaps: DEAD CODE IN SCALAR PIPELINE - never called by any scalar stepper/solver (grep-confirmed: no src consumer of StepControl). _second_stepsize FALLBACK UNTESTED. all-zero ArgumentError UNTESTED. length<2 ArgumentError UNTESTED. eps_rel!=eps_abs relative-mode branch UNTESTED. BigFloat/Complex coefs UNTESTED.

Ground truth needed: TI.jl stepsize on an odd-parity (sin-like) jet that triggers _second_stepsize - pin vs TaylorIntegration. A case with eps_rel*|c0|>eps_abs for relative mode. Cheap TI.jl-oracle additions.

**step_pade_root** (internal-algorithm) — `src/StepControl.jl:240`

FW 2011 3.1 pole-distance heuristic: roots of Polynomial(P.b) projected onto unit step direction t=Re((r-z)*conj(unit)); keep t>0; cap at min(forward roots, |target-z|). Constant denom -> full distance. Polynomials.jl roots - Float64/Complex{Float64} only (Arb deferred).

Regimes: constant denominator (no poles); target==z degenerate -> 0; real-axis forward pole; complex-conjugate poles projected to real axis (3+/-i -> 3); roots behind/sideways discarded; no forward roots -> full distance; Float64 & Complex{Float64} only

Current coverage: stepcontrol_test.jl 4.1.3 (1/(1-z/2) pole at 2 -> 2.0), 4.1.4 (1/(1+(z-3)^2) poles 3+/-i -> 3.0 real-projection). STRONG closed-form (Polynomials.roots oracle). Mutation C (real->imag) REDs both.

Gaps: DEAD CODE IN SCALAR PIPELINE - never called by any scalar stepper/solver (reused only conceptually by vector VectorWedgeStep, re-implemented). target==z degenerate->0 UNTESTED. complex z_current / off-axis direction UNTESTED (both tests use 0->5 real). no-forward-roots->full-distance UNTESTED (both cases HAVE a forward pole). constant-denominator early-return UNTESTED.

Ground truth needed: A complex-direction step (z and target off the real axis) with poles on either side vs hand-computed projection. A case with all poles behind to exercise full-distance return. A constant-denominator (nu=0) P.

**PadeStepperState / pade_step! / pade_step_with_pade!** (public-api) — `src/PadeStepper.jl:227 (state), :267 (pade_step!), :300 (pade_step_with_pade!)`

One Pade-Taylor step of u''=f(z,u,u') via FW 2011 5-step recipe: order Taylor jet -> rescale c~_k=h^k*c_k (FW 3.2 line 396) -> diagonal (order/2,order/2) robust_pade -> new_u=P_u(1), new_up=P_u'(1)/h (analytic quotient-rule derivative, NOT re-Pade). Complex-h aware (wedge steps).

Regimes: Float64 (Weierstrass & Painleve-I); ComplexF64 (wedge h=h*e^{i theta}); BigFloat-256 (transitive via solve_pade); near-pole large coefficients (|u| 100->400); odd order (m=n biased down); order<2 ArgumentError; denominator-vanishes-at-t=1 DomainError; step composition

Current coverage: padestepper_test.jl 5.1.1-5.1.4. STRONG: closed-form WeierstrassP + NDSolve WP=50 + mpmath.odefun oracle. 5.1.1 step h=0.5 rtol 1e-13; 5.1.2 composition->0.9 1e-12; 5.1.3 Painleve-I +z 1e-13; 5.1.4 near-pole 0.9->0.95 1e-11. Mutations A (eval at -1) and C (*h vs /h chain rule) RED all four. Also ComplexF64 in AS.1.6.

Gaps: order<2 ArgumentError UNTESTED at pade_step level. DomainError when local-Pade denominator vanishes exactly at t=1 UNTESTED (Rule-1 pole-hit branch). BigFloat pade_step! never DIRECTLY tested (only via solve_pade 6.1.6). Odd-order (order/2 down-bias) UNTESTED - all tests use even order 30/40. PadeStepperState Int->Float64 coercion untested directly.

Ground truth needed: A crafted step landing exactly on a denominator root to confirm DomainError fires vs silent NaN. An odd-order step vs even-order to confirm graceful +/-1 degree behaviour. Direct BigFloat-256 pade_step! vs WeierstrassP oracle.

**_evaluate_pade / _evaluate_pade_deriv** (oracle-helper) — `src/PadeStepper.jl:367 (eval), :397 (deriv)`

Rational + rational-derivative evaluation; the dense-output primitive reused by Problems.PadeTaylorSolution for u(z)=P_k(t) and u'(z)=P_k'(t)/h_k. Fail-fast pole detection (Rule 1).

Regimes: interior t in [0,1] dense eval; endpoint t=1 step eval; Float64/ComplexF64/BigFloat; D(z)=0 pole-hit DomainError

Current coverage: INDIRECT only - exercised inside every pade_step! (5.1.x) and solve_pade dense eval (6.1.x). The quotient-rule derivative is validated transitively by the up assertions (mutation C confirms /h chain rule). No standalone test.

Gaps: D(z)=0 DomainError branch in BOTH functions UNTESTED. No standalone unit test (always composed). The 4-sweep derivative is never checked against an independent symbolic derivative of a known rational - correctness rests entirely on WeierstrassP up-values matching.

Ground truth needed: Standalone: evaluate a known rational (1/(1-z/2)) and its analytic derivative at several z vs closed form; a z hitting a denominator root to confirm DomainError. Cheap closed-form pins currently absent.

**ffw_truncation_error / ffw_rescale_q / adaptive_pade_step!** (internal-algorithm) — `src/PadeStepper.jl:447 (T_h), :505 (q), :553 (adaptive_pade_step!)`

FFW 2017 2.1.2 adaptive Pade step: T(h)=|eps~_{n+1}/a_rescaled(1)| with eps~_{n+1}=c~_{n+1}+sum_{r=1..nu} b_r c~_{n+1-r}; q=(k*Tol/T_h)^(1/(n+1)); iterate T(h)>Tol -> h:=q*h until <=Tol or max_rescales (throws). Complex-h aware (q real in (0,1] preserves wedge direction). Controller memory via meta.h_used.

Regimes: Float64 & ComplexF64 (wedge); T_h=0 -> q=Inf defensive branch; max_rescales exhaustion -> ErrorException; single-rescale (h=2) vs multi-rescale (h=10); a_rescaled(1)=0 -> DomainError; order<2 ArgumentError; Tol/k/max_rescales validation throws; nu_eff<nu trimmed-denominator sum

Current coverage: adaptive_step_test.jl AS.1.1-AS.1.7. STRONG: AS.1.1 h^{n+1} scaling + prefactor pin (closed-form exp); AS.1.2 q-convergence <=5 iters; AS.1.3 adaptive vs fixed-h end-to-end <=1e-9 at z=30; AS.1.4 PIII transform agreement; AS.1.5 Tol-sweep h-monotonicity; AS.1.6 single-step contract; AS.1.7 multi-rescale at h=10. Mutations M1-M4 all bite.

Gaps: a_rescaled(1)=0 DomainError UNTESTED. max_rescales exhaustion ErrorException UNTESTED (no case forces it). ffw_rescale_q validation throws (Tol<=0,k<=0,order<1) and T_h=0->Inf branch UNTESTED. BigFloat adaptive step UNTESTED. nu_eff<nu trimmed-denominator path not explicitly targeted. order<2 ArgumentError in ffw_truncation_error untested.

Ground truth needed: A pathological RHS/h that drives T(h) un-shrinkable to force max_rescales (confirm fail-loud). Direct ffw_truncation_error vs an independently-derived FFW eps_{n+1}*h^{n+1}/a(h) (mpmath) on a non-exp function. BigFloat adaptive vs reference.

**PadeTaylorProblem** (problem-builder) — `src/Problems.jl:108 (struct), :115 (constructor)`

Analytic IVP container (ADR-0001 driver layer). Selects 1st/2nd-order branch by y0 type; promotes element type from zspan; order=FW 2011 5.1 default 30.

Regimes: 2nd-order (Tuple y0); 1st-order (scalar y0 - container only; solve_pade rejects it); Float64 / BigFloat / Complex zspan; order<2 ArgumentError; coincident zspan endpoints ArgumentError; mixed-type literal coercion

Current coverage: problems_test.jl 6.1.1-6.1.6 construct 2nd-order problems (Float64 + BigFloat-256). Widely reused downstream (painleve_*, path_network). Construction exercised end-to-end.

Gaps: order<2 ArgumentError UNTESTED. Coincident-endpoint (z_start==z_end) ArgumentError UNTESTED. 1st-order (scalar y0) construction path UNTESTED at this slice. promote_type with mismatched endpoint types untested.

Ground truth needed: Cheap throw-tests for order<2 and coincident endpoints; a scalar-y0 construction asserting the Y-type branch.

**solve_pade / PadeTaylorSolution** (public-api) — `src/Problems.jl:172 (solve_pade), :217 (callable), :153 (solution struct)`

Fixed-h_max Pade-Taylor IVP driver: iterates pade_step_with_pade! while state.z<z_end, stores per-segment (z,y,h,pade) for dense interpolation. The headline pole-bridge demo (one Pade spanning a lattice pole). Rejects 1st-order branch (v1).

Regimes: single-segment (h_max>=window - all current tests); multi-segment (state.z<z_end loop, >1 push); Float64 / BigFloat-256; h_max<=0 ArgumentError; 1st-order branch error; max_steps exhaustion error; dense eval interior + past-pole + out-of-range DomainError; multi-segment k-scan in callable

Current coverage: problems_test.jl 6.1.1-6.1.6. STRONG: 3-source pole-bridge oracle (closed-form Weierstrass == NDSolve == mpmath). 6.1.1 sol(0.5) 1e-13; 6.1.2 0.95 1e-9; 6.1.3 sol(1.05) bridges pole 1e-9 (headline); 6.1.4 1.4 1e-7; 6.1.5 Pade-vs-Taylor divergence; 6.1.6 BigFloat-256 order 40 1e-13. Mutations A and B documented.

Gaps: ALL tests SINGLE-SEGMENT (h_max=1.5>=window) -> the multi-segment while-loop, multi-push of h_vec/pade_vec, and the callable k-scan (while k<length(h)) NEVER exercised in this slice. h_max<=0 ArgumentError UNTESTED. 1st-order (scalar y0) rejection error UNTESTED. max_steps exhaustion error UNTESTED. Out-of-range DomainError (z<z_start, z>z_end) in callable UNTESTED. (Multi-segment IS exercised by PathNetwork/vector slices, not by solve_pade itself.)

Ground truth needed: A multi-segment solve (h_max<window, e.g. 0.5 over [0,1.5]) vs the Weierstrass oracle at several z, to exercise the segment loop + k-scan. Cheap throw-tests for h_max<=0, scalar-y0, max_steps=1 exhaustion, out-of-range dense eval.

**taylor_eval** (oracle-helper) — `src/Problems.jl:252`

Plain truncated-Taylor Horner evaluator - the foil for the headline Pade-vs-Taylor demonstration (test 6.1.5: same coefficients, Taylor diverges past the radius of convergence while Pade continues).

Regimes: Float64; h beyond radius of convergence (divergence demonstration)

Current coverage: problems_test.jl 6.1.5 (assert |u_taylor-u_ref|/|u_ref| > 0.1 at z=1.05 past the pole). Asserts DIVERGENCE, not a precise value.

Gaps: No assertion that taylor_eval matches the Taylor partial sum to a reference WITHIN the radius of convergence (only the divergence direction is checked). BigFloat/Complex untested. A bug making taylor_eval wrong-but-still-large would pass 6.1.5.

Ground truth needed: A within-radius assertion: taylor_eval(exp coefs, 0.5) ~ exp(0.5) partial sum to a closed-form value - pins correctness, not just divergence.

**Slice notes:** CROSS-CUTTING / SUSPECTED WEAK SPOTS (feeds a correctness-bug hunt):

1. DEAD-IN-SCALAR STEP CONTROL. StepControl.step_jorba_zou and step_pade_root (src/StepControl.jl) are NEVER called by the scalar pipeline. Confirmed by grep: no `using ..StepControl` consumer in src/; PadeStepper.jl and Problems.jl do not import StepControl at all. solve_pade uses FIXED h_max only; adaptive uses the separate FFW controller (PadeStepper.adaptive_pade_step!). Both StepControl functions are re-implemented from scratch in the vector pipeline (VectorStepControl / VectorWedgeStep). A correctness bug in step_pade_root's forward-projection or in step_jorba_zou's _second_stepsize fallback would NOT surface in any scalar integration test - only the 2 isolated unit tests in stepcontrol_test.jl guard them, and several branches (fallback, all-zero throw, no-forward-roots, constant-denom, complex direction) are unexercised.

2. BIGFLOAT 2ND-ORDER JET UNTESTED at coefficient level. taylor_coefficients_2nd is only ever directly oracle-checked at Float64 (3.1.3). The BigFloat 2nd-order recurrence is exercised only end-to-end (problems 6.1.6, value not coefficients). A precision-driven sign/off-by-one specific to BigFloat could hide. The up-resync invariant up[j-1]=j*u[j] is only indirectly observed via u.

3. solve_pade FAIL-FAST GUARDS + MULTI-SEGMENT UNTESTED. Every problems_test case is SINGLE-SEGMENT (h_max=1.5 >= window). The multi-segment while-loop, the callable's k-scan, the 1st-order-branch error, h_max<=0, max_steps exhaustion, and out-of-range DomainError are all unexercised for solve_pade itself.

4. POLE-HIT DomainError PATHS UNTESTED. _evaluate_pade, _evaluate_pade_deriv, and ffw_truncation_error all throw on a vanishing denominator/numerator; none of these Rule-1 branches is exercised. Compliance asserted in docstrings only.

5. robust_pade EDGE BRANCHES. Negative-m/n ArgumentErrors and the _trim_and_normalise error()/force-a=[0] branches are untested. default_tol is never directly unit-tested; the entire Float32 path (default_tol, dispatch, classical) is unexercised - Float32 appears nowhere in the corpus.

ORACLE INFRA FOUND (strong, reusable): test/_oracles.jl (Octave padeapprox.m commit 7574c77 - robust_pade :svd closed-form + reference mu/nu); test/_oracle_coefficients.jl (wolframscript Series + sympy/mpmath); test/_oracle_stepcontrol.jl (3-source TI.jl==mpmath==Mathematica 47-digit pin); test/_oracle_padestepper.jl (Mathematica WeierstrassP + NDSolve WP=50 + mpmath.odefun 40dps); test/_oracle_problems.jl (3-source pole-bridge). All marked do-NOT-hand-edit, regenerated by external/probes/*/verify.jl with cross-source agreement enforcement - strong substrate for adding BigFloat-coefficient and pole-hit oracles.

ARCHITECTURE (ADR-0001 four-layer): scalar core = LinAlg -> RobustPade -> Coefficients -> StepControl(orphaned in scalar) -> PadeStepper -> Problems. Public scalar exports from umbrella: PadeTaylorProblem, solve_pade, PadeTaylorSolution, taylor_eval, robust_pade, PadeApproximant, taylor_coefficients_1st, taylor_coefficients_2nd, PadeTaylorAlg. StepControl functions are NOT re-exported (reachable only via PadeTaylor.StepControl.*). ADR-0002 = BigFloat SVD via GenericLinearAlgebra one-sided Jacobi (relative-accuracy). ADR-0003 = Arb/CommonSolve weakdep extensions. ADR-0005 (referenced) = element-type-driven :classical/:svd default.</parameter>
</invoke>


### sharedpade

**build_square_cell (cell B)** (public-api) — `src/SharedPadeCellB.jl:109`

Wide-square Mano-Tsuda cell B at m_eff=ceil(m/d)*d, +1 window, degree-(m_eff-1) numerators; nothing->cell A on degeneracy; never throws. ADR-0028.

Regimes: d=2 pole; d=1->robust_pade(jet,m-1,m;:svd); d-not-div-m; d==m boundary; BigFloat/Complex; nothing on d>m/short/cnorm<=tol/post-QR/collapse

Current coverage: SPD.1 WEAK (1 d=2 rational, d=1 cross-check, mu_B<mu_A; guards 3a/3b). No mutation-proof/determinant/Calgo/AAA/d>=3/m=4/BigFloat/Rational

Gaps: CRITICAL: +1-window m_eff index no transcendental discriminator (off-by-one survives exact-rational SPD.1); d==m boundary (guard d>m) never tested; cnorm<=tol(126)/post-QR(145)/collapse(155) untested

Ground truth needed: determinant(m_eff)+AAA at d=2,3,4,m=4,BigFloat,Rational; transcendental d=1 discriminator; 4-mutation proof

### pathnet-pole

**path_network_solve** (public-api) — `src/PathNetwork.jl:375`

FW 2011 §3.1 whole-complex-plane path-network IVP solver. Stage 1: builds a tree of Pade-Taylor segments rooted at the IC, walking 5 candidate complex wedge directions per step and selecting min-|u| (heuristic 'furthest from next pole') to bridge poles by off-axis detours. Stage 2: dense-output fine-grid fill — for each grid cell, evaluate the nearest visited node's canonical real-h Pade at t=(z_f-z_v)/h_v inside |t|≤1. Ref: src/PathNetwork.jl:19-50; FW2011 md:155-166; ADR-0004.

Regimes: Float64 and BigFloat-256 (Complex element type forced); pole-crossing wedge steps vs regular steps; fixed-h (FW) vs :adaptive_ffw controller (FFW 2017 §2.1.2); uniform vs node_separation R(ζ) non-uniform Stage-1 nodes; branch_points refuse-mode vs cross_branch multi-sheet; enforce_real_axis_symmetry Schwarz reflection (real-coeff f + real IC); extrapolate=true skip-disc vs default fail-soft NaN; diagnose=true loop-closure certificate (needs DelaunayTriangulation)

Current coverage: test/pathnetwork_test.jl (PN.1.1-PN.7.1): closed-form ℘ oracle for Stage1/Stage2, published FW Table 5.1 numbers u(30) BF-256 ≤1e-13 and u(10⁴) F64 ≤5e-10, tree-structure invariants, enforce_real_axis_symmetry bit-exact + fail-fast, order-default. test/non_uniform_nodes_test.jl (NU.1.1-1.7): node_separation. test/adaptive_step_test.jl (AS.1.3-1.5): :adaptive_ffw end-to-end. test/path_network_branch_test.jl, sheet_aware_stage2_test.jl, extrapolate_test.jl, diagnose_test.jl. STRENGTH: closed-form + published-number on the ℘ problem; mutation-proofed (footers).

Gaps: :steepest_descent never combined with branch_points/node_separation/:adaptive_ffw/symmetry. enforce_real_axis_symmetry never combined with node_separation or :adaptive_ffw (though threaded through _solve_with_schwarz_reflection). grid_sheet+extrapolate not jointly tested. Only ℘ (u''=6u²) and PIII-transformed have any quantitative pin; actual Painlevé transcendents (PI/PII/PIV/PVI) only structural. R returning Complex (vs Real) not tested (would throw wrong error type). verbose path untested.

Ground truth needed: High-precision Painlevé pole/value tables from Mathematica (NDSolve + WeierstrassP) or mpmath for PI/PII/PIV/PVI at fixed z, to pin pole-field accuracy beyond the single ℘ closed-form. Closed-form already exists for ℘ (the lattice formula); extend the published-number coverage to a second equation.

**eval_at** (public-api) — `eval_at(sol::PathNetworkSolution{T}, z::Complex; extrapolate::Bool=false) -> (u, up)`

Sheet-blind per-point Stage-2 dense-output accessor: nearest visited node's local Pade evaluated at t=(z-z_v)/h_v; u'=deriv/h_v chain rule. Fail-soft NaN+NaN·im beyond disc unless extrapolate=true (FFW md:62). Ref: src/PathNetwork.jl:1101; ADR-0015.

Regimes: inside-disc vs beyond-disc (|t|>1); extrapolate true/false; lexicographic tiebreak via _nearest_visited; Float64 (tested); BigFloat untested

Current coverage: test/extrapolate_test.jl (SX.1.1-1.3, SX.1.6): default fail-soft NaN far cell, extrapolate=true finite, no-op inside disc, sparse-walker+dense-eval. Mutation X2 bites (6 RED). STRENGTH: smoke + structural (NaN/finite); value not pinned for extrapolated cells (intentionally).

Gaps: No BigFloat coverage. The _nearest_visited lexicographic tiebreak it relies on is never unit-tested in isolation. No value-accuracy assertion against a closed form for in-disc eval_at (only equality with default path).

Ground truth needed: Closed-form ℘ value at an off-node query point to pin eval_at accuracy directly (currently only the grid path is pinned against ℘).

**eval_at_sheet** (public-api) — `eval_at_sheet(sol::PathNetworkSolution{T}, z::Complex, sheet::AbstractVector{<:Integer}; extrapolate::Bool=false) -> (u, up)`

Sheet-restricted per-point Stage-2 accessor: nearest visited node whose visited_sheet[k]==sheet, then same t=(z-z_v)/h_v Pade eval. Returns NaN if no matching-sheet node within disc. For multi-sheet PVI heatmap rendering. Ref: src/PathNetwork.jl:1041; ADR-0013 A5, bead padetaylor-hed.

Regimes: branched vs unbranched (empty sheet) solutions; sheet match vs absent sheet (NaN); extrapolate true/false; shape-validation throws on wrong-length sheet

Current coverage: test/sheet_aware_stage2_test.jl (SA.1.5 matches vectorised grid_sheet path; SA.1.6 input validation). test/extrapolate_test.jl (SX.1.5 kwarg honoured, Mutation X3 bites). STRENGTH: cross-checks against grid_sheet path + analytic exp(z) when finite (SA.1.4).

Gaps: Value accuracy only checked against exp(z) for the entire (u'=u) problem, never a genuinely multi-sheet function. extrapolate+sheet-restricted not value-pinned. No BigFloat.

Ground truth needed: A genuinely multi-valued reference (e.g. analytic continuation of a log/sqrt across a known cut) with known per-sheet values to confirm sheet selection returns the correct branch value, not just any finite value.

**extract_poles (PathNetworkSolution & PadeTaylorSolution)** (public-api) — `src/PoleField.jl:269 and :308`

Pole-location extraction (FW Fig 4.7/4.8): poles of the solution = denominator roots t* of each node's local Pade, mapped to z=z_ctr+h·t*. Filters: far-root (z-plane |h·t*|≤radius_t·h_max, scale-covariant §S7), Froissart-doublet (residue |N(t*)/D'(t*)|≥min_residue), and cross-node support (≥min_support distinct nodes). Greedy cluster in increasing z-plane distance; representative = z-closest node. Ref: src/PoleField.jl:140-226 (_extract_poles_core); ADR-0026 §S7.

Regimes: path-network fan (min_support=3) vs single-trajectory chain (min_support=1); uniform-h vs varying/adaptive-h walk (§S7 scale-covariance); second-order pole double-root collapse; constant-denominator nodes (no poles, skipped); Froissart-doublet rejection; Float64 (all tests); BigFloat untested

Current coverage: test/polefield_test.jl (PF.1.1-PF.4.1): CLOSED-FORM ℘ lattice oracle (z(m,n)), spurious-pole rejection in-box ≤1e-6, completeness, conjugate symmetry, support-filter load-bearing (PF.2.2), degenerate constant-denom (PF.2.1), single-trajectory PadeTaylorSolution (PF.3.*), §S7 varying-h scale-covariance (PF.4.1). Mutations M1-M6 all bite. STRENGTH: closed-form oracle — strongest in slice.

Gaps: min_residue/Froissart filter NOT directly mutation-proofed (no constructed-doublet test confirming the residue drop branch fires). cluster_atol second-order-pole collapse described in docstring but never pinned with a 'double pole → exactly 1 cluster' test. No BigFloat extraction. Only ℘ closed-form; no Painlevé pole table. radius_t default 5.0 tuning untested across regimes.

Ground truth needed: A meromorphic test function with a KNOWN Froissart-prone Pade (or a hand-built PadeApproximant with a deliberate near-coincident pole-zero) to pin the min_residue drop. Mathematica/mpmath high-precision pole tables for a Painlevé transcendent to extend beyond the ℘ closed-form.

**laplacian_residual** (public-api) — `laplacian_residual(u_grid::AbstractMatrix{Complex{T}}, h::Real) -> Matrix{Complex{T}}  (boundary cells NaN+NaN·im)`

FW 2011 eq (3.3) 5-point cross-stencil Laplacian Δu≈([1;1,-4,1;1]·u)/h², O(h²) approx to ∇²u. Zero for harmonic (meromorphic-away-from-poles) u, large near poles. Ref: src/EdgeDetector.jl:172; FW2011 md:202-208.

Regimes: harmonic quadratic (exact-zero to roundoff); complex-analytic z² (both parts harmonic); near-pole 1/(z-z₀) large residual; boundary cells NaN; grid <3 / h≤0 fail-fast; Float64 (tested); element type T generic

Current coverage: test/edge_detector_test.jl (ED.1.1-1.3, ED.3.1): closed-form harmonic-zero (<1e-12), z² analytic, 1/(z-z₀) >100 near pole + corner ratio, fail-fast guards. Mutations A (−4→−3) and B (drop /h²) bite. STRENGTH: closed-form analytic (harmonicity is exact ground truth).

Gaps: No BigFloat residual test (T is generic but only Float64 exercised). NaN/Inf propagation from upstream u_grid documented as IEEE-passthrough but not asserted. Non-square grids only lightly exercised.

Ground truth needed: Already has exact ground truth (∇²(harmonic)=0). Could add a known-|∇²(1/z)|=2/|z|³ quantitative magnitude pin (currently only a ratio is checked, not the absolute 1/|z|³ scaling).

**pole_field_mask (+ _auto_level)** (public-api) — `src/EdgeDetector.jl:213, _auto_level:156`

Boolean pole-field classifier: mask[i,j] = log₁₀|Δu| > level (FW Fig 3.3 contour). Default :auto resolves to h-aware LEVEL0+2·log₁₀(min(h,H0)/H0), anchor (H0,LEVEL0)=(0.25,0.001), clamped at H0 — corrects the h² residual scaling so the flood-fill annulus survives on fine grids. Ref: src/EdgeDetector.jl:213-260; ADR-0009; bead padetaylor-f8l.

Regimes: smooth exp(z) all-false vs near-pole cluster; h-auto scaling across h∈{0.25,0.125,0.0625}; clamp at H0 for h≥0.25 (backward compat); 1-arg precomputed-Δu form (:auto unresolvable → throws); unknown Symbol fail-loud; BigFloat _auto_level threading

Current coverage: test/edge_detector_test.jl (ED.2.1-2.2, ED.4.1, ED.5.1-5.3): boundary-false, smooth/pole separation at FW level 0.001, two masks via precomputed residual, _auto_level anchor+clamp+sign+BigFloat, annulus-floor tracking across h, :auto propagation + fail-loud. Mutations M1-M4 bite (M3 cascades into EG.1.x proving clamp load-bearing). STRENGTH: closed-form + published FW level + calibration-anchor pins.

Gaps: The empirical calibration anchor (0.25,0.001) is taken as ground truth from a bead probe, not independently re-derived in a test. No test on a real path_network_solve output grid (only synthetic 1/z and exp grids). Auto-level at BigFloat only checks _auto_level, not full pole_field_mask at BigFloat.

Ground truth needed: Independent re-derivation of the h² residual scaling law on a known meromorphic function to confirm the 2·log₁₀(h/H0) slope (currently anchored on a single bead measurement).

**edge_gated_pole_field_solve** (public-api) — `src/EdgeGatedSolve.jl:273`

FW md:401 edge-gated IVP: confines the path-network to the pole field by region growing (seed from IC → dilate by grow_rings → solve → pole_field_mask → morphological-open → flood-fill connected to field → admit) to a fixpoint, then a final solve over field+1 ring. Prevents IVP error in unstable smooth sectors from blooming spurious poles (PI tritronquée pole-free sectors). Ref: src/EdgeGatedSolve.jl:273; FW2011 md:401; ADR-0004.

Regimes: PI tritronquée sector-confinement (the one tested problem); uniform isotropic lattice (Δx≈Δy enforced); coarse vs fine grid (resolution ≲1 required, caller's responsibility); grow_rings≥2, open_radius≥0, seed_radius default vs explicit; fixpoint vs max_iter exhaustion

Current coverage: test/edge_gated_solve_test.jl (EG.1.1-2.1): ONE problem (PI tritronquée FW eq 4.1) on ONE grid (33×33, [-12,12]², h=0.5). Sector-confinement counts (npop≥20, nfree≤5, nfree<0.1·npop), contrast vs plain solve (EG.1.3 leak), fail-fast validation. Mutations M1-M3 bite. Also test/fw_fig_41_test.jl:189 uses it. STRENGTH: structural/count-based (no closed-form pole oracle for PI tritronquée).

Gaps: Single problem, single resolution. open_radius=0 path, seed_radius explicit, max_iter exhaustion, bounded-smooth-band geometry all untested. The doomed-coarse-grid case (spacing≳1) documented as caller responsibility but never tested to confirm graceful degradation vs silent lie. No assertion that recovered pole POSITIONS are correct (only sector counts) — a systematic position bias would pass.

Ground truth needed: High-precision PI tritronquée pole positions (Mathematica NDSolve to high precision, or a published Boutroux-tritronquée pole table) to pin extracted pole accuracy, not just sector membership counts. A second problem with bounded smooth bands to exercise the non-tritronquée geometry.

**_dilate / _erode / _open / _flood_fill (morphology)** (internal-algorithm) — `src/EdgeGatedSolve.jl:140`

8-connected (Chebyshev) morphological operators powering the edge-gated region-growing gate: dilation grows the frontier, erosion+dilation (opening) removes false-positive specks and bridges thinner than 2r, flood-fill keeps only mask cells reachable from the seed (excludes stranded false positives). Ref: src/EdgeGatedSolve.jl:140-213.

Regimes: grid-edge windows (erode: a window hanging off the lattice erodes away); 8-connected vs 4-connected; r=0 (open disabled) vs r≥1; disconnected components (flood-fill excludes stranded); thin bridges (opening removes)

Current coverage: NONE direct. Only exercised transitively through edge_gated_pole_field_solve on the single tritronquée fixture, and referenced in mutation notes M2/M3 (which mutate _dilate). STRENGTH: smoke-only (no isolated unit test asserting a known input→output morphology result).

Gaps: ENTIRE module has no unit tests: erosion grid-edge handling (lines 172-174 treat off-lattice windows as false), the opening idempotence, flood-fill 8-connectivity, r=0 short-circuit, multi-component separation. A correctness bug in erosion boundary handling or a 4-vs-8-connectivity error would be invisible to the suite.

Ground truth needed: Closed-form/hand-built small BitMatrix fixtures with known morphological results (textbook erosion/dilation/opening/connected-component outputs) — trivial to construct, currently absent.

**_select_candidate** (internal-algorithm) — `src/PathNetwork.jl:987`

Wedge-direction selection. :min_u (FW 2011 default) = argmin|u_new| over 5 candidates. :steepest_descent (FW §5.4.1) = wedge edge closest on S¹ (circular geodesic via rem(θ_sd-off,2π,RoundNearest)) to θ_sd=arg(-u/u'). The S¹ vs linear distinction is load-bearing near the ±π branch cut. Ref: src/PathNetwork.jl:987.

Regimes: :min_u argmin (handles Inf-sentinel failed/forbidden candidates); :steepest_descent S¹ geodesic near ±π cut; up_cur=0 fallback to goal_dir; Float64 and BigFloat element types

Current coverage: test/pathnetwork_test.jl (PN.3.1-3.3): DIRECT unit test with INDEPENDENT brute-force S¹ geodesic oracle (circ_dist), RED-on-linear / GREEN-on-circular scenarios, wrong-fix witness (per-element wrap), BigFloat type-generality, :min_u⇄:steepest_descent parity. Mutation F (sign-flip θ_sd) bites. STRENGTH: independent oracle — strongest internal test in slice.

Gaps: :steepest_descent never tested with the Inf-sentinel filter active (branch_points refuse-mode nullifies u to Inf, but :steepest_descent selects purely on angle ignoring u — it can pick a forbidden candidate that only fails at the downstream isfinite guard). This interaction is a plausible latent correctness bug. :min_u's argmin-over-Inf behaviour when ALL candidates are Inf is only tested via the all-failed throw path.

Ground truth needed: A constructed evals fixture where a forbidden (Inf-u) candidate is the S¹-closest wedge edge, to confirm :steepest_descent + refuse-mode together either skip it or fail loud — not silently step through a cut.

**_nearest_visited / _nearest_visited_on_sheet** (internal-algorithm) — `src/PathNetwork.jl:882`

Nearest-visited-node lookup for Stage-1 walk start and Stage-2 dense-output, with deterministic lexicographic (Re then Im) tiebreak on exact distance ties. The on-sheet twin restricts to nodes whose visited_sheet==sheet, returns 0 (fail-soft) if none. Ref: src/PathNetwork.jl:882 and :905.

Regimes: exact-distance ties (lexicographic Re/Im tiebreak); sheet-restricted pool vs sheet-blind; no-match (returns 0 → NaN Stage-2); Float64 and BigFloat

Current coverage: _nearest_visited_on_sheet: DIRECT unit test test/sheet_aware_stage2_test.jl (SA.1.8) — hand-built 5-node fixture, sheet-match + tiebreak + not-found, Mutation S1 bites. _nearest_visited: only INDIRECT (deterministic-walk equality, ℘ Stage-2). STRENGTH: on-sheet closed/structural; sheet-blind smoke-only on tiebreak.

Gaps: Sheet-BLIND _nearest_visited lexicographic tiebreak has NO isolated unit test — a sign-flip or < vs > in the Re/Im tiebreak comparison would pass all current tests (exact ties are rare on ℘ fixtures). The tiebreak feeds eval_at, Stage-2 fill, and Stage-1 walk-start.

Ground truth needed: A hand-built visited_z fixture with deliberate exact-distance ties (e.g. ±0.5 equidistant from query) and known lexicographic winner — trivial, currently absent for the sheet-blind path.

**quality_diagnose (scalar PathNetworkSolution method)** (public-api) — `ext/PadeTaylorDiagnosticsExt.jl:215`

Loop-closure quality certificate (ADR-0016): Delaunay-triangulate sheet-0 visited nodes, extract non-tree edges, measure per-edge midpoint disagreement ΔP_rel=|P_A(M)-P_B(M)|/(|P_A|+|P_B|+ε), categorise (:well_closed/:noisy/:extrap_driven/:depth_driven via tol thresholds and extrap_max=max(|t_A|,|t_B|)>1), aggregate quantiles+worst-edges+bad-centroid. Tree distance via LCA on parent chains. Ref: src/Diagnostics.jl:170; ext/PadeTaylorDiagnosticsExt.jl:215; promoted from loop-closure-fig1 probe.

Regimes: sheet-0 via ζ-strip predicate (-2π<imag≤2π) when unbranched, or visited_sheet==[0] when branched; sheet≠0 throws (v1 deferral); empty/single-node walk → n_edges=0 NaN quantiles no-throw; Pade blow-up at midpoint → edge skipped; extension-not-loaded → MethodError rethrown as ArgumentError

Current coverage: test/diagnose_test.jl (DG.1-DG.8): diagnose=false leaves nothing, diagnose=true attaches report, tallies partition n_edges, quantile monotonicity, worst-edges sorted, post-hoc==eager, sheet≠0 throws citing padetaylor-8py, n_worst respected. Mutation A (well_closed→noisy) bites DG.3. STRENGTH: structural-invariant + self-consistency on small ℘ fixture; published probe trimodal distribution NOT re-pinned.

Gaps: _build_depths and _tree_path_distance (LCA) have NO isolated unit test on a hand-built tree with known LCA distances. :extrap_driven and :depth_driven category discrimination NEVER asserted to fire (benign fixture is all :well_closed) — only :well_closed is mutation-proofed. bad_centroid populated-case (>tol_bad edges) untested. The missing-extension ArgumentError branch is in-process untestable (documented). Ghost-edge filtering (negative indices) untested directly.

Ground truth needed: A hand-built parent-chain tree with known pairwise LCA distances to unit-test _tree_path_distance. A constructed walk with a deliberately bad in-disc loop closure to force :depth_driven, and one with midpoint past disc to force :extrap_driven, confirming the category split fires correctly.

**Slice notes:** Cross-cutting observations and suspected weak spots for the correctness-bug hunt:

ORACLE INFRA FOUND (strongest in repo for this slice):
- test/polefield_test.jl carries a CLOSED-FORM oracle: equianharmonic Weierstrass-℘ pole lattice z(m,n)=1+2Ω(m+n/2)+iΩ√3·n with Ω=Γ(1/3)³/(2^{13/6}π)≈1.3630340904278908 hard-coded (FW md:297). This is a genuine closed-form ground truth and the gold standard for the whole slice.
- test/pathnetwork_test.jl carries PUBLISHED reference numbers from FW Table 5.1: u(30)=1.0950982559597442 (BF-256 rel-err ≤1e-13) and u(10⁴)=21.02530339471055 (F64 rel-err ≤5e-10). Strong end-to-end pin on the wedge walker + canonical-Padé invariant.
- Every slice file has a documented mutation-proof footer; most load-bearing branches are mutation-verified. Treat these as evidence the named branch is tested, NOT that ALL branches are.

SUSPECTED WEAK SPOTS (where a bug could hide untested):
1. _nearest_visited (PathNetwork.jl:882) — the LEXICOGRAPHIC tiebreak (d==best && Re/Im compare) is exercised only INDIRECTLY (via deterministic-walk equality tests). The sheet-aware twin _nearest_visited_on_sheet HAS a direct unit test (SA.1.8) covering tiebreak, but the sheet-BLIND _nearest_visited tiebreak has no isolated unit test. A sign-flip or < vs > in the tiebreak would likely pass all current tests because exact ties are rare on the ℘ fixtures.
2. Morphological helpers in EdgeGatedSolve.jl (_dilate :140, _erode :163, _open :186, _flood_fill :191) have ZERO direct unit tests — only covered through the full edge_gated_pole_field_solve tritronquée integration and mutation notes. _erode's grid-edge handling (a window hanging off the lattice erodes away, :172-174), the 8-connected vs 4-connected choice, and the opening-nibbles-then-restores logic are all untested in isolation. A correctness bug in erosion boundary handling would be invisible.
3. EdgeGatedSolve is tested on EXACTLY ONE problem (PI tritronquée, FW eq 4.1) at ONE resolution (33×33 over [-12,12]², h=0.5). No coverage of: bounded smooth bands, the seed_radius=nothing default vs explicit, open_radius=0 path, max_iter exhaustion, the final-frontier single-ring solve correctness. The "Grid resolution matters" doomed-grid case (spacing≳1) is documented as caller's responsibility but never tested to confirm it degrades rather than silently lies.
4. extract_poles residue/Froissart filter (min_residue, PoleField.jl:196) is NOT directly mutation-proofed in the footer (M1-M3 cover root→z map, radius_t filter, support filter; M5/M6 cover §S7). The Froissart-doublet drop branch has no dedicated test that constructs a known doublet and confirms it is dropped — it is only implicitly relied upon on the ℘ fixture. cluster_atol behaviour (second-order-pole doublet collapse) is described in docstring but never pinned with a test that a double pole yields exactly ONE cluster.
5. eval_at / eval_at_sheet are tested only at Float64 (extrapolate_test.jl). No BigFloat coverage of the per-point accessors, and eval_at's lexicographic-tiebreak path inherits weak-spot #1.
6. :steepest_descent is NEVER tested in combination with branch_points/cross_branch, node_separation, or :adaptive_ffw — only :min_u is exercised on those paths. _select_candidate's :min_u argmin over Inf-sentinels (forbidden/failed candidates) is exercised; the :steepest_descent path's interaction with the Inf-sentinel filter (a forbidden candidate is NOT skipped by steepest_descent's angle-only logic — it can still pick a forbidden wedge index) is a plausible latent bug: _filter_forbidden_candidates nullifies u to Inf but :steepest_descent ignores u entirely and selects purely on θ_sd, so it could select a forbidden candidate and only fail at the downstream isfinite(u_new) guard. NOT tested.
7. enforce_real_axis_symmetry composes with node_separation and :adaptive_ffw in the _solve_with_schwarz_reflection signature (it threads them through) but the PN.6.* tests only exercise fixed-h/:min_u. The interaction is untested.
8. Diagnostics extension (PadeTaylorDiagnosticsExt.jl): _build_depths iterative relaxation and _tree_path_distance LCA are ported from the probe but have NO isolated unit test on a hand-built tree with known LCA distances — only validated via the aggregate DG.* assertions on a small ℘ fixture. The category thresholds (_categorise :207, the :extrap_driven vs :depth_driven split on extrap_max>1) are only mutation-proofed for :well_closed (Mutation A); the :extrap_driven / :depth_driven discrimination is never asserted to fire (the benign fixture lands everything in :well_closed). The "missing extension throws" branch is explicitly UNTESTABLE in-process (documented).
9. node_separation throw-guards (NU.1.6) cover negative/NaN/zero R but the docstring claims R must return positive FINITE — an Inf-returning R is tested via the isfinite check, OK. But R returning a Complex (the signature is R::Complex{T}->T, coerced via T(r)) — if R returns a complex value the T(r) coercion would throw an InexactError/MethodError, not the documented ArgumentError; not tested.
10. grid_sheet Stage-2 with extrapolate=true is not jointly tested (sheet-restricted + skip-disc-check). The two kwargs are tested separately.

GROUND-TRUTH GAPS: the whole slice leans on ONE closed-form (℘) and ONE published table (FW 5.1). The Painlevé transcendents themselves (PI/PII/PIV/PVI) have NO closed-form pole oracle here; pole-field correctness on actual Painlevé problems is only smoke/structural (sector-confinement counts), never against a high-precision reference pole table from Mathematica/mpmath. That is the biggest ground-truth hole for a correctness hunt.

### bvp-dispatch

**bvp_solve (2-arg RHS)** (public-api) — `src/BVP.jl:239`

Chebyshev-Newton spectral collocation for 2nd-order scalar analytic BVP u''=f(z,u) with Dirichlet BCs on a complex segment. D2=D1*D1, interior residual R=(D2 u)_int + bc_col - scale*f, diagonal analytic Jacobian J=D2_ii - scale*diag(df/du), scale=(z_b-z_a)^2/4 affine factor. FW 2011 sec 3.2 (FW2011 md:176-200) + Trefethen SMIM ch.6,13; references/bvp_recipe.md.

Regimes: Float64 LAPACK; BigFloat-256 GenericLinearAlgebra (BV.5.1); real segment; purely-imaginary segment (FW Fig 4.1); Complex{T} endpoints via complex scale; linear 1-step vs nonlinear multi-step Newton; N from 4 to 240

Current coverage: test/bvp_test.jl BV.1.2/1.3/1.4 (published-number: DMSUITE+mpmath oracle in _oracle_bvp.jl, atol 1e-10..1e-12), BV.5.1 (closed-form cosh BF-256), test/fw_fig_41_test.jl FF.1.* (published-number FW eq.4.1 to 1e-10), test/dispatcher_test.jl DP.1.2 cross-check. STRONG: oracle-pinned + mutation-proofed (A/B/C).

Gaps: Oblique complex segments (both Re and Im of z_b-z_a nonzero) never oracle-checked; all oracle tests use real or purely-imaginary segments. Singular-Jacobian (J\R failure) path untested. callable initial_guess returning Complex on oblique segment untested.

Ground truth needed: A closed-form complex-oblique BVP: u''=u with z_a=0, z_b=1+1im, BCs from cosh(z), assert node values vs cosh to spectral floor. Or mpmath/Mathematica PI reference on an oblique complex segment.

**bvp_solve (3-arg RHS overload, u'-dependent)** (public-api) — `src/BVP.jl:370`

Generalises the 2-arg path to RHS depending on u'. u'_int via (D1*u_nodes)/half_diff (chain rule). Jacobian J = D2_ii - scale*diag(df/du) - (scale/half_diff)*diag(df/dup)*D1_ii. Needed for P-tilde-III w''=(w')^2/w+... (FFW md:43). bead padetaylor-i76, ADR-0014.

Regimes: Float64 (hybrid is Float64-first); complex slices with purely-imaginary z difference; PIII RHS with 2*wp/w df/dup; per-node asymptotic-series initial_guess

Current coverage: ONLY indirectly via test/ivp_bvp_hybrid_test.jl IB.1.4/1.5 and test/ffw_fig_5_test.jl FF5.1.2-1.5 -- SELF-CONSISTENCY checks (BC recovered, |u| bounded, vs same asymptotic series) at GENEROUS thresholds (5e-2). Mutation M3/M1 bite it. SMOKE-TO-SELF-CONSISTENT, no closed-form oracle.

Gaps: The non-diagonal D1-coupled Jacobian term (scale/half_diff)*diag(df/dup)*D1_ii could be mis-scaled O(1) and survive -- no closed-form isolates it. bc_col_D1 = D1_ib*[u_b;u_a] boundary ordering untested in isolation. initial_guess_up silently ignored is a latent trap.

Ground truth needed: A closed-form 3-arg BVP: linear variable-coefficient u''=a(z)u'+b(z)u (exercises u'-Jacobian but exact), or nonlinear u''=(u')^2/u with solution u=A*e^{kz}. Assert node values to spectral floor.

**BVPSolution callable sol(z)->(u,u')** (public-api) — `src/BVP.jl:486`

Berrut-Trefethen 2004 sec5 Chebyshev-2 barycentric interpolation (w_j=(-1)^j halved at ends) + chain-rule derivative recovery via D1*u_nodes /half_diff. DomainError if |Re(t*)|>1+100eps.

Regimes: interior z; on-node z (coincident guard); endpoint z recovers BC; out-of-segment DomainError; BigFloat eps in guard; complex z with domain check on real(t*) only

Current coverage: test/bvp_test.jl BV.2.1 (chebint oracle, atol 1e-12), BV.3.1 (endpoint derivatives vs oracle, atol 1e-9), BV.4.1 (DomainError z=5.0). Mutation C (drop (-1)^j) bites BV.2.1. Published-number strong for real segments.

Gaps: Domain check uses real(t*) ONLY -- a complex z with |real(t*)|<=1 but large imag PASSES the guard and silently extrapolates off-segment; no test feeds an off-axis complex z to a real-segment solution. Derivative rebuilds D1 every call (perf).

Ground truth needed: Closed-form check whether sol(z) with complex off-axis z is correctly rejected or interpolated -- decide if the real(t*)-only guard is a bug for oblique-complex segments.

**_chebyshev_D1 / _barycentric_eval (BVP internals, duplicated 3x)** (internal-algorithm) — `src/BVP.jl:525`

Trefethen cheb.m / Weideman-Reddy chebdif.m first-derivative matrix (off-diag (c_i/c_j)(-1)^{i+j}/(t_i-t_j), diagonal by negative-row-sum D*1=0). Barycentric 2nd-form. Duplicated VERBATIM into VectorBVP.jl and Laplace2D.jl per ADR-0023.

Regimes: Float64; BigFloat; N=4..240; diagonal negative-sum identity

Current coverage: test/bvp_test.jl BV.1.1 (DMSUITE-pinned full D1/D2 at N=4,8 to 1e-11; partial diag/supdiag at N=16). Re-tested via VectorBVP VB.2.1/2.2 and Laplace2D L2D.8. STRONG published-number.

Gaps: THREE VERBATIM COPIES (BVP, VectorBVP, Laplace2D) with NO test asserting they are byte-identical -- a fix or bug in one will not propagate; they could silently diverge. Only BVP's copy is directly DMSUITE-oracled.

Ground truth needed: A cross-module test asserting BVP._chebyshev_D1 == VectorBVP._chebyshev_D1 == Laplace2D._cheb_D1 elementwise (catches copy drift).

**dispatch_solve (1D IVP<->BVP chain dispatch)** (public-api) — `src/Dispatcher.jl:212`

FW 2011 sec 4.4 (md:249-261) 1D ordered-chain composition of path_network_solve (IVP) and bvp_solve (BVP). Propagates (z,u,u') across junctions; IVP->BVP takes u_a from IVP terminus, BVP->IVP takes IC from barycentric (u,u') at BVP right end. Records junction_match (Du,Du') diagnostic (line 192, 1e-7/1e-8). Strict mode throws on mismatch.

Regimes: single IVP passthrough; single BVP passthrough; IVP->BVP; IVP->BVP->IVP; strict vs non-strict; wedge_angles forwarding; order defaulting to prob.order (padetaylor-9xf)

Current coverage: test/dispatcher_test.jl DP.1.1 (passthrough == direct path_network_solve), DP.1.2 (== direct bvp_solve), DP.1.3 (order default regression), DP.2.1/DP.3.1 (closed-form cosh + junction diagnostics, atol 1e-5..1e-12), DP.4.1 (fail-fast). Mutations A/B/C. Closed-form cosh STRONG for the linear spine.

Gaps: (1) bvp_tol-set arms of the 4-way kwarg fan (Dispatcher.jl:309-326) not directly tested. (2) BVP->IVP junction match HARD-CODED to (0,0) (line 380) -- a real Du' there is invisible by construction; only caught downstream via cosh(0.75) atol=1e-5. (3) No test with BVP segment FIRST (u_a from prob.y0[1]). (4) Nonlinear ODE through a full IVP->BVP->IVP chain never tested -- only linear u''=u and DP.1.1 weierstrass passthrough. (5) wedge_angles!==nothing path untested.

Ground truth needed: A nonlinear closed-form chain (weierstrass-P or PI with known values) through IVP->BVP->IVP. A test where the FIRST segment is a BVP. A fixture perturbing the BVP->IVP handoff derivative to confirm the downstream check catches what the zero-hardcoded junction cannot.

**lattice_dispatch_solve (2D lattice dispatch)** (public-api) — `src/LatticeDispatcher.jl:257`

FW 2011 sec4.4 + md:190 '161 BVPs one per grid line': per-row BVP fill of contiguous smooth runs flanked by pole-field cells. mask=nothing routes IVP through edge_gated_pole_field_solve (ADR-0017 / bead 0tj fix); mask supplied routes plain path_network_solve. strict=false swallows ONLY 'Newton did not converge' and tags :bvp_fail.

Regimes: edge-gated auto path (mask=nothing); manual mask path; mask all-false (no bridging); strict=true throw; strict=false fail-soft :bvp_fail; smooth runs touching boundary :ivp_only; anisotropic-grid rejection; edge_level :auto vs Real; NaN-carrying up_grid scatter

Current coverage: test/lattice_dispatcher_test.jl LD.1.2 (closed-form cosh on BVP-filled cells, atol 1e-10 -- STRONG but only MANUAL mask path), LD.1.1/LD.X.7 (auto path: SIZES + tag-membership + count<=5 only, NO value check), LD.1.3 (no-bridge), LD.2.* (fail-fast), LD.X.1-7 (strict-mode forced-divergence via bvp_tol=1e-300). Mutations E/F/X.

Gaps: (1) The edge-gated AUTO path is NEVER value-checked against a closed form -- all numeric BVP-fill correctness rides on the manual-mask path. (2) up_grid scatter (lines 338-345) assuming _dilate_one matches EdgeGatedSolve._solve_targets indexing has NO direct test that up_grid[i,j] is the right derivative. (3) _dilate_one 8-connectivity vs the gated solver's actual dilation is an untested cross-module assumption. (4) Multiple disjoint smooth runs per row not explicitly tested.

Ground truth needed: A closed-form (cosh) value check on the EDGE-GATED auto path. A test pinning up_grid[i,j] at field cells to the path_network grid_up to confirm scatter ordering. A two-disjoint-runs-per-row fixture.

**_dilate_one (LatticeDispatcher helper)** (oracle-helper) — `src/LatticeDispatcher.jl:447`

8-connected (Chebyshev) one-ring dilation lifting gated.field_mask to the one-ring frontier to index the IVP up_grid scatter, intended to match the final _solve_targets dilation inside edge_gated_pole_field_solve.

Regimes: interior cells; boundary cells clamped; 8-neighborhood

Current coverage: NONE directly. Only exercised inside the auto path of LD.1.1/LD.X.7 which do not value-check up_grid.

Gaps: Entirely unverified that its 8-connected dilation matches the connectivity used inside EdgeGatedSolve. If EdgeGatedSolve dilates 4-connected, the up_grid scatter misaligns at diagonal frontier cells and writes wrong derivative values silently.

Ground truth needed: A direct unit test of _dilate_one on a hand-built BitMatrix, plus a test asserting it reproduces EdgeGatedSolve's internal target ordering (or a refactor to share the helper).

**solve_pole_free_hybrid (pole-free-sector IVP+BVP hybrid)** (public-api) — `src/IVPBVPHybrid.jl:443`

FFW 2017 sec3 (md:203-247) four-step hybrid: (1) two PFS walks along sector boundaries, (2) harvest BCs by linear interp of pfs.grid_u/grid_up, (3) bvp_solve 3-arg on each Re zeta=const slice, (4) glue via sol(zeta) dispatch with glue_tol continuity. Cures kappa_r~(27/16)e^{2 Re zeta/3} ~ 157 IVP ill-conditioning on the pole-free sector. ADR-0014.

Regimes: PIII :transformed frame ONLY (others throw); degenerate_full_plane fast path == pure PFS; Float64 only v1; n_slices>=1; linear-interp between slices in Re zeta; _bvp_solve_on_slice hard-codes PIII analytic df/dw, df/dwp; per-node asymptotic initial_guess via z=exp(zeta/2)

Current coverage: test/ivp_bvp_hybrid_test.jl IB.1.3 (degenerate bit-exact to PFS -- STRONG), IB.1.4 (BC-match self-consistency), IB.1.5 (finite/|w| bounded/Newton-converged SELF-CONSISTENT, thresholds 5.0 not FFW 1e-7), IB.1.7 (fail-fast). test/ffw_fig_5_test.jl FF5.1.2-1.5 (|u| bounded + self-cross-check vs same series, threshold 5e-2). Mutations M1-M4. SELF-CONSISTENT/SMOKE -- no independent oracle on sector interior.

Gaps: (1) Sector-interior solution value NEVER checked against independent ground truth -- only against the SAME asymptotic series used as IC; FFW md:240 per-strip errors explicitly declared unreachable at v1. (2) _harvest_at_re linear interp + clamping (extrapolation r<=r_first / r>=r_last) untested for correctness. (3) A :transformed PV reaching _bvp_solve_on_slice's :III guard not tested (only :direct rejection in IB.1.7). (4) POSSIBLE DOC-VS-CODE DRIFT: docstring/callable promises a BVP-vs-PFS glue-continuity assertion that throws on violation, but the callable code (IVPBVPHybrid.jl:738) appears to only do inside/outside dispatch + slice interp -- no glue check found; the 'asserts agreement to glue_tol / Failure throws' may be stale.

Ground truth needed: An independent high-precision PIII tronquee reference on the sector (BF-256 + many-term series via Mathematica/mpmath) to replace the circular self-check. A direct test that the glue continuity assertion actually fires (feed a discontinuous slice stack, confirm throw) -- currently promised by docstring but possibly not implemented.

**pIII_asymptotic_ic (tronquee asymptotic IC builder)** (problem-builder) — `src/IVPBVPHybrid.jl:189`

FFW md:222/md:230 truncated tronquee asymptotic series u~z^{1/3}[1+sum a_n z^{-2n/3}]. Closed forms a_1=-beta/3, a_2=-(beta^2/9)(1+delta)/(2-delta) [=0 at delta=-1 Fig5 family], a_3+ deferred to zero. Derivative by termwise differentiation. bead padetaylor-qsj (corrected wrong a_2 ~ -0.222).

Regimes: principal sheet arg z in (-3pi/4, pi]; |z|>=1 else throw; delta=-alpha constraint else throw; (alpha,gamma)=(1,0) only; n_terms>=1; ComplexF64 coefficients; Float64-hardcoded coeff computation

Current coverage: test/ivp_bvp_hybrid_test.jl IB.1.1 (published-number FFW eq.6/md:243 u(z2)/u'(z2) at |z|=30, atol 1e-5/1e-6 -- STRONG independent oracle), IB.1.2 (leading u~z^{1/3} at z=1000,1e9), IB.1.7 (fail-fast). test/ffw_fig_5_test.jl FF5.1.1 (md:243 at z1,z2), FF5.1.8 (a_2=0 => n_terms 2..15 bit-identical). Mutations M2/M4. STRONG at |z|=30.

Gaps: (1) Only validated at |z|=30 and asymptotically large |z|; no mid-range |z|=3..10 where truncation error is larger and a_3 matters. (2) The arg z in (pi, 9pi/4) upper sub-range (sheet-unfolding) unreachable by direct call and untested. (3) _pIII_asymptotic_coeffs returns general-delta a_2 but the public entry forces delta=-alpha, so the a_2!=0 branch is never reached in production -- dead-ish with no test. (4) Float64-hardcoded coeffs vs BF-256 untested.

Ground truth needed: An mpmath/Mathematica PIII tronquee value at a mid-range |z| (say 5) to bound truncation. A test of the sheet-unfolding convention if/when the upper sub-range is used.

**IVPBVPSolution callable sol(zeta)->(w,w')** (public-api) — `src/IVPBVPHybrid.jl:738`

zeta-frame dispatch: nearest/bracketing Re zeta slice barycentric + linear-interp in Re. _eval_bvp_slice projects zeta to the slice's Re-line before calling the BVPSolution callable. Throws DomainError outside sector; ErrorException if bvp_slices empty (degenerate).

Regimes: inside sector; on boundary (glue_tol slack); outside DomainError; degenerate empty-slices ErrorException; below first / above last slice (clamp); between slices linear interp

Current coverage: test/ivp_bvp_hybrid_test.jl IB.1.4 (boundary query finite + |w|<1e6), IB.1.5 (finite at interior grid), IB.1.7 (DomainError outside). test/ffw_fig_5_test.jl FF5.1.2-5 (interior zeta values self-consistent). SMOKE-to-self-consistent.

Gaps: The _eval_bvp_slice Re-projection (ignores the query's actual Re, snaps to slice Re-line) introduces an interpolation error never bounded by a test. The docstring's promised boundary glue-assertion is not evidently implemented (see solve_pole_free_hybrid gap). Slice linear-interp accuracy in Re untested vs a known smooth Re-dependence.

Ground truth needed: Same independent sector reference as solve_pole_free_hybrid. A test isolating the Re-projection error magnitude.

**laplace2d_solve (2D Laplace harmonic fill)** (public-api) — `src/Laplace2D.jl:163`

Tensor-product Chebyshev spectral Dirichlet Laplace solve grad^2 phi=0 on a rectangle (Trefethen SMIM Prog.16/23). L=kron(I_y,D2x_ii)+kron(D2y_ii,I_x), D2 affine-scaled (2/(b-a))^2 and (2/(d-c))^2. Dirichlet absorbed into RHS. SINGLE linear solve (no Newton). Voter (2) of the KKG P_I^(2) triple-method harmonic fill (ADR-0024).

Regimes: Float64 LAPACK; BigFloat-192 (L2D.5); square [-1,1]^2; non-square rectangle (L2D.3 affine scale); callable vs vector boundary data; Nx!=Ny; polynomial-exact (deg<=2,3) vs transcendental geometric convergence

Current coverage: test/laplace2d_test.jl L2D.1-6: THREE closed-form harmonic oracles (x^2-y^2, exp(x)cos(y), x^3-3xy^2) to 1e-9..1e-12 (Float64) and 1e-20 (BF-192). L2D.7 fail-fast. L2D.8 standalone mutation-proof (drop affine scale on unequal-width rect; swap kron order on Nx!=Ny). STRONG closed-form oracle coverage.

Gaps: (1) Corner-consistency NOT enforced and NOT tested -- inconsistent/discontinuous boundary data has undefined behavior with no fail-fast guard; the real F3 use (conformal w=log x sector boundary) could hit a corner mismatch. (2) Complex-valued phi never tested (Re u / Im u are solved as separate real problems, likely fine by design). (3) barycentric tensor callable coincident-node guard eps(T) tested at one node only (L2D.4), not at BF corners.

Ground truth needed: A test feeding deliberately-inconsistent corner data to confirm either a fail-fast guard or to document the silent behavior. (Closed-form harmonic oracles are otherwise gold-standard.)

**vector_bvp_solve (first-order vector-system BVP)** (public-api) — `src/VectorBVP.jl:180`

Chebyshev collocation + Newton + tau-method for first-order vector ODE y'=f(z,y), y in C^d, under general linear two-point BC B_a*y(z_a)+B_b*y(z_b)=g. Stack D1 kron I_d, residual node-major, tau-method replaces node-0's d rows with BC. Node Jacobian by Taylor1 autodiff (coeff [1]) or caller-supplied. ADR-0023; for the P_I^(2) tritronquee SEPARATRIX where shooting cannot pin.

Regimes: d=1 reduction; d=2 harmonic-osc / PI companion; d=4 P_I^(2) companion; real and OBLIQUE-complex segments (VB.1.3); Float64 + BF-256 (VB.5.1/6.6); autodiff vs explicit Jacobian; linear vs nonlinear Newton; DMSUITE descending node order (B_a hits node-N=z_a, B_b hits node-0=z_b); singular-Jacobian rethrow

Current coverage: test/vector_bvp_test.jl VB.1.1-1.4 (closed-form e^z and (cos,-sin), incl. OBLIQUE complex segment -- STRONGER than scalar BVP here), VB.2.* (D1 pins), VB.3.1 (autodiff vs analytic Jacobian), VB.4.* (fail-fast), VB.5.1 (BF-256), VB.6.* (cross-validation vs oracle-pinned scalar BVP on nonlinear PI companion). test/vector_bvp_wirein_test.jl VW.* (analytic==autodiff==FD Jacobian triple-check). test/kkg_pi2_figure_test.jl PI2F.* (headline tritronquee, published KKG asymptotic). Mutations A/B/C/D + VW M1. STRONG.

Gaps: (1) GENERAL linear BC never tested -- every B_a/B_b is a Dirichlet selector ([1 0;0 0] etc.); a coupling/mixing or endpoint-coupling B (e.g. [1 1;0 0] or periodic B_a=I,B_b=-I) is the stated rationale of ADR-0023 yet has NO test. The B_a*Y_N+B_b*Y_0 assembly could be wrong for non-selector B and survive. (2) Singular-Jacobian path (_solve catch) untested -- no rank-deficient B fixture. (3) initial_guess wrong-length (_checked_guess) untested (VB.4.1 tests RHS wrong-length only). (4) d>4 / Noumi-Yamada A_n^(1) vector systems (claimed to generalise for free) untested.

Ground truth needed: A closed-form vector BVP with a GENUINELY GENERAL (non-Dirichlet, endpoint-coupling) linear BC and known solution -- e.g. periodic BC on y'=(y2,-y1), or a mixing BC. A rank-deficient-B fixture to exercise the singular-Jacobian throw.

**VectorBVPSolution callable + _residual + _autodiff_jacobian (VectorBVP internals)** (internal-algorithm) — `src/VectorBVP.jl:313`

Component-wise Berrut-Trefethen barycentric vector eval; the single canonical collocation+tau residual (Rule 6 de-dup, Newton + diagnostic share it); exact d x d Jacobian by Taylor1 (coeff [1], ADR-0020 element-type-generic RHS).

Regimes: d components; complex/BigFloat; Taylor1{CT} promotion; coeff [1] read (Mutation C/D [0] breaks it); shared residual prevents Newton/diagnostic drift

Current coverage: VB.3.1 (autodiff Jacobian real+complex vs analytic, atol 1e-13), VW.2 (analytic==autodiff==FD at m=1,m=2 random points), VB.1.* callable vs closed form. Mutations C/D bite autodiff. STRONG.

Gaps: _residual's BC-row line R[bc_rows]=Ba*Y_{nodeN..}+Bb*Y_{1..}-g only exercised with selector B (see general-BC gap). Callable's real(t*)-only domain check shares the BVP.jl oblique-complex concern. autodiff-specific _check_len only via one fixture (VB.4.1 f_bad).

Ground truth needed: Same general-BC fixture as vector_bvp_solve; would exercise _residual's BC assembly non-trivially.

**Slice notes:** CROSS-CUTTING for the correctness-bug hunt on the spectral-BVP / dispatch side.

ORACLE INFRA FOUND (strong): test/_oracle_bvp.jl pins DMSUITE/Octave 8.4.0 D1/D2 matrices (N=4,8,16) + mpmath-40dps PI BVP node values + chebint barycentric evals -- the strongest oracle in the slice, underpinning scalar BVP.jl. Closed-form harmonic oracles (x^2-y^2, exp(x)cos(y), x^3-3xy^2) anchor Laplace2D. cosh/cos/sin/e^z closed forms anchor the dispatchers and VectorBVP. KKG asymptotic +(6x)^{1/3} and FFW md:243 published u(z2) anchor the headline figures.

TOP SUSPECTED BUG-HIDING SPOTS (ranked):

1. THE 3-ARG bvp_solve OVERLOAD HAS NO CLOSED-FORM ORACLE (biggest gap). Its D1-coupled non-diagonal Jacobian term (scale/half_diff)*diag(df/dup)*D1_ii is exercised ONLY through hybrid PIII self-consistency at 5e-2 thresholds. The term could be mis-scaled O(1) and only M3/M1 mutations (drop BC / drop BVP) catch anything. A linear variable-coefficient or u=A*e^{kz} closed-form 3-arg test is missing.

2. POSSIBLE DOC-VS-CODE DRIFT in solve_pole_free_hybrid / IVPBVPSolution callable: the docstrings repeatedly promise a BVP-vs-PFS glue-continuity assertion that throws on violation (glue_tol), but the callable at IVPBVPHybrid.jl:738 appears to only do inside/outside dispatch + slice interpolation -- I found no glue/continuity check or throw in that code path. Either the assertion lives elsewhere or the promised fail-loud guard is stale/unimplemented. Worth a direct check.

3. COMPLEX-SEGMENT scalar bvp_solve is under-oracled: every scalar oracle uses a real or purely-imaginary segment; no OBLIQUE complex segment (Re and Im of z_b-z_a both nonzero) is oracle-checked. The domain guard in both BVP and VectorBVP callables checks real(t*) only, so an off-axis complex z could pass the guard and silently extrapolate. (VectorBVP VB.1.3 does test an oblique complex segment against closed form -- so the vector path is actually better covered here than the scalar path.)

4. LatticeDispatcher up_grid scatter (lines 338-345) + _dilate_one (8-connected) assume they reproduce EdgeGatedSolve._solve_targets internal column-major + dilation ordering, with NO direct test. If the gated solver dilates 4-connected, diagonal frontier derivative values are silently misplaced. The entire edge-gated AUTO path is never value-checked against a closed form -- all numeric correctness rides on the manual-mask path (LD.1.2).

5. VectorBVP only tested with Dirichlet SELECTOR B matrices ([1 0;0 0] etc.) despite ADR-0023's whole rationale being 'general linear two-point BC, not pointwise Dirichlet'. A bug in B_a*Y_N + B_b*Y_0 assembly that manifests only for coupling/mixing/periodic B would survive. Singular-Jacobian throw paths (both scalar and vector) are untested (no rank-deficient fixture).

6. Dispatcher BVP->IVP junction match is hard-coded to (0,0); a genuine derivative mismatch at that handoff is invisible at the junction and only caught downstream by a loose atol=1e-5 cosh check. No BVP-first-segment test; no nonlinear full-chain test; bvp_tol-set kwarg-fan arms untested.

7. THREE verbatim copies of _chebyshev_D1 (BVP/VectorBVP/Laplace2D) with no test asserting they are identical -- copy drift would be silent (only BVP's copy is DMSUITE-oracled directly).

8. Laplace2D does not enforce or test corner-consistency of boundary data (documented as caller's responsibility, no guard) -- a real risk once F3 feeds conformal w=log x sector boundary data.</parameter>
</invoke>


### multivalued-painleve

**pIII_z_to_zeta / pIII_zeta_to_z** (public-api) — `src/CoordTransforms.jl:120 (z_to_zeta), src/CoordTransforms.jl:132 (zeta_to_z)`

PIII fixed-branch-point-at-z=0 exponential coordinate transform (FFW 2017 sec2.1, FFW2017...md:43; derived worklog 017). Maps branch at z=0 to zeta=-inf so P~III is meromorphic. Strip width 4pi = one z-sheet. Principal-branch log; non-principal sheets at zeta+4pi*i*s.

Regimes: Complex z principal log; real positive z; non-principal sheet shift +4pi*i*s (documented, NOT tested); z near 0 branch point (guarded only at PainleveProblem layer); BigFloat/Arb (NEVER exercised)

Current coverage: coord_transforms_test CT.1.1 closed-form (zeta=0,w=1/4,wp=5/8 at z=1, exact + round-trip 1e-14). CT.1.3 end-to-end direct-vs-transformed single Pade step <=1e-10. ffw_fig_1 FF1.1.1 IC round-trip; ffw_fig_5 PIII via PainleveProblem. Strength: closed-form oracle + published FFW Fig1/5.

Gaps: Only z=1 (zeta=0) sample for the closed-form pin. No non-principal-sheet (4pi*i*s) round-trip despite documented API. No BigFloat. Round-trip only on principal branch.

Ground truth needed: Closed-form (algebraic map) at several z incl. complex and non-principal sheet -- complete oracle. BigFloat: same closed form at setprecision(256) to confirm type propagation + no precision loss.

**pIII_transformed_rhs** (public-api) — `src/CoordTransforms.jl:150`

Transformed P~III RHS (FFW 2017 sec2.1, md:43) -- meromorphic-in-w(zeta) equation for the pole-friendly solver. The /4 divisor and four parameter terms are the load-bearing algebra.

Regimes: w near 0 (1/w, d/w blow up -> IEEE Inf + downstream fail-loud per Rule 1); Complex/Taylor1 args; large Re zeta (pole density grows; FFW R(zeta) node adaptation NOT in v1; e^(2zeta) overflow); BigFloat (never tested)

Current coverage: coord_transforms_test CT.1.2 closed-form hand-pin at (0,1/4,5/8) FFW Fig1 params = 111/256 exactly (mutation /4->/3 and z_to_zeta sign-flip bite). CT.1.3 end-to-end. ffw_fig_1 FF1.1.3 sheet-0 absolute pin <=1e-5, FF1.1.4 conjugate-symmetry, FF1.1.5 >=30 poles, FF1.1.6 pole-density gradient. Strength: closed-form + multiple FFW Fig1 published-number pins.

Gaps: Single sample (zeta=0) for the hand-pin. w-near-0 fail-loud asserted only indirectly via stepper. No BigFloat. Large-Re-zeta tail degradation + e^(2zeta) overflow untested.

Ground truth needed: Closed-form RHS at several (zeta,w,wp) incl complex zeta + param sweep. FFW Fig1 captured tight-tol references exist; an independent mpmath/Mathematica PIII value would strengthen beyond self-consistent loose-vs-tight.

**pV_z_to_zeta / pV_zeta_to_z** (public-api) — `src/CoordTransforms.jl:169 (z_to_zeta), src/CoordTransforms.jl:181 (zeta_to_z)`

PV fixed-branch-point-at-z=0 exponential transform (FFW 2017 sec2.1, md:47). Strip width 2pi per z-sheet. Same map serves PVI zeta-plane (Painleve.jl _build_VI passes it for :transformed). u-component identity (no rescaling), unlike PIII.

Regimes: Complex z principal branch; real positive z; PV AND PVI zeta-plane usage (shared map, distinct equations); non-principal sheet +2pi*i*s (documented, untested); BigFloat (untested)

Current coverage: coord_transforms_test CT.1.4 closed-form (zeta=0,w=0.5,wp=0.3 at z=1)+round-trip 1e-14. CT.1.6 end-to-end PV. sheet_tracker_test ST.1.3 reuses for PVI end-to-end. ffw_fig_2 FF2.1.1, ffw_fig_4 FF4.1.1, ffw_fig_6 FF6.1.1, ffw_fig_3/7 round-trips. painleve_solution_test PS.3.2 (poles exp), PS.5.2/PS.5.3 (grid_values exp + up/z). Strength: closed-form + heaviest figure reuse in slice.

Gaps: Non-principal sheet shift untested. No BigFloat. Dual PV/PVI role means a bug corrupts BOTH; no test asserts the map is correct for PVI independently of its PV derivation.

Ground truth needed: Closed-form (algebraic) -- complete. For PVI-reuse: independent derivation check that u(z)=w(zeta), w'=z*u' is the correct PVI zeta-transform (cite md:146) rather than inheriting by assertion.

**pV_transformed_rhs** (public-api) — `src/CoordTransforms.jl:195`

Transformed P~V RHS (FFW 2017 sec2.1, md:47). Meromorphic in w(zeta). Singular at w=0 (b/w,1/(2w)) and w=1 (1/(w-1),(w+1)/(w-1)) -- the two fixed singular surfaces.

Regimes: w near 0 and w near 1 (two singular surfaces, IEEE Inf + downstream fail-loud); Complex/Taylor1 args; large Re zeta tail (FFW node deferral); BigFloat (untested)

Current coverage: coord_transforms_test CT.1.5 closed-form hand-pin (0,1/2,3/10)=-0.1675 exactly (mutation N (w+1)/(w-1) swap bites). CT.1.6 end-to-end. ffw_fig_4 FF4.1.2 sheet-0 <=3e-9, FF4.1.3 sheet-1, FF4.1.4 sheet-2 absolute pins (FFW md:236 3e-10/7e-7/1e-6). ffw_fig_6 FF6.1.3/4/5 three-sheet (FFW 3e-9/7e-6/2e-5). Strength: closed-form + strongest published-number multi-sheet pins in slice.

Gaps: Single hand-pin sample. w-near-0/w-near-1 fail-loud not asserted at this layer. No BigFloat. e^(2zeta) large-Re-zeta overflow untested.

Ground truth needed: Closed-form (done at one point; extend to complex zeta + param sweep). FFW Fig4/6 published per-sheet errors serve as the high-precision oracle.

**pVI_transformed_rhs** (public-api) — `src/SheetTracker.jl:155`

Transformed zeta-plane P~VI RHS (FFW 2017 sec2.2, md:144). Maps z=0 branch to zeta=-inf but leaves z=1 as lattice zeta=2pi*i*k (e^zeta=1). Fixed singular surfaces w in {0,1,e^zeta} and e^zeta=1. The equation that REQUIRES circumambulation -- the whole multi-sheet machinery exists for it.

Regimes: w in {0,1,e^zeta} singular surfaces; zeta=2pi*i*k branch lattice ((e^zeta-1)^2->0); Float64 rounding makes branch-point RHS ~1e30 huge-but-finite not Inf; Complex/Taylor1 args; multi-sheet via cross_branch routing; BigFloat (untested; e^zeta-1 cancellation near lattice is precision-sensitive)

Current coverage: sheet_tracker_test ST.1.1 closed-form a=b=g=d=0 (=-1 exact), ST.1.2 a=b=g=d=1 (=101/24 exact), ST.1.3 end-to-end. ST.1.7 branch-lattice blow-up (>1e10 at zeta=2pi*i, finite near). ffw_fig_2 FF2.1.4 zeta-refuse pin (w_ref=0.40810323599709286-0.05243643326242345im, FFW 6e-7), ffw_fig_3 multi-sheet phase, ffw_fig_7 zeta sheet [0,0,0]/[1,0,0]. Strength: TWO closed-form hand-pins + FFW Fig2/3/7 published numbers + multi-sheet eval_at_sheet.

Gaps: Hand-pins at one zeta=log2 point only. ST.1.7 asserts huge-magnitude sentinel but does NOT assert the stepper throws (fail-loud unverified end-to-end at lattice). The (e^zeta-1) cancellation near lattice is a catastrophic-cancellation hotspot with NO BigFloat coverage -- a cancellation bug invisible at Float64.

Ground truth needed: Closed-form at multiple (zeta,w,wp) incl complex zeta (done at one). For lattice fail-loud: end-to-end assertion that integration THROWS when a step lands on zeta=2pi*i*k. mpmath PVI for independent absolute pin beyond loose-vs-tight.

**pVI_eta_transformed_rhs** (public-api) — `src/SheetTracker.jl:206`

Double-exp eta-plane P~VI RHS (FFW 2017 sec2.2.1, md:154). Stacks zeta=e^eta atop z=e^zeta to push z=0 branch out of the finite plane, giving a COMPACT branch-free rectangle Re eta<log(2pi)~1.838 containing infinitely many zeta-sheets folded. The -1 chain-rule term is the documented easy-to-miss derivation mistake.

Regimes: Re eta<log(2pi) branch-free region; Re eta>log(2pi) lattice eta=log|2pi*k|+i*arg(2pi*i*k) for |k|>=1; E-1=exp(2pi*i)-1 ~O(eps) near lattice (precision hotspot); v in {0,1,E} singular; nested-exp overflow for large Re eta; BigFloat (untested)

Current coverage: eta_pvi_test ET.1.1 closed-form a=b=g=d=0 = (2/3)(1-log2) exact (mutation M1 drop -1 bites -- load-bearing chain-rule term), ET.1.1bis a=b=g=d=1. ET.1.3 end-to-end direct-vs-eta. ET.1.4 branch-free finiteness + lattice blow-up (>1e20). ffw_fig_2 FF2.1.3 (v_ref=0.4525871664876391-0.16405500482535595im, FFW 4e-7), ffw_fig_7 FF7.1.3 eta pin + FF7.1.5 Schwarz reflection. Strength: TWO closed-form hand-pins + FFW Fig2/7 published numbers.

Gaps: Hand-pins at one eta=log(log2) point. E-1 cancellation near lattice is a SEVERE double-exp precision hotspot with NO BigFloat coverage. Large-Re-eta nested-exp overflow untested. Lattice fail-loud (stepper throws) unverified. CRITICALLY: this RHS is tested only via raw PadeTaylorProblem, never via the PainleveSolution wrapper.

Ground truth needed: Closed-form at multiple (eta,v,vp). mpmath/Mathematica PVI on the eta-rectangle for independent absolute pin. BigFloat recompute of the E-1 cancellation to bound precision loss.

**pVI_z_to_eta / pVI_eta_to_z** (public-api) — `src/SheetTracker.jl:244 (z_to_eta), src/SheetTracker.jl:259 (eta_to_z)`

Composition of two exponential transforms for the eta-plane (FFW 2017 sec2.2.1, md:146/151). Dependent variable unchanged (v=u=w); chain rule vp=z*zeta*up. Singular at z=0 (zeta=-inf) and z=1 (zeta=0->eta=-inf).

Regimes: z>1 real (eta real, FFW Fig2 IC regime); complex z (log(log z) double-branch subtlety); z near 1 (zeta->0, eta->-inf); z near 0 (zeta->-inf); BigFloat (untested)

Current coverage: eta_pvi_test ET.1.2 round-trip at FFW Fig2 IC z=10 (eta=log(log10), v=u, vp=z*log(z)*up; round-trip <=1e-13/1e-15, symmetric eta->z->eta; mutation M2 sign-flip on vp bites). ffw_fig_2/7 round-trips. painleve_test ET.1.5 frame wiring. Strength: closed-form round-trip pin at the published FFW Fig2 IC.

Gaps: Only real z>1 (z=10) tested for round-trip. Complex z (where log(log z) has genuine double-branch ambiguity) NOT tested -- a real branch-cut hazard the principal-branch composition could get wrong. No BigFloat. z-near-1/z-near-0 limits (eta->-inf) untested.

Ground truth needed: Closed-form round-trip at COMPLEX z (e.g. z=10*exp(i*theta)) to expose any log-of-log principal-branch mismatch. Verified (eta,v,vp) at FFW Fig2 IC (done). BigFloat: extended-precision round-trip.

**winding_delta** (internal-algorithm) — `src/SheetTracker.jl:283`

Signed per-step angle change seen from a branch point, the atom of Riemann-sheet winding (FFW 2017 sec2.2.2, md:163-178). Used by accumulate_winding (caller-side) and step_sheet_update (walker-side BranchTracker).

Regimes: small step |dtheta|<pi (contracted-safe regime); step straddling theta=pi discontinuity (the normalisation's purpose); |dtheta|>=pi single step (SILENTLY loses a full revolution -- documented precondition, NOT enforced/tested); step through the branch point (degenerate)

Current coverage: sheet_tracker_test ST.1.4: +1->+i = +pi/2, reverse = -pi/2 (closed-form), near-cut hop normalised to ~0.02 (<0.05). Mutation P (drop normalisation) bites via downstream accumulate_winding. Used inside BT.1.5 and PNB.1.3 chain-walk. Strength: closed-form on canonical steps.

Gaps: The |dtheta|>=pi precondition-violation behaviour is UNTESTED -- no test confirms what happens (silent corruption vs detectable error) when a step encloses the branch. Highest-risk untested branch in the sheet machinery: a walker producing too-large steps near a branch would corrupt every downstream sheet index with NO test catching it. branch=non-origin only lightly via BT tests.

Ground truth needed: Closed-form (atan2 of complex differences) is the oracle. CRITICAL missing test: assert the documented failure when |dtheta|>=pi (either fail-loud per Rule 1, or an explicit pin of the known-wrong behaviour as a deferred limitation).

**accumulate_winding / sheet_index** (public-api) — `src/SheetTracker.jl:304 (accumulate_winding), src/SheetTracker.jl:320 (sheet_index)`

Caller-side after-the-fact sheet bookkeeping (FFW 2017 sec2.2.2, md:187-189): cumulative winding along a walked path -> integer Riemann-sheet index (+1 per CCW revolution, -1 per CW). Sheet 0 = principal in [-pi,pi] convention.

Regimes: closed CCW loop (+2pi->+1); closed CW loop (-2pi->-1); multi-revolution (+4pi->+2); non-enclosing path (~0->0); sub-2pi winding (rounds to nearest); any step violating |dtheta|<pi (inherits winding_delta corruption)

Current coverage: sheet_tracker_test ST.1.5 closed loops (CCW +2pi, CW -2pi, non-enclosing ~0 to 1e-12). ST.1.6 sheet_index (+2pi->+1, -2pi->-1, +4pi->+2, 0->0, pi/2->0, 3pi/2->+1; composed with accumulate_winding). Mutation Q (round->floor) bites ST.1.6 at 3pi/2. Strength: closed-form on canonical loops.

Gaps: Square-loop paths only (4-segment). Real walker paths (irregular many-step) only via PNB.1.3 chain-walk. No test of a path that legitimately winds multiple times around a branch with FINE steps (the realistic circumambulation). sheet_index round-half boundary at exactly (2k+1)pi only at 3pi/2; round-half-to-even quirk at e.g. pi untested.

Ground truth needed: Closed-form (sum of analytic angle differences) -- complete oracle. A realistic many-step circumambulation path with known winding number would test the accumulation under the actual walker step regime.

**segment_crosses_cut** (public-api) — `src/BranchTracker.jl:112`

Walker-side cut-crossing predicate (ADR-0013 Decision, FFW 2017 sec2.2.2 md:178). Half-line cut C(b,alpha)={b+s*e^(i*alpha):s>=0}; default alpha=pi matches Julia log convention. The geometric primitive enabling refuse/cross modes.

Regimes: cut at alpha=pi (log convention); cut at alpha=0 (steered away); branch off origin; endpoint-touching t=0/t=1 (excluded); branch-point-touching s=0 (excluded); parallel det=0 (false); backward-extension s<0 (false); oblique cut angles alpha not in {0,pi} (UNTESTED)

Current coverage: branch_tracker_test BT.1.1 (7 canonical cases on alpha=pi: crossing both directions, positive-side no-cross, above-axis no-cross, branch-point-touch excluded, far crossing, backward-extension excluded), BT.1.2 (parallel + endpoint t=0/t=1 excluded), BT.1.3 (branch off origin, alpha=0, left/right of b). Mutation B1 (drop s>0) bites 5 assertions. path_network_branch_test PNB.1.2 direct filter use. Strength: thorough axis-aligned closed-form geometry.

Gaps: ALL tested cuts axis-aligned (alpha in {0,pi}) and branches on real/imaginary axes. Oblique angles (pi/2, pi/4, irrational) and branches at general complex positions UNTESTED -- the imag(delta*conj(u))/imag(d*conj(u)) Cramer projection could carry a sign subtlety at oblique angles axis-aligned tests cannot expose. Near-parallel (small-but-nonzero det) numerical robustness untested.

Ground truth needed: Closed-form (2x2 linear intersection) -- complete oracle. Add oblique-angle cases (alpha=pi/4 with branch at 1+1im) hand-computed to expose any angle-dependent sign bug.

**any_cut_crossed / step_sheet_update** (public-api) — `src/BranchTracker.jl:130 (any_cut_crossed), src/BranchTracker.jl:156 (step_sheet_update)`

Multi-branch walker-side enforcement (ADR-0013): scan all cuts (refuse mode) and update per-branch sheet counters on deliberate crossing (cross mode). CCW crossing->+1, CW->-1, no-crossing->pass-through unchanged. Cut-crossing GATES the bump (small rotation not enclosing the branch leaves counter alone).

Regimes: single branch; multiple branches at-most-one crossed per step (tested); TWO cuts crossed by ONE step (UNTESTED); non-crossing pass-through; fresh-allocation non-aliasing; BigFloat element branches (untested)

Current coverage: branch_tracker_test BT.1.4 any_cut_crossed (3-branch lattice, crosses upper/origin/lower individually + no-cross). BT.1.5 step_sheet_update sign follows winding (CW->-1, CCW->+1, non-cross pass-through, non-aliasing). BT.1.6 per-branch isolation (only crossed branch updates). path_network_branch_test PNB.1.3 + ffw_fig_2 FF2.1.6 walker-level. Mutations B2 (sign flip) 9 RED, B3 (drop cut-gate) 52 RED. Strength: closed-form geometry + walker integration + figure pin.

Gaps: Every tested step crosses AT MOST ONE cut. A single step crossing TWO cuts simultaneously (two counters bumped in one step -- plausible near zeta=2pi*i lattice where adjacent branch cuts are close) is UNTESTED. Multi-branch winding explicitly deferred (SheetTracker docstring). No oblique-angle multi-branch case (inherits segment_crosses_cut gap).

Ground truth needed: Closed-form (composition of segment_crosses_cut + winding_delta) -- oracle exists. Add a hand-built step crossing two cuts to confirm both counters bump independently.

**resolve_cut_angles** (oracle-helper) — `src/BranchTracker.jl:186 (scalar), src/BranchTracker.jl:190 (collection)`

Kwarg normalisation for the public branch_cut_angles (ADR-0013): broadcast scalar or validate per-branch tuple. Fail-loud on length mismatch (CLAUDE.md Rule 1).

Regimes: scalar broadcast; matching tuple; matching vector (length-checked); mismatched length (throws); empty branches (); Float64 coercion of Int/Rational angles

Current coverage: branch_tracker_test BT.1.7: scalar->NTuple{3,Float64} all equal, matching tuple Float64-coerced, mismatched lengths (2 and 4 vs 3) throw ArgumentError, empty branches with scalar and with () both ->(). Strength: closed-form / exhaustive small-case.

Gaps: Minimal -- thin normalisation helper, well covered. No distinct Vector{Float64} (vs Tuple) angles input case, though the generic method handles both.

Ground truth needed: None beyond existing -- pure data normalisation, fully pinned by small-case tests.

**PainleveProblem builder (constructor + _build_I/II/III/IV/V/VI)** (problem-builder) — `src/Painleve.jl:202 (dispatch), src/Painleve.jl:241 (_build_transformed), src/Painleve.jl:276 (_build_VI)`

Per-equation problem-builder for all six Painleve equations (ADR-0006). For PIII/PV/PVI: validates, guards IC off fixed branch point(s) (z=0 for III/V; z in {0,1} for VI), maps IC+zspan endpoints into the zeta-frame (or eta for :transformed_eta), builds PadeTaylorProblem in the solve-frame, stores to_frame/from_frame + equation + params + name (always nothing here).

Regimes: the six equations + two PVI frames; complex vs real ICs/span (affects solve_pade applicability); IC on a fixed branch point (throws); missing/extra/unknown-equation kwargs (throws); PIV u=0 movable-zero (RHS throws DomainError); order override vs default 30

Current coverage: painleve_test PV.1.1-1.5 constructor correctness per equation (RHS sample-match against wrapped factory verbatim, frame, params, zeta-mapped IC, from_frame-of-to_frame=id). PV.2.1 fail-loud guards (unknown eq, missing param, extra param, IC on z=0/z=1 branch point). eta_pvi_test ET.1.5 :transformed_eta wiring (zspan[1]/y0=pVI_z_to_eta, unknown-frame throws). Mutations A/B (RHS perturb) + M3 (eta->zeta misdispatch) bite. Strength: closed-form RHS sample-match + exhaustive fail-loud.

Gaps: Transformed-frame IC mapping checked via pp.problem.y0 == z_to_zeta(...) equality but only at one IC per equation. Branch-point guard uses EXACT equality (any(b->z0==b, branch_points)) -- a z0 floating-point-close-to-but-not-equal a branch point passes the guard and the transform produces a huge-but-finite IC; this near-branch-point regime is untested and could silently produce garbage. zeta_end computed via z_to_zeta(zspan[2],u0,up0)[1] reusing the IC u0/up0 (valid since 1st component ignores u,up) -- this projection assumption is implicit and untested for the eta map at the builder level (ET.1.5 checks only zspan[1]).

Ground truth needed: Closed-form: builder output is deterministic assembly, fully pinnable (mostly done). Add a near-branch-point IC test (z0 = branch + 1e-12) to pin the exact-equality guard's blind spot.

**solve_pade / path_network_solve forwarding (::PainleveProblem)** (public-api) — `src/Painleve.jl:340 (solve_pade), src/Painleve.jl:365 (path_network_solve)`

Solver forwarding methods (ADR-0006/0007). :direct forwards transparently. :transformed: path_network_solve maps the caller's z-frame grid into the solve-frame, solves, returns a PainleveSolution presenting z-frame poles/grid. solve_pade serves a transformed problem only when the zeta-domain is real-typed.

Regimes: :direct transparent forward; :transformed complex-zeta via path_network_solve; :transformed real-zeta via solve_pade (PV real ICs); :transformed complex-zeta via solve_pade (THROWS); :transformed_eta (rides generic non-direct path -- UNTESTED at forwarding level)

Current coverage: painleve_test PV.3.1 (:direct forward byte-identical raw, :transformed path_network_solve returns PainleveSolution, solve_pade throws ArgumentError on complex-zeta III). PV.4.1 end-to-end PII pole field finite + grid order preserved. painleve_solution_test PS.2.3 (:transformed real-zeta PV solve_pade IC round-trip), PS.5.3 (:transformed path_network grid round-trips). Mutation C (invert guard) + M3 (grid mapping) bite. Strength: behaviour pinned for direct + transformed-zeta.

Gaps: NO test forwards a :transformed_eta problem through solve_pade or path_network_solve. The eta-frame grid mapping (z->eta via log-of-log) and the eta-frame point-value remap out are completely untested at the forwarding level -- code branches only on frame===:direct so :transformed_eta silently uses the generic path, but its from_frame (double-exp pVI_eta_to_z) is never exercised through grid_values/poles. The solve_pade real-zeta guard eltype(pp.problem.zspan)<:Complex -- a :transformed_eta with real eta-domain would attempt real-axis stepping; untested.

Ground truth needed: End-to-end: a :transformed_eta path_network_solve whose grid_values round-trip the z-frame grid (mirror of PS.5.3 for eta) and whose values match a direct PVI solve. Single biggest forwarding-coverage hole in the slice.

**PainleveSolution wrapper (struct + callable + poles + grid_values + accessors + show)** (public-api) — `src/PainleveSolution.jl:80 (struct), src/PainleveSolution.jl:132 (callable), src/PainleveSolution.jl:159 (poles), src/PainleveSolution.jl:194 (grid_values)`

Self-describing solve-output wrapper (ADR-0007): uniform z-frame access surface over six raw solve types. For :transformed, maps callable/poles/grid_values back to the natural z-frame at every boundary while the Pade store stays zeta-frame (documented honest leak). Fail-loud: grid-type raw not callable; unwired raw has no poles route.

Regimes: :direct (identity maps); :transformed (zeta<->z via from_frame); :transformed_eta (eta<->z via double-exp -- UNTESTED at access level); PadeTaylorSolution raw (callable); PathNetworkSolution raw (grid, not callable); unwired raw (String etc, fail-loud); name=:tritronquee/etc (only for :direct equations)

Current coverage: painleve_solution_test PS.1.1 provenance, PS.2.1 direct callable transparent, PS.2.2 grid-raw fail-loud, PS.2.3 transformed-PV callable round-trip, PS.3.1 poles direct identity (hand-built pole at 1.5), PS.3.2 poles transformed-PV maps exp, PS.4.1 show, PS.5.1 grid_values direct passthrough, PS.5.2 grid_values transformed-PV (exp grid, u unchanged, up/z), PS.5.3 transformed-PV grid round-trip, PS.6.1 unwired-raw fail-loud on every verb. Mutations M1-M4 bite. Strength: hand-built deterministic raw + closed-form PV map pins + fail-loud.

Gaps: (a) ALL :transformed access tests use the PV single-exp map (from_frame=pV_zeta_to_z, z=exp(zeta)). PIII (z=exp(zeta/2), u=w/z) and the eta double-exp map (z=exp(exp(eta)), up=vp/(z*zeta)) are NEVER exercised through callable/poles/grid_values. (b) poles() uses _coord (projection on dummy zeros) while grid_values uses the full from_frame triple -- DIFFERENT code paths; only PV-via-each is pinned, so a bug in the eta or PIII from_frame's up component (which _coord drops but grid_values keeps) is untested. (c) name never exercised as non-nothing for a :transformed solution (all named constructors are :direct). (d) No BigFloat raw.

Ground truth needed: Hand-built PathNetworkSolution wrapped as :transformed PIII and as :transformed_eta PVI, asserting poles map via exp(zeta/2)/exp(exp(eta)) and grid_values remap up correctly (mirror of PS.3.2/PS.5.2 for the other two maps). This isolates the wrapper's frame-mapping from the solver.

**Slice notes:** CROSS-CUTTING OBSERVATIONS for the bug hunt:

(1) ZERO BigFloat/Arb coverage across the ENTIRE slice. Every test in coord_transforms_test, sheet_tracker_test, branch_tracker_test, painleve_test, painleve_solution_test, eta_pvi_test, path_network_branch_test, sheet_aware_stage2_test, and ALL ffw_fig_*_test files uses only Float64/ComplexF64. Transform closures are generic and CLAUDE.md emphasises a BigFloat-256 / Arblib tier, but a precision-loss or type-instability bug specific to BigFloat/Complex{BigFloat} in any transformed RHS or coordinate map would be completely invisible. (painleve_hierarchy_test touches BigFloat but it is the PI^(m) hierarchy, OUT of slice.)

(2) STRONG published-number oracle infrastructure EXISTS and is well-used for the figure layer. ffw_fig_*_test.jl pin against FFW 2017's own per-sheet symmetry-based error estimates (Fig1 4e-10/4e-8; Fig2 md:195 4e-7/3e-6 eta, 6e-7 zeta-refuse, 6e-4/5e-4 zeta-cross; Fig4 md:236 3e-10/7e-7/1e-6; Fig6 3e-9/7e-6/2e-5) AND captured tight-tol reference complex numbers (Fig2 v_ref=0.4525871664876391-0.16405500482535595im, w_ref=0.40810323599709286-0.05243643326242345im). These are the rigorous oracles; the loose-vs-tight + absolute-pin two-pronged pattern (worklog 036) is the established methodology a new test should follow.

(3) :transformed_eta (eta-plane PVI, double-exp) is BUILDER-tested only. ET.1.5 verifies PainleveProblem(:VI; frame=:transformed_eta) constructs correctly. But NO test exercises the PainleveSolution ACCESS surface for an eta-frame solution: no solve_pade/path_network_solve forwarding, no poles(), no grid_values(), no callable. Forwarding branches only on frame===:direct vs not (Painleve.jl:341,366; PainleveSolution.jl:134,167,196), so eta silently rides the generic :transformed path. The eta from_frame pVI_eta_to_z has z=exp(exp(eta)) and up=vp/(z*zeta) -- these wrapper-level remaps are UNTESTED. ffw_fig_7 exercises the eta RHS only via raw PadeTaylorProblem, never the wrapper.

(4) poles() z-frame mapping tested ONLY for PV single-exp (PS.3.2: p_trans ~ exp.(p_dir)). PIII (exp(zeta/2)) and eta double-exp maps NOT covered by poles(). poles() uses _coord (projection on dummy zeros) while grid_values uses the FULL from_frame triple -- two different code paths; only PV-via-each is pinned.

(5) name field for transformed equations is ALWAYS nothing. Named constructors (tritronquee->I, hastings_mcleod->II, pii_rational/pii_airy->II, piv_entire->IV) are ALL :direct. No named constructor builds :III/:V/:VI, so PainleveSolution.name for a multivalued problem is exercised only as nothing.

(6) winding_delta single-step |dtheta|>=pi precondition is a documented landmine, NOT enforced. Silently loses a full revolution if a step encloses the branch. accumulate_winding/sheet_index/step_sheet_update all inherit this. No test asserts what happens (or that it fails loud) when violated. A walker producing too-large steps near a branch would corrupt sheet indices SILENTLY -- the most likely place a real sheet-tracking bug hides untested.

(7) segment_crosses_cut conventions well unit-tested (BT.1.1-1.3) but ONLY on cuts at alpha in {0,pi} and branches on real/imaginary axes. Oblique angles (pi/2, pi/4, irrational) and branches at general complex positions UNTESTED; the imag(delta*conj(u))/imag(d*conj(u)) Cramer projection could carry an angle-dependent sign subtlety axis-aligned tests miss.

(8) Multi-branch winding (around BOTH zeta=0 AND zeta=2pi*i simultaneously) explicitly deferred and loosely covered: BT.1.4/BT.1.6 use multiple branches but each step crosses at most ONE cut. A single step crossing TWO cuts at once (plausible near the lattice) updating two counters is untested.

(9) pVI zeta/eta branch-point blow-up tests (ST.1.7, ET.1.4) assert magnitude>1e10/1e20 as the float sentinel and rely on the downstream Pade stepper to fail-loud, but do NOT assert the stepper actually throws -- the fail-loud contract is assumed, not verified end-to-end at the branch lattice. The (e^zeta - 1) / (E-1) cancellation near the lattice is a catastrophic-cancellation precision hotspot with NO BigFloat coverage.

### vector-stack

**vector_taylor_coefficients (VectorCoefficients.jl:136)** (internal-algorithm) — ``

Vector Taylor jet y'=f(z,y) component-wise bootstrap (FW 2011 2.1.2b; ADR-0020)

Regimes: Float64; BigFloat-256; Complex jets; d=1 bit-identical; decoupled-coupled; entire-meromorphic

Current coverage: VC.1.1-1.6 exp/cos/sin closed-form, d=1<->scalar bit-identical, BF-256, mutation-proven

Gaps: Complex jets never closed-form-checked; large-d only smoke

Ground truth needed: Closed-form coeffs of complex vector ODE y'=i*y

**vector_step_jorba_zou (VectorStepControl.jl:158)** (internal-algorithm) — ``

Norm-based vector Jorba-Zou step (ADR-0021)

Regimes: d=1 bit-identical; 2-norm vs inf-norm; zero-norm skip; second-stepsize fallback; Complex jets untested

Current coverage: VSC.1.1-1.6 d=1<->scalar exact, 2/inf-norm, mutation-proven

Gaps: All real Float64; complex jets + fallback value untested

Ground truth needed: Both-last-zero fallback jet; complex-jet norm

**vector_pade_step_with_pade! (VectorStepper.jl:225)** (internal-algorithm) — ``

Shared-Q Pade step via shared_pade_select (ADR-0028 cell A/B), pole guard sqrt(eps)||Q|| (ADR-0019)

Regimes: d=1<->scalar; entire; meromorphic shared pole; step-on-pole; BF-256; Complex h; cell A vs B

Current coverage: VS.1.1-1.6 d=1<->scalar 1e-12, harmonic cell-B BF<1e-28, shared-pole bridge >=1e8

Gaps: Complex-h step no oracle; ADR-0028 cell choice not isolated; d>=3 coupled meromorphic no oracle

Ground truth needed: d>=2 shared-pole complex-h step to known z; hook for cell selected

**vector_solve_pade + Problem/Solution (VectorProblems.jl:240)** (public-api) — ``

Fixed/adaptive IVP driver, per-segment shared-Q, dense P_i/Q, :fixed|:jorba_zou

Regimes: :fixed vs :jorba_zou; d=1<->scalar; entire ~1e-8; meromorphic short of pole; dense+out-of-range; final clamp; BF; rejects complex h

Current coverage: VP.1.1-1.6+VSC.1.4 d=1<->scalar, harmonic BF<1e-16, meromorphic 1e-9, dense+clamp

Gaps: :jorba_zou only entire harmonic; entire atol=1e-8 hides ~1e-9; dense O(N) scan

Ground truth needed: :jorba_zou vs :fixed on meromorphic tighter atol; high-precision coupled trajectory

**_select_wedge/_candidate_pole_disc (VectorWedgeStep.jl:421)** (internal-algorithm) — ``

Wedge selector (ADR-0025 B2) max pole-free disc h_mag*min|t*|, CLEAR_CAP, Froissart filter

Regimes: multi/single pole; :all_candidates_failed; -Inf excluded; pole-free tie-break; Froissart filtered; :max_q_root vs :min_y

Current coverage: VPN.3.2 P_I^(2) max_q_root maximal disc; M4 16 RED

Gaps: CLEAR_CAP/tie-break implicit; Froissart filter only via mirror; no :all_candidates_failed fixture

Ground truth needed: Hand-built 5-candidate fixture closed-form discs; planted doublet not shrinking disc

**_adaptive_h (VectorWedgeStep.jl:569)** (internal-algorithm) — ``

Non-ratcheting step (ADR-0026 Am4) h=clamp(SAFETY*D_local), throws :step_collapse

Regimes: mid pole; far/constant Q->h_max; one-step recovery; near pole throw; numerators=nothing legacy vs supplied filtered

Current coverage: VPN.3.3 mid-pole 1e-12, one-step recovery LOAD-BEARING; VPN.3.4 throw; M5/M6/M8 4 RED

Gaps: Production numerators-Froissart-filter path only via live walk; no fixture asserts a planted doublet ignored so h does NOT collapse

Ground truth needed: Hand-built numerators+denominator one genuine + one doublet (residue<1e-4), assert no collapse

**vector_path_network_solve (VectorPathNetwork.jl:642)** (public-api) — ``

Stage-1 5-dir wedge walk, canonical shared-Q, Stage-2, VC-10 rng, resilient :skip, skip_covered (ADR-0026)

Regimes: adaptive vs fixed; :max_q_root vs :min_y; :throw vs :skip; rng VC-10; skip_covered; d=3/4/5; pole-crossing; 3 failures

Current coverage: VPN.1/3/4/5.x Riccati 1e-8, A_4 tree, P_I^(2) |x|>=15, VC-10 1e-8, resilient; M1-M9; VPO.2 u(30)

Gaps: :all_candidates_failed/:step_collapse never deterministically triggered (only :unreachable); skip_covered demo DISABLED

Ground truth needed: Pole-on-chord geometry to fire those reasons via driver; grid where skip_covered reduces nodes

**_stage2_fill/_validity_radius (VectorPathNetworkStage2.jl:386)** (internal-algorithm) — ``

Stage-2 Voronoi fill gated by B1 R_gate=min(s(tol)*h_JZ*h_v,0.5*h_v*min|t*|) (ADR-0025 Am1)

Regimes: covered vs NaN; extrapolate bypass; t=0; pole clamp vs truncation; cache vs fallback; tol-dep; degenerate throws

Current coverage: VPN.2.1-2.7 Riccati 1e-6, harmonic 1e-8, NaN, B1 gate s-factors/cache/clamp/throws, B4 honesty; M1-M7

Gaps: s(tol) clamp between calibration points not asserted; Voronoi tiebreak untested

Ground truth needed: High-precision fill equidistant between two nodes; pre-B1 empty-jets fill = cached fill

**extract_poles_shared_q (VectorPoleField.jl:245)** (public-api) — ``

Shared-Q pole extractor: z-dist filter radius_t*h_max (S7), Froissart min_residue (ADR-0028), cluster, min_support

Regimes: d>=2 shared pole; single node; cross-node cluster; scale vs legacy atol; radius_t S7; Froissart filtered; constant Q; empty sol; adaptive-h rep z-closest

Current coverage: VPN.1.1/1.3/1.3b/1.3c Riccati 1e-8 + hand-built z-map/cluster/S7/spurious; VPO.4 1e-5; NYF.1.3/1.4; M1/M3/M3b/M3c/M5

Gaps: NO dedicated test file; min_residue Froissart filter (line 294) NOT directly tested; rep z-closest only via M5; 1e-4 vs scalar 1e-8 untested

Ground truth needed: Hand-built sol one genuine + one doublet (assert filtered); two nodes diff h one pole (assert z-closest rep)

**noumi_yamada_rhs/NoumiYamadaProblem (NoumiYamada.jl:123)** (problem-builder) — ``

Even-parity A_{2n}^(1) cyclic system (NY1998 A2n); sum alpha=1; invariant sum f_j=t (ADR-0022)

Regimes: n=1 A_2=>PIV; n=2 A_4 d=5; n>=3 cyclic; Type C nonzero IC; sum f_j=t validation; failure modes

Current coverage: NY.1.1-1.7 RHS vs hand oracle+cyclic, A_4 Type C f_j=t/5 1e-9; piv_test V5b A_2<->scalar PIV

Gaps: Only n in {1,2,3}; Type C only solved-trajectory oracle; generic A_4 no closed form

Ground truth needed: High-precision generic A_4 trajectory; A_6 Type C cross-check

**noumi_yamada_rational/_backlund/_rotation (NoumiYamadaSymmetry.jl:142)** (oracle-helper) — ``

Affine-Weyl W(A_{2n}^(1)) Backlund + rational oracles (NY1998 WAl1/WAl2); :A/:C all n, :B n=2 only

Regimes: :A/:C all n + :B n=2 only; exact Rational ==; f_i=0 throws; Type C vector-solve; Backlund preserves variety; failure modes

Current coverage: NYS.1.1-1.7 rational satisfy ODE exactly, Type C 1e-10, all group relations, Backlund preserves; M1/M2/M3

Gaps: n in {1,2} only; n>=3 not exercised; no long Weyl word; Type C only solver-tie test

Ground truth needed: n=3 A_6 relation spot-checks; longer Weyl-word orbit to a known rational solution

### heun-closedform

**heun_g** (public-api) — `src/Heun.jl:288`

General Heun function HeunG(a,q,α,β,γ,δ;z), DLMF 31.2.1, four regular singular points {0,1,a,∞}; Fuchsian constraint ε=α+β+1-γ-δ. Frobenius series (DLMF 31.3.2-31.3.4) for |z|≤epsilon_start, else path-network continuation from z0=epsilon_start. Normalised HeunG(0)=1. ADR-0018.

Regimes: |z|≤epsilon_start (pure Frobenius, no walk); |z| inside disc min(1,|a|) via walk; past z=1 and z=a (multi-singularity, Schwarz-symmetric upper-half walk); off real axis (upper/lower half via Schwarz mirror); real params (enforce_real_axis_symmetry=true arm); complex params (branch_points=(0,1,a) arm — UNTESTED); Float64 vs BigFloat (BigFloat UNTESTED); degenerate-parameter identities (q=0,α=0 ⇒ ≡1)

Current coverage: test/heun_test.jl HG.1-HG.5, HX.1, HX.2 — published-number (wolframscript oracle, 30-digit) at rtol 1e-12 inside disc / 1e-8 past singularities; closed-form identity (q=0,α=0⇒1) at 1e-14. STRONG on the single tuple (a,q,α,β,γ,δ)=(2,1,1,2,3,4) but that is the ONLY parameter set pinned.

Gaps: Only one (a,q,α,β,γ,δ) tuple ever oracle-pinned; ε=-3 only value of the Fuchsian constraint exercised. Complex-parameter else-branch (branch_points/branch_cut_angles path, lines 341-347) entirely untested. BigFloat fallback (ADR-0018 decision #4) untested. a=0 guard untested. cross_branch=true / sheet kwarg untested. n_taylor=60 adequacy near |z|=0.05 unverified for stiff params. order/h kwargs not swept.

Ground truth needed: A multi-parameter wolframscript HeunG oracle table spanning a∈{complex, |a|<1, |a|>1}, several (q,α,β,γ,δ) including complex params and large |params|, plus a BigFloat high-precision Mathematica run to test the BigFloat-256 path. The existing 42-record oracles.txt should be more fully consumed.

**heun_c** (public-api) — `src/Heun.jl:365`

Confluent Heun function HeunC(q,α,γ,δ,ε;z), DLMF 31.12.1, three regular singular points {0,1} + irregular at ∞. Frobenius series via Motygin 2018 eqs 3-4 (DLMF does not publish HeunC recurrence) for |z|≤epsilon_start, else path-network from z0 with branch_points=(0,1). Normalised HeunC(0)=1. ADR-0018.

Regimes: |z|≤epsilon_start (Frobenius); |z|<1 inside disc (2-source pinned: Mathematica↔Motygin agree); z≥1 (different analytic sheets — Motygin & Mathematica disagree); off real axis; real params (symmetry arm) vs complex params (branch_points arm — UNTESTED); Float64 vs BigFloat (UNTESTED)

Current coverage: test/heun_test.jl HC.1-HC.4, HX.1, HX.2 — published-number 2-source-pinned (z<1) at 1e-12/1e-8; z≥1 (HC.4) pinned only at rtol 1e-3 against a hand-probed Mathematica upper-limit value HeunC[...,z+0.001i]; HC.1 leading-order Frobenius check at only 1e-3 (no high-precision oracle). Single tuple (q,α,γ,δ,ε)=(1,1,2,3,-1).

Gaps: z≥1 sheet is the weakest oracle in the whole slice (rtol 1e-3, single point, hand-probed) — prime hiding spot for a sheet-selection bug. Only one parameter tuple. HC.1 Frobenius regime has no high-precision oracle (only a 1e-3 leading-order sanity check). Complex-param branch_points=(0,1) arm untested. BigFloat untested. heun_c non-positive-integer-γ logarithmic case only tested at γ=0, not γ=-1/-2.

Ground truth needed: Connection-formula-resolved HeunC reference values for z≥1 on a single agreed sheet (the open follow-up bead padetaylor-?heun-monodromy), so HC.4 can tighten from 1e-3 to machine precision. A multi-parameter Motygin+Mathematica cross-val table for the Frobenius regime, plus BigFloat references.

**_heun_g_frobenius_at_0** (internal-algorithm) — `src/Heun.jl:111`

HeunG three-term coefficient recurrence DLMF 31.3.2-31.3.4: c_0=1, c_1=q/(aγ), c_{j+1}=((Q_j+q)c_j - P_j c_{j-1})/R_j with P_j=(j-1+α)(j-1+β), Q_j=j((j-1+γ)(1+a)+aδ+ε), R_j=a(j+1)(j+γ), ε=α+β+1-γ-δ. Horner eval of u and term-by-term u'. Throws on non-positive-integer γ (logarithmic 2nd solution) and a∈{0,1}.

Regimes: small |z| direct sum; z0=epsilon_start IC sampling for the walk; N=60 truncation; γ non-positive-integer throw; a=0 throw; a=1 throw; complex z / complex params

Current coverage: Indirectly via heun_g HG.1/HG.5 (q=0,α=0⇒1 identity exercises the P_j=0 vanishing branch) and HX.1 (γ=0,γ=-1,a=1 throws). No DIRECT unit test of the coefficient vector against an independent c_j computation; recurrence correctness is inferred only through the end-to-end value.

Gaps: No direct assertion on returned coeffs[] vs an independent recurrence (e.g. Mathematica SeriesCoefficient). The c_1=q/(aγ) seed and the j-indexed P_j/Q_j/R_j are only tested through one (a,q,...) tuple plus the degenerate q=0,α=0 case. a=0 throw untested. Off-by-one in the j-loop range (1:N-1) could be masked because high coeffs are tiny at |z|=0.05.

Ground truth needed: Independent Frobenius coefficients c_0..c_N from Mathematica SeriesCoefficient[HeunG[...],{z,0,k}] for several parameter sets including complex, asserted element-by-element — a closed-form coefficient oracle.

**_heun_c_frobenius_at_0** (internal-algorithm) — `src/Heun.jl:168`

HeunC three-term recurrence Motygin 2018 eqs 3-4: b_0=1, b_1=-q/γ, b_n=(Q_n b_{n-1}+R_n b_{n-2})/P_n with P_n=n(n+γ-1), Q_n=-q+(n-1)(γ+δ-ε+n-2), R_n=α+(n-2)ε. Throws on non-positive-integer γ.

Regimes: small |z| direct sum; IC sampling at z0; N=60 truncation; γ=0 / non-positive-integer throw; complex z/params

Current coverage: Indirectly via heun_c HC.1 (only a 1e-3 leading-order check b_1=-q/γ=-1/2) and HX.1 (γ=0 throw). The recurrence body (n≥2) is exercised only through end-to-end HC.2/HC.3 values at one tuple.

Gaps: Weakest-tested recurrence in the slice — only the b_1 seed is checked at 1e-3, no high-precision coefficient oracle, no direct coeff-vector assertion. R_n=α+(n-2)ε and Q_n's γ+δ-ε grouping only at (q,α,γ,δ,ε)=(1,1,2,3,-1). DLMF-vs-Motygin convention transcription (a known hallucination risk per ADR-0018) rests on a single end-to-end tuple.

Ground truth needed: Mathematica SeriesCoefficient[HeunC[...]] coefficient table for multiple parameter sets, asserted element-by-element, to independently validate the Motygin recurrence transcription (the ADR flags this as the convention DLMF does not publish).

**_heun_g_rhs / _heun_c_rhs** (internal-algorithm) — `src/Heun.jl:215`

Thread-safe RHS closures: HeunG u''=-((γ/z+δ/(z-1)+ε/(z-a))u' + (αβz-q)/(z(z-1)(z-a))u) with ε=α+β+1-γ-δ; HeunC u''=-((γ/z+δ/(z-1)+ε)u' + (αz-q)/(z(z-1))u). Feed path_network_solve. DLMF 31.2.1 / 31.12.1 rearranged.

Regimes: regular z away from {0,1,a}; near singular points (handled by walker avoidance); Taylor1 coefficient args (must capture no mutable state — thread-safe); complex z

Current coverage: Indirectly via every heun_g/heun_c walk test (HG.2-HG.4, HC.2-HC.4). No direct unit test of the closure value at a known (z,u,u') point. NOTE: _heun_g_rhs (line 216) recomputes ε=α+β+1-γ-δ independently of _heun_g_frobenius_at_0 (line 127) — a duplicated formula; a sign/transposition in one but not the other would only show as IC-vs-trajectory mismatch.

Gaps: No direct evaluation-point assertion. The duplicated ε formula across recurrence and RHS is not cross-checked for consistency by any test. Thread-safety / Taylor1-arg behaviour not directly asserted (only implicit in the parallel wedge evaluator). Only ε=-3 case.

Ground truth needed: Direct closure-value assertions: compute u'' at several (z,u,u') from an independent symbolic evaluation of DLMF 31.2.1/31.12.1, including a complex z and Taylor1 inputs, to pin the RHS independently of the walker.

**pii_rational** (problem-builder) — `src/PainleveClosedForm.jl:206`

PII rational solution u_n for α=n: u_1=-1/z, u_2=(4-2z³)/(4z+z⁴), u_3=3z²(160+8z³+z⁶)/(320-24z⁶-z⁹). FW2014 eq.6. IC derived from the exact closed form at zspan[1]; pole guards at z=0 (n∈{1,2}) and z=-∛4 (n=2). ADR-0010.

Regimes: n=1,2 (real-axis poles, default zspan avoids z=0); n=3 (regular at 0, IC (0,0)); zspan[1] at a pole (throws); u_3 pole near z≈1.508 (solver must Padé-bridge, but tests stop short of it); n∉{1,2,3} throw

Current coverage: test/painleve_closed_form_test.jl CF.1.1-CF.1.4, CF.4.1, CF.4.4 — closed-form oracle: exact-equality IC check (y0==formula, no tolerance) AND solver cross-check vs _u_pii_rational at rtol 1e-10. Mutation-proven (MC1, MC3). STRONG.

Gaps: Solver cross-check (CF.1.4) only integrates over short clear intervals (n=1,2 on [1,1.5]; n=3 on [0,1]) — the documented u_3 pole-crossing at z≈1.508 (the actual Padé-bridge stress case) is explicitly avoided, so the meromorphic-jet / pole-crossing regime is untested here. n=2 second pole at z=-∛4 guard is in code but the cbrt-tolerance branch (atol 1e-12, line 234) is not directly tested. No BigFloat.

Ground truth needed: Already has its exact closed form as oracle (the ideal). The one missing oracle is a pole-crossing acceptance test: integrate u_3 from (0,0) ACROSS z≈1.508 and assert the trajectory matches u_3(z) on the far side — needs the closed form (available) plus confidence the Padé bridge is engaged.

**pii_airy** (problem-builder) — `src/PainleveClosedForm.jl:266`

PII Airy solution u_{n+1/2}(z;θ) for α=n+1/2: u_{1/2}=-Φ, u_{3/2}=(2Φ³+zΦ-1)/(2Φ²+z), Φ=φ'/φ, φ=cos(θ/2)Ai(-z/2^⅓)+sin(θ/2)Bi(-z/2^⅓), Φ'=-z/2-Φ². FW2014 eqs.8-9. θ∈[0,2π) selects the one-parameter family. ADR-0010.

Regimes: n=0 (u_{1/2}); n=1 (u_{3/2}, quotient-rule derivative block); θ=0 pure-Ai; θ=π pure-Bi; generic θ (both Ai+Bi — UNTESTED); near Ai zero at z≈2.946 (default zspan stays clear); n∉{0,1} throw

Current coverage: test/painleve_closed_form_test.jl CF.2.1-CF.2.4, CF.4.2 — exact closed-form oracle via SpecialFunctions Ai/Bi: IC check at rtol 1e-12, solver cross-check (n=0 only) at 1e-9. Mutation-proven (MC4: φ' chain-rule sign). Reasonably strong but with specific holes.

Gaps: θ tested only at 0 and π — the two values where cos(θ/2) OR sin(θ/2) vanishes; a cos/sin swap in _airy_phi (line 110) survives both. Generic θ (e.g. π/2) where both Ai and Bi contribute is UNTESTED. pii_airy(1)/u_{3/2} DERIVATIVE has NO closed-form assertion anywhere: CF.2.2 punts the z=0 derivative (line 135-137) and the solver cross-check CF.2.3 only runs n=0 — so the entire n==1 quotient-rule block _u_pii_airy lines 144-148 (Np,Dp) is never validated against ground truth. No BigFloat.

Ground truth needed: (1) A pii_airy(0; θ=π/2) IC + solver cross-check (closed form already available — just add the θ=π/2 case). (2) A pii_airy(1) solver cross-check OR a direct _u_pii_airy(1,θ,z) derivative assertion against a finite-difference-of-closed-form or a Mathematica PII-Airy reference, to cover the u_{3/2} derivative quotient rule.

**piv_entire** (problem-builder) — `src/PainleveClosedForm.jl:302`

PIV entire solutions u(z)=-2z at (α,β)=(0,-2) and u(z)=-(2/3)z at (α,β)=(0,-2/9). RF2014 md:91; (α,β) derived algebraically from substituting into u''=(u')²/(2u)+(3/2)u³+4zu²+2(z²-α)u+β/u. Linear ⇒ Padé-Taylor exact to roundoff. ADR-0010.

Regimes: kind=:minus_2z; kind=:minus_two_thirds_z; zspan[1]=0 (u=0 ⇒ β/u singular ⇒ throw); unknown kind throw; linear ⇒ machine-precision solver match

Current coverage: test/painleve_closed_form_test.jl CF.3.1-CF.3.3, CF.4.3 — exact IC equality + solver cross-check at atol 1e-10 (linear, so exact). Mutation-proven (MC2: (α,β) swap). ALSO cross-validated in test/noumi_yamada_piv_test.jl NYPIV.3 (line 278): the A_2^(1) Noumi-Yamada α=(0,1,0) map reproduces (α,β)=(0,-2) and the solve matches u=-2z at 1e-9 — a genuine independent cross-module check. STRONG.

Gaps: The (α,β) pairs are DERIVED algebraically in the docstring (PainleveClosedForm.jl:160-173) but the derivation (esp. the -2/9 case stated as 'analogous' without shown algebra) is only validated through the solver matching the linear formula — a tautology if the same wrong β were used to both pose and check. NYPIV.3 mitigates this for :minus_2z only via the independent NY parameter map; :minus_two_thirds_z has no independent (α,β) cross-check. No BigFloat.

Ground truth needed: An independent derivation/oracle for the (0,-2/9) pair — e.g. plug u=-(2/3)z into the PIV equation symbolically (Mathematica) and confirm (α,β)=(0,-2/9) is the unique pair, mirroring what NYPIV.3 does for :minus_2z.

**_u_pii_rational / _u_pii_airy / _u_piv_entire** (oracle-helper) — `src/PainleveClosedForm.jl:76`

The exact closed-form (value, derivative) oracles backing the three closed-form families — quotient-rule derivatives in closed form for rationals, Φ'=-z/2-Φ² folded derivatives for Airy, constant derivatives for the linear PIV. They are BOTH the IC source and the test oracle (ADR-0010 'the formula IS the oracle'). Exposed via PadeTaylor.Painleve._u_*.

Regimes: u_3 at z=0 special-cased to (0,0) avoiding 0/0 (line 90-91); Airy Φ near a zero of φ (Φ→∞ — not tested); complex z; n/kind out of range throw

Current coverage: Used as the oracle in CF.1.4/CF.2.3/CF.3.3 solver cross-checks; their own correctness is checked by direct hand-computed values in CF.1.1/CF.1.2 (u_2(1)=2/5, u_2'(1)=-46/25 spelled out) and CF.2.1 (Ai-based). Self-consistency: since they serve as oracle they are mostly trusted, with anchor points hand-verified.

Gaps: _u_pii_airy n==1 derivative (Np/Dp block) has no independent hand-computed anchor (only the trusted-oracle role). _u_pii_rational u_3 derivative (lines 93-97) has no hand-computed anchor point like u_2 got. The z=0 special-case for u_3 returns (0,0) — only checked via the IC, not as a limit of the general formula. Φ→∞ near φ-zeros untested (could NaN).

Ground truth needed: A handful of hand-computed or Mathematica-computed (u,u') anchor points for u_3 and for _u_pii_airy(1,θ,z) at interior z, asserted directly on the helper (not only through the solver), to remove the oracle's self-reference for those two branches.

**tritronquee** (problem-builder) — `src/PainleveNamed.jl:92`

PI tritronquée — the u''=6u²+z solution with one pole-free sector. Carries FW2011 §4.1 IC u(0)≈-0.1875543083404949, u'(0)≈0.3049055602612289 (32-digit Maple BVP over [-20i,20i]). The IC IS the definition. ADR-0008.

Regimes: equation=:I only (others throw); zspan[1] must equal IC point 0.0 (else throw); name tag propagation through solve_pade/path_network_solve

Current coverage: test/painleve_named_test.jl NT.1.1, NT.2.1, NT.3.1 — IC asserted verbatim against the literature literal (tautological by construction but mutation-proven M1); name propagation into solutionname(sol); fail-loud guards (:II/:V throw, zspan[1]≠0 throws). Closed-form-oracle: NONE (BVP-computed, no closed form).

Gaps: No independent verification that integrating THIS IC actually yields a tritronquée (pole-free sector, correct pole locations, asymptotic u~-√(z/6)). The IC literal is both subject and oracle. phase9_tritronquee_test.jl / fw_fig_41_test.jl carry the same IC as a raw literal but do NOT call the tritronquee() constructor, so the constructor's solved behaviour is unvalidated against pole-field ground truth.

Ground truth needed: An independent oracle for the tritronquée trajectory: FW2011 published pole locations (the first pole on the negative real axis / the pole-field figure) OR a high-precision asymptotic match u(z)~-sqrt(-z/6) on the pole-free sector, integrated FROM the constructor's IC — to confirm the IC truly selects the tritronquée, not just that the digits were copied correctly.

**hastings_mcleod** (problem-builder) — `src/PainleveNamed.jl:132`

Hastings-McLeod solution of PII at α=0 (u''=2u³+zu) — pole-free, non-oscillatory on ℝ, the PII behind Tracy-Widom. FW2014 IC (±0.3670615515480784, ∓0.2953721054475501). branch selects the u→-u sign-symmetric copy. ADR-0008.

Regimes: branch=:positive (u(0)>0); branch=:negative (u(0)<0); branch invalid throw; zspan[1]≠0 throw; PII α=0 symmetry u→-u (genuine invariant); decay u~Ai(x)→0 on x>0

Current coverage: test/painleve_named_test.jl NT.1.2, NT.2.1-NT.2.3, NT.3.1 — IC verbatim (mutation-proven M3); name propagation; GENUINE cross-checks: NT.2.2 the two branches solve to exact negatives at 1e-10 (real invariant, not tautology), NT.2.3 decays 0<u(3)<u(0). Stronger than tritronquee thanks to the symmetry check.

Gaps: NT.2.3 decay check is weak (only 0<u3<u0, no quantitative Ai-asymptotic match). The IC literal is otherwise self-referential. No check that the solution stays pole-free across a longer span, nor a match to the known u(0) Tracy-Widom constant from an independent source. General-α HM is out of scope (asymptotic-BC, ADR-0008).

Ground truth needed: An independent reference for the α=0 HM connection constant (the IC is widely tabulated — e.g. u(0) relates to the Hastings-McLeod constant; an independent literature value other than FW2014 would de-tautologise the IC) and a quantitative u(x)~Ai(x) asymptotic match on x→+∞ from the constructor's solve.

**painleve_hierarchy** (problem-builder) — `src/PainleveHierarchy.jl:231`

Companion-form first-order RHS for the m-th PI-hierarchy member. m=1: y'=(y2, 6y1²+x) (PI, same normalisation as v0.1). m=2: y'=(y2,y3,y4, -10y2²-20y1y3-40(y1³-6t·y1+6x)) (P_I^(2), KKG2015 eq.(p12) u_xxxx+10u_x²+20u·u_xx+40(u³-6tu+6x)=0). m≥3 throws (Lenard recursion deferred).

Regimes: m=1 (2-vector); m=2 (4-vector), t∈ℂ parametric; Taylor1 / BigFloat element-type RHS args; equation≠:I throw; m<1 throw; m≥3 throw; wrong state-vector length throw

Current coverage: test/painleve_hierarchy_test.jl PH.1.1 (m=2 RHS == hand-written companion, t∈{0,0.7,-1.3}, plus per-component pinned values), PH.1.2 (m=1 trajectory == v0.1 PainleveProblem(:I) at 1e-9 — genuine consistency anchor), PH.1.3 (P_I^(2) ODE residual via finite-diff of solved y4, NOT tautological f[4]), PH.1.5 (all throws). Mutation-proven (M1, M2). STRONG.

Gaps: m=2 RHS only spot-checked at 3 (x,y) points per t (PH.1.1) — fine for an exact algebraic identity. The only real gap: no BigFloat/Taylor1 element-type RHS evaluation directly asserted here (relied upon by the vector solver but not unit-tested in this file). t complex (not just real ±) untested.

Ground truth needed: Already has a hand-written verbatim companion-RHS oracle (the ideal for an algebraic transcription). Marginal addition: a complex-t evaluation and a Taylor1-arg evaluation to pin the generic-T contract.

**painleve_hierarchy_jacobian** (internal-algorithm) — `src/PainleveHierarchy.jl:291`

Exact hand-derived companion-form Jacobian. m=1: [[0,1],[12y1,0]]. m=2: super-diagonal identity chain plus last row [-20y3-120y1²+240t, -20y2, -20y1, 0]. The override path so the BVP Newton solve skips per-node Taylor1 autodiff. Generic element type via zero/one(eltype(y)).

Regimes: m=1 (2×2); m=2 (4×4), t parametric; Taylor1/BigFloat element type; equation≠:I / m≥3 / m<1 throw; wrong state length throw

Current coverage: test/vector_bvp_wirein_test.jl VW.* — cross-validated against the Taylor1 autodiff Jacobian (the independent oracle: hand-derived == autodiff), tested for m=1 and m=2 across several t, with VW.3 failure-mode guards and a mutation record (line 231). STRONG (autodiff is a genuine independent oracle).

Gaps: Validated against autodiff which is itself in-repo machinery — if both shared a bug (unlikely, different code paths) it would pass. The 240t term (∂/∂y1) is only exercised at the t values used in VW.*; t complex untested. Did not read VW.* in full (out-of-slice file) — confirmed existence and cross-validation by grep only.

Ground truth needed: The autodiff cross-check is a strong oracle. A symbolic (Mathematica/SymPy) ∂f/∂y for P_I^(2) would be a fully independent third source if wanted; otherwise coverage is adequate.

**pI2_tritronquee_ic** (internal-algorithm) — `src/PainleveHierarchy.jl:515`

Tritronquée asymptotic IC for P_I^(2) (KKG2015 §7). Leading Y=-∛6·|x0|^{1/3}·e^{iθ/3} with sheet angle θ∈(2π,4π] (neg real axis at θ=3π); branch-free rational derivatives y2=Y/(3x0), y3=-2Y/(9x0²), y4=10Y/(27x0³). n_terms=2 adds the c_6=1 correction Y^{-6}=x0^{-2}/36 (t=0 only). Negative-real-axis byte-identical legacy r-power path. Fail-fast: |x0|<1, n_terms∉{1,2}, n_terms=2 & t≠0, all throw.

Regimes: negative real axis (legacy r-power path, byte-identical); off-axis V_0-sheet rays φ∈{2π/3,π,4π/3}; n_terms=1 (leading) vs n_terms=2 (c_6 correction, t=0 only); Float64 / BigFloat / Complex{BigFloat} element type; |x0|<1 throw; t≠0 with n_terms=2 throw; sign correctness (POSITIVE on x<0 — bead 0ln.31)

Current coverage: test/painleve_hierarchy_test.jl PH.1.4/PH.1.7/PH.1.8 — GROUND-TRUTH ODE residual (u'''' from closed form, NOT tautological companion f[4]), residual-decay bounds 6/|x| (n=1) and 20/|x|^4 (n=2) with res2<res1, sign check u>0, off-axis ray reconstruction at 1e-12, negative-axis continuity at 1e-12, BigFloat + Complex{BigFloat} element-type propagation. Mutation-proven (M3,M4,M5). ALSO kkg_pi2_figure_test PI2F.1.2 pins u(-15) sign+magnitude against the BVP-solved figure. The GOLD STANDARD of the slice.

Gaps: Very few. n_terms=2 c_6 correction is t=0-only by design (t≠0 throws, tested). The c_7+ corrections are not implemented (documented v0.2 corner). Sheet-angle θ boundary cases (φ0 exactly 0, or the (2π,4π] half-open boundary at the pole wedge) not explicitly tested. Only φ∈{2π/3,π,4π/3} rays sampled — the full pole-free sector edges (arg x near 3π±6π/7) untested.

Ground truth needed: Already excellent (genuine ODE-residual ground truth). Marginal: a θ-boundary test (φ0→0+ vs 0-) to pin the half-open sheet window, and a ray nearer the pole-free-sector boundary to confirm the seed degrades gracefully where KKG's asymptotics weaken.

**PainleveHierarchyProblem** (problem-builder) — `src/PainleveHierarchy.jl:380`

Self-describing wrapper around a VectorPadeTaylorProblem for a PI-hierarchy member; validates length(y0)==2m and non-degenerate xspan; forwards to vector_solve_pade (IVP) and vector_bvp_solve (BVP, defaulting jacobian to the exact analytic painleve_hierarchy_jacobian). ADR-0022.

Regimes: m=1 (2-vector); m=2 (4-vector); IVP forwarding; BVP forwarding with analytic-Jacobian default / override / nothing(autodiff); wrong y0 length throw; degenerate xspan throw; m∉{1,2} throw

Current coverage: test/painleve_hierarchy_test.jl PH.1.2 (m=1 IVP forwarding consistency), PH.1.3 (m=2 IVP), PH.1.5 (length/xspan/m throws); test/vector_bvp_wirein_test.jl VW.4 (BVP convenience, analytic-Jacobian default). Solid wiring coverage.

Gaps: BVP z_a/z_b override path (sub-segment pinning, PainleveHierarchy.jl:445-453) lightly covered — VW.4 grep shows construction but full assertion depth not read. jacobian=nothing (force autodiff) path on the wrapper not confirmed tested. t stored-as-given (Tt type param) for complex t untested.

Ground truth needed: Coverage is structural (wiring) and adequate; the underlying solver oracles live in the vector-pipeline slice. Marginal: a BVP solve with z_a/z_b overridden to a sub-segment, asserting the result matches the full-window solution on the overlap.

**Slice notes:** Cross-cutting observations for the correctness-bug hunt:

(1) HEUN PARAMETER-SET MONOCULTURE. Every single oracle-pinned Heun assertion in test/heun_test.jl uses exactly ONE parameter tuple per function: heun_g(2,1,1,2,3,4;...) and heun_c(1,1,2,3,-1;...). There is NO second parameter family pinned to an oracle. The oracle file external/probes/heun-oracle/oracles.txt is described as 42 records but the test only consumes a handful of (a=2,q=1,...) points. A coefficient-recurrence bug that happens to be inert for THIS (a,q,α,β,γ,δ) but fires for other parameters (e.g. a term that vanishes because q=1 or a=2) would pass the whole suite. The Fuchsian-constraint ε = α+β+1-γ-δ in _heun_g_rhs (src/Heun.jl:216) and _heun_g_frobenius_at_0 (line 127) is only exercised at ε = 1+2+1-3-4 = -3; an arithmetic transposition that still yields -3 on this tuple is untested.

(2) HEUN HAS NO BigFloat / COMPLEX-PARAMETER COVERAGE AT ALL. ADR-0018 locked decision #4 ("element type follows the user; BigFloat-256 fallback by passing BigFloat-typed parameters") and the promote_type machinery (src/Heun.jl:298-301, 372-374) are entirely untested — every test call passes Float64/Int real params. The `real_params` branch selector (line 318-320, 385-387) that switches between enforce_real_axis_symmetry and branch_points/branch_cut_angles is therefore only ever exercised on its `true` arm; the complex-parameter `else` arm (the branch_points=(0,1,a) walk) has ZERO test coverage.

(3) HEUNC z>=1 SHEET IS DELIBERATELY LOOSE. HC.4 (z=3) is pinned at only rtol=1e-3 against a hand-probed Mathematica upper-limit value (z+0.001i), per the CROSSVAL.md caveat that Motygin and Mathematica disagree on sheet past z=1. This is the weakest oracle in the slice and the most likely place for a sheet-selection regression to hide undetected.

(4) HEUNG EXISTENCE GUARDS PARTIALLY TESTED. HX.1 tests γ=0, γ=-1 (heun_g), a=1, and γ=0 (heun_c). NOT tested: the a=0 guard (src/Heun.jl:120-122), the heun_c non-positive-integer-γ logarithmic case for γ=-2 etc., and the n_taylor truncation adequacy (N=60 hardcoded default — no test that 60 terms suffice near |z|=epsilon_start=0.05 for stiff parameters).

(5) CLOSED-FORM FAMILIES ARE THE BEST-TESTED PART OF THE SLICE. pii_rational / pii_airy / piv_entire each have an exact-equality IC check (no tolerance — the formula IS the oracle) AND a solver cross-check; mutation-proven (4 mutations, MC1-MC4). This is genuinely strong. Gaps: pii_airy θ is only tested at θ=0 and θ=π (the two arms where one of cos/sin vanishes) — the GENERIC θ where BOTH Ai and Bi contribute (e.g. θ=π/2) is untested, and that is exactly where a cos(θ/2)/sin(θ/2) swap in _airy_phi (src/PainleveClosedForm.jl:110) would hide. pii_airy(1) derivative at z=0 is explicitly NOT asserted (CF.2.2 comment line 135-137 punts it to the solver cross-check, but CF.2.3 only cross-checks n=0, not n=1) — so the u_{3/2} derivative quotient-rule block (_u_pii_airy n==1, lines 144-148) has NO direct closed-form derivative assertion at any point.

(6) NAMED POINT-IC TRANSCENDENTS: ICs are self-referential oracles. tritronquee/hastings_mcleod assert the stored IC == the literature literal verbatim (the IC is both the thing-under-test and the oracle — a tautology by construction, but mutation-proven via M1). The genuine cross-checks are NT.2.2 (HM branches are exact negatives — real invariant) and NT.2.3 (HM decays on x>0 — weak, only checks 0 < u3 < u0). There is NO independent verification that the tritronquée IC actually PRODUCES a pole-free-sector / tritronquée solution when integrated (no pole-location check, no asymptotic match). phase9_tritronquee_test.jl and fw_fig_41_test.jl carry the same IC as a literal but (per the grep) do not call tritronquee() the constructor.

(7) PI HIERARCHY tritronquée IC (pI2_tritronquee_ic) is the MOST rigorously tested capability in the slice: genuine ODE-residual ground truth (u'''' from closed form, NOT the tautological companion f[4]), residual-decay bounds (6/|x| for n=1, 20/|x|^4 for n=2), off-axis V_0-sheet rays, BigFloat + Complex{BigFloat} element-type propagation, continuity at the negative real axis, and 5 mutations (M1-M5). This is the gold standard the rest of the slice should aspire to. The one untested corner: the n_terms=2 c_6 correction is only validated at t=0 (correctly, since t≠0 throws); the general c_n(t) series is out of scope.

(8) painleve_hierarchy_jacobian m=1/m=2 analytic Jacobian: tested in vector_bvp_wirein_test.jl VW.* against the Taylor1 autodiff Jacobian, with failure-mode guards (VW.3). Did not read that file fully (out-of-slice-file) but the grep confirms cross-validation against autodiff exists and is mutation-proven (line 231).

Oracle infra found: external/probes/heun-oracle/oracles.txt (42-record wolframscript, underutilised), external/probes/heun-oracle/CROSSVAL.md (Motygin MATLAB cross-val), and the in-source _u_pii_rational/_u_pii_airy/_u_piv_entire closed forms which double as test oracles (exposed via PadeTaylor.Painleve._u_*).

### test-audit

**padeapprox.m / Calgo-766 oracle pins (_oracles.jl, _oracle_shared_pade.jl)** (oracle-helper) — `test/_oracles.jl:1-39; test/_oracle_shared_pade.jl:1-75`

Pinned double-precision GROUND TRUTH for the scalar robust-Pade (GGT 2013 Algorithm 2 + Chebfun reweighting) AND the type-II shared-denominator Hermite-Pade. _oracles.jl: Chebfun padeapprox.m outputs (Octave 8.4.0, commit 7574c77) for 6 cases: exp(2,2), exp(20,20), log(1+z/2), tan(z^4) pole-set+coefs, 1/(1+z^2), noisy geometric. _oracle_shared_pade.jl: ACM TOMS Algo 766 (Cabay-Jones-Labahn 1997, Beckermann-Labahn striped-Sylvester) for SP_1_1 (honest deg-2 Q, coefficient-exact) and SP_1_2 (deg-5 Q with deliberate spurious common factor; compare by VALUES + pole-subset only).

Regimes: Float64 only (no BigFloat pins captured); real-coefficient AND ComplexF64 (tan(z^4) pole set); entire (exp) vs meromorphic (1/(1+z^2), tan); honest-degree Q (SP_1_1) vs spurious-common-factor Q (SP_1_2); near-tol degenerate blocks (noisy_geom tol=1e-5, mu/nu reduction)

Current coverage: robustpade_test.jl (includes _oracles.jl; rtol~1e-12, two orders looser than Chebfun internal 1e-14 for QR sign non-uniqueness); shared_pade_test.jl (includes _oracle_shared_pade.jl; SP.2/SPM testsets; published-number strength, cross-checked LIVE vs determinant + AAA oracles). Strength: published-number / independent-reference-impl.

Gaps: All padeapprox/Calgo pins are Float64 ONLY -- no BigFloat-256 padeapprox oracle exists, so the BF path of robust_pade is pinned only against the closed-form Pade(2,2) of exp and the _oracle_coefficients exp(80dps), NOT against padeapprox at high precision. tan(z^4) is the only Complex pole-set pin. No oracle for the QR-reweighting branch (padeapprox.m:278-280) in isolation -- it is only exercised implicitly.

Ground truth needed: A BigFloat/mpmath re-run of padeapprox-equivalent at 256-bit for the same 6 cases (or Chebfun in variable precision) to pin the high-precision robust_pade path against an external reference rather than only closed-form exp.

**Coefficients oracle (_oracle_coefficients.jl)** (oracle-helper) — `test/_oracle_coefficients.jl:1-99`

GROUND TRUTH for Taylor-jet generation of scalar 1st-order ODEs. Four cases cross-validated wolframscript(Mathematica) vs sympy/mpmath: exp(z) order-10 (c_k=1/k!); t*J_{3/4}/J_{-1/4} Bessel closed-form of y'=z^2+y^2 (FW 2011 sec 2.2.1) order-14; Weierstrass u''=6u^2 order-30 with g_3=2 cross-check AND independent WeierstrassP-series route; 1/k! at 80 decimal digits k=0..60 for the BigFloat-256 path.

Regimes: Float64 (order 10/14/30) AND BigFloat-256 (exp 80-dps string literals); entire (exp) vs meromorphic (Weierstrass, Bessel-ratio); two independent oracle routes for the Weierstrass case (AsymptoticDSolveValue vs WeierstrassP series)

Current coverage: coefficients_test.jl (includes _oracle_coefficients.jl). Strength: published/closed-form (1/k!, g_3=2) + cross-oracle (wolfram == sympy/mpmath, asserted by verify.jl before write). 80-dps string round-trip for BF-256.

Gaps: Only scalar 1st-order jets are oracle-pinned here. No vector-jet (VectorCoefficients) oracle file exists -- vector_coefficients_test.jl uses only closed-form exp/cos/sin inline and the d=1 reduction, so coupled-system high-order jets have no external pin. The Weierstrass jet is pinned only to order 30 at Float64; no BF jet pin for a meromorphic ODE.

Ground truth needed: A vector-jet oracle: mpmath/Mathematica Taylor coefficients for a genuinely coupled d>=2 system (e.g. the harmonic or a meromorphic Riccati) to high order, cross-validated, mirroring _oracle_coefficients.jl for the vector layer.

**StepControl oracle (_oracle_stepcontrol.jl)** (oracle-helper) — `test/_oracle_stepcontrol.jl:1-35`

GROUND TRUTH for the Jorba-Zou adaptive step selector and the pole-distance step heuristic. Three-way cross-validated: TaylorIntegration.stepsize (canonical Julia) == mpmath/sympy == Mathematica. Pins: exp Taylor coefs k=0..30 + h=4.5012... at eps=1e-12; and Polynomials.roots-based pole distances for 1/(1-z/2) (pole z=2) and 1/(1+(z-3)^2) (poles 3+-i).

Regimes: Float64 (and verify.jl claims 47-digit high-precision pin); entire (exp) vs meromorphic (real pole z=2, complex conjugate poles 3+-i); single-real-pole vs complex-conjugate-pair pole distance

Current coverage: stepcontrol_test.jl (includes _oracle_stepcontrol.jl). Strength: three-way independent-tool agreement (TI.jl primary), enforced by verify.jl. The h=4.50 value is a canonical TaylorIntegration output -- strong.

Gaps: Only Float64 step values pinned despite the header claiming a 47-digit pin (the pinned literals shown are Float64). No vector-norm step oracle (VectorStepControl uses vnorm(c_k)); vector_step_control_test.jl reduces to the scalar selector and otherwise uses hand-computed closed-form -- no external 3-way pin for the vector eps/norm selector.

Ground truth needed: A vector StepControl oracle: cross-validated min-over-k (eps/||c_k||)^(1/k) for a known coupled jet (2-norm AND inf-norm) from an independent tool, to pin vector_step_jorba_zou beyond the d=1 reduction.

**PadeStepper oracle (_oracle_padestepper.jl)** (oracle-helper) — `test/_oracle_padestepper.jl:1-64`

GROUND TRUTH for one scalar Pade-Taylor step on u''=6u^2 (Weierstrass) and u''=6u^2+z (Painleve I). Cross-validated Mathematica closed-form WeierstrassP + NDSolve(WorkingPrecision=50) == mpmath.odefun(40dps). Pins one step h=0.5, continuation h=0.4, PI step, and a near-pole step h=0.05 toward pole at z=1.

Regimes: Float64 (NDSolve WP=50 / mpmath 40dps reference); Weierstrass (no x-term) vs Painleve I (+z RHS); regular step vs near-pole step (z=0.9->0.95, pole at z=1); pole-crossing not covered here (stops at 0.95, before the pole)

Current coverage: padestepper_test.jl (includes _oracle_padestepper.jl). Strength: closed-form + dual-tool cross-check (WeierstrassP == NDSolve == mpmath). The PI step is NDSolve-only-cross-checked (no closed form), still 2-source.

Gaps: Near-pole step stops at z=0.95 -- the actual pole-CROSSING step (Pade bridging through z=1) is pinned only in _oracle_problems.jl at z=1.05 via closed-form, NOT here. No BigFloat step oracle (the BF stepper path has no pin in this file). No vector stepper oracle -- vector_stepper_test.jl uses only inline closed-form harmonic + Riccati.

Ground truth needed: A vector PadeStepper oracle: one shared-Q step on a d>=2 meromorphic system pinned to mpmath/NDSolve, and a high-precision BF step pin for the scalar stepper, to remove the all-inline-closed-form-only situation for the vector and BF step paths.

**Problems oracle (_oracle_problems.jl) -- FW Table 5.1 + pole-bridge** (oracle-helper) — `test/_oracle_problems.jl:1-87`

The single STRONGEST external pin in the suite: FW 2011 Table 5.1 PUBLISHED reference values for the Weierstrass test problem u''=6u^2, u(z)=P(z-1;0,2). Pins u(30)=1.095098255959744 (FW md:301), u(10^4)=21.02530339471055, u(28.261)=9.876953517025014e6 (high on the pole wall), each cross-validated closed-form == FW-published; plus pole-bridge demo z in {0.5,0.95,1.05,1.4} 3-source (closed-form == NDSolve == mpmath.odefun at rtol 1e-13/1e-11), with z=1.05 PAST the pole (closed-form sole primary -- mpmath cannot cross). 80-dps strings for u(30) and u(1.05) for the BF-256 round-trip.

Regimes: Float64 AND BigFloat-256 (80-dps strings for u(30), u(1.05)); short-range (z=0.5) / near-pole (0.95) / past-pole (1.05) / medium (z=30) / long-range (z=10^4) / pole-wall (z=28.261); pole-crossing IS covered (z=1.05 closed-form bridges pole at z=1); meromorphic (Weierstrass lattice poles)

Current coverage: Included by 9 test files via joinpath: problems_test.jl, adaptive_step_test.jl, diagnose_test.jl, dispatcher_test.jl, ext_commonsolve_test.jl, non_uniform_nodes_test.jl, pathnetwork_test.jl, polefield_test.jl, vector_pipeline_oracle_test.jl. Strength: PUBLISHED-NUMBER (FW Table 5.1) + multi-source closed-form/NDSolve/mpmath. The flagship oracle.

Gaps: Entirely confined to ONE ODE (equianharmonic Weierstrass u''=6u^2). The FW Table 5.1 row deferred-to-v2 note in the header is now wired. But: no FW Table 5.1 pin for any OTHER FW test problem; no Painleve-I (u''=6u^2+z) long-range published pin (PI has different asymptotics and no closed form). The vector pipeline reuses THIS SAME Weierstrass oracle (vector_pipeline_oracle_test) -- so the vector stack's only published-number oracle is a re-encoding of the scalar one, not independent.

Ground truth needed: Independent published or high-precision oracles for additional FW/FFW test problems beyond Weierstrass: e.g. a Painleve-I medium-range value table, or a second meromorphic ODE with a published reference, so the long-range/pole-bridge claims are not all resting on a single closed-form lattice.

**BVP oracle (_oracle_bvp.jl) -- DMSUITE Chebyshev** (oracle-helper) — `test/_oracle_bvp.jl:1-147`

GROUND TRUTH for the Chebyshev-collocation Newton BVP solver. Octave 8.4.0 via DMSUITE chebdif.m/chebint.m. Pins: D1/D2 differentiation matrices (N=4,8,16) bit-level; linear BVP u''=u with exact cosh(t)/cosh(1) (N=8 err 7e-10 = genuine Cheb truncation, N=16 err 6e-15); nonlinear PI BVP u''=6u^2+z on [-18,-14] (asymptote sqrt(-z/6) BCs, residual 9.7e-13) and on [0.1,0.5] (BCs from PI mpmath 40dps, endpoint deriv agrees with mpmath up(0.5)=16.204); barycentric chebint eval.

Regimes: Float64 only; linear (u''=u, exact cosh) vs nonlinear (PI u''=6u^2+z); real segment vs the [0.1,0.5] segment with mpmath BCs; Newton convergence pinned (iters + residual floor per cond(J)); D-matrix primitives pinned independent of any solve

Current coverage: bvp_test.jl (includes _oracle_bvp.jl), and referenced as cross-check in vector_bvp_test.jl Group 3. Strength: exact closed-form (cosh) for linear; mpmath-pinned BCs + endpoint-derivative cross-check for nonlinear PI; DMSUITE bit-level for D-matrices. PUBLISHED-method (DMSUITE) + closed-form.

Gaps: Float64 ONLY -- no BigFloat BVP oracle (the BF BVP path is unpinned externally). The nonlinear PI BVP is pinned by residual + endpoint-deriv, NOT by a full interior closed-form solution (PI has no closed form), so interior accuracy rests on the self-residual + mpmath BC consistency. The COMPLEX-segment BVP (genuinely off-real-axis, the actual FW Fig 4.1 use case) is NOT in this oracle -- only real segments and the misleadingly-named '_complex_N24' which is actually real [0.1,0.5].

Ground truth needed: A genuinely complex-segment BVP oracle (e.g. PI on a vertical imaginary-axis segment cross-validated against mpmath.odefun shot from a known IC), and a BigFloat-256 BVP pin, to cover the off-real-axis and high-precision BVP regimes that the figures actually exercise.

**Heun oracle (external/probes/heun-oracle, inlined in heun_test.jl)** (oracle-helper) — `external/probes/heun-oracle/README.md:1-40; test/heun_test.jl:1-45`

GROUND TRUTH for heun_g/heun_c evaluators via Mathematica HeunG/HeunC built-ins at 30 significant digits (WorkingPrecision 50, 42 records), cross-validated against the Motygin MATLAB/Octave implementation (CROSSVAL.md). The 30-digit literals are inlined into heun_test.jl as @pinned assertions.

Regimes: Float64 evaluation pinned against 30-digit Mathematica reference; Fuchsian HeunG (z<1 and z>1 sheets) and confluent HeunC; complex z (off-real-axis anchors); HeunC z>=1 deliberately on a different analytic sheet (pinned vs Mathematica HeunC[...,z+0.001i], looser tol)

Current coverage: heun_test.jl: 16 @pinned assertions, tolerances 1e-15..1e-3 (1e-3 for the sheet-ambiguous z>=1 HeunC), each with a 2s wall-time gate. Strength: published-built-in (Mathematica) + independent MATLAB (Motygin) cross-val. Strong for the pinned points.

Gaps: The z>=1 HeunC sheet is pinned at only 1e-3 (the two oracles disagree across the branch -- connection-formula work deferred). Only 16 of the 42 captured records are asserted (the rest of oracles.txt is unused). No BigFloat Heun pin despite 30-digit oracle availability -- the high-precision Heun path is unexercised.

Ground truth needed: Connection-formula reference values to pin the HeunC z>=1 second sheet tighter than 1e-3; and BigFloat assertions using the already-captured 30-digit oracle to exercise the high-precision Heun evaluation path.

**noumi_yamada_rational closed-form oracles (NoumiYamadaSymmetry)** (oracle-helper) — `test/noumi_yamada_symmetry_test.jl:1-26 (NYS.* testsets)`

SELF-VALIDATING closed-form rational solutions of the even-parity A_{2n}^{(1)} Noumi-Yamada systems (kind in {:A,:B,:C}); each is an exact rational solution that doubles as its own oracle (the constructor's solution satisfies the RHS exactly by Rational algebra). Also the affine-Weyl W(A_{2n}^{(1)}) Backlund symmetry group relations (s_i^2=id, braid, rotation pi^{2n+1}=id).

Regimes: Rational{BigInt} exact algebra (NYS.1.1 ODE satisfaction); Float64 vector-solve cross-check (NYS.1.2 Type C A_4 f_j=t/5 to ~1e-10); n>=2 even-parity systems (2n+1 components); group-theoretic (Weyl relations) -- algebraic not numerical

Current coverage: noumi_yamada_symmetry_test.jl (NYS.1.1-1.6). Strength: closed-form self-oracle (rational sol satisfies ODE exactly) + ONE vector-solve numerical cross-check (Type C, ~1e-10). The Weyl relations are exact algebra.

Gaps: The rational solutions are SELF-REFERENTIAL oracles -- they validate that the solver reproduces a solution the SAME codebase constructed, not an externally published Noumi-Yamada value. Only Type C is numerically cross-checked via vector_solve; Types A and B have NO numerical solver cross-check (only the exact-algebra ODE-satisfaction). No external (NY1998/NY1999 published) numerical value is pinned anywhere.

Ground truth needed: An independent numerical oracle for a Noumi-Yamada transcendent (non-rational) solution -- e.g. mpmath/Mathematica integration of the A_2^(1) or A_4^(1) system from a fixed IC to a downstream t -- so the vector pole-field and walk are validated against something external, not only against in-house rational solutions and the flow invariant Sigma f_j = t.

**Calogero-Moser exact-oracle test** (oracle-helper) — `test/calogero_moser_test.jl:1-55`

GENUINE exact closed-form oracle for the v0.2 vector stack: the repulsive rational N=2 Calogero-Moser pole system r''=8/r^3 (Krichever 1980 degeneration / Airault-McKean-Moser 1977). First-integral r^2 = C - 8/r^2 gives a smooth hyperbola; the EOM constant 4 is independently confirmed by RK4 integration in the test.

Regimes: Float64 vector solve along a smooth real-axis hyperbola; repulsive rational pole system (poles repel -- no collision); first-integral conserved quantity as the oracle

Current coverage: calogero_moser_test.jl. Strength: closed-form (hyperbola from first integral) + independent RK4 cross-check of the EOM constant. A real exact oracle for the vector stepper/driver/step-controller (V2->V3c).

Gaps: Smooth real-axis stretch only -- does NOT exercise pole-crossing (the poles repel and stay separated, so the shared-Q never bridges a pole here). The header itself calls this a 'smoke test' needing a smooth oracle. No complex-plane or pole-bridging Calogero-Moser regime. N=2 only.

Ground truth needed: A Calogero-Moser (or other) vector oracle where the trajectory passes NEAR/THROUGH a shared pole in the complex plane, to validate the shared-Q pole-bridging on a system with a known closed-form pole location (the v0.2 keystone claim).

**Painleve closed-form families (pii_rational, pii_airy, piv_entire)** (oracle-helper) — `test/painleve_closed_form_test.jl:1-130 (CF.* testsets)`

SELF-VALIDATING parametrised Painleve families with EXACT closed-form solutions: PII rational u_n (FW 2014 eq.6), PII Airy via the Airy phi (FW 2014 sec 2.3, cross-checked against SpecialFunctions airyai/airyaiprime/airybi), PIV entire (RF 2014, linear so exact). The constructor IC is derived from the formula, the solver output is cross-checked against the formula downstream.

Regimes: Float64 (rtol 1e-10 solver cross-check, 1e-12 IC pin); PII rational (poles at 0, -cbrt(4)) vs PII Airy (entire phi) vs PIV entire (linear); Airy branch cross-checked against SpecialFunctions (genuine external special-function oracle)

Current coverage: painleve_closed_form_test.jl (CF.1-CF.4). Strength: closed-form self-oracle for rational/entire; pii_airy IC is pinned against SpecialFunctions.airyai (EXTERNAL special-function oracle) at 1e-12 -- the strongest external pin among the Painleve transcendents. Solver cross-check at 1e-10.

Gaps: These are special (rational/Airy/entire) Painleve solutions -- the GENERIC transcendent (e.g. tritronquee PI, Hastings-McLeod PII as actual non-closed-form transcendents) is NOT pinned against any external numerical value. The solver-vs-formula cross-check at rtol 1e-10 is loose relative to the 1e-12/1e-14 used elsewhere. Only the Airy IC (not the downstream solve) touches SpecialFunctions.

Ground truth needed: A high-precision external pin (mpmath/Mathematica NDSolve at WP=50) for a GENERIC Painleve transcendent solve (tritronquee PI value at a downstream z, or Hastings-McLeod at x>0) so the actual transcendent integration -- not just rational/Airy special cases -- is validated against ground truth.

**PoleField extract_poles vs Weierstrass lattice (polefield_test)** (oracle-helper) — `test/polefield_test.jl:1-40 (PF.* testsets)`

GROUND TRUTH for pole extraction from a solved path-network: the equianharmonic Weierstrass-P test problem (FW 2011 sec 5.1.1) has ANALYTICALLY KNOWN second-order poles on a rhombic lattice z(m,n)=1+2omega(m+n/2)+i*omega*sqrt3*n with omega=Gamma(1/3)^3/(2^{13/6}pi). Extracted poles checked against this exact closed form -- no pole-finder oracle needed.

Regimes: Float64; 2D grid pole field (PF.1.2: no spurious poles, all in-region lattice poles, conjugate symmetry); single-node nearest-pole (PF.1.1) vs full field; cross-node-support filter (PF.2.2 load-bearing); degenerate network (PF.2.1)

Current coverage: polefield_test.jl (34 @test; tolerances 4.57e-7..1e-3). Strength: CLOSED-FORM lattice (exact Weierstrass pole positions). The vector analogue extract_poles_shared_q is pinned in vector_path_network_test VPN.1.1 against a d=3 Riccati exact pole. Strong closed-form oracle.

Gaps: Pinned only to the Weierstrass lattice (one ODE) and the Riccati exact pole (vector). No pole-extraction oracle for a system with KNOWN higher-order or off-lattice poles, or for the Painleve transcendents whose pole positions FW pins in Table 5.1 (those pole positions are NOT cross-checked here). The cross-node-support filter and degenerate cases are checked structurally.

Ground truth needed: A second closed-form or published pole-position oracle (e.g. FW-tabulated tritronquee pole locations, or a meromorphic system with known higher-order poles) to confirm extract_poles is not overfit to the rhombic Weierstrass lattice.

**FFW figure quantitative pins (ffw_fig_1..7) -- loose-vs-tight self-oracle** (oracle-helper) — `test/ffw_fig_1_test.jl:1-40; ffw_fig_4_test.jl:1-45; ffw_fig_6_test.jl:1-40`

Reproduction pins for the FFW 2017 Riemann-surface figures (PIII/PIV/PV/PVI on multiple sheets). The PRIMARY oracle is SELF-REFERENTIAL: a tighter adaptive_tol run (1e-12/1e-13) is treated as the reference for the figure's looser 1e-10 run; the difference is an honest tol-induced error. Secondary: comparison against FFW's PUBLISHED symmetry-method per-sheet error estimates (Table 2 / figure captions), but deliberately accepted 'within one order of magnitude' looser. Plus conjugate/Schwarz-reflection symmetry pins (FFW md:122 estimator).

Regimes: Float64 only; multi-sheet (sheets 0, +-1, +-2) -- the hardest regime; PIII three-sheet spiral, PIV/PV three-sheet, PVI eta/zeta/z planes; cross-mode vs refuse-mode branch tracking; pole-spike-adjacent Stage-2 lookups (sheet +-1 admitted as poorly resolved)

Current coverage: ffw_fig_1..7_test.jl. Strength: SELF-CROSS-CHECK (loose-vs-tight) + published-number comparison ONE ORDER looser than FFW + symmetry pins. This is the WEAKEST oracle class for a correctness-bug hunt: the absolute pins are against in-house tight-tol runs (self-referential), and the published-number acceptance is deliberately loose (sheet +-1 thresholds as loose as 8e-1).

Gaps: Sheet +-1/+-2 values are essentially UNPINNED against external ground truth (thresholds far looser than FFW's own numbers; FFW admit large sheet +-1 error). The loose-vs-tight oracle cannot catch a bug that affects BOTH the loose and tight runs identically (a shared systematic error). No external high-precision multi-sheet oracle. FF5 documents the n_terms=2 helper hard-codes a_3..a_N=0 (|Delta|~7e-3) -- a known accuracy gap accepted without external pin.

Ground truth needed: An independent high-precision multi-sheet oracle: mpmath/Mathematica integration of PIII/PV/PVI along the SAME winding path to a sheet-1 sample point, pinned to FFW's stated per-sheet error. This is the highest-value corpus-search target -- the figure tests' self-referential loose-vs-tight oracle is where a systematic multi-sheet bug could hide undetected.

**KKG/Noumi-Yamada figure pins (kkg_pi2, noumi_yamada_a4) -- structural-invariant only** (oracle-helper) — `test/kkg_pi2_figure_test.jl:1-45; test/noumi_yamada_a4_figure_test.jl:1-45`

Acceptance pins for the P_I^(2) tritronquee pole-field figure and the A_4^(1) Noumi-Yamada pole-field figure. EXPLICITLY no published machine-readable reference exists (KKG 2015 plotted surfaces only; no A_4^(1) pole-field reference). Tests assert STRUCTURAL invariants instead of pinned values: BVP convergence (residual ~1.5e-11), tritronquee SIGN+magnitude vs asymptote (6*15)^{1/3}, P_I^(2) ODE residual at interior nodes (FD vs equation), flow invariant Sigma f_j=t, pole genuineness (Q near-vanishes at extracted root), conjugate symmetry, reproducibility.

Regimes: Float64; P_I^(2) 4th-order (4-vector companion) tritronquee; A_4^(1) 5-component vector system; pole-field march (guarded by march_ok honesty flag); BVP anchor + path-network + extract_poles_shared_q pipeline

Current coverage: kkg_pi2_figure_test.jl (PI2F.1.1-1.6), noumi_yamada_a4_figure_test.jl (NYF.1.1-1.6). Strength: STRUCTURAL invariants + asymptote-magnitude check + self-residual. NO external numerical value pinned (none exists). The tritronquee sign check (PI2F.1.2) guards a real prior sign bug (bead 0ln.31).

Gaps: These are the LEAST externally-anchored capability tests in the suite -- by necessity (no published data). Pole positions are validated only by self-consistency (Q vanishes there) + conjugate symmetry, NOT against any external pole location. The P_I^(2) interior ODE residual is a self-check (FD of the same trajectory). A bug that produces a self-consistent but WRONG pole field would pass. The march is guarded by march_ok, so a march failure does not RED the suite.

Ground truth needed: A high-precision independent integration of the P_I^(2) tritronquee (e.g. mpmath Taylor-method along the real axis from a verified asymptotic IC) to pin at least one downstream u(x) value, and an independent A_4^(1) integration -- to convert these from structural-only to value-pinned. This is the top corpus-search target for the vector/higher-order pole-field claims.

**eta-plane PVI (FFW published u(10)) usage** (oracle-helper) — `test/eta_pvi_test.jl:1-26, 108-138, 246-269`

FFW 2017 PUBLISHED PVI reference: (alpha,beta,gamma,delta)=(4,-4,8,-8), u(10)=0.429534600325223, u'(10)=-1.61713114374804e-3 (FFW md:195). The eta-plane PVI RHS (FFW eq.5) and z<->eta IC helpers are pinned, plus the branch-point-free region Re eta < log(2pi).

Regimes: Float64 / Complex; PVI eta-plane vs zeta-plane (:transformed) frames; branch-point-free region check (finite inside Re eta<log2pi, blows up outside); RHS hand-pin at vanishing-parameter point (closed-form inline)

Current coverage: eta_pvi_test.jl (ET.1.1-1.5). Strength: RHS hand-pin closed-form (1e-14); IC round-trip machine precision; direct-vs-eta-transformed step agreement (1e-10/1e-9). The FFW published u(10) value is used at lines 112/249.

Gaps: CRITICAL: the FFW published u(10)=0.4295... value is used ONLY as an IC SEED and a frame-wiring check (lines 112, 249) -- it is NEVER asserted against an independent end-to-end PVI SOLVE. So the one available external published PVI number is not exercised as a solver oracle. ET.1.3 only checks direct-vs-transformed CONSISTENCY (self-referential), not absolute accuracy.

Ground truth needed: An end-to-end PVI solve from a DIFFERENT IC integrated TO z=10 (or eta) asserting u(10)=0.429534600325223 against the FFW published value -- converting the available published number from a seed into a genuine solver oracle. mpmath/NDSolve could also provide a second PVI point.

**Cross-validating verify.jl harness (per-oracle probes)** (oracle-helper) — `test/_oracle_coefficients.jl:1-10; test/_oracle_problems.jl:1-18; shared_pade_test.jl:14-40`

The META-oracle infrastructure that ENFORCES multi-source agreement BEFORE pinning: each of coefficients/problems/padestepper/stepcontrol oracles has a verify.jl that asserts wolfram==python(==julia) within stated rtol and only then emits the pinned _oracle_*.jl. Plus the LIVE cross-oracles in shared_pade_test (determinant via Mano-Tsuda 2017 bordered block-Toeplitz, exact on Rational{BigInt}; AAA via Nakatsukasa-Sete-Trefethen 2018).

Regimes: multi-tool agreement enforcement (Mathematica/mpmath/sympy/TaylorIntegration); LIVE vs PINNED (determinant+AAA run live in-suite; Calgo/wolfram pinned as literals); Rational{BigInt} bit-exact route (determinant oracle)

Current coverage: verify.jl scripts are offline capture-time gates (not run in the suite). shared_pade_test runs the determinant + AAA oracles LIVE during testing -- the only test that executes independent reference implementations rather than comparing to frozen literals. Strength: the strongest verification discipline in the repo.

Gaps: verify.jl agreement is a CAPTURE-TIME check -- if an oracle file is hand-edited or drifts, the suite would not re-detect disagreement (only the pinned literals are tested). The LIVE cross-oracle pattern (determinant/AAA) exists ONLY for shared_pade -- no other capability runs an independent impl in-suite. Several probes (pathnetwork-long-range oracles.txt) are EMPTY (capture not run) -- so PN.2.3 long-range BF claim rests on worklog numbers, not a regenerated oracle.

Ground truth needed: Extend the live-cross-oracle pattern (as in shared_pade) to the vector stepper and BVP -- run an independent in-suite reference (e.g. a small mpmath-equivalent or RK reference) rather than only comparing to frozen Float64 literals; and regenerate the empty pathnetwork-long-range BF-256 oracle so the long-range claim has a captured pin.

**Slice notes:** ORACLE INVENTORY (what is pinned, at what precision, from which source):\n\nSTRONG external/published oracles (a bug would be caught):\n- _oracle_problems.jl: FW 2011 Table 5.1 PUBLISHED numbers u(30), u(10^4), u(28.261) + closed-form Weierstrass-P, 3-source (closed-form/NDSolve/mpmath), Float64 AND BF-256 (80-dps strings). THE flagship oracle. But confined to ONE ODE (u''=6u^2).\n- _oracle_coefficients.jl: wolfram==sympy/mpmath cross-validated Taylor jets, Float64 + BF-256 (exp 80dps). Scalar only.\n- _oracle_stepcontrol.jl: 3-way (TaylorIntegration==mpmath==Mathematica). Scalar only.\n- _oracle_padestepper.jl: closed-form WeierstrassP==NDSolve==mpmath. Scalar, stops BEFORE pole-cross.\n- _oracle_bvp.jl: DMSUITE Octave + closed-form cosh + mpmath PI BCs. Float64, REAL segments only.\n- _oracle_shared_pade.jl: Calgo-766 FORTRAN (pinned) + LIVE determinant (Mano-Tsuda, exact on Rational) + LIVE AAA. The ONLY capability with independent impls run IN-SUITE.\n- _oracles.jl: Chebfun padeapprox.m (Octave). Float64 only.\n- heun-oracle: Mathematica HeunG/HeunC 30-digit + Motygin MATLAB cross-val. 16 of 42 records asserted.\n- polefield_test: closed-form Weierstrass rhombic pole LATTICE. Strong.\n- calogero_moser: closed-form hyperbola (first integral) + RK4. Smooth real-axis only, NO pole-cross.\n- painleve_closed_form: SpecialFunctions.airy external pin for the pii_airy IC (1e-12).\n\nWEAK / SELF-REFERENTIAL oracles (where bugs can hide -- TOP corpus-search targets):\n1. FFW figure tests (ffw_fig_1..7): PRIMARY oracle is loose-vs-tight SELF-CROSS-CHECK (a tighter in-house run is the 'reference'); published-number comparison deliberately ONE ORDER looser than FFW; sheet +-1/+-2 thresholds as loose as 8e-1. A systematic multi-sheet bug present in BOTH loose and tight runs would pass undetected. NO external high-precision multi-sheet oracle anywhere.\n2. KKG P_I^(2) + Noumi-Yamada A_4^(1) figure tests: STRUCTURAL-INVARIANT ONLY (no published data exists). Pole positions validated only by self-consistency (Q vanishes there) + symmetry; a self-consistent-but-wrong pole field passes. Pole-field march guarded by march_ok so a march failure does not RED the suite.\n3. eta_pvi: the ONE available FFW PUBLISHED PVI number u(10)=0.429534600325223 is used ONLY as an IC seed / frame-wiring check -- NEVER asserted against an end-to-end solve. Wasted external oracle.\n4. noumi_yamada_rational: SELF-VALIDATING (in-house rational sol vs solver); only Type C has a numerical cross-check; no external NY-published numerical value anywhere.\n\nSTRUCTURAL GAP -- the VECTOR (v0.2) stack: its ONLY published-number oracle (vector_pipeline_oracle_test) is a RE-ENCODING of the scalar Weierstrass oracle, not independent. No vector-jet oracle file, no vector-stepper external pin (inline closed-form harmonic/Riccati only), no vector-StepControl 3-way pin, no complex/pole-bridging Calogero-Moser. The shared-Q pole-BRIDGING keystone claim (the v0.2 headline) is exercised against closed-form poles (Riccati d=3 in VPN.1.1, Weierstrass) but NOT against any externally-published vector transcendent value.\n\nBigFloat-256 gap: external BF pins exist ONLY for exp coefficients (80dps) and the Weierstrass problem (u(30)/u(1.05) 80dps). NO BF oracle for padeapprox, BVP, the stepper step, Heun, or any vector capability -- the high-precision paths of most modules are validated only by self-consistency or Float64 literals.\n\nPole-CROSSING gap: the padestepper oracle stops at z=0.95 (before pole z=1); the actual pole-bridging step is pinned only via closed-form at z=1.05 in _oracle_problems.jl. No external oracle exercises a step that lands a Pade approximant straddling a pole for a NON-Weierstrass system.\n\nverify.jl harnesses enforce wolfram==python(==julia) at CAPTURE time only -- not re-run in-suite; the live-independent-impl pattern (determinant/AAA) exists ONLY for shared_pade. pathnetwork-long-range/oracles.txt is EMPTY (BF-256 capture never run) -- the PN.2.3 long-range BF claim rests on worklog numbers, not a regenerated oracle.\n\nINCLUDE-PATTERN NOTE: only 4 oracle helpers (_oracles, _oracle_coefficients, _oracle_padestepper, _oracle_stepcontrol) are included by the simple include(\"_oracle_*.jl\") form; _oracle_problems (9 consumers), _oracle_bvp, _oracle_shared_pade are included via joinpath(@__DIR__,...). heun oracle is inlined as literals, not a helper file.
