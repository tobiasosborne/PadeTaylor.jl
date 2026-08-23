using Test
using PadeTaylor

@testset "PadeTaylor.jl" begin
    @testset "umbrella loads" begin
        # Confirms the umbrella module + its sub-modules (count asserted below) load cleanly.
        @test isdefined(PadeTaylor, :LinAlg)
        @test isdefined(PadeTaylor, :RobustPade)
        @test isdefined(PadeTaylor, :SharedPade)
        @test isdefined(PadeTaylor, :Coefficients)
        @test isdefined(PadeTaylor, :VectorCoefficients)
        @test isdefined(PadeTaylor, :StepControl)
        @test isdefined(PadeTaylor, :PadeStepper)
        @test isdefined(PadeTaylor, :VectorStepper)
        @test isdefined(PadeTaylor, :VectorStepControl)
        @test isdefined(PadeTaylor, :VectorProblems)
        @test isdefined(PadeTaylor, :VectorBVP)
        @test isdefined(PadeTaylor, :NoumiYamada)
        @test isdefined(PadeTaylor, :NoumiYamadaSymmetry)
        @test isdefined(PadeTaylor, :PainleveHierarchy)
        @test isdefined(PadeTaylor, :VectorPathNetwork)
        @test isdefined(PadeTaylor, :VectorPoleField)
        @test isdefined(PadeTaylor, :Problems)
        @test isdefined(PadeTaylor, :PathNetwork)
        @test isdefined(PadeTaylor, :BVP)
        @test isdefined(PadeTaylor, :Dispatcher)
        @test isdefined(PadeTaylor, :EdgeDetector)
        @test isdefined(PadeTaylor, :LatticeDispatcher)
        @test isdefined(PadeTaylor, :CoordTransforms)
        @test isdefined(PadeTaylor, :Laplace2D)
        @test isdefined(PadeTaylor, :SheetTracker)
        @test isdefined(PadeTaylor, :BranchTracker)
        @test isdefined(PadeTaylor, :Painleve)
        @test isdefined(PadeTaylor, :IVPBVPHybrid)
    end

    # Per-module test files. Each file is a self-contained @testset.
    # Phases 1–6 + Phase 10 PathNetwork + Phase 11 BVP + Phase 12
    # Dispatcher fill in these.
    include("linalg_test.jl")
    include("robustpade_test.jl")
    include("shared_pade_test.jl")
    include("shared_pade_dispatch_test.jl")
    include("classical_pade_test.jl")
    include("coefficients_test.jl")
    include("vector_coefficients_test.jl")
    include("stepcontrol_test.jl")
    include("padestepper_test.jl")
    include("vector_stepper_test.jl")
    include("vector_step_control_test.jl")
    include("vector_problems_test.jl")
    include("noumi_yamada_test.jl")
    include("noumi_yamada_symmetry_test.jl")
    # Metamorphic-relation test layer (oracle-free invariants), bead
    # padetaylor-krgy.3.  Symmetry MRs (conjugate, PII α-negation, PIV parity)
    # + consistency MRs (step additivity, forward∘reverse).  The canonical
    # algebraic NY Bäcklund MR lives in noumi_yamada_symmetry_test.jl (NYS.1.6).
    include("metamorphic_symmetry_test.jl")
    include("metamorphic_consistency_test.jl")
    include("noumi_yamada_piv_test.jl")
    include("painleve_hierarchy_test.jl")
    include("vector_path_network_test.jl")
    include("vector_path_network_stage2_test.jl")
    include("calogero_moser_test.jl")
    include("problems_test.jl")
    include("corpus_scalar_pole_bridge_test.jl")
    include("corpus_robust_pade_test.jl")
    include("corpus_taylor_jet_test.jl")
    include("corpus_vector_test.jl")
    include("corpus_vector_polefield_test.jl")
    include("corpus_bvp_test.jl")
    include("corpus_bvp_hybrid_test.jl")
    include("corpus_morphology_test.jl")
    include("corpus_multisheet_test.jl")
    include("corpus_winding_test.jl")
    include("corpus_painleve_rational_test.jl")
    include("corpus_noumi_yamada_test.jl")
    include("pathnetwork_test.jl")
    include("branch_tracker_test.jl")
    include("path_network_branch_test.jl")
    include("sheet_aware_stage2_test.jl")
    include("extrapolate_test.jl")
    include("adaptive_step_test.jl")
    include("non_uniform_nodes_test.jl")
    include("polefield_test.jl")
    include("bvp_test.jl")
    include("dispatcher_test.jl")
    include("edge_detector_test.jl")
    include("laplace2d_test.jl")
    include("phase9_tritronquee_test.jl")
    include("lattice_dispatcher_test.jl")
    include("edge_gated_solve_test.jl")
    # Seed-invariance gate for the FW 2011 Fig 4.7 "seam" (bead padetaylor-sny7,
    # P0 epic padetaylor-vwgl).  Guards that windowed_path_network_solve
    # (src/WindowedComposite.jl) produces a SEED-INVARIANT pole field where the
    # monolithic path_network_solve does not.  FSEAM.1 pins NON-TRIVIALITY (the
    # two global seeds genuinely re-randomised every window — the anti-gaming
    # guard) so the FSEAM.2 pole-set seed-invariance (≥0.90 both ways; measured
    # mono 77.6% RED vs composite 97.2% GREEN) certifies real path-independence,
    # not a frozen pipeline.  See docs/worklog/078-fig47-seam-cure-scoping.md.
    include("field_seam_test.jl")
    # Public-API kwarg deprecation shims (beads xds, 0xn; api-review §3(a)).
    # Asserts each renamed kwarg's deprecated alias still WORKS (same result)
    # AND WARNS.  The WARNS half is meaningful under Pkg.test's `--depwarn=yes`.
    include("api_kwarg_deprecation_test.jl")
    include("ext_commonsolve_test.jl")
    include("ext_arblib_test.jl")
    include("ext_makie_test.jl")
    # Certified ball oracles (Arblib midpoint-radius), bead padetaylor-krgy.5.
    # Tier-2 §2.2: certifies the closed-form Weierstrass-℘ ORACLE as a
    # provably-tiny Arb ball (radius < ε_cert), THEN pins the package's
    # solve_pade value to that ball's midpoint at the package tolerance
    # (tol_pkg ≫ ε_cert).  Certifies the oracle side, NOT the SVD output
    # (Arblib has no SVD; ADR-0002).  See docs/adr/0029-certified-ball-oracles.md.
    include("certified_oracle_test.jl")
    # Live differential / back-to-back testing (bead padetaylor-krgy.14).
    # Tier-2 §2.3: integrates manufactured smooth (pole-free) IVPs with BOTH
    # solve_pade AND an in-file classical RK4 sharing no code with src/, then
    # asserts the two LIVE codes agree in value AND derivative (Hatton:
    # independent codes disagree more than expected).  Distinct from the
    # certified-ball oracle above (which BOUNDS the true value) — here the
    # invariant is AGREEMENT with a second, different-method code path.  No
    # new dependency (RK4 is hand-rolled), so this runs standalone too.
    include("differential_test.jl")
    # Convergence-order / MMS / observed-order-of-accuracy V&V (bead
    # padetaylor-krgy.6).  Tier-2 §2.4: CODE verification (Roy 2005) — refine
    # the discretisation and confirm the error decays at the THEORETICAL rate,
    # distinct from value pinning.  CV.1 asserts the IVP single-step observed
    # order tracks the diagonal-Padé law 2·(order÷2)+1 (NOT the textbook
    # h^(N+1) — the stepper evaluates an (m,m) Padé, m=order÷2, not the order-N
    # Taylor polynomial); CV.2 asserts exponential (spectral) convergence of
    # the Chebyshev BVP until the conditioning floor.  Catches order-
    # degradation bugs a single-resolution value test passes.
    include("convergence_test.jl")
    # Accuracy-regression ledger (bead padetaylor-krgy.8).  Tier-2 §2.5: stores
    # each representative case's BEST ACHIEVED relERR in a committed,
    # human-reviewable baseline (test/accuracy_baseline.toml) and asserts
    # achieved ≤ baseline·(1+slack).  Catches SILENT accuracy DRIFT a fixed
    # absolute-tol value pin passes (the mutation-proof shows a degradation the
    # corpus z=1.4 rtol-1e-7 pin misses but the ledger catches), and surfaces
    # IMPROVEMENT (advisory to tighten + re-approve the baseline).  Complements,
    # never replaces, the absolute-tol pins.  Local no-CI baseline; the slack
    # band absorbs cross-platform float jitter.  Regeneration is explicit +
    # reviewed (ties to krgy.9 snapshot/approval).
    include("accuracy_ledger_test.jl")
    # Oracle PROVENANCE contract (bead padetaylor-krgy.9).  Tier §2.6: asserts
    # every test/_oracle*.jl data file carries a non-trivial header declaring
    # BOTH a SOURCE (ground-truth class) and a REGENERATION reference (capture
    # script / regen command) — the executable companion to the snapshot-
    # approval runbook docs/test_corpus/04_snapshot_approval_workflow.md.  A
    # golden whose source/regeneration is undocumented turns the suite RED.
    include("oracle_provenance_test.jl")
    # Static-analysis + lint quality gate (Aqua/JET/ExplicitImports), bead
    # padetaylor-krgy.2.  Infra-tier: asserts package hygiene + abstract-
    # interpretation cleanliness, not a numerical invariant.
    include("quality_test.jl")
    # Allocation-budget + type-stability gate (AllocCheck/JET report_opt/@inferred),
    # bead padetaylor-krgy.12.  EFFICIENCY-tier: asserts the per-step inner loop is
    # type-stable (no runtime dispatch / Any) and within an O(1)-per-step allocation
    # budget — the upstream, deterministic causes of the wall-time the krgy.11 perf
    # gate measures downstream.  Shares the JET setup with quality_test.jl (krgy.2).
    include("type_stability_test.jl")
    # Property-based invariant gate (Supposition.jl), bead padetaylor-krgy.1.
    # Tier-1 §1.2: generates hundreds of inputs per property and shrinks
    # failures to a minimal counterexample; pins INVARIANTS where the corpus
    # pins VALUES.
    include("property_test.jl")
    include("diagnose_test.jl")
    include("fw_fig_41_test.jl")
    include("coord_transforms_test.jl")
    include("sheet_tracker_test.jl")
    include("eta_pvi_test.jl")
    include("painleve_test.jl")
    include("painleve_solution_test.jl")
    include("painleve_named_test.jl")
    include("painleve_closed_form_test.jl")
    include("ffw_fig_6_test.jl")
    include("ffw_fig_1_test.jl")
    include("ffw_fig_4_test.jl")
    include("ivp_bvp_hybrid_test.jl")
    include("ffw_fig_5_test.jl")
    include("ffw_fig_2_test.jl")
    include("ffw_fig_3_test.jl")
    include("ffw_fig_7_test.jl")
    include("noumi_yamada_a4_figure_test.jl")
    include("kkg_pi2_figure_test.jl")
    include("vector_bvp_test.jl")
    include("vector_bvp_wirein_test.jl")
    include("corpus_vector_bvp_test.jl")
    include("vector_pipeline_oracle_test.jl")
    include("heun_test.jl")
    include("corpus_heun_test.jl")
    # Stage 4c of the P1 Heun epic (bead padetaylor-4lh): Schwarzschild QNM
    # FREQUENCY via Leaver's radial continued fraction, pinned to the
    # cross-verified Leaver 1985 targets in docs/teukolsky_heun_mapping.md §3.2.1.
    include("teukolsky_qnm_test.jl")
    include("corpus_special_fn_test.jl")
    include("corpus_riccati_rational_test.jl")
    include("corpus_riccati_special_test.jl")
    include("corpus_periodic_pole_test.jl")
    include("corpus_orthopoly_bvp_test.jl")
    include("corpus_special_fn_bvp_test.jl")
    include("corpus_oblique_guard_test.jl")
    include("corpus_higher_order_pole_test.jl")
    include("corpus_elliptic_lattice_test.jl")
    include("corpus_elementary_branch_test.jl")
    include("corpus_pathnet_walls_rows_test.jl")
    include("corpus_pathnet_lattice_sectors_test.jl")
    include("corpus_pathnet_winding_test.jl")
    include("corpus_algebraic_pvi_test.jl")
    include("corpus_out_of_class_test.jl")
end
