# DESIGN.md §6 — Path-network and beyond (Tier 2–5 phase table)

> Drop-in fragment for `DESIGN.md §4`.  Every algorithmic claim cites
> `references/markdown/<paper>/<file>.md:<lines>` or `docs/<doc>.md:§`.

| phase | name | status | file(s) | LOC | key invariant | test plan summary | dependencies |
|---|---|---|---|---|---|---|---|
| **10** | `PathNetwork` | Planned | `src/PathNetwork.jl` (+ `PathNetworkStage2.jl` if split) | ≤200 | Stage-1 visited map covers all grid nodes or throws; Stage-2 is read-only Padé eval — no new Taylor jets | 5 testsets (PN.1.1–PN.3.1); mutation-proven; PN.2.2 = FW Table 5.1 z=30 to ≤1e-13 (`FW2011...md:385–391`) | `PadeStepper`, `Problems.PadeTaylorProblem`; ADR-0004 |
| **11** | `BVPSmooth` | Planned (TBD interface) | `src/BVPSmooth.jl` | ≤200 | IVP/BVP derivative match ≤1e-7 abs at junction (FW 2011 §4.4, `FW2011...md:249–261`) | FW Fig 4.6 acceptance: BVP z∈[-18,-14] vs `√(-z/6)` asymptotic ≤1e-8 abs at endpoints | Phase 10, `EdgeDetector` (`padetaylor-c2p`); bead `padetaylor-804` |
| **12** | `Dispatcher` | Planned | `src/Dispatcher.jl` | ≤100 | No discontinuities at IVP/BVP boundary; edge threshold 0.001 user-tunable (`docs/unified_path_network_spec.md §12` gap 4) | Integration test: PI tritronquée full figure matches FW Fig 4.1 (`figure_catalogue.md §1`) | Phases 10 + 11; bead TBD |
| **13** | `CoordTransforms` | Planned | `src/CoordTransforms.jl` | ≤100 | `z = e^(ζ/2)` (PIII) and `z = e^ζ` (PV) round-trip exact; multi-sheet partition by 4π/2π strip (FFW 2017 §2.1–§2.2, `FFW2017...md:40–48`) | FFW 2017 Fig 1, 4, 5, 6 pole-pattern matches FFW (`figure_catalogue.md §5`) | Phase 10; wraps `PathNetwork` externally |
| **14** | `SheetTracker` | Planned | `src/SheetTracker.jl` | ≤150 | Sheet index `s ∈ ℤ`, initial s=0; cross-sheet eval forbidden without `cross_branch=true`; branch-cut increments logged (FFW 2017 §2.3, `FFW2017...md:163–189`) | FFW 2017 Fig 2, 3, 7 multi-sheet results match FFW (`figure_catalogue.md §5`) | Phase 13; wraps `PathNetwork` + `CoordTransforms` |

**Bead `padetaylor-1jf`** (Phase 10) is the unblocking P0; all Tier-3–5
phases depend on it.  Phases 11–14 are Tier 3–5 work items with P2
priority pending Phase 10 GREEN.

**Tier-3 interface point** in `src/PathNetwork.jl`:
```julia
# TIER-3 INTERFACE: BVP dispatcher reads PathNetworkSolution.{grid_z,
# grid_u, grid_up}. No code here until Phase 11 (bead padetaylor-804).
```

**No parallel Julia agents** across Phases 10–14 (CLAUDE.md Rule 7).
