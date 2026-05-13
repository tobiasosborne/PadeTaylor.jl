# docs/make.jl — Documenter build script for PadeTaylor.jl.
#
# Quality gates run locally per CLAUDE.md Rule 11; no `deploydocs` call.
# To build:
#
#     julia --project=docs docs/make.jl
#
# Self-bootstraps a develop-edge of PadeTaylor.jl into the docs env on first
# run; subsequent runs are no-ops on the develop entry.

using Pkg
Pkg.activate(@__DIR__)

let parent_path = joinpath(@__DIR__, ".."),
    padetaylor_uuid = Base.UUID("00000000-0000-0000-0000-000000000001"),
    project = Pkg.project()

    # Develop PadeTaylor into the docs environment if it isn't already.
    # On first run the manifest is empty; `Pkg.develop` is idempotent so
    # we always call it once before `instantiate`.
    if !haskey(project.dependencies, "PadeTaylor")
        Pkg.develop(path = parent_path)
    end
    Pkg.instantiate()
end

using Documenter
using PadeTaylor

DocMeta.setdocmeta!(PadeTaylor, :DocTestSetup, :(using PadeTaylor); recursive = true)

# All sub-modules whose docstrings should appear in the API reference.
# Order matches the four-layer architecture (ADR-0001), then the higher
# composition tiers, then the helpers.
const MODULES = [
    PadeTaylor,
    PadeTaylor.LinAlg,
    PadeTaylor.RobustPade,
    PadeTaylor.Coefficients,
    PadeTaylor.StepControl,
    PadeTaylor.PadeStepper,
    PadeTaylor.Problems,
    PadeTaylor.PathNetwork,
    PadeTaylor.BVP,
    PadeTaylor.Dispatcher,
    PadeTaylor.EdgeDetector,
    PadeTaylor.LatticeDispatcher,
    PadeTaylor.CoordTransforms,
    PadeTaylor.SheetTracker,
]

makedocs(
    sitename = "PadeTaylor.jl",
    authors  = "Tobias Osborne",
    modules  = MODULES,
    format   = Documenter.HTML(
        prettyurls = false,
        sidebar_sitename = false,
        canonical = nothing,
        edit_link = nothing,
        repolink = nothing,
    ),
    pages = [
        "Home"         => "index.md",
        "Architecture" => "architecture.md",
        "API"          => "api.md",
        "Figures"      => "figures.md",
    ],
    checkdocs = :none,
    warnonly = [:missing_docs, :cross_references, :docs_block],
)
