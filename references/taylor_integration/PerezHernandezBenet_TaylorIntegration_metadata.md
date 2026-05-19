# TaylorIntegration.jl — software metadata

- **Repository**: https://github.com/PerezHz/TaylorIntegration.jl
  (Note: the scope doc listed two candidate URLs; `PerezHernandez93/TaylorIntegration.jl`
  returns HTTP 404 and does not exist; `JuliaDiff/TaylorIntegration.jl` returns HTTP 404.
  The canonical repository is `PerezHz/TaylorIntegration.jl`, confirmed via Zenodo
  record DOI 10.5281/zenodo.2562353 which links to
  https://github.com/PerezHz/TaylorIntegration.jl/tree/v0.4.1 — verified 2026-05-19.)
- **Authors**:
  - Jorge A. Pérez Hernández (PerezHz), Minor Planet Center, Center for Astrophysics,
    Harvard & Smithsonian
  - Luis Benet, Instituto de Ciencias Físicas, Universidad Nacional Autónoma de México
    (UNAM)
  - Additional contributors: Marcelo Forets, Christopher Rackauckas, Blas Kolic,
    and others (visible in Zenodo v0.13.0 record)
- **Citation DOI** (Zenodo concept DOI): 10.5281/zenodo.2562352
  (This is the concept DOI covering all versions; the version DOI for v0.4.1 cited
  in the scope doc is 10.5281/zenodo.2562353. The current most-cited Zenodo record
  appears to be the concept DOI `2562352` — README badge links to this.)
- **Current version** (as of 2026-05-19): **v0.18.4**, released 2026-05-02
  (confirmed via GitHub API: latest release tag `v0.18.4`, published_at
  2026-05-02T04:22:35Z; last push 2026-05-17T02:05:49Z)
- **License**: MIT ("Expat" license; confirmed in repository README and GitHub API)
  (Note: GitHub API returns `NOASSERTION` for `spdx_id` — this is a known GitHub
  artefact when the license file is non-standard; the README explicitly states MIT.)
- **Last commit date** (main branch, GitHub API 2026-05-19): 2026-05-17
  (commit sha 68334059, message "Update README.md")
- **Stars** (GitHub API 2026-05-19): 139
- **Relevance to v0.2** (≤ 3 sentences): PRD line 311 explicitly names
  `jetcoeffs!` and `Vector{Taylor1{T}}` as the substrate for the v0.2 vector /
  multi-component Taylor jet infrastructure; TaylorIntegration.jl is the package
  that ships this interface and is the Julia reference implementation of the
  Jorba–Zou Taylor-method integrator that PadeTaylor already uses in v0.1. The
  v0.2 extension to vector Painlevé-type systems requires generating Taylor jets
  for systems of ODEs (not just scalars), and TaylorIntegration.jl's `jetcoeffs!`
  mechanism is the intended substrate per the PRD. No standalone arXiv paper exists;
  the Zenodo DOI is the citable form; the JuliaCon 2017 slides (linked from README)
  give the best overview of the package's design philosophy.
- **Version range PadeTaylor v0.2 targets**: TBD (set during ADR work). Current
  release is v0.18.4; v0.1 of PadeTaylor targets an earlier version — the ADR
  should pin the minimum version that ships the `jetcoeffs!` interface for
  `AbstractArray` initial conditions, which was introduced in the v0.6–v0.8 range.
- **Coverage limitations** relevant to PadeTaylor v0.2:
  - TaylorIntegration.jl integrates systems of ODEs using Taylor's method to
    arbitrary order; it is not specialised to Painlevé equations.
  - The `jetcoeffs!` interface requires the ODE RHS to be written in a form that
    TaylorSeries.jl can differentiate symbolically — this may require wrapping
    Noumi–Yamada RHS terms (which involve quotients) carefully to avoid division-
    by-Taylor-series issues at poles.
  - Performance for high-order Taylor jets (order 30+, the v0.1 default) scales
    linearly in the system dimension for vector ODEs; no known bottleneck for
    n ≤ 10 (the maximum Noumi–Yamada dimension at A_9^{(1)}).

> Stub placeholder for code-side citations.
> Cite as: Pérez Hernández, J.A. and Benet, L., *TaylorIntegration.jl*,
> Zenodo, DOI 10.5281/zenodo.2562352, https://github.com/PerezHz/TaylorIntegration.jl.
