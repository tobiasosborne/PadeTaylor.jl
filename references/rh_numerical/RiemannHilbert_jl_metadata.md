# RiemannHilbert.jl — software metadata

- **Repository**: https://github.com/JuliaHolomorphic/RiemannHilbert.jl
  (Note: the scope doc and PRD line 408–410 list `JuliaApproximation` as the GitHub
  organisation; the current canonical repository is under `JuliaHolomorphic`. The
  `JuliaApproximation/RiemannHilbert.jl` path does not exist (GitHub API returns the
  `JuliaHolomorphic` org entry when queried). The README's Gitter badge links to
  `JuliaApproximation/ApproxFun.jl`, which explains the org confusion — the package
  is in the `JuliaHolomorphic` org but closely associated with the `JuliaApproximation`
  ecosystem. Verified 2026-05-19.)
- **Authors**: Sheehan Olver (primary author; package cited as implementing the
  numerical method of Olver 2011 and Olver 2012); affiliated contributors unlisted
  in README.
- **Citation DOI** (Zenodo / JOSS): none found. The package has no JOSS paper or
  Zenodo DOI as of 2026-05-19. Cite via GitHub URL and the companion papers:
  - Olver, S. (2011). "Numerical solution of Riemann–Hilbert problems: Painlevé II."
    *Found. Comput. Math.* 11: 153–179. DOI 10.1007/s10208-010-9079-8
  - Olver, S. (2012). "A general framework for solving Riemann–Hilbert problems
    numerically." *Numer. Math.* 122: 305–340. DOI 10.1007/s00211-012-0459-7
  - Trogdon, T. & Olver, S. (2016). *Riemann–Hilbert Problems, Their Numerical
    Solution and the Computation of Nonlinear Special Functions*, SIAM OT146.
    (See `references/rh_numerical/TrogdonOlver2016_RH_OT146_metadata.md`)
- **Current version** (as of 2026-05-19): **v0.1.0** (the only tagged release;
  released 2019-12-09, confirmed via GitHub API `pushed_at` field and WebFetch of
  the repository page)
- **License**: MIT (confirmed via GitHub API `license.spdx_id: MIT`)
- **Last commit date** (master branch, GitHub API 2026-05-19): 2024-12-19
  (commit sha 8f08c0a, message "Update .gitignore"; this is a maintenance-only
  commit with no functional change)
- **Last push to remote** (GitHub API `pushed_at`): 2025-08-15
  (Note: `pushed_at` 2025-08-15 but latest commit date 2024-12-19 — this discrepancy
  may reflect a force-push or branch operation that did not change the main commit
  history; functionally the package has been dormant since the v0.1.0 release in
  December 2019)
- **Stars** (GitHub API 2026-05-19): 5
- **Language composition**: Wolfram Language 58.8%, Julia 41.2%
  (the Wolfram Language component is likely the RHPackage predecessor code)
- **Open issues**: 14; **Forks**: 3 (GitHub API 2026-05-19)
- **Relevance to v0.2** (≤ 3 sentences): This package is the primary RH-numerical
  competitor that PRD lines 408–410 and the scope-doc Pillar E narrative reference;
  it implements the Olver 2011/2012 collocation method for Riemann–Hilbert problems
  and ships a worked example computing the Hastings–McLeod solution to Painlevé II
  (confirmed in README), which is the canonical PII-only example in the Trogdon–Olver
  book. The package's PII scope, dormancy since v0.1.0 (December 2019), and the
  absence of any multi-component or higher-hierarchy example confirm PRD's claim that
  "no released software exists for higher-hierarchy or multi-component systems."
  PadeTaylor v0.2 should explicitly cite this package in ADR documentation as
  the state of the art for RH-numerical PII, positioning PadeTaylor's Padé–Taylor
  approach as the complementary method that extends to higher hierarchy and vector systems.
- **Version range PadeTaylor v0.2 targets**: TBD. PadeTaylor v0.2 does not depend on
  RiemannHilbert.jl directly; this entry is for cross-validation reference only.
  If a benchmark comparison is implemented, target v0.1.0 (the only released version).
- **Coverage limitations** relevant to PadeTaylor v0.2:
  - **PII only in the nonlinear-special-functions examples**: the README's worked
    example is the Hastings–McLeod solution to PII; no examples for PI, PIV, PV,
    PVI, Noumi–Yamada hierarchies, or Garnier systems are present.
  - **Dormant since December 2019**: v0.1.0 is the only release; the package
    has not been updated to track Julia version changes (e.g., the Travis CI badge
    is dead; modern Julia CI uses GitHub Actions). May require dependency updates
    to run on Julia 1.10+.
  - **Wolfram Language dependency**: the 58.8% Wolfram Language composition suggests
    the package may call Mathematica internally for symbolic steps, which would make
    it non-portable without a Mathematica licence.
  - **13 open pull requests** (WebFetch of the repo page, 2026-05-19): indicates
    pending changes that have not been merged, reinforcing the dormant-development
    classification.

> Stub placeholder for code-side citations.
> Cite as: Olver, S., *RiemannHilbert.jl*, Julia package,
> https://github.com/JuliaHolomorphic/RiemannHilbert.jl, v0.1.0 (2019).
