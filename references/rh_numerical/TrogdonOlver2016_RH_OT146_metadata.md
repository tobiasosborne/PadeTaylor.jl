# Trogdon, Olver 2016: Riemann–Hilbert Problems, Their Numerical Solution and the Computation of Nonlinear Special Functions

- **Canonical citation**: Thomas Trogdon and Sheehan Olver,
  *Riemann–Hilbert Problems, Their Numerical Solution, and the Computation of
  Nonlinear Special Functions*, SIAM, Other Titles in Applied Mathematics vol. 146,
  Philadelphia, 2016 (copyright page: 2015). ISBN 978-1-61197-419-5 (print),
  eISBN 978-1-61197-420-1.
- **Publisher URL**: https://epubs.siam.org/doi/book/10.1137/1.9781611974201
- **DOI**: 10.1137/1.9781611974201
- **Pages**: xviii + 356 pp.
- **Open sample chapter**: The front matter, back matter, and Chapter 1 (plus the
  table of contents) are openly available at
  https://www.math.uci.edu/~ttrogdon/publications/TrogdonOlver-Sample.pdf
  Confirmed accessible 2026-05-19; file fetched (701.7 KB PDF).
- **Sample chapter local path**: `references/rh_numerical/TrogdonOlver2016_RH_OT146_sample.pdf`
  (sample PDF has been fetched and saved to this path by the Batch 4+5 agent,
  2026-05-19; confirmed pdftotext extraction yields complete TOC and Chapters 1–2 text)
- **Acquisition status**: book_paywalled for body chapters. Routes:
  - Sample (front matter + Ch. 1 + back matter): open at URL above
  - Full book: SIAM member access at epubs.siam.org; institutional library subscription;
    or SIAM bookstore purchase https://bookstore.siam.org/ot146/

- **Table of contents** (fetched from open sample PDF via pdftotext 2026-05-19;
  source: `TrogdonOlver-Sample.pdf` p. i–iii):

  **Part I — Riemann–Hilbert Problems**
  - Chapter 1: Classical applications of Riemann–Hilbert problems
    - §1.1 Error function
    - §1.2 Elliptic integrals
    - §1.3 Airy function
    - §1.4 Monodromy
    - §1.5 Jacobi operators and orthogonal polynomials
    - §1.6 Spectral analysis of Schrödinger operators
  - Chapter 2: Riemann–Hilbert Problems (precise statement, Hölder / Hardy theory,
    singular integral equations)
  - Chapter 3: Inverse Scattering and Nonlinear Steepest Descent

  **Part II — Numerical Solution of Riemann–Hilbert Problems**
  - Chapter 4: Approximating Functions (DFT, Chebyshev, mapped series, vanishing bases)
  - Chapter 5: Numerical Cauchy Transforms (convergence, unit circle, unit interval,
    square root singularities, endpoint approximation)
  - Chapter 6: Numerical Solution of RH Problems (projection and collocation methods;
    case studies: Airy equation, monodromy of ODE with three singular points)
  - Chapter 7: Uniform Approximation Theory (numerical RH framework, collocation
    realisation)

  **Part III — The Computation of Nonlinear Special Functions**
  - Chapter 8: The KdV and mKdV Equations (uniform approximation, small-dispersion)
  - Chapter 9: The Focusing and Defocusing NLS Equations (inverse scattering,
    Robin boundary conditions, singular solutions)
  - Chapter 10: The Painlevé II Transcendents (positive/negative x regimes,
    Stokes-data parametrisation, numerical results)
  - Chapter 11: The Finite-Genus Solutions of the Korteweg–de Vries Equation
    (Riemann surfaces, regularisation, numerical computation)
  - Chapter 12: The Dressing Method and Nonlinear Superposition

  **Part IV — Appendices**
  - Appendix A: Function Spaces and Functional Analysis
  - Appendix B: Fourier and Chebyshev Series
  - Appendix C: Complex Analysis (inferred analyticity)
  - Appendix D: Rational Approximation (bounded contours, Lipschitz graphs)
  - Appendix E: Additional KdV Results (comparison with existing methods, g-function)

  Bibliography (pp. 379–384), Notation and Abbreviations (pp. 385–389), Index (pp. 391–).

- **Chapters most relevant to v0.2**:
  - Chapter 10 (Painlevé II Transcendents) — the canonical RH-numerical treatment
    of PII that defines the "competing methodology" v0.2 docs must honestly compare
    against; confirms the Hastings–McLeod computation RiemannHilbert.jl implements.
  - Chapters 6–7 (Numerical Solution of RH Problems / Uniform Approximation Theory)
    — the algorithmic foundation for Olver's numerical RH framework; explains why
    this approach generalises poorly to higher-hierarchy / multi-component without a
    new RH problem formulation per system.
  - Chapter 3 (Inverse Scattering and Nonlinear Steepest Descent) — background for
    the Painlevé-VI / Garnier isomonodromy side of the v0.2 story.
  - Appendix D (Rational Approximation) — confirms the book's awareness of rational
    approximation as an adjacent methodology; reinforces PadeTaylor's positioning.

- **Why this book matters for v0.2** (≤ 3 sentences): The canonical book reference
  for the Olver/Trogdon numerical RH approach that PRD lines 408–412 names as the
  primary competing methodology to PadeTaylor's Taylor–Padé approach; Chapter 10
  treats Painlevé II specifically and is the intellectual basis for RiemannHilbert.jl.
  The openly available sample chapter and TOC confirm that the book covers PII only
  in the nonlinear-special-functions part, validating PRD's claim that "no released
  software exists for higher hierarchy or multi-component systems." Appendix D's
  rational approximation section is directly relevant as cross-methodology context.

- **Blocking status**: non-blocking — sample chapter + the open paper arXiv:1210.2199
  (Olver 2012 FoCM "Numerical solution of Riemann–Hilbert problems: Painlevé II")
  cover what v0.2 ADRs need for the "competing methodology" narrative. Full book
  useful for a thorough v0.3 survey of the RH-numerical landscape.

> Stub placeholder — future commits cite this book via
> `references/rh_numerical/TrogdonOlver2016_RH_OT146_metadata.md`.
> The open sample PDF (TOC + Ch. 1 + front/back matter) is at
> `references/rh_numerical/TrogdonOlver2016_RH_OT146_sample.pdf`.
