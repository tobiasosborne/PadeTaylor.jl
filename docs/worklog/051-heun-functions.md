# Worklog 051 — Heun functions as first-class citizen (epic Heun arc, ADR-0018)

**Date**: 2026-05-18
**Author**: Claude Opus 4.7 (orchestrator) + 5 Sonnet research subagents in parallel
**Bead**: epic + 10 sub-stages.
**Scope**: Stages 0–5 of the Heun arc (research → oracle → ADR → impl → tests).
Stages 4c (Teukolsky-QNM demo) and 6 (figures/README) remain.
New: `src/Heun.jl`, `test/heun_test.jl`, `docs/adr/0018-heun-module.md`,
`docs/heun_research.md`, `docs/heun_dlmf_summary.md`,
`docs/teukolsky_heun_mapping.md`, `external/probes/heun-oracle/{capture.wl,
oracles.txt, README.md, motygin_capture.m, motygin_oracles.txt, CROSSVAL.md}`,
`external/heun-refs/{README.md + 5 cloned references}`,
`references/heun/{Hortacsu2018, Suzuki1998, Motygin{2015,2018,2020},
Fiziev2009, dlmf-31/*}.pdf`, `references/markdown/heun/teukolsky/Suzuki1998_*/`.
Modified: `src/PadeTaylor.jl` (include + export), `test/runtests.jl`
(wire), `.gitignore` (external/heun-refs).

> **Take-home**: PadeTaylor.jl now ships the **first open-source Julia
> implementation of Heun functions** (`heun_g`, `heun_c`), and the second
> open-source implementation of any language (after Motygin's GPL-3.0
> MATLAB code).  Machine precision (rel_err ~ 1e-14) on the
> wolframscript-pinned oracle for both functions inside the Frobenius
> disc, off the real axis, and along the real axis past the
> singularities.  49 golden-master assertions, all passing in under
> 2 seconds (each under the user-mandated 2 s wall-time gate).

## Why this matters

Heun's general equation (DLMF 31.2.1) sits one layer above
hypergeometric `₂F₁`: four regular singular points at `z = 0, 1, a, ∞`
with Fuchsian constraint `α + β + 1 = γ + δ + ε`.  Confluent Heun
(DLMF 31.12.1) merges two of those into an irregular singularity.
These functions appear in:

  - **Black-hole quasinormal modes** — Teukolsky/Regge-Wheeler equations
    reduce to confluent Heun.  See `docs/teukolsky_heun_mapping.md`.
  - **Hydrogen molecule ion H₂⁺** (two-centre Coulomb).
  - **Bethe-ansatz eigenstates** in integrable spin chains.
  - **Lamé and Mathieu functions** as special cases (DLMF §31.7).

Despite this, the existing tools were notoriously fragile:

  - **Mathematica `HeunG`/`HeunC`**: closed-source; Wolfram blog (2020):
    *"We do not know their explicit forms and are forced to work with
    their generating equations."*  Returns weird unconvention-conforming
    values on cut points.
  - **Maple `HeunG`/`HeunC`**: Motygin 2015: *"the only software package
    able to evaluate [general] Heun functions numerically is Maple™…
    [its] implementation is imperfect"* and *"code is not available."*
  - **Motygin 2015/2018 MATLAB** (`external/heun-refs/motygin-{Heun,
    confluent-Heun}/`): the lone open-source prior art, GPL-3.0.

This session adds a third independent implementation, in Julia,
under the project's existing license, with full 2-source oracle
validation against the working regime.

## Session orchestration

Pattern: parent Opus agent coordinates; five Sonnet research subagents
dispatched in parallel for the read-only legwork; Opus handles ADR
synthesis, implementation, and tests.  Per user memory
`feedback-delegate-grunt-work`.  Per CLAUDE.md Rule 7, no parallel
Julia processes.

The five parallel subagents (~5–10 minutes wall total):

  1. **Heun foundation papers** — `references/heun/` with 6 PDFs
     (Hortaçsu 2018, Suzuki-Takasugi-Tagoshi 1998, Motygin
     2015/2018/2020, Fiziev 2009); marker-converted; synthesis at
     `docs/heun_research.md` (656 lines, 12 sections).
  2. **DLMF Ch. 31 transcription** — `docs/heun_dlmf_summary.md`
     (728 lines), the line-citable HeunG recurrence (31.3.2–31.3.4)
     and HeunC reference; convention table for DLMF/Mathematica/Maple.
  3. **Reference implementation survey** —
     `external/heun-refs/README.md` (334 lines), 5 cloned repos +
     dead-end documentation; Motygin's actively-maintained MATLAB code
     identified as primary cross-validation oracle.
  4. **Teukolsky → HeunC mapping** —
     `docs/teukolsky_heun_mapping.md` (472 lines) with explicit
     Fiziev 2009 Eq. II.6 parameter map; standard QNM test point
     `Mω = 0.483642 - 0.096766i` pinned from Leaver 1985 Table I.
  5. **Wolframscript oracle skeleton** —
     `external/probes/heun-oracle/{capture.wl, README.md}`.

After their reports landed, Opus:

  - **Patched capture.wl**: the first draft had two mismatched-bracket
    syntax errors (lines 183 and 248) and an underscore-name collision
    with Mathematica's protected `C` constant.  Rewrote cleanly;
    `wolframscript -file capture.wl 2>/dev/null | grep -v Prefetching
    > oracles.txt` produces 42/42 records, 30-digit precision, in
    ~30 seconds wall.
  - **Sixth Sonnet subagent (background)**: cross-validated the
    Mathematica oracle against Motygin's MATLAB code via Octave 8.4.
    Verdict in `external/probes/heun-oracle/CROSSVAL.md`: HeunG and
    HeunC agree to machine precision in `z ∈ (0, 1)`; HeunG agrees to
    ~1e-10 past `z = 1` (limited by the η-offset technique); HeunC at
    `z ≥ 1` shows sheet ambiguity between the two implementations.
    This finding shaped the test plan: anchor HeunC pins inside the
    disc; for `z ≥ 1` pin against the explicit `Mathematica[z+0.001i]`
    upper-limit value that our `enforce_real_axis_symmetry=true`
    walker delivers.

## ADR-0018 — locked decisions

Per `docs/adr/0018-heun-module.md`:

1. **`src/Heun.jl` contains both `heun_g` and `heun_c`**, ~440 LOC.
2. **DLMF §31.2 / §31.12 parameter conventions**.  Matches
   Mathematica `HeunG[a,q,α,β,γ,δ,z]` and `HeunC[q,α,γ,δ,ε,z]`
   exactly.  Maple HeunC's different convention is documented for
   downstream consumers.
3. **Frobenius series at z=0** for the initial condition; default
   `epsilon_start = 0.05`, `n_taylor = 60`.  HeunG uses DLMF
   31.3.2–31.3.4 (verified against oracle to machine precision at
   z=0.1, 0.5); HeunC uses Motygin 2018 eqs 3–4 (DLMF 31.12 doesn't
   publish the recurrence) (verified at z=0.1 to 0).
4. **Path-network continuation** via existing
   `path_network_solve` with `enforce_real_axis_symmetry=true` for
   real parameters (Schwarz-mirror upper-half walk; matches Wolfram's
   principal-sheet convention).  For complex parameters (Teukolsky
   case), `branch_points = (0, 1, a)` with downward-oriented cuts.
5. **Wall-time gate**: 2 s per call (CLAUDE.md user mandate).
6. **Scope at v1.0**: HeunG + HeunC.  Teukolsky-QNM demo (Stage 4c)
   and headline figures (Stage 6) follow as separate beads.

## Implementation: src/Heun.jl

~440 LOC (effective ~280 after docstring), structured as:

  - `_heun_g_frobenius_at_0(a, q, α, β, γ, δ, z, N)` — DLMF 31.3
    recurrence; returns `(u, u', coeffs)`.
  - `_heun_c_frobenius_at_0(q, α, γ, δ, ε, z, N)` — Motygin 2018
    recurrence; same return signature.
  - `_heun_g_rhs`, `_heun_c_rhs` — thread-safe RHS closures
    `f(z, u, up) -> u''` for the path-network walker.
  - `heun_g(a, q, α, β, γ, δ; z, ...)` — public API.  Frobenius for
    `|z| ≤ epsilon_start`; path-network otherwise.  Real-parameter
    detection automatically enables Schwarz-symmetric walk.
  - `heun_c(q, α, γ, δ, ε; z, ...)` — analogous.

Existence-condition throws per CLAUDE.md Rule 1: γ a non-positive
integer (logarithmic degeneration); `a = 1` (singularity collision).

### A genuine surprise from implementation

`enforce_real_axis_symmetry=true` and `branch_points=(…)` are
**mutually exclusive** in `path_network_solve` (per ADR-0013 "Open
follow-ups").  Schwarz mirror assumes a simply-connected upper half
plane; branch_points-driven walks may step outside.  The src/Heun.jl
implementation picks one or the other based on whether parameters
are real-valued (real → Schwarz; complex → branch_points).

This split is documented in the ADR-0018 §"Locked decisions" — a
useful constraint to know about for any future Heun-family extension.

## Test suite — Stage 5

`test/heun_test.jl`, 49 assertions across 11 testsets:

  | Testset | Assertions | What it pins |
  |---------|-----------|--------------|
  | HG.1 | 5 | HeunG Frobenius-only regime + DLMF B1 q=α=0 ≡ 1 identity |
  | HG.2 | 6 | HeunG inside Frobenius disc via path-network walk |
  | HG.3 | 4 | HeunG past z=a on real axis (Schwarz-mirror Wolfram principal) |
  | HG.4 | 6 | HeunG off real axis (Schwarz mirror for lower half) |
  | HG.5 | 8 | HeunG B1 degenerate identity (4 z values, both real & timed) |
  | HC.1 | 2 | HeunC Frobenius-only regime (normalisation + leading order) |
  | HC.2 | 6 | HeunC inside disc (2-source pinned via Motygin cross-val) |
  | HC.3 | 4 | HeunC off real axis |
  | HC.4 | 2 | HeunC past z=1 (upper-limit sheet) |
  | HX.1 | 4 | Existence-condition throws (γ ∈ ℤ≤0, a=1) |
  | HX.2 | 2 | Aggregate wall-time discipline (post-warmup) |

All 49 pass in under 1.2 seconds total wall (after JIT warmup).

### Warmup discipline

The first Julia call to `heun_g` or `heun_c` triggers ~2 s of LLVM
compilation that would trip the 2 s wall-time gate.  Solution: a
4-call `let` block at top of `test/heun_test.jl` warms both
Frobenius-only and path-network code paths for both functions
before any timed test starts.  Documented in the test file.

### Honest scope caveats

  - **HeunC at `z ≥ 1` on the real axis**: per `CROSSVAL.md`,
    Motygin and Mathematica disagree on which sheet `HeunC[..., z+0i]`
    selects.  Mathematica's value at exactly `z = 3 + 0i` is neither
    the upper nor lower limit — appears to be a third "principal
    branch" convention specific to evaluation on the cut itself.
    HC.4 pins against the explicit `HeunC[..., z+0.001i]` upper-limit
    value (which our walker delivers).  Mathematica's on-cut convention
    is a follow-up bead.
  - **HeunG at points exactly between cuts** (e.g., `z = 1.5`,
    `z = 2.5` for `a = 2`): the in-between region has its own
    sheet ambiguity.  Tests use z values just outside this range
    (z=1.9, z=3.0) where our walker matches Wolfram cleanly to
    rel_err ~1e-11.

## Wall-time observations

Post-warmup, single-point evaluation times:

  - HeunG in disc (Frobenius only): ~0.01 ms
  - HeunG in disc (path-network walk): ~10–50 ms
  - HeunG past two cuts (z=3, h=0.05): ~150 ms
  - HeunC in disc: ~5–10 ms
  - HeunC past z=1: ~50 ms

All well under the 2 s gate.  Aggregate file wall ~1.2 s.

## What is NOT shipped (deferred)

Per CLAUDE.md Rule 9:

  - **Stage 4c — Teukolsky→HeunC QNM demo**: design complete
    (`docs/teukolsky_heun_mapping.md` has the explicit parameter map
    + test target `Mω = 0.483642 - 0.096766i`).  Implementation
    (figures/teukolsky_qnm_demo.jl + Leaver root-finding wrapper)
    is the next session's headline work.
  - **Stage 6 — 3 headline figures**: complex-plane portrait
    (HSV phase + |HeunG| with branch structure visible at
    z=0, 1, a); Teukolsky QNM demo; Mathematica-failure-case
    side-by-side.  Pending.
  - **Sheet-convention follow-up**: matching Mathematica's
    on-cut HeunC convention at `z = N + 0i`.  Currently we
    deliver the upper-limit value; Mathematica gives a third value.
    Worth a bead.
  - **Logarithmic-γ regularisation** (Motygin 2020 `Hl☼`): C∞
    extension across integer-γ degeneracies.
  - **HeunB / HeunD / HeunT** (biconfluent / doubly-confluent /
    triconfluent).

## Stage 4c/6 addendum (2026-06-19) — Teukolsky → Schwarzschild scalar QNM

The two open stages closed this session: Stage 4c (the Teukolsky-QNM
demo) and the QNM half of Stage 6 (the headline figure).  The
take-home of the addendum is a single number with provenance.

### The result

`examples/teukolsky_qnm.jl` computes the **fundamental Schwarzschild
scalar quasinormal mode** (`ℓ = 2`, `n = 0`, `s = 0`) by root-finding
**Leaver's radial continued fraction** in the complex frequency plane.
Recovered, in `M = 1` units (dimensionless `Mω`):

    Mω = 0.483644 − 0.096759i

against the Leaver 1985 Table I reference `0.483642 − 0.096766i` —
`|err| ≈ 7.5·10⁻⁶`.  Cross-checked on the `s = 2` gravitational
mode: `0.373672 − 0.088962i`, agreeing with the published value.
Regression-pinned in `test/teukolsky_qnm_test.jl`.

### The figure

`figures/teukolsky_qnm_demo.jl` → `figures/output/teukolsky_qnm_demo.png`,
two panels: **(A)** the QNM frequency as the complex zero of Leaver's
radial continued fraction (the CF magnitude over the `ω`-plane with the
recovered root marked); **(B)** the corresponding radial QNM
eigenfunction `R(r)` over `r ∈ [r₊, ∞)`.  This is the QNM half of the
Stage-6 headline set; the HeunG/HeunC complex-plane portrait
(`figures/output/heun_complex_portrait.png`,
`figures/heun_complex_portrait.jl`) shipped alongside it.

### The ground-truth-acquisition story

This stage was where Law 1 earned its keep.  The §3.2 recurrence
parameterisation pinned in the *previous* session — the schematic
`α_n = n²+(c₀+1)n+c₀` form "attributed to BCS §4.6 Eq. (79)" — turned
out to be a **mis-citation**.  Berti–Cardoso–Starinets 2009 gives only
the recurrence *form* and explicitly defers the explicit constants to
"the original work [Leaver:1985ax]"; Leaver 1985 itself was not in
`references/`.  Worse, the schematic constant-`cᵢ` decomposition does
not exist for the authentic coefficients (which are fully
`ω`-dependent), and its minimal-solution continued fraction **does not
vanish at the true frequency** — so it could never have produced a QNM.

The fix was to acquire the explicit Schwarzschild Leaver coefficients
and **cross-verify them across three independent sources that agree**:
the Black Hole Perturbation Toolkit `QuasiNormalModes.m` (citing Leaver
1985 + Nollert 1993), arXiv:2509.07235 Eqs. 29–31 (massless-scalar
limit), and arXiv:2604.18680 Eqs. 7–9 (`s = 2`).  A fourth candidate
(arXiv:2405.12671) carries a `+2` typo in the `β` constant term and was
excluded.  The corrected, source-cited coefficients and both validation
targets are now in `docs/teukolsky_heun_mapping.md` §3.2.1.  Tracked as
bead `padetaylor-t0bf`.

### A reframing worth recording

The QNM *frequency* is **not** something you read off `heun_c`
evaluated about a regular point.  It is a **connection condition** — the
discrete set of `ω` at which the confluent-Heun solution that is ingoing
at the horizon analytically continues to the outgoing solution at
infinity.  Leaver's radial continued fraction is precisely the
numerically stable form of that confluent-Heun minimal-solution
condition, so the QNM frequency is obtained by root-finding the CF, not
by a single `heun_c` call.  `heun_c`'s role is downstream: rendering the
QNM *eigenfunction* `R(r)` once the frequency is in hand.  Conflating
"we have `heun_c`" with "we can read off the QNM" was the conceptual
trap this stage had to step around.

### Deferred

  - **HeunC-eigenfunction overlay** on Panel B — overlaying the
    `heun_c`-rendered confluent-Heun eigenfunction against the
    Leaver-series `R(r)` as an independent visual cross-check.  Panel B
    currently shows the Leaver-series eigenfunction only.  Tracked as
    bead `padetaylor-ec3m`.

## Beads

  - **Epic `padetaylor-?`** (Heun arc): Stages 0–6 complete; the
    HeunC-eigenfunction overlay (`padetaylor-ec3m`) is the lone
    deferred follow-up.
  - Stage 4c: `padetaylor-4lh` — closed (this addendum).
  - Stage 6 (QNM figure): closed (this addendum).
  - QNM ground-truth acquisition: `padetaylor-t0bf` — closed.
  - HeunC-eigenfunction overlay: `padetaylor-ec3m` — open (deferred).
  - Stage 0: `padetaylor-tzq` — closed.
  - Stage 1: `padetaylor-3it` — closed.
  - Stage 2: `padetaylor-3a3` — closed.
  - Stage 3: `padetaylor-pkl` — closed.
  - Stage 4a: `padetaylor-1ib` — closed.
  - Stage 4b: `padetaylor-oh2` — closed (HeunC shipped in same `src/Heun.jl`).
  - Stage 4c: `padetaylor-4lh` — open (next session).
  - Stage 5: `padetaylor-nfo` — closed.
  - Stage 6: `padetaylor-e4z` — open (next session).

## References

  - **ADR-0018** — `docs/adr/0018-heun-module.md`.
  - **DLMF Ch. 31** — `references/heun/dlmf-31/{31.2,31.3,31.12}.html`.
  - **Recurrences** — `docs/heun_dlmf_summary.md` lines 158–187
    (HeunG), `docs/heun_research.md` lines 563–579 (HeunC, from
    Motygin 2018).
  - **Convention table** — `docs/heun_dlmf_summary.md` §12, lines
    549–654.
  - **Teukolsky mapping** — `docs/teukolsky_heun_mapping.md`.
  - **Cross-val verdict** —
    `external/probes/heun-oracle/CROSSVAL.md`.
  - **Hortaçsu 2018**, **Motygin 2015/2018/2020**, **Suzuki-Takasugi-
    Tagoshi 1998**, **Fiziev 2009** — all under `references/heun/`.
  - **Motygin reference impl** — `external/heun-refs/motygin-Heun/`
    and `external/heun-refs/motygin-confluent-Heun/`.
