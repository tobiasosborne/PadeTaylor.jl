# Bug sweep D4 — Heun module (DLMF 31.2 general + 31.12 confluent)

Date: 2026-06-01. Auditor area: `src/Heun.jl` equation coefficients +
local Frobenius series. Read-only audit (no Julia run, no edits).

## Area

`src/Heun.jl` — the general Heun equation (DLMF 31.2.1) Frobenius
recurrence and RHS, the confluent Heun equation (DLMF 31.12.1) Frobenius
recurrence and RHS, the Fuchsian-constraint computation, the accessory
parameter `q` placement, the initial conditions `c_1`/`b_1`, and the
real-vs-complex parameter dispatch into `path_network_solve`.

Headline result: the Heun **equation coefficients and local series are
faithfully transcribed and verified correct** — both against the
authoritative references (DLMF 31.3 HTML; Motygin 2018 TeX source) and
numerically against the Mathematica oracle to full double precision. No
sign error, off-by-one, factorial confusion, or `q`-placement error was
found in the coefficient/series core. The only Heun-specific
intermittency mechanism is a *dispatch* discontinuity (the algorithm
silently switches strategy at the real/complex parameter boundary and
between the two cut conventions); that is a design/branch-cut concern
documented as a deliberate v1 choice, not a mis-transcribed equation,
and the sheet-selection consequences belong to the path-network audit.

## References checked

- DLMF 31.3 recurrence, raw HTML alttext extracted from
  `references/heun/dlmf-31/31.3.html`:
  - initial condition: `a\gamma c_{1}-qc_{0}=0` (so `c_1 = q/(aγ)`)
  - recurrence: `R_{j}c_{j+1}-(Q_{j}+q)c_{j}+P_{j}c_{j-1}=0`
  - `P_{j}=(j-1+\alpha)(j-1+\beta)`
  - `Q_{j}=j((j-1+\gamma)(1+a)+a\delta+\epsilon)`
  - `R_{j}=a(j+1)(j+\gamma)`
- Motygin 2018 confluent recurrence, TeX source
  `references/tex/heun/Motygin2018_confluent_Heun_numerical_1804.01007/motygin.tex:140-207`:
  - ODE (`:140-141`): `HeunC'' + (γ/z + δ/(z-1) + ε)HeunC' + (αz-q)/(z(z-1))HeunC = 0`
  - independence of ε, α (`:148`)
  - recurrence (`:195`): `P_n b_n = Q_n b_{n-1} + R_n b_{n-2}`
  - coefficients (`:201-203`): `P_n=n(γ-1+n)`, `Q_n=-q+(n-1)(γ+δ-ε+n-2)`, `R_n=(n-2)ε+α`
  - ICs (`:207`): `b_{-1}=0`, `b_0=1`, `HeunCl'(...;0)=-q/γ` (so `b_1=-q/γ`)
- Motygin 2015 general-Heun recurrence (Motygin's own index convention),
  TeX `references/tex/heun/Motygin2015_Heun_evaluation_1506.03848/motygin.tex:250-266`:
  `P_n=a n(γ-1+n)`, `Q_n=q+(n-1)[(a+1)(γ+n-2)+ε+aδ]`, `R_n=-(n-2+α)(n-2+β)`,
  `b_{-1}=0`, `b_0=1`; radius `R_0=min{1,|a|}` (`:269-270`). (Algebraically
  the DLMF form after reindexing; used as an independent cross-check.)
- `docs/heun_dlmf_summary.md:158-187` (HeunG recurrence), `:382-475`
  (HeunC ODE + recurrence cautions), `:548-654` (convention tables).
- `docs/heun_research.md:106-159` (HeunG recurrence), `:280-321` (HeunC
  recurrence Motygin 2018), `:543-585` (recurrence summary both forms).
- `docs/adr/0018-heun-module.md:111-155` (recurrence transcription),
  `:178-211` (RHS factories + branch handling), `:38-42` (declared
  `branch_cut_angles = π`).
- `docs/worklog/051-heun-functions.md:113-122` (Frobenius verified to
  machine precision; downward-oriented cuts for complex params decided
  deliberately).
- `external/probes/heun-oracle/oracles.txt` (Mathematica 14.3 oracle,
  50-digit working precision) — used for end-to-end numeric replication.
- `src/Problems.jl:98-132` (`PadeTaylorProblem` IC convention `y0=(u,u')`).
- `src/PathNetwork.jl:339-348, 375-404` (`enforce_real_axis_symmetry`
  semantics; mutual exclusivity with `branch_points`).

## Findings

### [LOW] Strategy/branch dispatch can jump value at the real/complex parameter boundary and across the chosen cuts

Location: `src/Heun.jl:318-347` (HeunG) and `src/Heun.jl:385-405` (HeunC).

Ground truth (cited): the principal-sheet convention this module targets
is Wolfram's, with cuts `[1,+∞)` and `[a, e^{i arg a}·∞)`
(`docs/heun_research.md:226-241`, `:533-539`). ADR-0018 declares
`branch_cut_angles = π per branch (Wolfram convention)`
(`docs/adr/0018-heun-module.md:38-42`). The actual code instead uses
`branch_cut_angles = (-π/2, -π/2, -π/2)` (`src/Heun.jl:345, 403`), a
deliberate "downward cuts" choice (`docs/worklog/051-heun-functions.md:122`).

Code behavior: when all six (HeunG) / five (HeunC) parameters are detected
real to within `< eps()` (`src/Heun.jl:318-320, 385-387`), the walk uses
`enforce_real_axis_symmetry=true` with NO `branch_points` (Schwarz mirror
in `Im z ≥ 0`); otherwise it uses `branch_points=(0,1,a)` (or `(0,1)`)
with cuts at angle `-π/2`. These are two different continuation
strategies with two different cut placements, so the same target `z`
(off the positive real axis, or on the far side of a singular point) can
be returned on different sheets depending on whether the parameters are
exactly real.

Mechanism for intermittent discontinuity: a parameter sweep that crosses
from real to slightly-complex parameters (e.g. a Teukolsky frequency
sweep where `ω` acquires a small imaginary part) flips the dispatch from
the Schwarz branch to the `-π/2`-cut branch, which can land on a
different Riemann sheet and produce a jump in the returned value at the
crossing. Likewise, because the chosen cuts (`-π/2`) differ from the
declared/Wolfram cuts, continuation paths that cross the real axis pick
a different sheet than a downstream consumer expecting Wolfram's
`[1,∞)`/`[a,∞)` cuts. This is data-/parameter-dependent, hence
intermittent.

Intermittent? Yes (parameter-locus dependent).

Confidence: 0.4. This is a real mechanism, but it is a *dispatch / sheet
selection* concern, not a mis-transcribed Heun coefficient. The cut
choice is documented as deliberate; the actual sheet-correctness of
`path_network_solve` is outside this area (path-network audit owns it).
The coefficient transcription itself — my mandate — is correct (below).

### [LOW] Docstring/ADR `h` default disagrees with code default (cosmetic; no numeric effect)

Location: docstrings `src/Heun.jl:249, 274, 352`, ADR
`docs/adr/0018-heun-module.md:25` all say default `h=0.1`; the actual
keyword defaults are `h::Real = 0.05` (`src/Heun.jl:290, 367`).

Ground truth (cited): Law 2 (docs in lockstep with code),
`CLAUDE.md`. The numerical value used is well-defined (`0.05`); only the
prose is stale.

Code behavior: callers who don't pass `h` get `0.05`, not the documented
`0.1`. No sign/coefficient impact.

Mechanism: none that produces a discontinuity — a smaller fixed step is
strictly safer. Listed only because Law 2 requires flagging doc drift.

Intermittent? No.

Confidence: 0.95 (the mismatch is a literal fact in the source).

## Areas verified correct

- **HeunG Frobenius recurrence coefficients** (`src/Heun.jl:133-136`):
  `P_j=(j-1+α)(j-1+β)`, `Q_j=j((j-1+γ)(1+a)+aδ+ε)`, `R_j=a(j+1)(j+γ)`,
  and `c_{j+1}=((Q_j+q)c_j − P_j c_{j-1})/R_j`. Byte-for-byte match to
  DLMF 31.3.3/31.3.4 (extracted alttext from
  `references/heun/dlmf-31/31.3.html`) and to
  `docs/heun_dlmf_summary.md:170-173`. The sign in front of `P_j` is the
  correct minus (the DLMF relation `R c_{j+1} − (Q+q)c_j + P c_{j-1}=0`
  rearranges to a minus on the `P` term). Loop `for j in 1:(N-1)` fills
  `c[3..N+1]` = `c_2..c_N`; with `c[1]=c_0`, `c[2]=c_1` all `N+1`
  coefficients are produced — no off-by-one.

- **HeunG initial conditions** (`src/Heun.jl:130-131`): `c_0=1`,
  `c_1=q/(aγ)` (sign `+`). Matches DLMF `aγ c_1 − q c_0 = 0`
  (`references/heun/dlmf-31/31.3.html` alttext) and independently
  re-derived from DLMF 31.2.1 (coefficient of `z^0` after clearing
  denominators gives `γa c_1 − q = 0`).

- **HeunG Fuchsian constraint** (`src/Heun.jl:127` and `:216`):
  `ε = α + β + 1 − γ − δ`, the correct solution of `α+β+1 = γ+δ+ε`
  (`docs/heun_dlmf_summary.md:34`, `docs/heun_research.md:26`). Same `ε`
  used in both the recurrence `Q_j` and the RHS drift term.

- **HeunC Frobenius recurrence coefficients** (`src/Heun.jl:182-185`):
  `P_n=n(n+γ-1)`, `Q_n=-q+(n-1)(γ+δ-ε+n-2)`, `R_n=α+(n-2)ε`, and
  `b_n=(Q_n b_{n-1}+R_n b_{n-2})/P_n`. Byte-for-byte match to Motygin
  2018 TeX `motygin.tex:201-203` and `docs/heun_research.md:307-310`.
  `ε` here is correctly the *independent* confluent parameter (Motygin
  TeX `:148`), NOT the Fuchsian-derived `ε` — `_heun_c_frobenius_at_0`
  takes `ε::T` as a direct argument, so there is no contamination from
  the HeunG constraint. Loop `for n in 2:N` fills `b[3..N+1]` = `b_2..b_N`;
  `b[1]=b_0`, `b[2]=b_1` — all `N+1` coefficients produced, no off-by-one.

- **HeunC initial condition** (`src/Heun.jl:180`): `b_1 = -q/γ` (sign
  `-`, deliberately opposite to HeunG's `+q/(aγ)`). Matches Motygin TeX
  `:207` (`HeunCl'(...;0) = -q/γ`) and independently re-derived from DLMF
  31.12.1 (coefficient of `z^0` after clearing `z(z-1)` gives
  `-γ b_1 - q = 0`). The sign genuinely differs between the two functions;
  it is NOT a copy-paste error — both match their own reference.

- **HeunG RHS factory** (`src/Heun.jl:215-222`):
  `u'' = -((γ/z + δ/(z-1) + ε/(z-a))u' + (αβz-q)/(z(z-1)(z-a))u)`,
  exactly DLMF 31.2.1 rearranged (`docs/heun_dlmf_summary.md:26-29`).
  The accessory parameter `q` sits in the numerator `αβz − q` with the
  correct sign and is not in any drift term.

- **HeunC RHS factory** (`src/Heun.jl:233-239`):
  `u'' = -((γ/z + δ/(z-1) + ε)u' + (αz-q)/(z(z-1))u)`, exactly DLMF
  31.12.1 (`docs/heun_dlmf_summary.md:393-395`, Motygin TeX `:140-141`).
  The `ε` term is correctly a *constant* (the confluent-equation
  signature), not `ε/(z-a)`; `q` is in the numerator `αz − q`.

- **Horner evaluation of series + derivative** (`src/Heun.jl:139-145`
  HeunG, `:188-193` HeunC): `u(z)=Σ_{k=0}^N c_k z^k` and
  `u'(z)=Σ_{k=1}^N k c_k z^{k-1}` both evaluated by correct Horner
  schemes; the derivative guard `k >= 1` correctly omits the `k=0`
  (constant) term. No factorial/coefficient confusion — these are Taylor
  *coefficients* `c_k` (not derivatives), and `u'` multiplies by `k`
  exactly as the power rule requires.

- **IC handoff to the walker** (`src/Heun.jl:323-326, 390-393`): passes
  `(u0, up0)` = `(u(z0), u'(z0))` as `y0`; `PadeTaylorProblem`'s
  2nd-order branch expects `y0=(u,u')` (`src/Problems.jl:100-101`), and
  the RHS closure takes `f(z,u,up)` with `up=u'`. No sign flip or
  rescaling of the derivative at the interface.

- **Existence-condition guards** (`src/Heun.jl:113-125, 170-175`):
  HeunG throws on non-positive-integer `γ` (DLMF 31.3.1 fails), on `a=0`
  and on `a=1` (singularity collisions); HeunC throws on
  non-positive-integer `γ` (`P_1=γ` vanishes). These fail-loud per
  CLAUDE.md Rule 1; they do not silently return a wrong value.

- **End-to-end numeric replication vs Mathematica oracle** (read-only
  Python re-implementation of the exact coded recurrences,
  `external/probes/heun-oracle/oracles.txt`): HeunG(a=2,q=1,α=1,β=2,γ=3,
  δ=4) at z=0.1 → 1.016981839745486 (oracle 1.01698183974548859...), at
  z=0.5 → 1.0908683419235652 (oracle 1.0908683419235652...); HeunC(q=1,
  α=1,γ=2,δ=3,ε=-1) at z=0.1 → 0.9471955542198833 (oracle
  0.94719555421988328...), at z=0.5 → 0.6121246449119607 (oracle
  0.61212464491196072...). Agreement to ~15-16 digits inside the
  convergence disc confirms the coefficient transcription is correct.
  (z=0.9 deviates at the 3rd digit, but that is expected near-boundary
  truncation of a slowly-converging series — Choun 2020,
  `docs/heun_research.md:262-268` — not a recurrence error; it is exactly
  why the module uses `epsilon_start=0.05` and the path-network for
  larger `|z|`.)
