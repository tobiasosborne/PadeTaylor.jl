# Bug sweep 2026-06-01 — Area E5: comparison-operator & tolerance bugs across `src/`

## Area

WHOLE `src/`: comparison-operator and tolerance bugs — strict-vs-nonstrict at
decision boundaries (step acceptance, rank cutoff, dense-output segment
ownership, walk termination), abs-vs-rel tolerance, float equality / `≈`,
NaN gates (NaN slipping through a `<` gate), and the tolerance constants
(`1e-14`, `1e-15`, `√ε`, `τ = tol·‖c‖`) versus their reference values.

Method: read every comparison/tolerance site in the Padé/Taylor/step core
(`Coefficients`, `VectorCoefficients`, `RobustPade`, `SharedPade`, `LinAlg`,
`StepControl`, `VectorStepControl`, `PadeStepper`, `VectorStepper`,
`VectorWedgeStep`) and the path-network solvers (`PathNetwork`,
`VectorPathNetwork`, `VectorPathNetworkStage2`), plus the branch-cut /
sheet sites (`CoordTransforms`, `SheetTracker`) and the pole extractors
(`PoleField`, `VectorPoleField`). Each comparison was diffed against the
canonical references below.

## References checked (path:line cites)

- `external/chebfun/padeapprox.m:36` — default `tol = 1e-14` (relative).
- `external/chebfun/padeapprox.m:49-50` — `tc = 1e-15*norm(f)`, discard `abs(f) < tc` (strict `<`).
- `external/chebfun/padeapprox.m:66` — absolute SV tolerance `ts = tol*norm(c)`.
- `external/chebfun/padeapprox.m:69` — special-case `r ≡ 0`: `norm(c(1:m+1),inf) <= tol*norm(c,inf)` (`<=`).
- `external/chebfun/padeapprox.m:93` — rank test `rho = sum(svd(C) > ts)` (strict `>`).
- `external/chebfun/padeapprox.m:95-101` — break on `rho == n`; else `m -= (n-rho); n = rho`. (For the `n×(n+1)` `C`, `rho ≤ n`, so `rho==n` ⇔ `rho≥n`.)
- `external/chebfun/padeapprox.m:112-113` — `D = diag(abs(b)+sqrt(eps))`, QR of `(C*D)'` (adjoint; documents the historical `.'` transpose bug).
- `external/chebfun/padeapprox.m:123,130,134` — trim with `abs(b)>tol` (leading/trailing-b), `abs(a)>ts` (trailing-a). `tol` vs `ts` distinction.
- `references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md:76` — `εₖ = cₖ + Σ_{r=1..ν} b_r c_{k-r}` (sign `+`, index `k-r`).
- `.../FFW2017...md:80` — `T(h) ≈ |ε_{n+1} h^{n+1} / p(h)|` (`p` = numerator).
- `.../FFW2017...md:82,93` — accept on `T(h) ≤ Tol`; rescale while `T(h) > Tol`.
- `.../FFW2017...md:43,47` — P̃_III / P̃_V transformed equations.
- `references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/JorbaZou2005_taylor_IVP_package_ExpMath14.md:617,629,633` — `pₙ = -½ln εₙ + 1`, `ρₙ = min{ρ̄_{pₙ-1}, ρ̄_{pₙ}, ρ̂_{pₙ-1}, ρ̂_{pₙ}}`, `hₙ = ±ρₙ/e²`.
- `docs/v0p2_pillarA_hermite_pade_findings.md:210-215,222-227,240-243,422-427` — the shared-Q rank rule: A_full is `dn×(m+1)`; for d≥2 the system is over-determined and "consistent only generically"; the implementation **must check σ₂/σ₁** (isolation), not σ₁ alone; "If ρ < m: reduce m to ρ".

## Findings

### [HIGH] `shared_denominator_pade` accepts a full-column-rank A_full (no null space) as if it had an isolated 1-D null space — break uses `ρ ≥ m_cur` instead of `ρ == m_cur`, and the isolation guard is dead code for d ≥ 2

- **Location:** `src/SharedPade.jl:196` (`ρ = count(s -> s > τ, S)`), `src/SharedPade.jl:208` (`if ρ ≥ m_cur`), `src/SharedPade.jl:213` (`m_cur = ρ`), `src/SharedPade.jl:219` (`b_init = Vector{T}(Vt[end, :])`), `src/SharedPade.jl:224-230` (Failure-mode-4 guard `n_near = count(s -> s ≤ τ, S); n_near > 1 && throw`).
- **Ground truth (cited):**
  - `docs/v0p2_pillarA_hermite_pade_findings.md:210-215`: "A_full is dn×(m+1) … The relevant comparison is whether **σ_min is small relative to σ₂** … a large gap σ₂/σ_min ≫ 1 indicates a genuine one-dimensional null space (unique shared denominator up to scale)."
  - `:222-226`: "**Subtlety not present in scalar case:** for d ≥ 2 … the system is overdetermined at n = m and consistent only generically."
  - `:240-243` (Conclusion): "the implementation must check the ratio σ₂/σ₁ (smallest two singular values) rather than just σ₁ alone."
  - `:422-424` (Step 4): "Let ρ = count(S .> τ). If ρ < m: reduce m to ρ … re-run." plus `:426`: "If σ₂/σ₁ < 10 (null space poorly isolated): warn but proceed."
  - Scalar oracle `external/chebfun/padeapprox.m:95` breaks on `rho == n`, which for the `n×(n+1)` matrix can never overshoot (`rho ≤ n`).
- **Code behavior:** For a `d`-component vector ODE the matrix `A_full = vcat(blocks)` is `(d·m_cur) × (m_cur+1)`. For `d ≥ 2` and `m_cur ≥ 1` it has `m_cur+1` rows-worth of singular values (`min(d·m_cur, m_cur+1) = m_cur+1`), so `ρ = count(s > τ)` ranges `0 … m_cur+1`. The break `if ρ ≥ m_cur` (line 208) therefore fires for **two** distinct cases:
  1. `ρ == m_cur` — exactly one SV at/below τ, an isolated 1-D null space → `Vt[end,:]` is a genuine annihilating denominator. Correct.
  2. `ρ == m_cur+1` — **all** SVs above τ, A_full is full column rank, there is **no** null space. `b_init = Vt[end,:]` (line 219) is then the right-singular vector of `S[end] > τ`, i.e. `A_full·b_init = S[end]·u_min` with norm `> τ ≠ 0`. The Hermite–Padé system is inconsistent at this degree, yet a "denominator" is silently produced from the least-dominant (non-null) direction.
  The Failure-mode-4 guard at line 224-225 cannot catch case 2: after the break `ρ ≥ m_cur`, so `n_near = count(s ≤ τ) = (m_cur+1) − ρ ≤ 1`, hence `n_near > 1` is **never** true post-loop — the guard is dead code for the over-determined d≥2 matrix. (The comment at `:222` — "if ρ > m_cur the smallest σ is degenerate (multiple near-equal σ at the bottom)" — is inverted: ρ > m_cur means `n_near = 0`, i.e. *no* null vector, not "multiple near-equal σ.") No σ₂/σ₁ isolation check exists anywhere in the function, contrary to the cited Conclusion.
  For `d == 1` the matrix is `m_cur × (m_cur+1)` with only `m_cur` SVs, so `ρ ≤ m_cur` and `ρ ≥ m_cur` collapses to `ρ == m_cur` — the d=1 reduction stays bit-identical to the scalar `padeapprox.m`, which is why the documented d=1 oracle (SP.1.1) does not expose this.
- **Mechanism (intermittent discontinuity):** `shared_denominator_pade` is the heart of the *vector* solver (`VectorStepper.vector_pade_step_with_pade!` → `VectorWedgeStep` / `VectorPathNetwork` Stage-1 walk → every `d ≥ 2` solve, including the Painlevé systems carried as first-order systems). Whether the over-determined A_full's smallest singular value sits just below or just above `τ = tol·‖c‖` is purely data-dependent and drifts continuously as the walk steps along the trajectory and the local jet changes. While `σ_min/‖c‖ < tol` the code takes the isolated-null-space branch (a true annihilating Q); the instant the ratio drifts to `σ_min/‖c‖ > tol` the same code silently switches to extracting Q from a *non-null* singular direction — a different, discontinuous denominator polynomial. The shared Q is then evaluated at `t = 1` to advance the state (and its roots drive `_adaptive_h`, `_candidate_pole_disc`, and the Stage-2 `R_gate`), so the computed `y` and the chosen step jump at the threshold crossing. This is exactly an intermittent, threshold-straddling discontinuity in the computed solution, off in the scalar code and on only for d≥2 and only near the tolerance boundary.
- **Intermittent?** Yes — fires only when the over-determined A_full's smallest SV straddles `τ = tol·‖c‖` (data/step-dependent), and only for `d ≥ 2`.
- **Confidence:** 0.78. The break-condition divergence from the cited ground truth (`ρ ≥ m_cur` vs the reference's `ρ == m_cur` / "If ρ < m reduce" + mandatory σ₂/σ₁ isolation check) and the dead-code guard are demonstrated against `docs/v0p2_pillarA_hermite_pade_findings.md:210-243,422-427`. Not 0.9+ because I could not run code to construct a concrete jet that lands `ρ = m_cur+1` (read-only), so the in-practice firing frequency is inferred, not measured.

### [LOW] Adaptive controller: a NaN truncation-error estimate slips through the `T(h) > Tol` rescale gate and is accepted as a converged step

- **Location:** `src/PadeStepper.jl:569` (`while Th > adaptive_tol`), with `Th` from `ffw_truncation_error` (`src/PadeStepper.jl:486`, `return abs(eps_nplus1 / a_at_one)`). Same shape in `ffw_rescale_q` (`:514` `T_h == 0 && return Inf`).
- **Ground truth (cited):** `.../FFW2017...md:82,93` — "Suppose T(h) is greater than some specified tolerance Tol, then we must rescale … This is repeated until T(h) ≤ Tol." CLAUDE.md Rule 1 — never accept a NaN/Inf as a truthful value.
- **Code behavior:** `Th = abs(eps_nplus1 / a_at_one)`. If the RHS evaluation produced a non-finite Taylor coefficient (so `eps_nplus1` is `NaN`/`Inf`) and `a_at_one` is finite-nonzero, `Th = NaN`. The loop guard `Th > adaptive_tol` is then `false` (any comparison with NaN is false), so the `while` exits immediately and the step is **accepted** with `meta.T_h = NaN` — a too-large/garbage step passed off as converged. The `a_at_one == 0` case is guarded (`:481`), but a NaN/Inf `eps_nplus1` over a finite denominator is not.
- **Mechanism:** A pixel/segment where the local jet briefly goes non-finite (near a pole, or a transient overflow at large `|z|`) yields `Th = NaN`, the gate passes, and the accepted step writes a non-finite or wrong state that then propagates — a localized, data-dependent discontinuity. Requires upstream non-finiteness, so it is a defensive gap rather than an always-on error.
- **Intermittent?** Yes (only when an upstream coefficient goes non-finite).
- **Confidence:** 0.4. Real NaN-through-`<`-gate hazard, but contingent on upstream non-finite input that other fail-loud checks may catch first.

### [LOW] Stage-2 vector fill divides by `q_t` with no near-zero denominator guard (inconsistent with both the scalar and the Stage-1 vector path)

- **Location:** `src/VectorPathNetworkStage2.jl:437-438` (`q_t = _eval_poly(visited_denominator[idx_v], t); grid_y[i] = C[_eval_poly(num, t) / q_t for num in numerators]`).
- **Ground truth (cited):** scalar oracle `src/PadeStepper.jl:376` (`iszero(den) && throw(DomainError…)`) and the Stage-1 vector gate `src/VectorStepper.jl:258` (`abs(q_at_one) > eps_T * norm(denominator) || throw(DomainError…)`). CLAUDE.md Rule 1.
- **Code behavior:** Three different denominator-zero policies coexist — scalar stepper guards the exact `iszero` case, the Stage-1 vector stepper guards the *relative* `|Q(1)| ≤ √ε·‖Q‖` case, and the Stage-2 fill guards **nothing**, relying solely on the geometric B1 gate `abs(z_f − z_v) > R_gate` (`:429`) with `R_gate ≤ 0.5·h_v·min|t*|` (`:328`) to keep grid points off poles. A grid point inside `R_gate` but where `q_t` is nonetheless tiny would produce a silent `Inf`/`NaN` pixel.
- **Mechanism:** In practice the `0.5·h_v·min|t*|` clamp keeps `t` at ≤ half the distance to the nearest root, so `q_t` is bounded away from zero; the inconsistency only bites if the B1 radius is mis-estimated. Genuine but low-probability discontinuity source.
- **Intermittent?** Yes (grid-point/data-dependent), but suppressed by the geometric gate in normal operation.
- **Confidence:** 0.3.

### [LOW] Pole-field clustering sort is not order-stable, so equal-`|t*|` candidates pick a non-deterministic cluster representative

- **Location:** `src/PoleField.jl:156` (`sort!(candidates; by = c -> c[2])`) feeding `:159-167` (greedy first-in-cluster representative).
- **Ground truth (cited):** CLAUDE.md determinism contract (RobustPade.jl:77-84 "bit-identical given the platform fingerprint"); the path-network `_nearest_visited` tie-break (`PathNetwork.jl:880-894`) shows the intended lexicographic-tie-break discipline this site omits.
- **Code behavior:** `sort!` defaults to an unstable algorithm; two candidates with exactly equal `|t*|` may swap order between runs, and since the greedy loop makes the *first* candidate to land in a cluster its representative (`:161-163`), the reported pole location can shift. Affects the diagnostic pole-field output, not the computed `(u, u')` solution.
- **Mechanism:** Diagnostic figures (pole scatter) could show a non-reproducible jitter at exact-tie roots; not a discontinuity in the solution itself.
- **Intermittent?** Yes (run-to-run on ties), but confined to a diagnostic accessor.
- **Confidence:** 0.3.

### [LOW] `step_pade_root` discards a pole that projects to exactly `t == 0` (strict `t > 0`)

- **Location:** `src/StepControl.jl:255` (`if t > 0`).
- **Ground truth (cited):** FW 2011 §3.1 pole-distance heuristic (StepControl.jl:84-100 docstring); the dense-output gates use the inclusive convention `abs(z−z_v) > h_v` (`PathNetwork.jl:1053,1095`), i.e. the disc boundary is *included*.
- **Code behavior:** A denominator root whose forward projection is exactly `0` (a pole sitting on the current point, perpendicular component aside) is treated as "behind/sideways" and contributes no cap, so the step is not shortened for it. This is a measure-zero boundary and the scalar stepper's own `iszero(den)` guard would fire on a true on-point pole, so the practical exposure is negligible.
- **Mechanism:** Could permit one over-long step exactly when a pole projects to the origin of the step direction; vanishingly rare.
- **Intermittent?** Yes, but measure-zero.
- **Confidence:** 0.2.

## Areas verified correct

- **`RobustPade.robust_pade` SVD path — comparison/tolerance operators match `padeapprox.m` exactly.** Special case `≤` (`RobustPade.jl:397` vs `padeapprox.m:69` `<=`); rank test strict `>` (`:418` `s > ts_typed` vs `:93` `> ts`); break `ρ == n` (`:420` vs `:95`); leading/trailing-b trim uses `tol` (`:463,479` vs `:123,130`); trailing-a trim uses `ts` (`:484` vs `:134` — the `tol`-vs-`ts` distinction is preserved); `ts = tol*norm(c)` absolute tolerance (`:385` vs `:66`); `default_tol(Float64)=1e-14` (`:154` vs `:36`); QR-reweight `sqrt(eps)` and `adjoint(C*D)` (`:441-443` vs `:112-113`, including the documented adjoint-not-transpose fix).
- **`Coefficients` Taylor recurrences — factorial/coefficient relations correct.** 1st order `y_k = f_{k-1}/k` (`Coefficients.jl:129`); 2nd order `u_j = f_{j-2}/(j(j-1))` and `up_{j-1} = j·u_j` (`:198-199`), matching FW 2011 §2.1.2 method (b); the `up` resync maintains the read-invariant for `f_t[j-2]`.
- **`VectorCoefficients.vector_taylor_coefficients`** — component-wise lift `y[i][k] = f_t[i][k-1]/k` (`VectorCoefficients.jl:166`) identical to the verified scalar relation.
- **`PadeStepper` adaptive controller error formula** — `ε_{n+1} = c̃_{n+1} + Σ b_r c̃_{n+1-r}` (`PadeStepper.jl:472-477`) matches FFW md:76 sign (`+`) and index (`k-r`); `T(h)=|ε/a(1)|` over the numerator `P_u.a` (`:480,486`) matches md:80 (`p(h)` = numerator); accept `T(h) ≤ Tol` ⇔ loop `Th > Tol` (`:569`) matches md:82,93; `q = (k·Tol/T(h))^{1/(n+1)}`, `k = 1e-3` (`:515` and the `k_conservative` default `:522`) matches md:88-89.
- **`StepControl` / `VectorStepControl` Jorba–Zou** — fixed-order TI.jl reduction `h = min_{k∈{p-1,p}} (ε/|c_k|)^{1/k}`, ε resolution `eps_abs ≥ eps_rel·|c₀|`, zero-skip, and the `_second_stepsize` fallback over `1≤j≤p-2` are faithfully ported (`StepControl.jl:177-204`; `VectorStepControl.jl:187-214`); the documented `/e²` absorption into the ε substitution is consistent with JorbaZou md:617-633; the vector module reduces bit-identically to the scalar at d=1.
- **Dense-output / segment-ownership gates are nearest-node Voronoi with a consistent inclusive disc.** `PathNetwork.eval_at` / `eval_at_sheet` (`:1053,1095`) and `VectorPathNetworkStage2._stage2_fill` (`:429`) all use `abs(z − z_v) > R/h_v` (strict `>`, boundary included), and a single nearest owner is chosen (`_nearest_visited` with a lexicographic Re/Im tie-break, `PathNetwork.jl:880-894`, `VectorPathNetwork.jl:921`), so there is no two-segments-claim-one-point join discontinuity from an operator mismatch.
- **Walk termination & coverage skip** — `while abs(z_cur − target) > real(visited_h[parent])` (`VectorPathNetwork.jl:798`) with `parent` correctly advanced to the freshly pushed node (`:836`); `skip_covered` uses `≤ COVER_FRAC·real(visited_h[i])` (`:765`) and the default exact-coincidence skip `≤ 10·eps(T)` (`:769`) — consistent, no NaN exposure (distances are finite).
- **`VectorWedgeStep` selection determinism** — `@threads` writes `evals[k]` by exactly one thread per k (no aliasing), results ordered by k; lexicographic `(primary, −|wedge_angle|)` max with `score > best_score` and a `best_k == 0` prime (`VectorWedgeStep.jl:432-443`); `_candidate_pole_disc` / `_adaptive_h` use the same deterministic `minimum(abs, roots)` measure. No RNG/order intermittency in candidate selection.
- **`CoordTransforms` transformed RHS and branch convention** — P̃_III (`:154`) and P̃_V (`:199-202`) match FFW md:43,47 term-by-term (β e^ζ, δ e^{2ζ}/w; the PV `(w−1)` / `(w+1)/(w−1)` structure); principal-branch `log` used consistently with documented sheet offsets.
- **`SheetTracker` winding** — `winding_delta` normalises to `(−π, π]` with consistent `Δθ ≤ −π` / `Δθ > π` boundaries (`:286-289`); `sheet_index = round(Int, total/2π)` (`:320`) matches FFW md:187-189; internally consistent, no off-by-one at the wrap.
- **`SharedPade` d=1 reduction and the rest of the tolerance arithmetic** — `τ = tol·‖vcat(jets)‖` (`SharedPade.jl:182-183`) matches the scalar `ts = tol·‖c‖`; the trailing-trim `findlast(|b|>tol)` / `findlast(|a|>τ)` (`:254,263`) mirror `padeapprox.m:130,134`; the d=1 stack is bit-identical to the scalar `C̃`. (The d≥2 rank-break is the HIGH finding above; everything else in this module checks out.)
