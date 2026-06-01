# Bug sweep 2026-06-01 — Dispatcher + LatticeDispatcher (IVP↔BVP composition)

## Area

`src/Dispatcher.jl` (1D ordered-chain IVP↔BVP composition, FW 2011 §4.4)
and `src/LatticeDispatcher.jl` (2D per-row BVP fill flanked by IVP cells,
FW 2011 §4.4 + md:190). Special focus, per the assignment: the
IVP↔BVP endpoint handoff continuity, the smooth/pole classification
threshold that decides IVP-cell vs BVP-row, the per-row BVP fill flanked
by IVP cells, and strict-mode NaN behaviour.

The headline suspicion was a units/index/sign mismatch at the IVP↔BVP
handoff producing a literal discontinuity at every segment boundary, or
an in-place buffer-aliasing / index-shift error that produces an
*intermittent* (grid-geometry-dependent) discontinuity in the stitched
2D output.

Bottom line: **no demonstrated mis-transcription was found in these two
files.** The handoff orientation, the BVP boundary-value mapping, the
chain-rule derivative units, and the `up_grid` re-scatter index
alignment are all internally consistent and consistent with the cited
references. Two low-severity *documentation* inconsistencies were found
(they do not change behaviour). One genuinely fragile invariant (the
`_dilate_one` re-scatter alignment) was audited closely and found
**correct**, but it is the single place in this area where a future
edit could silently introduce exactly the intermittent-discontinuity
symptom; it is called out below as a risk note, not a present bug.

## References checked

- `references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:185`
  — eq. (3.2): `D₂u = (1/4)(z_b − z_a)²(6u² + …)`; the BVP affine scale
  factor `(z_b − z_a)²/4`.
- Same file, line 188 (post-eq.3.2 prose) — "enforce the BCs, `u_0 = u_b`
  and `u_N = u_a`"; the endpoint↔node orientation of the BVP.
- Same file, line 190 — "161 separate BVP solutions; one for each grid
  line" (the per-row BVP fill); and "≤ six iterations".
- Same file, line 192 — "The pole field solver … also provides
  derivative values `u_a'` and `u_b'`… derivatives at the endpoints can
  be estimated via `D₁u`… matched with those of the pole field solver.
  If the agreement is not to within some tolerance (typically set to
  10⁻⁷ or 10⁻⁸), an increase in N … is indicated." (the junction
  derivative-match diagnostic.)
- Same file, lines 249-261 (§4.1, Fig 4.1 caption) — the IVP↔BVP
  composition workflow: BVP over `[−20i, 20i]` gives `u, u'` at `z=0`
  and `z=20i`; "use this data to separately run out the two pole
  fields"; "fill in the area between pole field edges by further BVP
  solutions".
- Same file, line 401 (§5.5) — "smooth regions are unstable regions for
  any IVP solver, and the associated loss of accuracy ought not be
  carried back into the pole field… force the path selection algorithm
  to complete pole fields before stepping into smooth regions…"
- Same file, lines 202-208 (§3.2.2) — the 5-point Laplacian edge
  detector and the `log₁₀|Δu|` contour-at-0.001 classifier.
- `docs/adr/0017-lattice-dispatcher-strict-mode.md:44-79` — the v3
  edge-gated default + `strict::Bool` design; explicitly states
  `mask_used = gated.field_mask` is the *un-dilated* mask.
- `src/BVP.jl:268-297, 408-409, 486-505, 560-585` — BVP node/endpoint
  orientation, BC enforcement, `scale = (z_b−z_a)²/4`, the
  barycentric callable and its chain-rule derivative.
- `src/PathNetwork.jl:497-650, 665-693, 882-894` — Stage-1 walk
  termination (within `term_dist = h` of every target) and Stage-2
  input-order-preserving output (`grid_z = collect(grid)`), guaranteeing
  the terminus index and the dense-grid index alignment the Dispatcher
  relies on.
- `src/EdgeGatedSolve.jl:140-158, 217-233, 346-350` — `_dilate`,
  `_solve_targets` scatter order (`for j in 1:ny, i in 1:nx`), and the
  final solve over `_dilate(field, 1)`.
- `src/EdgeDetector.jl:166-195, 245-260` — interior-only stencil,
  boundary cells set NaN→mask=false.
- `src/Problems.jl:97-134` — `PadeTaylorProblem(f, y0, zspan)` with
  `y0 = (u, u')` tuple ordering.
- `test/dispatcher_test.jl`, `test/lattice_dispatcher_test.jl` — the
  pinned cosh closed-form cross-checks and the documented mutation-proof
  procedures (Mutations A/B/C, E/F, X).

## Findings

### [LOW] LatticeSolution docstring says `mask` is "dilated by one ring"; code uses the un-dilated `field_mask`

- **Location**: `src/LatticeDispatcher.jl:197-201` (docstring) vs
  `src/LatticeDispatcher.jl:354` (code `mask_used = copy(gated.field_mask)`).
- **Ground truth (cited)**: `docs/adr/0017-lattice-dispatcher-strict-mode.md:50-52`
  — "sets the partition `mask_used = gated.field_mask` (the **un-dilated**
  edge-confirmed mask — only edge-confirmed pole-field cells serve as BVP
  flanks; one-ring frontier cells are smooth-side)." Also
  `LatticeDispatcher.jl:346-354` (the inline comment block) correctly
  states the partition is the un-dilated `field_mask`, NOT the dilation.
- **Code behavior**: The code at line 354 stores the *un-dilated*
  `field_mask` as `mask_used`. The one-ring dilation (`_dilate_one`) is
  used only for the `up_grid` re-scatter (line 339), exactly as the ADR
  intends. So the code is correct per the ADR; only the `LatticeSolution`
  struct docstring (line 197-201) wrongly says "dilated by one ring".
- **Mechanism**: None — purely a documentation mismatch. A reader who
  trusts the struct docstring over the ADR could later "fix" the code to
  dilate `mask_used`, which the ADR explicitly rejected
  (`0017…md:95-98`, "One-ring-dilate field_mask as the partition.
  Rejected") because frontier cells were not edge-confirmed and using
  them as BVP flanks re-introduces stealth wrongness. That *would* be a
  real bug, so the stale docstring is a latent trap.
- **Intermittent?**: No.
- **Confidence**: 0.95 (direct doc-vs-ADR-vs-code contradiction, read
  in all three places).

### [LOW] Per-row loop comment claims "boundary rows have mask = false everywhere"; the grown `field_mask` can be true on boundary rows/cols

- **Location**: `src/LatticeDispatcher.jl:376-378`.
- **Ground truth (cited)**: `src/EdgeDetector.jl:96-100` and
  `:245-254` — only `laplacian_residual`/`pole_field_mask` *output* is
  forced false on boundary cells. But under the default path
  `mask_used = gated.field_mask` is the *grown region* `field`
  (`src/EdgeGatedSolve.jl:309-350`), which is seeded from cells within
  `seed_r` of the IC (`EdgeGatedSolve.jl:309-314`) and persists through
  growth (`mask .|= field`, `EdgeGatedSolve.jl:334`). Nothing forces
  `field` false on the grid boundary.
- **Code behavior**: If the IC seed or the grown field reaches `j=1` or
  `j=ny` (IC near the grid edge), `mask_used` is `true` on a boundary
  row. The per-row loop iterates only `j in 2:(ny-1)`
  (`LatticeDispatcher.jl:378`), so boundary-row cells are never bridged;
  pole-field boundary cells keep `:ivp`, smooth boundary cells keep
  `:ivp_only`. Behaviour is correct; the *comment*'s stated premise is
  not always true.
- **Mechanism**: None at present (the loop bounds make the premise
  irrelevant to behaviour). It is only a misleading comment.
- **Intermittent?**: No.
- **Confidence**: 0.7 (the premise's falsity is provable; that it never
  bites under current loop bounds is also provable).

## Areas verified correct

- **IVP→BVP left-BC handoff orientation is correct.**
  `Dispatcher.jl:310-326` calls `bvp_solve(f_bvp, ∂f_∂u_bvp, cur_z,
  CT(seg.z_end), cur_u, seg.u_b; …)`, i.e. `z_a=cur_z, z_b=seg.z_end,
  u_a=cur_u, u_b=seg.u_b`. `BVP.jl:209-210, 290-291` enforce
  `u(z_a)=u_a` at node `N+1` (`t=−1`) and `u(z_b)=u_b` at node `1`
  (`t=+1`), with the affine map `BVP.jl:274-277` giving `t=−1↔z_a,
  t=+1↔z_b`. So the BVP's left boundary value is exactly `cur_u` (the
  IVP terminus). At the junction `bvp_sol(cur_z)` hits `t*=−1`, the
  coincident-node guard (`BVP.jl:568-569`) returns `u_nodes[N+1]=u_a`
  exactly, so `Δu = |cur_u − u_a| = 0` by construction — matching the
  Dispatcher docstring (`Dispatcher.jl:30-32`) and FW md:188's
  "`u_N = u_a`". Mutation B in `test/dispatcher_test.jl:258-268` and
  Mutation E in `test/lattice_dispatcher_test.jl:451-456` both pin the
  no-swap orientation.

- **Junction derivative-match `Δu'` is unit-consistent (both are
  `du/dz`).** `Dispatcher.jl:335-337` compares `cur_up` (the IVP's
  `du/dz`, produced by `PathNetwork.jl:691` as
  `_evaluate_pade_deriv(...)/h_v`) against `up_bvp_left` from
  `bvp_sol(cur_z)`, which is `dudt_at_t / half_diff` with
  `half_diff = (z_b−z_a)/2` (`BVP.jl:488, 502`) — i.e. the chain-rule
  `du/dz = (du/dt)/((z_b−z_a)/2)`. Both quantities are derivatives with
  respect to `z`, so `Δu' = |cur_up − up_bvp_left|` is apples-to-apples.
  This is the FW md:192 diagnostic. No `1/h` or `2/(z_b−z_a)` factor is
  missing or doubled.

- **BVP→IVP handoff feeds the next IVP its IC from the BVP terminus.**
  `Dispatcher.jl:362-365` sets `new_cur_u, new_cur_up = bvp_sol(seg.z_end)`.
  At `z=z_b` the barycentric guard returns `u_nodes[1]=u_b=seg.u_b`
  exactly, and `up_at_end` is the spectral `du/dz` at `t*=+1`. The next
  IVP problem is then built `PadeTaylorProblem(prob_ivp.f, (cur_u,
  cur_up), (cur_z, seg.z_end))` (`Dispatcher.jl:277`) with the
  `(u, u')` tuple in the order `Problems.jl:125-128` expects. The
  recorded BVP→IVP junction match of `(0, 0)` (`Dispatcher.jl:378-382`)
  is *genuinely* zero — the IVP IC is literally the BVP terminus, so no
  discontinuity is introduced and none is masked. This is the FW Fig 4.1
  workflow (md:249-261: "use this data to separately run out the two
  pole fields"). Mutation A (`test/dispatcher_test.jl:249-256`) pins the
  `up_at_end` propagation.

- **IVP-segment terminal extraction index is correct and never NaNs
  under fixed `h`.** `Dispatcher.jl:282` sends
  `targets = vcat(seg.dense_grid, [seg.z_end])`; `path_network_solve`
  preserves input order in its outputs (`PathNetwork.jl:665`,
  `grid_z = collect(CT, grid)`), so `ivp_sol.grid_*[n_dense+1]` is the
  terminus at `seg.z_end`, and `grid_*[1:n_dense]` are the dense-grid
  cells (`Dispatcher.jl:293-304`). The Stage-1 walk terminates within
  `term_dist = h_T` of every target (`PathNetwork.jl:521-523`), so the
  nearest visited node to `seg.z_end` is ≤ `h_T` away and Stage-2's
  disc-radius check (`PathNetwork.jl:685`) passes — the terminus does
  *not* receive `NaN` under the default `:fixed` policy. Mutation C
  (`test/dispatcher_test.jl:270-278`) pins the dense-grid
  concatenation.

- **BVP affine scale factor matches FW eq. (3.2).** `BVP.jl:294`
  `scale = (z_b_CT − z_a_CT)^2 / 4` matches
  `FW2011…md:185` `D₂u = (1/4)(z_b − z_a)²(…)`. (BVP internals are
  another auditor's primary scope; verified here only because the
  Dispatcher composes `bvp_solve` and the diagnostic depends on it.)

- **`up_grid` re-scatter index alignment is correct (the fragile
  invariant).** `LatticeDispatcher.jl:339-345` walks
  `final_targets = _dilate_one(gated.field_mask)` in `for j in 1:ny,
  i in 1:nx` order, incrementing `k` and assigning
  `up_grid[i,j] = pn_sol.grid_up[k]`. For this to be correct,
  `_dilate_one(field_mask)` must reproduce *exactly* the mask
  `_dilate(field, 1)` that `EdgeGatedSolve` used in its final
  `_solve_targets` call (`EdgeGatedSolve.jl:347-348`), and the scatter
  order must match `_solve_targets`'s build order
  (`EdgeGatedSolve.jl:222`, also `for j in 1:ny, i in 1:nx`). Both
  conditions hold: `_dilate_one` (`LatticeDispatcher.jl:447-462`) and
  `_dilate(mask, 1)` (`EdgeGatedSolve.jl:140-158`) are byte-for-byte
  identical single-ring 8-connected dilations, including identical
  boundary clamping `1 ≤ ii ≤ nx && 1 ≤ jj ≤ ny`, and both read
  neighbours from the original snapshot. `gated.field_mask === field`
  is returned by `EdgeGatedSolve.jl:350`. Therefore the `k` counter
  walks exactly `length(pn_sol.grid_up)` cells with perfect physical-cell
  alignment, and `u_grid` (copied directly from `gated.u_grid`,
  `LatticeDispatcher.jl:331`) is scattered by the *same* mask/order, so
  `u_grid[i,j]` and `up_grid[i,j]` refer to the same cell. RISK NOTE
  (not a present bug): this alignment depends on the two private
  `_dilate*` helpers staying identical and on the column-major scatter
  order matching. If a future edit changed `_dilate_one`'s connectivity,
  ring count, or iteration order — or `EdgeGatedSolve._solve_targets`'s
  build order — the `k` index would desynchronise and *every*
  `up_grid` cell past the first divergence would carry a neighbour's
  derivative: a silent, grid-geometry-dependent (hence intermittent)
  discontinuity in `u'`. The code is correct today; this is the single
  highest-leverage place in the area to guard with a regression test
  that asserts `up_grid` agrees with a direct
  `path_network_solve` re-scatter on the auto path.

- **Per-row smooth-run scan: no off-by-one, no infinite loop, no flank
  aliasing.** `LatticeDispatcher.jl:378-437`. For fixed interior row
  `j`, `i` advances from 2; a mask=true cell advances `i` (lines
  381-384); a smooth run `[run_start, run_end]` is found with the inner
  `while` always advancing `i` at least once (entry condition
  `!mask_used[i,j]` is true), so `i` strictly increases — no infinite
  loop. After the inner loop `i = run_end+1 = right_flank`, which is
  mask=true, so the next outer iteration advances past it — no double
  processing. Flanks `left_flank=run_start−1`, `right_flank=run_end+1`
  with the `2 ≤ … ≤ nx−1` interior guard (line 394) feed
  `z_a/z_b/u_a/u_b` in the same left=smaller-x, right=larger-x
  orientation `bvp_solve` expects (lines 397-400). Flank cells are
  mask=true and the BVP only overwrites mask=false cells `run_start:run_end`
  (lines 427-433), so a prior run's BVP fill can never corrupt a later
  run's flank BC — no in-place aliasing. Fixing `j` and varying `i`
  (real axis) is a *horizontal* grid line, matching FW md:190 "one BVP
  per grid line".

- **Smooth/pole classification threshold is the documented FW detector,
  consumed self-consistently.** Under the default path the partition is
  the edge-gated `field_mask` (`EdgeGatedSolve` region-growing over the
  `pole_field_mask` of FW eq. 3.3, `EdgeDetector.jl:213-260`); under the
  manual path it is the caller-supplied `mask`. The
  `grid_mat[i,j]=xs[i]+im·ys[j]` convention (rows=x, cols=y,
  `LatticeDispatcher.jl:304`) is the transpose of `EdgeDetector`'s
  *stated* convention (`EdgeDetector.jl:166`, rows=y, cols=x), but the
  5-point Laplacian is isotropic and every matrix inside the auto path
  (`u_grid`, `field_mask`, `up_grid`, the per-row scan) uses the *same*
  rows=x/cols=y convention, so the transpose is self-consistent and does
  not change which cells are classified pole-field.

- **Strict / fail-soft NaN behaviour is correct.** Default
  `strict=true` rethrows the `bvp_solve` non-convergence
  `ErrorException` verbatim (`LatticeDispatcher.jl:410-423`,
  `ADR-0017:59-79`); `strict=false` catches *only* the message
  `"bvp_solve: Newton did not converge"` (matching `BVP.jl:332-338`) and
  tags the run `:bvp_fail` while leaving those cells at their input
  value (NaN under the edge-gated default — the honest Rule-1 signal).
  Converged-`:bvp` cells are always finite (tested,
  `lattice_dispatcher_test.jl:289-293`). Any other exception rethrows
  (`LatticeDispatcher.jl:420-421`). Mutation X
  (`lattice_dispatcher_test.jl:473-493`) pins the gate polarity. The
  Dispatcher's own strict mode (`Dispatcher.jl:340-349`) throws only on
  `Δu > tol` or `Δu' > tol`; the BVP→IVP `(0,0)` junctions never trip it
  (`Dispatcher.jl:381`).
