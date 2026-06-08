# PadeTaylor.jl Ground-Truth Corpus — Extension Plan v2 (pole-zoo + path-network)

This is the **Phase-1 catalogue** for the *second* ground-truth corpus epic
(`padetaylor-p1v2`), extending the v1 corpus (`padetaylor-p1v0`,
`docs/test_corpus/00_capability_map.md` + `01_corpus_catalogue.md`).

The v1 corpus proved the numerical core clean across ~691 assertions but is
**overwhelmingly 2nd-order scalar + small vector systems in the FW/Painlevé
world**. It has dense coverage of *isolated* poles but almost no coverage of
**pole *morphology variety*** (order, arrangement, periodicity, lattice
geometry) and no first-order equations with movable poles framed as such.
This extension fills that.

## Provenance

Produced by **five read-only derivation/verification subagents** (Rule 7: no
Julia invoked — Python/sympy/mpmath/scipy only, the project's sanctioned oracle
stack). **Every exact solution, RHS recast, pole/branch location, and residue
below was independently verified** (sympy symbolic residual ≡ 0 and/or mpmath
`odefun` integration matching the closed form to ≤1e-40 at dps=50–80). Where a
subagent's verification *contradicted* the design brief, the correction is
recorded — skepticism (CLAUDE.md Rule 3) caught five such cases plus one new
library bug.

## The organizing principle: the blind-spot matrix

The corpus lives in a 3-axis space **{pole order} × {pole arrangement} ×
{residue/Laurent known}**. v1 occupies a narrow corner. This extension
systematically sweeps the rest:

| Pole order | v1 coverage | This extension |
|---|---|---|
| 1 (simple) | tan, coth, 1/(1−z), sn (isolated) | **real rows** (Airy/Bessel/Hermite/Chebyshev zeros), **periodic vertical row** (logistic) |
| 2 (double) | equianharmonic ℘ | **lemniscatic ℘** (square lattice), Jacobi cn/dn |
| 3,4,5+ | **none** | order-k linear `(1−z)^{−k}`, ℘′ (triple), mixed-order vector |
| branch pt | only via PIII/PV/PVI transforms | **elementary** √(C−z²) finite cut, log, atanh; algebraic PVI |
| essential | **none** | `e^{1/z}` fail-loud target (**found a bug**) |

| Arrangement | v1 | This extension |
|---|---|---|
| 2D lattice | equianharmonic ℘ only | **lemniscatic** (square) + lattice-vs-lattice routing contrast |
| horizontal periodic | tan/coth | — |
| **vertical periodic** | none | **logistic** (2πi-spaced) |
| **finite/semi-infinite real row** | none | **Airy/Bessel/Hermite/Chebyshev-zero rows** |
| wedge-confined | tritronquée (qualitative) | **quantitative** far-target value + zero-pole pins |

---

## Key discoveries & corrections (read before building)

1. **NEW BUG — essential-singularity silent bridge (Rule-1 violation).** `solve_pade`
   on `u''=u(1+2z)/z⁴` (exact `u=e^{1/z}`) integrates toward and *across* the
   essential singularity at z=0 returning finite, plausible values (relerr
   1e-60 → 1e-11; it even "bridges" z=0 to 11 digits) with **no throw, no
   NaN**. The meromorphic-only contract (Rule 1, fail-loud) has **no guard**:
   nothing inspects the local-jet super-geometric coefficient growth
   (`|c_k|^{1/k}` drift) that signals non-meromorphy. → bug bead
   `essential-singularity-silent-bridge` (P3, out-of-class input; maintainer
   may downgrade to documented-limitation). CFail.1 ships an `@test_broken`
   that flips green when a diagnostic lands.

2. **CRic.5 has no real poles (brief was wrong).** The recessive solution of
   `w''=z²w` is the parabolic-cylinder `U(0,√2 z)`, **non-oscillatory** → zero
   real zeros. `u'=u²−z²` has poles only as **complex-conjugate pairs** near
   `arg≈±3π/4`, residue −1. The real-row analog is the *oscillatory* sister
   **`u'=u²+z²`** (`w=√z·J_{1/4}(z²/2)`, real row on the positive axis) — but
   that overlaps the existing CPB.5; keep CRic.5 as a complex-pole record and
   defer the sister (CRic.6) pending a dedup decision with CPB.5.

3. **`extract_poles` reports locations only — no multiplicity/order field**
   (`src/PoleField.jl:229-314`, `src/VectorPoleField.jl:245`). An order-k pole
   plants **exactly k near-coincident Q-roots** (measured cluster diameter:
   3.4e-8 at k=2 → 3.5e-3 at k=5; always ≪ the ~0.9 gap to the next root) that
   `cluster_atol` collapses to one location. A "sees order k" assertion must
   read the **pre-collapse root count**, not a returned order. → follow-on API
   bead `extract-poles-cluster-multiplicity` (P3).

4. **Pure-power `(1−z)^{−k}` makes the classical Float64 Toeplitz solve
   singular** (the exact `[k/0]` Padé *is* the function). CHop.1 must route
   through the SVD/RobustPade path with the pole interior (`z0=0, h=1.5`).

5. **℘ oracle precision.** The package has no native ℘; the v1 oracle is
   Mathematica `WeierstrassP[z,{g2,g3}]` strings. A validated mpmath
   ℘-via-Jacobi-sn recipe exists (`℘ = e₃+(e₁−e₃)/sn²(z√(e₁−e₃),m)`), but the
   **equianharmonic** (complex-root) case needs **dps≥80** for 45-digit
   agreement; the **lemniscatic** (real-root) case is dps=50-clean. Capture
   equianharmonic ℘/℘′ from Mathematica or at dps≥80.

6. **Algebraic PVI branch-at-zero limitation.** `1/√z` solves PVI with
   `(α,β,γ,δ)=(0,0,1/8,−5/8)` and `(√z+1)/2` with `(1/2,0,1/8,−5/8)` (both
   symbolic residual ≡ 0) — the first **non-self-referential** PVI oracles in
   the suite (v1 PVI values are all self-pinned, gap #9). But their branch
   point z=0 coincides with PVI's own fixed singularity, so they do **not**
   test circumambulation of a *generic* (movable, off-{0,1,∞}) branch. A
   Dubrovin–Mazzocco / Lisovyy–Tykhyy algebraic solution with a tabulated
   named-sheet value at an off-singular z₀ is **not in-repo** → reference-
   acquisition bead (P3, blocked).

---

## New test files & bead map

All beads are children of epic **`padetaylor-p1v2`**. Discipline unchanged:
one Julia process at a time (Rule 7), every assertion pinned to independent
ground truth (Rule 5), every load-bearing assertion mutation-proven (Rule 4),
`src/` never touched. Oracle-capture scripts go under
`external/probes/corpus-oracles/<area>/capture.py` (mpmath/sympy, dps=50 unless
noted, residual==0 gate before any value is emitted).

| # | Test file | IDs | Pri | Catches bug | Oracle |
|---|---|---|---|---|---|
| 1 | `corpus_riccati_rational_test.jl` | CRic.3,4 | P1 | — | exact rational (sympy) |
| 2 | `corpus_riccati_special_test.jl` | CRic.1,2,5 | P2 | — | mpmath airyaizero/besseljzero/PCF-root |
| 3 | `corpus_higher_order_pole_test.jl` | CHop.1,2,3 | P1 | — | exact rational + Mathematica ℘′ |
| 4 | `corpus_elliptic_lattice_test.jl` | CEl.1,2,3 | P1 | — | mpmath ellipfun + ℘ (dps≥80/MMA) |
| 5 | `corpus_periodic_pole_test.jl` | CVrow.1 | P1 | — | closed form + mpmath dps=50 |
| 6 | `corpus_elementary_branch_test.jl` | CBr.1,2,3 | P1 | **61um** | sympy/mpmath dps=50 |
| 7 | `corpus_algebraic_pvi_test.jl` | CPvi.1,2 | P2 | — | mpmath dps=50 (±sheet) |
| 8 | `corpus_out_of_class_test.jl` | CFail.1 | P2 | **new** | mpmath e^{1/z} dps=50 |
| 9 | `corpus_orthopoly_bvp_test.jl` | CBvx.1 | P1 | — | exact rational (sympy) |
| 10 | `corpus_special_fn_bvp_test.jl` | CBvx.2,3 | P2 | — | scipy mathieu + mpmath whitm |
| 11 | `corpus_oblique_guard_test.jl` | CBvx.4 | P2 | **53tu** | mpmath cosh dps=30 |
| 12 | `corpus_pathnet_walls_rows_test.jl` | CPN.1,2,4,8 | P1 | 53tu (inc.) | mpmath airy/logistic/rational |
| 13 | `corpus_pathnet_lattice_sectors_test.jl` | CPN.3,5,6 | P1 | 53tu (inc.) | odefun ray-integration |
| 14 | `corpus_pathnet_winding_test.jl` | CPN.7 | P1 | **61um** | mpmath sqrt/log sheet values |

Follow-on (not test files): `essential-singularity-silent-bridge` (bug, P3),
`extract-poles-cluster-multiplicity` (API, P3), `algebraic-pvi-generic-branch-
oracle` (P3, blocked on reference acquisition),
`vector-bvp-coupling-bc-oracle` (CVB.4, P3, residual gap #10).

---

## Family A — Riccati = log-derivative of a linear ODE (simple poles at special-function zeros)

For any `w″+p·w′+q·w=0`, `u=w′/w` solves `u′=−u²−p·u−q` with a **simple pole at
every zero of `w`, residue +1** (or −1 for `u=−w′/w`). Special-function zeros
are tabulated to arbitrary precision → exact pole locations *and* residues. This
is the corpus's first **rows of simple poles**. First-order RHS recast to
2nd-order for `solve_pade` (`u''=∂R+(∂R/∂u)·R`); off-axis / negative-axis rows
use `path_network_solve`+`eval_at`; `extract_poles` pinned to the tabulated
zeros.

- **CRic.1 Airy** — `u'=u²−z`, recast `u''=2u³−2uz−1`, `u=−Ai′/Ai`. Poles at
  Airy zeros `a_k` (semi-infinite **negative** real row), residue −1. Neg-axis
  is z-decreasing → `solve_pade` can't reach it; use `path_network_solve`.
  Oracle `mp.airyaizero(k)`. Tol: pos-axis value 1e-11, neg-axis bridge 1e-10,
  pole-extract 1e-7. *Distinct from PII-Airy (CPR.3 is scaled, +3uz+3/2) and
  CSF.1 (`Ai` itself).*
- **CRic.2 Bessel** — `u'=−u²−u/z−1`, recast `u''=2u³+3u²/z+2u+2u/z²+1/z`,
  `u=−J₁/J₀`. Poles at `j_{0,k}` (**positive** real row), residue +1; **1/z RHS
  coefficient-singularity** → seed at z₀=1. `solve_pade` reaches it directly.
  Oracle `mp.besseljzero(0,k)`. Tol value 1e-9, pole 1e-7. *Distinct from CSF.2
  (J₀ itself) and CPB.5 (scaled Bessel).*
- **CRic.3 Hermite (rational)** — `u'=−u²+2zu−2n`, recast
  `u''=2u³−6u²z+4uz²+2u+4nu−4nz`, `u=2n·H_{n−1}/H_n`, **n=4** (odd n puts a zero
  at origin). Machine-precision rational oracle (sympy Hermite roots). Tol value
  **1e-13**, pole 1e-9. *No Hermite-Riccati exists in corpus.*
- **CRic.4 Chebyshev (rational)** — first-order `u'=−u²+(zu−n²)/(1−z²)`, recast
  derived & verified, `u=T_n′/T_n`, **n=4**. Poles `cos((2k−1)π/2n)`; **±1 RHS
  coefficient-singularities** → stay strictly inside (−1,1). Tol value 1e-13,
  pole 1e-9.
- **CRic.5 parabolic cylinder (complex poles)** — `u'=u²−z²`, recast
  `u''=2u³−2uz²−2z`, `u=−U′(0,√2z)/U`. **Complex-conjugate pole pairs** near
  `arg≈±3π/4`, residue −1; oracle harder (PCF root-find, no closed-form zero
  list). `path_network_solve`+`eval_at`. Tol bridge 1e-7, pole 1e-5.
- **CRic.6 (deferred)** — real-row sister `u'=u²+z²` (`w=√z·J_{1/4}(z²/2)`);
  P3, **fold-or-promote decision vs CPB.5** which already pins one such pole.

**Beads 1–2.** File 1 (rational, P1) ships first — it validates the new
row-extraction machinery against machine-precision numbers before the
transcendental cases add their own floors. Oracle dir
`external/probes/corpus-oracles/riccati-logderiv/`.

---

## Family B — Higher-order poles (orders 3–5)

- **CHop.1 order-k fixed pole** — recast `u''=k(k+1)u/(1−z)²`, `u=(1−z)^{−k}`,
  IC `(1,k)`, sweep **k=2,3,4,5**. A **fixed** (non-movable) pole — contrast vs
  the movable Painlevé/Riccati poles. Residue 0 (pure order-k). Past-pole exact
  rationals: `u(1.05)=(−0.05)^{−k}` ∈ {400,−8000,160000,−3200000}. **Must use
  the SVD path** (classical Toeplitz singular, see discovery #4). Multiplicity =
  pre-collapse cluster count (discovery #3): per-k location tol 1e-7 (k=2) →
  5e-3 (k=5), justified by the measured split width. Value tol 1e-13 (F64),
  1e-30 (BF-256).
- **CHop.2 ℘′ triple pole** — first-order d=2 vector `f=[y₂, 6y₁²−g₂/2]`,
  `(y₁,y₂)=(℘,℘′)`, invariant `(℘′)²=4℘³−g₂℘−g₃`. ℘′ has a **triple** pole,
  Laurent `−2/w³`, residue 0. (No 2nd-order *scalar* form — the scalar ℘′ ODE is
  3rd order.) Oracle Mathematica `WeierstrassPPrime[z,{g2,g3}]` (mpmath has no
  ℘′; mpmath ℘-via-sn at dps≥80). `vector_path_network_solve`+
  `extract_poles_shared_q`. Tol past-pole 1e-7.
- **CHop.3 mixed-order shared pole (vector)** — first-order d=2 `f=[2y₂, 3y₁²]`,
  IC `(1,1)`, `y=((1−z)^{−2},(1−z)^{−3})` sharing z=1 with **different orders**
  over one cubic Q. Extends CVE.1 (orders 1,2) to max order 3. Exact rationals;
  tol 1e-11.

**Bead 3** (`corpus_higher_order_pole_test.jl`, P1). Oracle dir
`external/probes/corpus-oracles/higher-order-poles/` (only ℘′ needs an external
pin).

---

## Family C — New elliptic lattices + Jacobi vector triple

- **CEl.1 lemniscatic ℘ (square lattice)** — recast `u''=6u²−1/2` (g₂=1,g₃=0),
  `u=℘(z−1;1,0)`, IC `u(0)=1.050839791…`, `u'(0)=1.894935232…`. **Square**
  lattice (ω'/ω=1, ω=Γ(1/4)²/(4√π)) vs equianharmonic's 60° rhombus — isolates
  lattice geometry at fixed (double) pole order. Discriminating signal: at
  z=1.05 the lattice correction above the bare 400 is 1.25e-4 (lemniscatic) vs
  4.5e-7 (equianharmonic). Oracle Mathematica `WeierstrassP[z,{1,0}]` or mpmath
  sn-recipe (real roots, dps=50-clean). Tol past-pole 1e-9.
- **CEl.2 Jacobi cn, dn** — recasts `u''=(2m−1)u−2m·u³` (cn) and
  `u''=(2−m)u−2u³` (dn), IC `(1,0)`, m=k². Simple poles at `iK'`+lattice,
  residues −i/k (cn), −i (dn). Second invariant `dn²+m·sn²=1` cross-check.
  Oracle `mp.ellipfun('cn'/'dn', z, m)`. `path_network_solve`+`eval_at`. Tol
  complex bridge 1e-5.
- **CEl.3 Jacobi triple (the clean shared-Q oracle)** — first-order **d=3**
  vector `sn'=cn·dn, cn'=−sn·dn, dn'=−m·sn·cn`, IC `(0,1,1)`. Three components,
  **one identical pole lattice** — the textbook external oracle for the shared-Q
  keystone (vs VPN.1.1 which has *distinct* per-component poles). Use **m=1/2**
  (self-dual K=K', square lattice — links to CEl.1). Far-side targets are
  **exact algebraic** at half-periods: past `iK'`, `(sn,cn,dn)=(−i/√k,
  −√((1+k)/k), −√(1+k))`. Oracle `mp.ellipfun(...,m=1/2)`.
  `vector_path_network_solve`+`extract_poles_shared_q`. Tol 1e-5.

**Bead 4** (`corpus_elliptic_lattice_test.jl`, P1). Oracle dir
`external/probes/corpus-oracles/elliptic-lattices/`.

---

## Family D — Periodic vertical pole row (logistic / Bernoulli)

- **CVrow.1 logistic** — `u'=u(1−u)`, recast `u''=u(1−u)(1−2u)`,
  `u=1/(1+c·e^{−z})`. Poles at `z=log(c)+iπ+2πi·k` — an **infinite vertical row
  spaced 2πi**, all simple, residue **+1** (contour-integral confirmed). A pole
  *geometry* absent from the corpus (CPB.6 tanh uses only the *single* nearest
  pole, not the periodicity). Oracle = closed form + mpmath spot values dps=50.
  Tol residual/spot 1e-13, end-to-end approach-pole run 1e-9.

**Bead 5** (`corpus_periodic_pole_test.jl`, P1). Self-contained oracle.

---

## Family E — Elementary branch points (feed the sheet machinery against trivial ground truth)

v1 exercises sheet machinery **only via PVI/PIII transforms or the entire ODE
`u''=u'`** (branches as pure routing constraints). These bind the primitives to
genuinely multivalued closed forms.

- **CBr.1 √ finite cut** — `y'=−z/y`, recast `y''=−C/y³`, `y=√(C−z²)`, **C=2**.
  Two branch points ±√2, a **finite cut**. `path_network_solve(branch_points=
  (−√2,√2))`+`eval_at_sheet` vs `±√(2−z²)`. Tol 1e-9.
- **CBr.2 log single branch** — `u'=1/z`, recast `u''=−1/z²`, `u=log z`. Branch
  at 0. Winding test: CCW unit circle → total 2π, `sheet_index=+1`, `u` jumps
  exactly 2πi/loop (binds CWD.1's angle arithmetic to log's monodromy).
- **CBr.3 atanh two-branch (61um catcher)** — `u'=2/(1−z²)`, recast
  `u''=4z/(1−z²)²`, `u=log((1+z)/(1−z))=2·atanh(z)`. Branch points ±1.
  **Designed to catch `padetaylor-61um`** under a *realistic two-branch walk*:
  a single step subtending `|Δθ|≈1.1π` CCW about z=1 (e.g.
  `z_old=1+0.1e^{−0.45πi}`→`z_new=1+0.1e^{0.65πi}`). `winding_delta` returns the
  wrapped −0.9π (loses a revolution); assert the wrong value (pin current) +
  `@test_broken` the correct +1.1π. **Distinct from CWD.5** (branch at 0,
  endpoints on the real axis): here the branch is at z=1, endpoints off-axis at
  nonzero arg, inside a real ODE walk. *Pin the derivative identity & the ±1
  branch points, NOT the raw log-vs-atanh string (sympy shows a cosmetic
  principal-branch relabel).* Tol angles 1e-12, values 1e-9.

**Bead 6** (`corpus_elementary_branch_test.jl`, P1, references `61um`).
Self-contained sympy/mpmath oracle.

---

## Family F — Algebraic PVI (external multi-sheet oracle)

- **CPvi.1** `u=1/√z` solves PVI `(α,β,γ,δ)=(0,0,1/8,−5/8)` (symbolic residual
  ≡ 0). Two sheets ±. Oracle dps=50: z=4 → ±1/2; z=2+3i →
  `0.4643254526…−0.2484994409…i` (other sheet = negation). Pin
  `pVI_transformed_rhs(0,0,1/8,−5/8)` against the mpmath FFW-formula at the
  algebraic state — **breaks the v1 self-reference** (gap #9).
- **CPvi.2** `u=(√z+1)/2` solves PVI `(1/2,0,1/8,−5/8)` (residual ≡ 0),
  companion.
- **Limitation (discovery #6):** branch at z=0 = PVI's fixed singularity → no
  generic-branch circumambulation. Generic-branch oracle → reference-acquisition
  bead (P3, blocked).

**Bead 7** (`corpus_algebraic_pvi_test.jl`, P2). Oracle
`_oracle_corpus_algebraic_pvi.jl` (mpmath dps=50).

---

## Family G — Out-of-class / fail-loud

- **CFail.1 essential singularity** — `u''=u(1+2z)/z⁴`, `u=e^{1/z}`. **The
  solver silently bridges it** (discovery #1). Pin the error-inflation curve
  (relerr at the z=0.5 step <1e-12; at the z=0.1 step >1e-4) + `@test_broken`
  "should detect out-of-class / refuse" (flips green when a diagnostic lands).
  Oracle mpmath `e^{1/z}` dps=50.
- **CFail.2 Chazy (deferred)** — `u‴=2uu″−3(u′)²`, 3-vector
  `[y₂,y₃,2y₁y₃−3y₂²]`. Movable **natural boundary** (radius ~2.18, a wall of
  essential singularities) → **no far-side oracle** → not a value-pinned corpus
  test. Document as a deferred bead / optional finite-radius diagnostic only.

**Bead 8** (`corpus_out_of_class_test.jl`, P2) + bug bead
`essential-singularity-silent-bridge` (P3).

---

## Family H — BVP exact oracles (gaps #6, #10)

- **CBvx.1 orthogonal-polynomial BVPs (exact polynomial → machine precision)** —
  all 3-arg (u'-dependent) → isolate the D1-coupled Jacobian (gap #6):
  - Hermite `u''=2zu'−2nu` (regular), n=5 on [−1,1]; oracle exact ints
    (`u(−1)=8, u(1)=−8, u(1/4)=881/32`).
  - Legendre `u''=(2zu'−n(n+1)u)/(1−z²)` (singular ±1 → segment inside), n=5 on
    [−4/5,4/5]; exact rationals.
  - Chebyshev `u''=(zu'−n²u)/(1−z²)`, n=5 on [−9/10,9/10]; exact rationals.
  Node-varying `∂f/∂u'` (`2z/(1−z²)`, `z/(1−z²)`) is a stronger Jacobian probe
  than v1's CBV.1–3. Tol **1e-13**.
- **CBvx.2 Mathieu** — `u''=−(a−2q·cos2z)u`, fix `a=a₂(1)=4.371300982735086`
  (`scipy.special.mathieu_a(2,1.0)`), q=1, on [0,π/2], BC = `ce_2`. 2-arg
  periodic var-coeff with an externally-pinned eigenvalue (new class). Tol 1e-10
  (scipy double only).
- **CBvx.3 Kummer/Whittaker** — Kummer `z·u''+(b−z)u'−au=0` with a=−2,b=5/2 →
  `M=4z²/35−4z/5+1` (degree-2 polynomial, exact, 3-arg, segment away from 0,
  [1/2,2]), tol 1e-13; optional Whittaker `whitm(2,1.5,z)` on [1,4] (mpmath
  dps=50), tol 1e-12.

**Beads 9–10** (`corpus_orthopoly_bvp_test.jl` P1;
`corpus_special_fn_bvp_test.jl` P2).

- **CBvx.4 oblique-complex guard (53tu catcher)** — `u''=u`, `u=cosh`, on
  [0,1+i], N=24, but with **interior-Re(t*)** off-segment query (`t*=0.5+3i`,
  i.e. z*=−0.75+2.25i) — the insidious case CBV.7 (perpendicular `t*=10i`)
  misses. Guard never throws; barycentric extrapolates >0.1 off `cosh`. Assert
  on-segment <1e-12, off-segment error >0.1, `@test_broken` should-throw; mirror
  for `vector_bvp_solve`. Oracle mpmath cosh dps=30.

**Bead 11** (`corpus_oblique_guard_test.jl`, P2, references `53tu`). *Gap-#10
residual:* a genuine off-diagonal **coupling** vector BC (e.g.
`B_a=[[1,1],[0,1]]` on `(sin,cos)`) is still untested → follow-on bead CVB.4.

---

## Family I — Path-network test targets

Background ODEs in 2nd-order form (`path_network_solve` needs `y0=(u,u')`); each
has an exact closed form so both **values** and the **chosen route / extracted
poles / sheet index** are pinnable.

- **CPN.1 vertical-wall threading (logistic)** — `f=u(1−u)(1−2u)`, c=2, vertical
  pole wall at Re z=log2. Targets straddling the wall `{2±2i, −1±2i}`; assert no
  visited node within ε of `log2+iπ(2k+1)`. Oracle mpmath dps=50.
- **CPN.2 semi-infinite real row, zero-adjacent-to-pole (Airy-Riccati)** —
  `f=2u·up−1`, `u=−Ai′/Ai`, IC z₀=0.5. The adversarial case for `:min_u`: a
  small-|u| valley (Airy-function zero of `u`) sits next to a pole spike.
  Targets march past `a_k`; the `−2+0.5i` target also exercises **53tu**.
- **CPN.3 lattice-vs-lattice routing (equianharmonic vs lemniscatic ℘)** — same
  target grid, different lattice → spanning trees must differ; extracted poles
  on the rhombic vs square lattice. Oracle odefun ray-integration + lattice
  formulae.
- **CPN.4 near-coalescent pair δ-sweep** — `f=−2u·up+g'(z)`,
  `u=1/(z−a)+1/(z−a−δ)`, a=1, δ∈{0.4,0.1,0.02,0.005}. Couples routing to the
  Froissart filter: resolved as 2 poles for δ≥4·(cluster floor), merged below;
  no spurious doublet between them. Exact rationals.
- **CPN.5 tritronquée empty-sector quantitative pins** — `f=6u²+z`, tritronquée
  IC. Far targets in the pole-free sector `{−30,−50,−30+10i,−40−15i}`: (a) value
  vs `−√(−z/6)` (or tight odefun ray); (b) **`extract_poles` returns ZERO**
  there — the naïve-solver failure mode. Upgrades v1's qualitative-only
  `phase9_tritronquee_test.jl`.
- **CPN.6 long-range multi-direction** — ℘ and logistic to large |z| along
  several rays (0°,±30°,±60°), not just FW Table 5.1's single real axis. Oracle
  = independent mpmath odefun ray-integration. Tol 1e-13 (F64).
- **CPN.7 spiral/winding (61um headline catch, real walk)** — `y=√(1−z²)`,
  recast `f=(−1−yp²)/y`, `branch_points=(1,)`, `cross_branch=true`, targets
  circling z=1. Pin `sheet_index`/`visited_sheet` vs the exact winding; after one
  CCW loop the sqrt sheet flips sign. A single radius-0.45 wedge step subtends
  −2.948 rad (>π) → `winding_delta` loses a revolution → corrupts the returned
  sheet. `@test_broken` the 61um-corrupted sheet on a **fine** ring; pin correct
  agreement on a **coarse** ring (proving it's step-size-triggered). This is the
  first time 61um is exercised under a real walk (CWD.5 is only a unit fixture).
- **CPN.8 off-node dense-eval** — reuse logistic; `eval_at` strictly *between*
  tree nodes (offset ~h/3), pinned to the closed form; assert NaN sentinel
  beyond every disc (default), and `extrapolate=true` recovers. Exercises Stage-2
  fill + per-node Padé store.

**Beads 12–14:** `corpus_pathnet_walls_rows_test.jl` (CPN.1,2,4,8 — P1);
`corpus_pathnet_lattice_sectors_test.jl` (CPN.3,5,6 — P1);
`corpus_pathnet_winding_test.jl` (CPN.7 — P1, references `61um`). Oracle dirs
`external/probes/corpus-oracles/pathnet-{routing,lattice,winding}/`.

---

## Sequencing

1. **Oracle-capture first** per area (residual==0 gate before emitting any
   value), so every test file pins against frozen, independently-derived
   numbers — matching the v1 convention.
2. **P1 rational/exact families first** (CRic.3/4, CHop.1/3, CBvx.1, CVrow.1):
   they validate new machinery (row extraction, mixed-order shared-Q, D1-coupled
   Jacobian) against **machine-precision** numbers before transcendental floors
   enter.
3. **Then P1 transcendental + path-network** (CEl, CBr, CPN.*).
4. **Then P2** (CRic special, CPvi, CFail, CBvx.2/3/4).
5. **Bug + follow-on beads** filed now, worked as the `@test_broken`s demand.

Build **serially** (Rule 7), mutation-prove every load-bearing assertion (Rule
4), full `Pkg.test()` only at end-of-phase (the new `@test_broken`s for 61um /
53tu / essential-singularity auto-flip green when those guards are fixed).

---

## Build status — DELIVERED (2026-06-08)

All 14 new test files + the CVB.4 addition shipped GREEN, mutation-proven (Rule
4), every assertion pinned to independent ground truth (Rule 5), `src/`
untouched, built **serially** by one Opus coding agent each (Rule 7). Per-file
standalone results:

| Bead | File | Result | Notes |
|---|---|---|---|
| bnqq | corpus_riccati_rational | 19/19 | Hermite/Chebyshev exact-rational rows |
| 5rks | corpus_periodic_pole | 12/12 | logistic vertical row; extract_poles resolved it at default min_support=3 |
| 40q1 | corpus_orthopoly_bvp | 12/12 | gap #6 D1-Jacobian **clean**; BC-sign + tol-metric errata |
| gf8i | corpus_higher_order_pole | 43/43 | orders 3–5 incl ℘′ triple (not deferred); SVD-routing premise erratum |
| dctk | corpus_elliptic_lattice | 30/30 | Jacobi-triple shared-Q keystone pinned by independent oracle (closes gap #5) |
| g9lg | corpus_elementary_branch | 29 + 1 broken | CBr.3 = 61um auto-flip marker; CBr.1 cross-sheet deferred |
| vt02 | corpus_pathnet_walls_rows | 36/36 | δ-sweep 2→1 pole merge at cluster floor; no spurious doublet |
| 3jfw | corpus_pathnet_lattice_sectors | 21/21 | tritronquée sector z=-50 + zero-pole; off-axis sector-edge deferred |
| nkw8 | corpus_pathnet_winding | 21 + 1 broken | **first real-walk 61um catch**; min_u-into-branch finding documented |
| tyef | corpus_riccati_special | 19/19 | Airy/Bessel real rows + PCF complex pole (not deferred) |
| by02 | corpus_algebraic_pvi | 29/29 | PVI ζ-RHS wiring verified vs external algebraic solution (closes gap #9) |
| nzcj | corpus_out_of_class | 12 + 1 broken | **v1ub CONFIRMED** (silent essential-singularity bridge); Chazy deferred |
| evlo | corpus_special_fn_bvp | 10/10 | Kummer exact-poly + Mathieu pinned eigenvalue + Whittaker |
| p0yv | corpus_oblique_guard | 47 + 2 broken | 53tu confirmed in **both** BVP + VectorBVP; cosh-sign erratum |
| t6z8 | corpus_vector_bvp (CVB.4) | 112/112 file | off-diagonal coupling τ-assembly verified **correct** (closes gap #10) |

**Library bugs found: 1** — `padetaylor-v1ub` (essential-singularity silent bridge,
CONFIRMED with the real solver; relERR curve annotated on the bead). The 3
prior-known guard bugs (`q0yq`/`53tu`/`61um`) now have **5 in-suite `@test_broken`
auto-flip markers** (CBr.3, CPN.7, CFail.1, and the two CBvx.4 scalar+vector).
**No wrong computed value anywhere** — the numerical core remains clean.

**Brief errata: 5** (`docs/test_corpus/ERRATA.md`, corpus-extension section) — all
brief/recipe slips caught by independent verification, package correct in every
case.

**Clean-verdict checks (no bug, gap closed):** gap #5 (shared-Q keystone, dctk),
gap #6 (D1-coupled BVP Jacobian, 40q1 + evlo), gap #9 (non-self-referential PVI,
by02), gap #10 (coupling vector BC τ-assembly, t6z8).

**Deferred / follow-on (open beads):** `padetaylor-v1ub` (the fix needs a careful
out-of-class diagnostic that must NOT misfire on genuine pole-bridging — its own
design effort + ADR), `padetaylor-90oh` (extract_poles per-cluster multiplicity
API — additive enhancement + ADR), `padetaylor-vimm` (generic-branch algebraic
PVI oracle — blocked on reference acquisition). In-test deferrals: CBr.1
cross-sheet eval, CFail.2 Chazy, CRic.6 real-row sister, the off-axis tritronquée
sector-edge targets.
