# Bug sweep D3 — NoumiYamada + NoumiYamadaSymmetry

**Date**: 2026-06-01
**Auditor area**: `A_{2n}^{(1)}` cyclic RHS and affine-Weyl Bäcklund transforms
**Files audited**: `src/NoumiYamada.jl`, `src/NoumiYamadaSymmetry.jl`

## Area

The even-parity Noumi–Yamada `A_{2n}^{(1)}` system builder (`noumi_yamada_rhs`,
`NoumiYamadaProblem`) and its affine-Weyl symmetry layer
(`noumi_yamada_rational`, `noumi_yamada_backlund` = reflection `s_i`,
`noumi_yamada_rotation` = diagram rotation `π`). Special focus per assignment:
cyclic index arithmetic / off-by-one in the modular wrap; the alternating-sign
neighbour sum in the RHS; the affine-Weyl reflection and rotation formulas; the
`π s_i = s_{i+1} π` self-test; the rational-solution oracle.

## References checked

- `references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex:85-88`
  — equation `\eqref{A2n}`, the verbatim even-parity RHS
  `f_j' = f_j ( Σ_{1≤r≤n} f_{j+2r-1} − Σ_{1≤r≤n} f_{j+2r} ) + α_j`, `0 ≤ j ≤ 2n`.
- `…/NoumiYamada1998…/main.tex:107` — "Formulas for the other `f_j'` are obtained
  … simply by the rotation of indices."
- `…/NoumiYamada1998…/main.tex:137-151` — `\eqref{WAl1}` (reflection `s_i`) and
  `\eqref{WAl2}` (rotation `π`); line 147 fixes indices in `ℤ/(l+1)ℤ`.
- `…/NoumiYamada1998…/main.tex:158-165` — group relations `s_i²=1`, braid,
  `π^{l+1}=1`, `π s_i = s_{i+1} π`.
- `…/NoumiYamada1998…/main.tex:1198-1199` — Appendix small-`l` `f_0'`:
  `A_2: f_0(f_1−f_2)+α_0`, `A_4: f_0(f_1−f_2+f_3−f_4)+α_0`.
- `references/tex/noumi_yamada/Matsuda2012_rational_A4_NoumiYamada_JMP53/main.tex:223-228`
  — explicit `A_4^{(1)}` 5-component system.
- `…/Matsuda2012…/main.tex:264-300` — the explicit `A_4` Bäcklund table
  (`s_0..s_4`, `π` action on `f_0..f_4` and `α_0..α_4`) — the independent oracle.
- `…/Matsuda2012…/main.tex:323-338` — rational-solution Theorem 1, Types A/B/C.
- `…/Matsuda2012…/main.tex:304-310` — `f_i ≡ 0 ⇒ α_i = 0`, `s_i = identity`.
- `references/tex/noumi_yamada/NoumiYamada1999_PIV_symmetries_okamoto_NagoyaMathJ153/main.tex:689-694`
  — `π(f_0)=f_1,…,π(f_2)=f_0`, `π(α_0)=α_1,…` (A_2 rotation, independent check).
- `references/tex/noumi_yamada/Clarkson_etal_2020_cyclic_maya_higher_painleve_StudApplMath144/a4sol_revised.tex:266,1489`
  — `f_0' + f_0(f_1−f_2) = α_0` and `f_0' + f_0(f_1−f_2+f_3−f_4) = α_0` (third
  independent confirmation of the A_2 / A_4 RHS).
- `docs/v0p2_pillarB_noumi_yamada_findings.md` §1, §4, §5; `docs/adr/0022-vector-problem-types.md`.

## Findings

No bugs found. Every formula in both files was diffed line-by-line against at
least two independent canonical sources and matched exactly. Details of the
verification are recorded under "Areas verified correct" below. There is no
sign error, no off-by-one in the cyclic wrap, no factorial/coefficient slip, no
branch-cut issue, no aliasing, no RNG/order dependence, and no
transpose-vs-adjoint hazard in either file (there are no matrix operations,
`adjoint`, `transpose`, or `conj` calls anywhere in these two modules — the
entire arithmetic is `+ − * /` on scalars/`Taylor1`).

## Areas verified correct

### RHS cyclic structure and alternating-sign neighbour sum — CORRECT

`src/NoumiYamada.jl:142-147`:
```julia
slot(j) = mod(j, d) + 1
bracket = sum(fvec[slot(j + 2r - 1)] - fvec[slot(j + 2r)] for r in 1:n)
fvec[slot(j)] * bracket + α[slot(j)]
```
Ground truth `…/NoumiYamada1998…/main.tex:85-88`:
`f_j' = f_j ( Σ_{1≤r≤n} f_{j+2r-1} − Σ_{1≤r≤n} f_{j+2r} ) + α_j`.
The per-term form `Σ_r (f_{j+2r-1} − f_{j+2r})` equals
`Σ f_{j+2r-1} − Σ f_{j+2r}` exactly. For `n=2`, `d=5` this expands to
`f_{j+1} − f_{j+2} + f_{j+3} − f_{j+4}`, matching Matsuda2012:223-227 and
Clarkson2020 a4sol:1489. For `n=1` it is `f_{j+1} − f_{j+2}` (= PIV/A_2),
matching main.tex:1198 and a4sol:266.

The 0-based→1-based mapping `slot(j) = mod(j, d) + 1` is the classic spot for an
off-by-one. It is correct: in Julia `mod(x, d)` returns a value in `0..d-1` for
**all** integer `x` when `d > 0` (sign of the divisor), so `slot` lands in
`1..d` even for negative arguments. The code uses `mod`, not `rem` (`rem(-1,5) =
-1` would have broken the wrap). The maximum bracket argument is
`(2n) + 2n = 4n`, still wrapped safely. No index overflow, no off-by-one.

The cyclic-rotation property of main.tex:107 is enforced by construction (the
comprehension `for j in 0:(d-1)` builds every component from the same `slot`
expression), and is independently cross-checked by test NY.1.3
(`test/noumi_yamada_test.jl:106-124`).

### Bäcklund reflection `s_i` — CORRECT (diffed against Matsuda2012 table)

`src/NoumiYamadaSymmetry.jl:228-238`:
```julia
ip, im = mod(i + 1, d), mod(i - 1, d)
ratio = α[slot(i)] / fi
α′[slot(i)]  = -α[slot(i)]
α′[slot(ip)] = α[slot(ip)] + α[slot(i)]   # α_{i+1} += α_i
α′[slot(im)] = α[slot(im)] + α[slot(i)]   # α_{i-1} += α_i
f′[slot(ip)] = f[slot(ip)] + ratio        # f_{i+1} += α_i/f_i
f′[slot(im)] = f[slot(im)] - ratio        # f_{i-1} -= α_i/f_i
```
Ground truth — NY1998 `\eqref{WAl1}` main.tex:137-146:
`s_i(α_i)=−α_i`, `s_i(α_{i±1})=α_{i±1}+α_i`, `s_i(f_i)=f_i`,
`s_i(f_{i±1})=f_{i±1} ± α_i/f_i` (+ for `i+1`, − for `i−1`).
Verified column-by-column against the explicit Matsuda2012 A_4 table
(main.tex:264-300):
- `s_0`: `f_1 += α_0/f_0` (line 272), `f_4 −= α_0/f_0` (281), `f_0` fixed (269),
  `α_0 ↦ −α_0` (285), `α_1 += α_0` (285), `α_4 += α_0` (297) — all match.
- `s_1`: `f_2 += α_1/f_1` (275), `f_0 −= α_1/f_1` (269), `α_0 += α_1` (285),
  `α_2 += α_1` (291) — all match.

The neighbour-sign convention (`+` for the `i+1` neighbour, `−` for `i−1`) is
the load-bearing detail; it is correct and is the only thing test NYS.1.6 (the
ODE-residual / solution-variety check, `test/noumi_yamada_symmetry_test.jl:232-257`)
catches — the mutation note at lines 287-293 confirms a sign swap bites there
and only there. No aliasing: all RHS reads are from the original `f`/`α`, all
writes to fresh copies (`copy(collect(...))`), so the three α-assignments are
order-independent. For the smallest `d=3`, `i+1 ≠ i−1` (mod 3), so the two
neighbour writes never collide.

The `iszero(fi)` guard (line 223) throws rather than emitting `±Inf` — this is
the CLAUDE.md Rule-1 fail-loud choice, intentional and documented; it diverges
from Matsuda's "treat `s_i` as identity when `f_i ≡ 0`" convention
(main.tex:304-310) but in a strictly safer direction (the caller's Type A/B
zero-component tuples are flagged, not silently mishandled). Not a defect.

### Diagram rotation `π` — CORRECT

`src/NoumiYamadaSymmetry.jl:270-271`:
`f′ = [f[mod(j+1,d)+1] for j in 0:(d-1)]` (and same for `α`), i.e. new
component `j` = old component `j+1`. Ground truth NY1998 `\eqref{WAl2}`
main.tex:149-151 `π(f_j)=f_{j+1}, π(α_j)=α_{j+1}`; independently NY1999
main.tex:689-694 (`π(f_0)=f_1,…,π(f_2)=f_0`) and Matsuda2012 π-column
(main.tex:269-298, `π(f_0)=f_1,…,π(f_4)=f_0`). All match.

### `π s_i = s_{i+1} π` self-test — CORRECT (and the contravariance note is sound)

The test `test/noumi_yamada_symmetry_test.jl:201-210` checks the operator
identity as the value-action `vs_i ∘ vπ == vπ ∘ vs_{i+1}`. I hand-traced the
`d=3, i=0` case symbolically against the actual code semantics: both sides
reduce to `f = (f1, f2+α1/f1, f0−α1/f1)`, `α = (−α1, α2+α1, α0+α1)`. They are
equal. The contravariance reasoning in the test comment (lines 191-200) is
correct and the implementation is genuinely consistent (not merely
self-referentially passing). `π^{2n+1}=id` is checked at lines 178-184.

### Rational-solution oracle — CORRECT

`src/NoumiYamadaSymmetry.jl:148-176`. Type A `(t,0,…,0)`/`α=(1,0,…,0)`, Type C
`f_j=t/d`/`α_j=1/d`, Type B (n=2 only) `(t/3,t/3,t/3,0,0)`/`α=(1/3,1/3,1/3,0,0)`
all match Matsuda2012 Theorem 1 (main.tex:323-338) verbatim. The `t/d` and `1/d`
arithmetic uses `/` and `//` (no integer truncation): with `t::Rational` the
oracle stays exact, as test NYS.1.1 relies on (`t = 7//3`, exact `==`). The
`:B`-for-`n≠2` throw is a documented Rule-9 scope guard, not a bug. Substituting
each closed form back into `noumi_yamada_rhs` reproduces `f_j'` exactly (test
NYS.1.1, lines 62-86).

### Parameter / constraint guards — CORRECT, and cannot cause intermittency

`_check_alpha_sum` (NoumiYamada.jl:160-168) and the `Σf0 = tspan[1]` check
(lines 240-248) run only at **construction** of `noumi_yamada_rhs` /
`NoumiYamadaProblem` — never inside the integration loop — so they cannot
introduce a mid-trajectory discontinuity. `real(eltype)`/`float` handle complex
parameters correctly (NY1998 allows complex `α`). The RHS closure allocates a
fresh result vector per call via comprehension: no shared buffer, no in-place
mutation, no RNG, no element-order dependence — none of the recognised
intermittency mechanisms is present in this area.
