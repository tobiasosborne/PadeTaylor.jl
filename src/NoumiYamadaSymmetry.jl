"""
    PadeTaylor.NoumiYamadaSymmetry

The affine-Weyl-group symmetry layer for the even-parity **Noumi–Yamada**
`A_{2n}^{(1)}` higher-order Painlevé systems, and its companion
rational-solution oracles.  Where `NoumiYamada` builds the ODE itself,
this module supplies the *structure on the solution variety*: the exact
rational solutions and the Bäcklund transformations that permute them.

## The affine-Weyl group `W(A_{2n}^{(1)})`

Noumi & Yamada 1998 show the differential system `A_{2n}^{(1)}` admits a
representation of the extended affine Weyl group
`W̃ = ⟨s_0, …, s_{2n}, π⟩` of type `A_{2n}^{(1)}` as a group of Bäcklund
transformations — automorphisms of the rational-function field
`ℂ(α; f)` that commute with the derivation `d/dt` (NY1998 Theorem
`\\ref{thmA}`, main.tex:174–178).

The generators, transcribed **verbatim** from
`references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex`
lines 137–151 (`\\eqref{WAl1}`, `\\eqref{WAl2}`):

**Reflection `s_i`** (for each `i = 0, …, l`, indices mod `l+1`):

    s_i(α_i)   = −α_i
    s_i(α_j)   = α_j + α_i        (j = i ± 1)
    s_i(α_j)   = α_j              (j ≠ i, i±1)
    s_i(f_i)   = f_i
    s_i(f_j)   = f_j ± α_i/f_i    (j = i ± 1 — `+` for j=i+1, `−` for j=i−1)
    s_i(f_j)   = f_j              (j ≠ i, i±1)

**Diagram rotation `π`** (cyclic shift by 1):

    π(α_j) = α_{j+1},   π(f_j) = f_{j+1}

They satisfy (NY1998 main.tex:158–165, `\\eqref{rel-s}`):

    s_i² = 1
    s_i s_j = s_j s_i            (j ≠ i, i±1 — commutativity)
    s_i s_j s_i = s_j s_i s_j    (j = i±1 — braid relation)
    π^{l+1} = 1
    π s_i = s_{i+1} π

These purely algebraic identities are the **machine-checkable self-test**
of pillar B §4.3: no ODE integration is required — the group relations
are exact rational arithmetic on the `(α; f)` vectors.

## Representation choice

`noumi_yamada_backlund` and `noumi_yamada_rotation` act on a *value
tuple*: `f` and `α` are length-`2n+1` vectors of numbers (typically
`Rational`, so the group relations hold by exact `==`, but any number
type works).  This is the field-element representation of NY1998 §1
specialised to a point of the solution variety.  The reflection's
`α_i/f_i` term forces `f_i ≠ 0`; a zero `f_i` is a genuine breakdown of
the Bäcklund action and is reported loudly (Rule 1) rather than emitting
`±Inf`.

## Rational-solution taxonomy

`noumi_yamada_rational(n, kind)` returns an exact closed-form rational
solution of `A_{2n}^{(1)}` — a function `f(t)` plus its parameter vector
`α` — which doubles as its own oracle (the `PainleveClosedForm` pattern).
Per pillar B §5 (Matsuda 2012 Theorem 1, lines 316–370, for `A_4`; NY1999
lines 580–588 for `A_2`):

  - **`:A`** — `f = (t, 0, …, 0)`,  `α = (1, 0, …, 0)`.  Valid for every
    `n`: the bracket of every component vanishes (all-but-one component
    is zero), so `f_0′ = 1 = α_0` and `f_j′ = 0 = α_j` for `j ≥ 1`.

  - **`:C`** — `f_j = t/(2n+1)` for all `j`,  `α_j = 1/(2n+1)` for all
    `j`.  Valid for every `n`: each bracket has `n` odd-offset and `n`
    even-offset terms all equal to `t/(2n+1)`, so the bracket is zero and
    `f_j′ = 1/(2n+1) = α_j`.  All components nonzero — the **only** type
    safe for a `vector_solve_pade` cross-check (see the caveat below).

  - **`:B`** — `f = (t/3, t/3, t/3, 0, 0)`, `α = (1/3, 1/3, 1/3, 0, 0)`
    — the intermediate `A_4` solution (pillar B §5.2; Matsuda 2012).
    Defined for `A_4` (`n = 2`) **only**: its cancellation structure
    (`f_1 = f_2`, `f_3 = f_4`) relies on the four-term bracket of `A_4`.
    The naive `n`-generalisation fails the ODE for `A_2` (the single
    difference `f_1 − f_2` no longer cancels against a partner); rather
    than ship a wrong formula, `kind = :B` throws for `n ≠ 2`.  A general
    intermediate family would require the full Matsuda classification per
    `n` — deferred (Rule 9): the follow-up condition is "a v0.2+ target
    needs a non-trivial intermediate rational solution for `A_{2n}`,
    `n ≠ 2`".

## Zero-component caveat (bead V5b)

The vector solver's shared-`Q` Padé throws `Q(0) ≈ 0` when a component
is **identically zero**.  Type A `(t, 0, …, 0)` and Type B
`(t/3, t/3, t/3, 0, 0)` have zero components — they **cannot** be
vector-solved directly.  Validate Type A/B by RHS-residual only
(substitute the closed form into `noumi_yamada_rhs`, check `f_j′ = RHS_j`
exactly); use Type C for the `vector_solve_pade` cross-check (all
components `t/(2n+1) ≠ 0`).

## References

  - `references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex`
    lines 136–166 — `\\eqref{WAl1}`, `\\eqref{WAl2}`, the group relations.
  - `docs/v0p2_pillarB_noumi_yamada_findings.md` §4 (Bäcklund generators,
    group relations, §4.3 self-test), §5 (rational-solution taxonomy).
  - `src/NoumiYamada.jl` — the `A_{2n}^{(1)}` RHS / problem builder.
  - `src/PainleveClosedForm.jl` — the v0.1 closed-form-family pattern
    (exact solution doubling as its own oracle) mirrored here.
"""
module NoumiYamadaSymmetry

export noumi_yamada_rational, noumi_yamada_backlund, noumi_yamada_rotation

# -----------------------------------------------------------------------------
# Rational-solution oracles.
# -----------------------------------------------------------------------------

"""
    noumi_yamada_rational(n::Integer, kind::Symbol) -> (α, f)

Exact rational solution of the even-parity Noumi–Yamada system
`A_{2n}^{(1)}` (`d = 2n+1` components, `k = 1` normalisation).

Returns `(α, f)`:

  - `α::Vector{Rational{Int}}` — the `(2n+1)` parameters (`Σα = 1`);
  - `f::Function` — a closure `t -> Vector` giving the exact closed-form
    solution `(f_0(t), …, f_{2n}(t))`.  Each `f_j` is linear in `t`, so
    `f` doubles as its own oracle: substituting `f(t)` into
    `noumi_yamada_rhs` reproduces `f_j′` exactly.

`kind ∈ {:A, :B, :C}` selects the Bäcklund-orbit type (pillar B §5):

  - `:A` — `(t, 0, …, 0)`,           `α = (1, 0, …, 0)`        (all `n`);
  - `:C` — `(t/d, …, t/d)`,          `α = (1/d, …, 1/d)`       (all `n`);
  - `:B` — `(t/3, t/3, t/3, 0, 0)`,  `α = (1/3, 1/3, 1/3, 0, 0)`
           — the intermediate `A_4` solution, `n = 2` **only**.

Throws `ArgumentError` (Rule 1) for `n < 1`, an unknown `kind`, or
`kind = :B` with `n ≠ 2` (see the module docstring for the `:B` scope
rationale and the deferred-work condition).
"""
function noumi_yamada_rational(n::Integer, kind::Symbol)
    n ≥ 1 || throw(ArgumentError(
        "noumi_yamada_rational: n must be ≥ 1 (got $n); A_{2n}^{(1)} is " *
        "defined for n = 1, 2, … (n = 1 ⇒ A_2, n = 2 ⇒ A_4). " *
        "Suggestion: pass n ≥ 1."))
    d = 2n + 1
    if kind === :A
        # (t, 0, …, 0), α = (1, 0, …, 0).  Every bracket vanishes (only one
        # component is nonzero), so f_0′ = 1 = α_0 and f_j′ = 0 = α_j.
        α = [j == 0 ? 1 // 1 : 0 // 1 for j in 0:(d - 1)]
        f = t -> [j == 0 ? t : zero(t) for j in 0:(d - 1)]
        return (α, f)
    elseif kind === :C
        # f_j = t/d for all j, α_j = 1/d.  Each bracket has n odd-offset and
        # n even-offset equal terms ⇒ bracket = 0, f_j′ = 1/d = α_j.
        α = [1 // d for _ in 0:(d - 1)]
        f = t -> [t / d for _ in 0:(d - 1)]
        return (α, f)
    elseif kind === :B
        n == 2 || throw(ArgumentError(
            "noumi_yamada_rational: kind = :B (the intermediate rational " *
            "solution) is defined for A_4 (n = 2) only — its cancellation " *
            "structure relies on the four-term A_4 bracket and the naive " *
            "n-generalisation fails the A_$(2n) ODE for n ≠ 2 (got n = $n). " *
            "Suggestion: use :A or :C for general n, or n = 2 for :B."))
        # A_4 only: (t/3, t/3, t/3, 0, 0), α = (1/3, 1/3, 1/3, 0, 0).
        α = [1 // 3, 1 // 3, 1 // 3, 0 // 1, 0 // 1]
        f = t -> [t / 3, t / 3, t / 3, zero(t), zero(t)]
        return (α, f)
    end
    throw(ArgumentError(
        "noumi_yamada_rational: kind must be :A, :B, or :C (got " *
        "$(repr(kind))).  Suggestion: :A = (t,0,…,0); :C = (t/d,…,t/d); " *
        ":B = the A_4 intermediate (t/3,t/3,t/3,0,0)."))
end

# -----------------------------------------------------------------------------
# Bäcklund transformations — the W(A_{2n}^{(1)}) generators.
# -----------------------------------------------------------------------------

"""
    noumi_yamada_backlund(i::Integer, f, α) -> (f′, α′)

Apply the affine-Weyl-group reflection generator `s_i` of
`W(A_{2n}^{(1)})` to a solution tuple `(f, α)`.  `f` and `α` are
length-`2n+1` vectors of numbers (values at a fixed `t`); the action is
exact arithmetic, so with `Rational` inputs the group relations
(`s_i² = 1`, braid, …) hold by exact equality.

The transform (NY1998 main.tex:137–146, `\\eqref{WAl1}`; indices mod
`2n+1`):

    α_i′   = −α_i
    α_{i±1}′ = α_{i±1} + α_i
    f_i′   = f_i
    f_{i+1}′ = f_{i+1} + α_i/f_i
    f_{i−1}′ = f_{i−1} − α_i/f_i

with all other components unchanged.

Throws `ArgumentError` (Rule 1) for `i ∉ {0, …, 2n}`, `f`/`α` of
mismatched length, or `f_i = 0` (the reflection's `α_i/f_i` term is then
undefined — a genuine breakdown of the Bäcklund action, reported rather
than yielding `±Inf`).
"""
function noumi_yamada_backlund(i::Integer, f, α)
    d = length(f)
    length(α) == d || throw(ArgumentError(
        "noumi_yamada_backlund: f has length $d but α has length " *
        "$(length(α)); both must be the (2n+1)-vector of an A_{2n}^{(1)} " *
        "tuple.  Suggestion: pass f and α of equal length."))
    (d ≥ 3 && isodd(d)) || throw(ArgumentError(
        "noumi_yamada_backlund: an A_{2n}^{(1)} tuple has 2n+1 components " *
        "(odd, ≥ 3), but f has length $d. " *
        "Suggestion: pass a (2n+1)-vector."))
    0 ≤ i ≤ d - 1 || throw(ArgumentError(
        "noumi_yamada_backlund: generator index i = $i is out of range; " *
        "W(A_{2n}^{(1)}) has generators s_0, …, s_$(d - 1). " *
        "Suggestion: pass i ∈ {0, …, $(d - 1)}."))
    slot(j) = mod(j, d) + 1            # 0-based component → 1-based Julia slot
    fi = f[slot(i)]
    iszero(fi) && throw(ArgumentError(
        "noumi_yamada_backlund: reflection s_$i needs f_$i ≠ 0 (the action " *
        "f_{i±1} ↦ f_{i±1} ± α_i/f_i divides by f_i), but f_$i = $fi. " *
        "Suggestion: apply s_$i only where component $i is nonzero — Type A " *
        "/ Type B rational solutions have zero components (pillar B §5)."))
    ip, im = mod(i + 1, d), mod(i - 1, d)
    ratio = α[slot(i)] / fi
    f′ = copy(collect(f))
    α′ = copy(collect(α))
    # α: reflect — α_i ↦ −α_i, α_{i±1} ↦ α_{i±1} + α_i.
    α′[slot(i)]  = -α[slot(i)]
    α′[slot(ip)] = α[slot(ip)] + α[slot(i)]
    α′[slot(im)] = α[slot(im)] + α[slot(i)]
    # f: f_i fixed; f_{i+1} ↦ f_{i+1} + α_i/f_i; f_{i−1} ↦ f_{i−1} − α_i/f_i.
    f′[slot(ip)] = f[slot(ip)] + ratio
    f′[slot(im)] = f[slot(im)] - ratio
    return (f′, α′)
end

"""
    noumi_yamada_rotation(f, α) -> (f′, α′)

Apply the diagram-rotation generator `π` of the extended affine-Weyl
group `W̃(A_{2n}^{(1)})` to a solution tuple `(f, α)`.  `f` and `α` are
length-`2n+1` vectors.

The transform (NY1998 main.tex:149–151, `\\eqref{WAl2}`):

    π(α_j) = α_{j+1},   π(f_j) = f_{j+1}    (indices mod 2n+1)

i.e. a cyclic left shift by one.  `π^{2n+1} = 1` and `π s_i = s_{i+1} π`
(NY1998 main.tex:163–164).

Throws `ArgumentError` (Rule 1) for `f`/`α` of mismatched length, or a
length that is not a valid `2n+1`.
"""
function noumi_yamada_rotation(f, α)
    d = length(f)
    length(α) == d || throw(ArgumentError(
        "noumi_yamada_rotation: f has length $d but α has length " *
        "$(length(α)); both must be the (2n+1)-vector of an A_{2n}^{(1)} " *
        "tuple.  Suggestion: pass f and α of equal length."))
    (d ≥ 3 && isodd(d)) || throw(ArgumentError(
        "noumi_yamada_rotation: an A_{2n}^{(1)} tuple has 2n+1 components " *
        "(odd, ≥ 3), but f has length $d. " *
        "Suggestion: pass a (2n+1)-vector."))
    # π(g_j) = g_{j+1}: the new component j is the old component j+1.
    f′ = [f[mod(j + 1, d) + 1] for j in 0:(d - 1)]
    α′ = [α[mod(j + 1, d) + 1] for j in 0:(d - 1)]
    return (f′, α′)
end

end # module NoumiYamadaSymmetry
