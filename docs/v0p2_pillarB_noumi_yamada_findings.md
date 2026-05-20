# Pillar B — Noumi–Yamada A_n^{(1)} Systems: Literature Findings

**Session**: 2026-05-20  
**Sources read**: all LaTeX sources under `references/tex/noumi_yamada/`  
**Law 1 citations**: every claim below cites file + line numbers.

---

## 1. Explicit A_n^{(1)} System

### 1.1 General form

The Noumi–Yamada system is an ODE for `(l+1)` unknowns
`f_0, ..., f_l` and parameters `α_0, ..., α_l` (all indices mod l+1).
The form depends on the parity of `l = n` (the subscript of `A_l^{(1)}`).

**Even case** `l = 2n` (n = 1, 2, ...):

```
A^{(1)}_{2n}:   f_j' = f_j ( Σ_{1≤r≤n} f_{j+2r-1}  -  Σ_{1≤r≤n} f_{j+2r} ) + α_j
```

for `j = 0, 1, ..., 2n` (indices mod 2n+1).

Source: `references/tex/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41/main.tex`, lines 85–88, equation `\eqref{A2n}`.

**Odd case** `l = 2n+1` (n = 1, 2, ...):

```
A^{(1)}_{2n+1}:  f_j' = f_j ( Σ_{1≤r≤s≤n} f_{j+2r-1} f_{j+2s}
                              - Σ_{1≤r≤s≤n} f_{j+2r} f_{j+2s+1} )
                       + (k/2 - Σ_{1≤r≤n} α_{j+2r}) f_j
                       + α_j ( Σ_{1≤r≤n} f_{j+2r} )
```

for `j = 0, 1, ..., 2n+1` (indices mod 2n+2), where `k = Σ α_j`.

Source: NY1998 main.tex lines 92–102, equation `\eqref{A2n+1}`.

### 1.2 Normalization / periodicity constraints

For all `l`, the parameters satisfy:

```
α_0 + α_1 + ... + α_l = k
```

For the applications in this project, `k = 1` (unit normalization).
The dependent variable sum satisfies:

- **Even** `l = 2n`:  `f_0 + f_1 + ... + f_{2n} = t`  (for k=1)
- **Odd** `l = 2n+1`: `f_0 + f_2 + ... + f_{2n} = f_1 + f_3 + ... + f_{2n+1} = t/2`  (for k=1)

Sources:
- Matsuda A_4 paper: `references/tex/noumi_yamada/Matsuda2012_rational_A4_NoumiYamada_JMP53/main.tex`, lines 228–233 (explicit `f_0+f_1+f_2+f_3+f_4=t`, `Σα_j=1`).
- Matsuda A_5 paper: `references/tex/noumi_yamada/Matsuda2007_rational_A5_NoumiYamada_0708.2960/main.tex`, lines 368–377 (explicit `f_0+f_2+f_4 = f_1+f_3+f_5 = t`, `Σα_j=1`).
- NY1998 lines 392–401 (Poisson radical: the invariant g = Σ f_j for even, (g_0, g_1) for odd).

### 1.3 Dependent-variable dimension

For `A_l^{(1)}`, there are `l+1` unknowns `f_0, ..., f_l`.
One constraint `Σ f_j = t` (or two constraints for odd l) reduces the
true degrees of freedom:

| System      | Variables | Constraint(s) | Free d.o.f. | Phase-space dim |
|-------------|-----------|---------------|-------------|-----------------|
| A_2^{(1)}  | 3         | 1             | 2           | 2n=2 (scalar PIV)|
| A_3^{(1)}  | 4         | 2 (g_0,g_1)   | 2           | 2n=2 (scalar PV) |
| A_4^{(1)}  | 5         | 1             | 4           | 2n=4 (4th order) |
| A_5^{(1)}  | 6         | 2             | 4           | 2n=4 (4th order) |
| A_{2n}^{(1)} | 2n+1   | 1             | 2n          | 2n (n pairs)    |
| A_{2n+1}^{(1)}| 2n+2  | 2             | 2n          | 2n (n pairs)    |

Source: NY1998 Theorem 3.1 (Hamiltonian structure): A_{2n} gives n canonical pairs (q,p); A_{2n+1} also gives n canonical pairs. Lines 604–793.

### 1.4 Explicit f_j' formulas for small l (Appendix, NY1998)

From NY1998 Appendix, lines 1196–1204:

```
A^{(1)}_2:   f_0' = f_0(f_1 - f_2) + α_0

A^{(1)}_3:   f_0' = f_0(f_1 f_2 - f_2 f_3) + (k/2 - α_2) f_0 + α_0 f_2

A^{(1)}_4:   f_0' = f_0(f_1 - f_2 + f_3 - f_4) + α_0

A^{(1)}_5:   f_0' = f_0(f_1 f_2 + f_1 f_4 + f_3 f_4 - f_2 f_3 - f_2 f_5 - f_4 f_5)
                   + (k/2 - α_2 - α_4) f_0 + α_0(f_2 + f_4)
```

All other components are obtained by cyclic index rotation `j → j+1`.

---

## 2. A_2^{(1)} ⟹ Scalar PIV Reduction

### 2.1 The A_2^{(1)} system (k=1 normalization)

```
f_0' = f_0(f_1 - f_2) + α_0
f_1' = f_1(f_2 - f_0) + α_1
f_2' = f_2(f_0 - f_1) + α_2
f_0 + f_1 + f_2 = t,   α_0 + α_1 + α_2 = 1
```

Source: Matsuda A_4 paper lines 244–251 (explicit subsystem after setting f_3=f_4=0);
NY1999 Theorem 1.1 lines 521–538 (Theorem `\ref{symmform}` with normalization `f_0+f_1+f_2=3x`, rescaled).

### 2.2 Explicit substitution

NY1999 uses normalization `f_0+f_1+f_2 = 3x` with `k=3`. The v0.1 (FW/standard) normalization uses `f_0+f_1+f_2 = t` with `k=1`. These are related by `t = 3x`.

From NY1999, lines 501–505:
```
f_0 = c(q - p + 2t),   f_1 = -c q,   f_2 = c p
```
with `x = -t/c`, and `c = sqrt(-3/2)`.

The scalar PIV equation is obtained by differentiating the `f_1` equation:

```
f_1'' - (1/2)(f_1')^2/f_1 - (3/2)f_1^3 + 6x f_1^2 + (-(9/2)x^2 + (α_0 - α_2))f_1 + α_1^2/(2 f_1) = 0
```

Source: NY1999 main.tex lines 568–575.

This is the standard PIV:  `y'' = (1/2y)(y')^2 + (3/2)y^3 + 4ty^2 + 2(t^2-a)y + β/y`
under `y = -f_1/c`, `t → -x/c`, with
```
a = 1 + 3v_3,    β = -2(v_1 - v_2)^2
```
and parameter relations `α_0 = 3(1-v_1+v_3)`, `α_1 = 3(v_1-v_2)`, `α_2 = 3(v_2-v_3)`.

Source: NY1999 lines 450–457, 515–519.

### 2.3 Backwards-compatibility test anchor

For the v0.2 bit-identical backward-compat self-test:
- Set `n=2`, initialize `(f_0, f_1, f_2)` with `f_3 = f_4 = 0` (the Matsuda A_4 paper explicitly demonstrates this reduction at lines 240–252).
- The v0.1 scalar PIV integration of `-c * f_1(t/c)` must agree with v0.2's `NoumiYamadaProblem(n=2)`.
- Normalization: for direct comparison, use `f_0+f_1+f_2 = t`, `Σα_j = 1` (not the NY1999 k=3 convention).

---

## 3. A_3^{(1)} ⟹ PV and Explicit Higher Systems

### 3.1 A_3^{(1)} ⟹ PV

NY1998 lines 108–109: "the differential systems for l=2 and l=3 are equivalent to the fourth and the fifth Painlevé equations, respectively (see [NY5])."

The A_3^{(1)} system (odd case, n=1, l=3) is:
```
f_0' = f_0(f_1 f_2 - f_2 f_3) + (k/2 - α_2) f_0 + α_0 f_2
```
(plus cyclic rotations for f_1', f_2', f_3'), with `f_0 + f_2 = f_1 + f_3 = t`.

Source: NY1998 main.tex lines 1200–1201, Appendix.

The PV system has 4 components but 2 constraints, giving a 2-dimensional phase space = scalar PV.
The Hamiltonian for A_3^{(1)} from NY1998 Appendix lines 1259–1266 is:
```
H(A_3^{(1)}) = (x_0^2 - q_1) p_1 q_1 (x_1/x_0 - p_1)
              + (α_1 + α_3) q_1 p_1 - α_1 (x_1/x_0) q_1
              + α_2 x_0^2 p_1
              + (1/4)(α_1 - 2α_2 - α_3) x_0 x_1
              + (1/4)(α_1 + α_3)^2
```
with auxiliary variables `x_0 = f_0+f_2`, `x_1 = f_1+f_3`.

### 3.2 A_4^{(1)} explicit system (v0.2 primary target)

From Matsuda 2012, lines 220–233:
```
A^{(1)}_4 (α_j)_{0≤j≤4}:
  f_0' = f_0(f_1 - f_2 + f_3 - f_4) + α_0
  f_1' = f_1(f_2 - f_3 + f_4 - f_0) + α_1
  f_2' = f_2(f_3 - f_4 + f_0 - f_1) + α_2
  f_3' = f_3(f_4 - f_0 + f_1 - f_2) + α_3
  f_4' = f_4(f_0 - f_1 + f_2 - f_3) + α_4
  f_0 + f_1 + f_2 + f_3 + f_4 = t,   Σα_j = 1
```
Dimension: 5 variables, 1 constraint, 4-dimensional phase space (2 canonical pairs).

Also given in NY1998 Appendix lines 1199, 1213, 1255–1258 (Hamiltonian):
```
H(A_4^{(1)}) = (t - q_1 - q_2)(q_1 p_1 + q_2 p_2) - q_1 p_1^2 - q_2 p_2^2
             - q_1(p_1 - p_2)q_2
             - α_1 q_1 - (α_1 + α_3) q_2 + α_2 p_1 + α_4 p_2
             + (1/5)(2α_1 - α_2 + α_3 - 2α_4) t
```
with `dx/dt = 1` (i.e., `x = t + const`).

### 3.3 A_5^{(1)} explicit system

From Matsuda 2007, lines 253–379 (the full system `(*)`):
```
A^{(1)}_5 (α_j)_{0≤j≤5}:
  (t/2) f_0' = f_0(f_1 f_2 + f_1 f_4 + f_3 f_4 - f_2 f_3 - f_2 f_5 - f_4 f_5)
              + (1/2 - α_2 - α_4) f_0 + α_0(f_2 + f_4)
  (t/2) f_1' = f_1(f_2 f_3 + f_2 f_5 + f_4 f_5 - f_3 f_4 - f_3 f_0 - f_5 f_0)
              + (1/2 - α_3 - α_5) f_1 + α_1(f_3 + f_5)
  (t/2) f_2' = f_2(f_3 f_4 + f_3 f_0 + f_5 f_0 - f_4 f_5 - f_4 f_1 - f_0 f_1)
              + (1/2 - α_4 - α_0) f_2 + α_2(f_4 + f_0)
  [+ cyclic for f_3', f_4', f_5']
  f_0 + f_2 + f_4 = f_1 + f_3 + f_5 = t,   Σα_j = 1
```
Dimension: 6 variables, 2 constraints, 4-dimensional phase space.

Note the A_5 system uses independent variable t with a (t/2) prefactor on the left,
corresponding to the logarithmic independent variable in standard PV form.

NY1998 Appendix gives `f_0'` (without the t/2 prefactor) at line 1202–1203; the Matsuda 2007
form with the t/2 is the correct formulation matching PV's singularity at t=0.

---

## 4. Affine-Weyl Symmetry W(A_n^{(1)})

### 4.1 Bäcklund generators s_0, ..., s_l and π

From NY1998 main.tex lines 136–152, equations `\eqref{WAl1}` and `\eqref{WAl2}`:

**Reflections** (for each i = 0, ..., l, indices mod l+1):
```
s_i(α_i) = -α_i
s_i(α_j) = α_j + α_i       for j = i ± 1
s_i(α_j) = α_j             for j ≠ i, i±1

s_i(f_i) = f_i
s_i(f_j) = f_j + α_i/f_i   for j = i+1  (the +1 neighbor)
s_i(f_j) = f_j - α_i/f_i   for j = i-1  (the -1 neighbor)
s_i(f_j) = f_j             for j ≠ i, i±1
```

**Diagram rotation** (cyclic shift by 1):
```
π(α_j) = α_{j+1},   π(f_j) = f_{j+1}
```

Source: NY1998 lines 136–152. The action on α is the standard reflection in the generalized
Cartan matrix A_{ij} = 2δ_{ij} - δ_{j,i±1} of type A_l^{(1)} (NY1998 eq. `\eqref{GCM}`,
lines 188–191).

The explicit table for A_4^{(1)} is at Matsuda 2012 lines 264–300 and for A_5^{(1)} at
Matsuda 2007 lines 386–430.

### 4.2 Group relations

From NY1998 lines 158–165:
```
s_i^2 = 1
s_i s_j = s_j s_i        for j ≠ i, i±1  (commutativity)
s_i s_j s_i = s_j s_i s_j  for j = i±1  (braid relation)
π^{l+1} = 1
π s_i = s_{i+1} π
```

### 4.3 Concrete identity for the W(A_n^{(1)}) self-test

The cleanest machine-checkable identity (no integration required) is:

**s_i^2 = identity**: for each generator i, applying s_i twice must return the original (f, α).

For A_4^{(1)}, the concrete check (using Matsuda 2012 lines 262–300):
```
s_0: f_4 ↦ f_4 - α_0/f_0,  then applying s_0 again: (f_4 - α_0/f_0) + α_0/f_0 = f_4
     f_1 ↦ f_1 + α_0/f_0,  then: (f_1 + α_0/f_0) - α_0/f_0 = f_1
     α_0 ↦ -α_0,           then: -(-α_0) = α_0  ✓
```

**Commutativity**: `s_0 ∘ s_2 = s_2 ∘ s_0` (indices differ by 2, hence commute in A_4).

**Braid relation**: `s_0 s_1 s_0 = s_1 s_0 s_1` (adjacent indices).

**π^5 = identity** for A_4^{(1)}.

The recommended PRD self-test is:
> For a generic solution (f_0,...,f_4,α_0,...,α_4), verify that s_i composed twice
> returns the original tuple (both f and α components), and that π^5 returns the original.
> These are purely algebraic: no ODE integration needed.

Source: NY1998 Theorem 1.1 (lines 152–166) + Matsuda 2012 table (lines 264–300).

---

## 5. Rational Solutions and Closed-Form Oracles

### 5.1 A_2^{(1)} (= PIV): seed rational solutions

From NY1999 lines 580–588:
```
Type A: (α_0,α_1,α_2; f_0,f_1,f_2) = (1,1,1; t,t,t)      [note: normalization f_0+f_1+f_2=3t,
                                                               so with k=1 normalization: (1/3,1/3,1/3; t/3,t/3,t/3)]
Type B: (α_0,α_1,α_2; f_0,f_1,f_2) = (3,0,0; 3t,0,0)      [k=3 normalization;
                                                               k=1: (1,0,0; t,0,0)]
```

These are the two seed solutions; all rational solutions of PIV are obtained from them via
Bäcklund transformations (NY1999 lines 589–591, citing Murata).

### 5.2 A_4^{(1)}: complete classification of rational solutions

From Matsuda 2012 Theorem 1 (lines 316–370), there are exactly 3 B̈acklund-orbit types:

```
Type A: (f_0,f_1,f_2,f_3,f_4) = (t, 0, 0, 0, 0),
        (α_0,α_1,α_2,α_3,α_4) = (1, 0, 0, 0, 0)

Type B: (f_0,f_1,f_2,f_3,f_4) = (t/3, t/3, t/3, 0, 0),
        (α_0,α_1,α_2,α_3,α_4) = (1/3, 1/3, 1/3, 0, 0)

Type C: (f_0,f_1,f_2,f_3,f_4) = (t/5, t/5, t/5, t/5, t/5),
        (α_0,α_1,α_2,α_3,α_4) = (1/5, 1/5, 1/5, 1/5, 1/5)
```

These are trivial (constant-proportionality) rational solutions, valid for all t.
Any rational solution is equivalent (via Bäcklund transformations) to a shifted/permuted
version of these three types. The conditions on parameters for which rational solutions
exist (in general, not just at the seed) are:
- Type A: all α_j ∈ Z
- Type B: some cyclic permutation makes (α_i,...,α_{i+4}) ≡ ±(1/3)(1,1,1,0,0) mod Z
- Type C: some cyclic permutation makes parameters ≡ j/5(1,1,1,1,1) mod Z or j/5(1,2,1,3,3) mod Z

Source: Matsuda 2012 lines 316–370 (Theorem 1.1 with full conditions).

**Wronskian-determinant form**: Clarkson et al. (2020) give explicit Wronskian expressions for
Type A solutions (signature class (5)) as ratios of Wronskians of Hermite polynomials.
Example from Clarkson2020 lines 1589–1611:
```
f_0(z) = d/dz [log H_{M_2}(c_1 z) - log H_{M_0}(c_1 z)],  α_0 = 1,
f_1(z) = z + d/dz [log H_{M_3}(c_1 z) - log H_{M_1}(c_1 z)],  α_1 = -2,
...
```
where H_{M_i} are Wronskian determinants of Hermite polynomials, and
c_1^2 = -1/2.

### 5.3 A_5^{(1)}: rational solutions (Matsuda 2007)

From Matsuda 2007 Theorem (lines 443–486), the rational solutions up to Bäcklund transformations are:

```
(a-1): (f_0,...,f_5) = (t,t,0,0,0,0),   (α_j) = (α_0, 1-α_0, 0, 0, 0, 0)

(a-2): (f_0,...,f_5) = (t,0,0,t,0,0),   (α_j) = (α_0, 0, 0, 1-α_0, 0, 0)

(a-3): (f_0,...,f_5) = (t,t,0,0,0,0) or (0,t,t,0,0,0) or (0,t,0,0,t,0) or (t,t,t,0,-t,0),
       (α_j) = (0, 1, 0, 0, 0, 0)  [one free parameter α_0 absorbed]

(b): (f_0,...,f_5) = (t/2, t/2, t/2, t/2, 0, 0),
     (α_j) = (α_0, -α_0+1/2, α_0, -α_0+1/2, 0, 0)

(c): (f_0,...,f_5) = (t/3, t/3, t/3, t/3, t/3, t/3),
     (α_j) = (α_4, -α_4+1/3, α_4, -α_4+1/3, α_4, -α_4+1/3)
```

Source: Matsuda 2007 lines 443–486.

### 5.4 Best oracle for unit-test: A_4^{(1)} Type C seed solution

**Recommended oracle** (simplest exact rational, no Hermite Wronskians needed):

For A_4^{(1)} with parameters
```
(α_0, α_1, α_2, α_3, α_4) = (1/5, 1/5, 1/5, 1/5, 1/5)
```
the exact rational solution is:
```
f_j(t) = t/5   for all j = 0, 1, 2, 3, 4
```

**Verification**: Substitute into the ODE. For `f_j = t/5`:
- LHS: `f_j' = 1/5`
- RHS: `f_j(f_1 - f_2 + f_3 - f_4) + α_j = (t/5)(t/5 - t/5 + t/5 - t/5) + 1/5 = 0 + 1/5 = 1/5` ✓

Source: Matsuda 2012 lines 322–338 (Type C); Clarkson 2020 lines 1499–1502.

Similarly for A_5^{(1)}, the Type C oracle uses `f_j = t/3` for all j with 1-parameter family
`(α_0,...,α_5) = (α_4, 1/3-α_4, α_4, 1/3-α_4, α_4, 1/3-α_4)`.

---

## 6. Pole Structure for n ≥ 4

### 6.1 Known results for PIV (A_2^{(1)}, scalar reference)

Masoero & Roffelsen (2018) study the distribution of poles of rational solutions of PIV
expressed via roots of generalized Hermite polynomials H_{m,n} and generalized Okamoto
polynomials Q_{m,n}.

Key results (MasoeroRoffelsen sigma18-002.tex, lines 89, 222–238):
- All poles and zeros of rational PIV solutions are roots of H_{m,n} or Q_{m,n}.
- The roots of H_{m,n} have an approximately rectangular lattice pattern in C.
- Asymptotically (m→∞, n fixed): roots scale as √(2m+n) and are distributed in a
  region with elliptic boundary; the bulk approximation is accurate to O(m^{-1/3}).

Source: MasoeroRoffelsen lines 286–346.

### 6.2 A_4^{(1)} pole patterns (Clarkson et al. 2020)

The zeros of the special polynomials (Wronskian determinants of Hermite polynomials)
that appear in rational solutions of A_4^{(1)} show similar regular patterns.

Clarkson 2020 lines 1999–2046:
> "The zeros of Okamoto and generalized Hermite polynomials that appear in the rational
> solutions to PIV are known to form very regular patterns in the complex plane.
> In this section we show the equivalent patterns for their A_4 counterparts,
> which are also very regular but show a richer structure."

The paper includes pole-field figures for several classes of A_4 solutions (Figures in
the `figures/` subdirectory: `Ok3_8-16-*.pdf`, `Ok5_*.pdf`, `2H_*.pdf`).

Relevant: for A_4^{(1)} the pole distribution depends on the signature class:
- Signature (5): genus-2 generalization of generalized Hermite; rectangular-ish pattern.
- Signature (1,1,1,1,1): all simple interlacing; most complex pattern.

Source: Clarkson 2020 lines 1999–2052.

### 6.3 Practical recommendation for pole-field figure (n ≥ 4)

For n=4 (A_4^{(1)}), the simplest non-trivial pole-field to plot is obtained by applying
Bäcklund transformations to the Type A seed `(f_0,f_1,f_2,f_3,f_4) = (t,0,0,0,0)`:

1. **Window**: Since rational solutions of A_4 have poles that scale like Hermite polynomial
   roots (roughly √(2m)), for modest integer parameters m the poles spread over a region
   of radius ~√(2|parameter|). For m~10, the window |z| ≤ 10 captures most poles.

2. **For generic (non-rational) transcendent solutions**: no explicit pole-field theory for
   A_4 is available analogous to FW's PIV figure. The FW-style figure would use a large
   domain in t, e.g., Re(t) ∈ [0,10], Im(t) ∈ [-5,5], with generic initial conditions
   at t=1, e.g., `(f_0,...,f_4) = (0.2, 0.2, 0.2, 0.2, 0.2) + small perturbation` to
   break the Type C symmetry, and let poles accumulate naturally.

3. **No "tronquée" analog known for A_4**: unlike PIV where Boutroux asymptotics give
   tritronquée solutions with pole-free sectors, no analogous rigorous theory exists for
   A_4^{(1)} at present. The pole structure is expected (but not proven) to show lattice-type
   patterns from the Wronskian/Maya diagram theory.

---

## 7. Concrete Recommendations for v0.2

### 7.1 Recommended target for pole-field figure (n ≥ 4)

**Use A_4^{(1)}** (the simplest case with n ≥ 4, i.e., l = 4):

```julia
# v0.2 pole-field figure: A_4^{(1)} with all-equal parameters, pertturbed IC
α = (1/5, 1/5, 1/5, 1/5, 1/5)   # parameters summing to 1
f0 = [0.3, 0.25, 0.15, 0.2, 0.1]  # IC at t0=1; NOT the symmetric equilibrium t/5
# (perturb from (t/5)^5 so transcendent behavior emerges)
t_span = (1.0, 20.0)
# Plot poles of f_0 (or Σ f_j - t = residual) in complex t-plane
```

**Why A_4^{(1)}**: 
- First NY system outside the classical Painlevé hierarchy (l=4 > l=3).
- Fully characterized by Matsuda 2012 (rational case) and Clarkson 2020 (explicit Wronskians).
- Has known rational solutions for oracle testing (Type C: f_j = t/5).
- Phase-space dimension 4 is already non-trivial but manageable.

Source: Matsuda 2012 abstract (lines 95–99); Clarkson 2020 lines 317–319.

### 7.2 Recommended oracle for `noumi_yamada_rational` unit test

**Best oracle: A_4^{(1)} Type C seed solution**

```julia
# Exact rational solution: f_j(t) = t/5 for all j
# Parameters: α = (1/5, 1/5, 1/5, 1/5, 1/5)
# Satisfies: f_j' = 1/5 = f_j * 0 + α_j ✓  (all pairwise differences cancel)

function noumi_yamada_type_c_a4(t)
    return fill(t/5, 5)  # exact solution
end
```

The next level of oracle (non-trivial rational, uses Bäcklund shift by 1):

Apply s_0 to the Type A seed (f_0 = t, f_1=...=f_4=0, α_0=1,rest=0):
```
s_0(f_1) = f_1 + α_0/f_0 = 0 + 1/t = 1/t
s_0(f_4) = f_4 - α_0/f_0 = 0 - 1/t = -1/t
s_0(f_0) = f_0 = t  (fixed)
s_0(α_0) = -1,  s_0(α_1) = 1,  s_0(α_4) = 1
```
This gives the solution `(f_0,f_1,f_2,f_3,f_4) = (t, 1/t, 0, 0, -1/t)` with parameters
`(α_0,...,α_4) = (-1, 1, 0, 0, 1)` — a rational function of t with simple poles at t=0.

Source: Matsuda 2012 table at lines 264–300 (explicit s_0 action).

### 7.3 Recommended W(A_n^{(1)}) self-test phrasing

The self-test should verify:

```julia
# For any (f, α) = (f_0,...,f_n, α_0,...,α_n) satisfying the ODE:

# Test 1: s_i^2 = identity (purely algebraic)
for i in 0:n
    fi_new, alpha_new = apply_s(i, f, α)
    fi_back, alpha_back = apply_s(i, fi_new, alpha_new)
    @test fi_back ≈ f && alpha_back ≈ α
end

# Test 2: π^{n+1} = identity
fp, αp = f, α
for _ in 1:(n+1)
    fp, αp = apply_pi(fp, αp)
end
@test fp ≈ f && αp ≈ α

# Test 3: s_i commutes with ODE flow (harder; requires integration)
# For the Type C rational solution f_j = t/5:
#   s_0(f_j = t/5) should satisfy the ODE with shifted parameters
```

For Test 1 (purely algebraic), use any rational solution as input, e.g., Type C.
For Test 3, the Bäcklund transformation of f_j=t/5 under s_0 gives a known rational function
(the 1-shifted solution from Matsuda 2012), which can be verified against direct integration.

Source: NY1998 Theorem 1.2 (lines 173–178: "action commutes with derivation").

---

## Summary of Key File Paths and Line Numbers

| Claim | File | Lines |
|-------|------|-------|
| General A_{2n} system | NY1998 main.tex | 85–88 |
| General A_{2n+1} system | NY1998 main.tex | 92–102 |
| A_2 Appendix formula | NY1998 main.tex | 1198 |
| A_3 Appendix formula | NY1998 main.tex | 1200–1201 |
| A_4 Appendix formula | NY1998 main.tex | 1199 |
| A_5 Appendix formula | NY1998 main.tex | 1202–1203 |
| A_4 full system | Matsuda2012 main.tex | 220–233 |
| A_5 full system | Matsuda2007 main.tex | 253–379 |
| PIV reduction from A_2 | NY1999 main.tex | 501–576 |
| A_2 → PIV substitution | NY1999 main.tex | 501–505 |
| A_2 scalar PIV equation | NY1999 main.tex | 568–575 |
| Bäcklund generators s_i | NY1998 main.tex | 136–152 |
| Group relations W(A_l^{(1)}) | NY1998 main.tex | 158–165 |
| A_4 BT table | Matsuda2012 main.tex | 264–300 |
| A_5 BT table | Matsuda2007 main.tex | 386–430 |
| A_4 rational classification | Matsuda2012 main.tex | 316–370 |
| A_4 Type C seed | Matsuda2012 main.tex | 332–338 |
| A_5 rational classification | Matsuda2007 main.tex | 443–486 |
| A_4 seed solutions | Clarkson2020 a4sol_revised.tex | 1497–1502 |
| A_4 Wronskian construction | Clarkson2020 a4sol_revised.tex | 1586–1612 |
| A_4 pole patterns | Clarkson2020 a4sol_revised.tex | 1999–2046 |
| PIV pole asymptotics | MasoeroRoffelsen sigma18-002.tex | 89, 286–346 |
| A_4 Hamiltonian | NY1998 main.tex | 1255–1258 |
| A_5 Hamiltonian | NY1998 main.tex | 1267–1281 |
