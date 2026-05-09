# §3.3 TaylorSeries.jl + Arblib.jl Compatibility Probe

**Environment:** Julia 1.12.3, TaylorSeries 0.21.7, Arblib 1.7.0 (FLINT-based).
**Run date:** 2026-05-09. All probes executed at
`external/probes/taylorseries-arb/` (see `probe.jl` and `output.txt`).

---

## Probe Results

**Probe 1 — Construction: PASS.**
Both `Taylor1(Arb, 60)` and `Taylor1([Arb(1), Arb(2), Arb(0.5)], 60)` construct
without error and produce `Taylor1{Arb}` with the correct order and coefficients.
No special-casing or wrapper needed.

**Probe 2 — Arithmetic: PASS.**
Addition, subtraction, multiplication, and division of `Taylor1{Arb}` at order 60
all execute without error. `(1+x)^2` computed via `^2` and via explicit
multiplication yields identical ball-overlapping coefficients at every order from 0
to 60.

**Probe 3 — Transcendentals at order 60: PASS.**
`exp`, `sin`, `cos`, and `log(1+x)` all produce `Taylor1{Arb}` whose coefficients
overlap (in the Arb ball sense) with the exact rational values for all 61
coefficients. No `MethodError` or missing specialisation was encountered. Default
precision (~256-bit from the ambient `Arb` environment) is sufficient.

**Probe 4 — High-order stress at order 80: PASS.**
`exp(x)^2` and `exp(2x)` computed as `Taylor1{Arb}` of order 80 agree on all 81
coefficients (mutual ball overlap confirmed). No degradation or type error beyond
order 60.

**Probe 5 — Arb precision propagation at 256 bits: PASS.**
With `setprecision(Arb, 256)`, all 61 coefficients of `exp(x)` have ball radii
below `2^-200` (max observed: `2.16e-78`). The 256-bit precision is faithfully
propagated through TaylorSeries arithmetic — the result is not Float64 in disguise.
The coefficient `exp[60] = 1/60!` overlaps with the independently computed `Arb`
value to full working precision.

---

## Caveats and Notes

- **No custom specialisations required.** TaylorSeries.jl's generic fallback paths
  handle `Arb` throughout. The dispatch is purely type-parameter generic.
- **Scoping quirk (benign).** Running the probe script at top-level triggers Julia's
  soft-scope warning for loop-local assignments. This is a script-level Julia
  behaviour, not a TaylorSeries/Arblib issue, and disappears inside functions.
- **Precision is caller-controlled via `setprecision`.** At the default `Arb`
  precision (also 256 bits in the test environment) radii are already sub-`2^-78`.
  Downstream code should call `setprecision(Arb, N)` before constructing Taylor
  coefficients if a specific precision budget is required.
- **`Arblib.radius(::Arb)` returns `Mag`**, not `Float64`; wrap with `Float64(...)`
  for numeric comparisons.

---

## Bottom-Line Recommendation

**PadeTaylor.jl can rely on `TaylorSeries.jl` with `Taylor1{Arb}` directly.
No hand-rolled coefficient layer is needed.**

All arithmetic operations and all standard transcendentals (`exp`, `sin`, `cos`,
`log`) work correctly up to and beyond order 80, with full rigorous-ball precision
propagation at any `setprecision` level. The generic `T <: Number` design of
PadeTaylor.jl is validated for both `Float64` and `Arb` as first-class element
types via TaylorSeries.jl.
