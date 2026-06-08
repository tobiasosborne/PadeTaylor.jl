# test/_oracle_corpus_higher_order_pole.jl
#
# Pinned mpmath oracle constants for the HIGHER-ORDER-POLE corpus — the
# ℘ / ℘′ (Weierstrass) transcendental case CHop.2 only (bead
# padetaylor-gf8i, epic padetaylor-25og).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/higher-order-poles/capture.py  (sympy 1.14
#   + mpmath 1.3, mp.dps = 50).  Every constant below is a verbatim copy of
#   that script's labelled output.  Before any value is emitted the script
#   GATES on (a) the elliptic invariant (℘′)² = 4℘³ − g₂℘ − g₃ to ≤1e-40 at
#   four test points, and (b) the triple-pole Laurent leading coefficient
#   z³·℘′(z) → −2 as z→0.  To update: re-run `python3 capture.py` and re-pin.
#
# Ground-truth recipe: docs/test_corpus/02_corpus_extension_plan.md:187-192
# (Family B, CHop.2) + DLMF 23.x (Weierstrass ℘).  LEMNISCATIC lattice
# g₂=1, g₃=0 — the real-root case (e = {1/2, 0, −1/2}) where the validated
# ℘-via-Jacobi-sn recipe ℘ = e₃ + (e₁−e₃)/sn²(z√(e₁−e₃), m) is dps=50-clean
# (avoids the equianharmonic dps≥80 requirement; plan discovery #5).
#
# CHop.2 vector system (first-order, d=2):  f = [y₂, 6 y₁² − g₂/2],
# (y₁, y₂) = (℘, ℘′).  ℘′ carries a TRIPLE pole (Laurent leading −2/w³),
# residue 0.  IC at z₀ = 0.9 (real, away from the lattice), bridge across
# the nearest real pole at z = 2ω ≈ 3.708 to z_tgt = 2ω + 0.5 ≈ 4.208.
#
# Exact rationals for CHop.1 (u(1.05), u(0.5), u′(1.05)) and CHop.3
# (y(2.0), y(2.5)) are written INLINE in the test (sympy gives exact
# integers/rationals) — only the transcendental ℘/℘′ pins live here.

# ---------------------------------------------------------------------------
# Nearest real pole of the lemniscatic ℘ (z = 2ω, ω = K(m)/√(e₁−e₃)).
# ---------------------------------------------------------------------------
const CHOP2_POLE_REAL = 3.708149354602743836867700694390520092435

# ---------------------------------------------------------------------------
# IC point z₀ = 0.9 (real) and the past-pole evaluation target z_tgt = 2ω+0.5.
# ℘ and ℘′ are REAL on the real axis between lattice points.
# ---------------------------------------------------------------------------
const CHOP2_Z0      = 0.9
const CHOP2_P_Z0    = 1.275513014682938691249293372997689270956
const CHOP2_PP_Z0   = -2.650506771638350317732865123061542108155

const CHOP2_ZTGT    = 4.208149354602743836867700694390520092435
const CHOP2_P_ZTGT  = 4.012513027096227403740473529263381151415   # ℘ past pole
const CHOP2_PP_ZTGT = -15.94984362471908463655459027872105851262  # ℘′ past pole
