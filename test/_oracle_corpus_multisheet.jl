# test/_oracle_corpus_multisheet.jl
#
# Pinned mpmath / closed-form oracle constants for the multi-sheet corpus
# (bead padetaylor-5qqr, epic padetaylor-p1v0).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/multisheet/capture.py  (mpmath 1.3).
#   Every value below is a verbatim copy of that script's labelled output.
#   To update: re-run `python3 capture.py` and re-pin.  Ground-truth recipes
#   are from docs/test_corpus/01_corpus_catalogue.md "Full candidate records":
#   pIII-tilde-sheet-oracle, pV-tilde-sheet-oracle, pIII-tilde-conjugate-
#   symmetry, pVI-eta-subexpr-rhs, pVI-eta-cancellation.
#
# NON-SELF-REFERENCE (gap #9): the four *_RHS_INTERIOR constants are mpmath
# evaluations of the FFW transformed equations (eq. III~, V~, 3, 5) -- the
# Julia test asserts the PACKAGE's RHS closure equals them.  Two independent
# typings of the same formula (mpmath vs Julia) ⇒ a real cross-check, not a
# self-pin.  The *_W_SHEET0 constants are an independent mpmath-RK4
# integration; they reproduce the catalogue's published sheet-0 numbers to
# ~15 digits (capture.py n=4000, dps=40), cross-validating the catalogue.
#
# Closed-form constants (winding deltas, Cramer t/s, residue ±2, the
# constant-solution residual) are written inline in the test, not here.

# ---------------------------------------------------------------------------
# P~III / P~V sheet-0 solution values + conjugate-symmetry residual
# (FFW Fig 1 / Fig 6 params; eval point ζ = -1 + 0.5i).  catalogue:790,807,909
# ---------------------------------------------------------------------------
const PIII_W_SHEET0     = -0.1221859635339962347241257071123534871532 +
                           0.004121037592403508230004318702999635368441im
const PIII_CONJ_SYM_ERR = 0.0
const PV_W_SHEET0       = 0.5415947315160408946982238572572480041575 -
                          1.218296755243758334511186919029554201611im
const PV_CONJ_SYM_ERR   = 0.0

# ---------------------------------------------------------------------------
# Transformed-RHS interior-point values (independent FFW-formula mpmath evals).
# The package RHS closure must reproduce these to Float64 round-off.
#   PIII : a=b=-1/2, g=1, d=-1   @ ζ=0.3+0.7i, w=0.5+0.2i,  w'=0.1-0.05i
#   PV   : a=1,b=-1,g=1,d=-1/2   @ ζ=0.4+0.6i, w=1.5+0.3i,  w'=-0.2+0.1i
#   PVI-η: (4,-4,8,-8)           @ η=log(log 10)+0.5i, v=…, v'=…  (FFW Fig 2)
#   PVI-ζ: (4,-4,8,-8)           @ ζ=0.7+0.3i, w=0.42+0.1i, w'=-0.05-0.02i
# ---------------------------------------------------------------------------
const PIII_RHS_INTERIOR     = -0.5760979817725261688047596701270620025896 -
                               0.8411858602516896528777439853774983097074im
const PV_RHS_INTERIOR       = -2.497301897027132675942574910248933077609 -
                               4.230771820642255466724423795909275061652im
const PVI_ETA_RHS_INTERIOR  = -2.290288413874789514772160274802685159044 +
                               4.497933815598983914480194314636747393389im
const PVI_ZETA_RHS_INTERIOR = -1.271815735114251313891787606724562557571 +
                              13.06289317998171010128068096274005548781im

# ---------------------------------------------------------------------------
# PVI-η double-exp cancellation reference (gap #3): E-1 = e^{e^η}-1 at
# η = log(2π) + iπ/2 + 1e-10.  mpmath dps=60; the BigFloat-256 test recovers
# the REAL part (Float64 gets it wrong by a factor ~2.25e3).  catalogue:841
# Stored as strings so the BigFloat parse keeps all digits.
# ---------------------------------------------------------------------------
const ETA_EM1_REAL_STR = "-1.97392088041526381173526051678469494204e-19"
const ETA_EM1_IMAG_STR =
    "6.283185307493745741881321043451487727874e-10"
# Float64-naive real part at this point (the cancellation casualty), for the
# documented-wrong pin: cmath.exp(cmath.exp(η))-1 → real ≈ +4.44e-16.
const ETA_EM1_REAL_FLOAT64_WRONG = 4.440892098500626e-16
