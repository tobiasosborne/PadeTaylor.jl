# test/_oracle_corpus_pathnet_lattice_sectors.jl
#
# Pinned mpmath oracle constants for the PATH-NETWORK lattice/sectors corpus
# (bead padetaylor-3jfw, epic padetaylor-25og; plan Family I,
# docs/test_corpus/02_corpus_extension_plan.md:345-385).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/pathnet-lattice/capture.py  (mpmath 1.3,
#   mp.dps = 30).  Every constant below is a verbatim copy of that script's
#   labelled output.  Before any ray value is emitted the script GATES:
#     * CPN.6: the FW ℘ IC satisfies the first integral (u')²=4u³-2 to <1e-13;
#     * CPN.5: each deep-fan odefun ray value agrees with the leading
#       asymptote -√(-z/6) to <1e-4 — the empirical proof the target really
#       sits in the POLE-FREE sector (a target in the pole sector diverges
#       from the asymptote by O(1)).
#   To update: re-run `python3 capture.py` and re-pin.
#
# Ground-truth recipe: plan Family I; FW 2011 §5 Table 5.1 (the ℘ long-range
# problem) + FW 2011 §4.1 eq. 4.1 (the tritronquée IC).  The INDEPENDENT
# oracle is mpmath.odefun RAY-integration of the SAME ODE — a different
# integrator from the package's Taylor–Padé walk, so the agreement is a real
# cross-validation, not a tautology.

# ===========================================================================
# CPN.5  PI tritronquée empty-sector: odefun ray oracle u(target) along
# z(t) = t·dir from the tritronquée IC.  Tight pin (~1e-13 vs package).  The
# leading asymptote -√(-z/6) agrees to ~1e-5 (computed in-test as the loose
# sanity cross-check, per plan md:366).
# ===========================================================================
const CPN5_U_M30    = complex(-2.236091119782981116195, 0.0)
const CPN5_U_M40    = complex(-2.582001916697002318646, 0.0)
const CPN5_U_M50    = complex(-2.886759678692253780191, 0.0)
const CPN5_U_M30P5I = complex(-2.243785562569020419521,  0.1856927067757596218981)
const CPN5_U_M30M5I = complex(-2.243785562569020419521, -0.1856927067757596218981)

# ===========================================================================
# CPN.6  equianharmonic ℘ multi-direction long-range: u(R·e^{iθ}), R = 11,
# odefun ray oracle.  Rays ±30°/±60° stay clear of the ℘ lattice (clearance
# ≥ 0.5).  The 0° (real-axis) value at z = 30 reuses the FW Table 5.1 paper
# number 1.095098255959744 (pinned inline in the test).
# ===========================================================================
const CPN6_U_P30 = complex(-0.9212534224805953301393, -1.062942011066967926564)
const CPN6_U_M30 = complex(-0.9212534224805953301393,  1.062942011066967926564)
const CPN6_U_P60 = complex( 1.134491235223783072278,   0.1687260812812803758035)
const CPN6_U_M60 = complex( 1.134491235223783072278,  -0.1687260812812803758035)
const CPN6_RAY_R = 11.0

# ===========================================================================
# CPN.3  lattice-vs-lattice constants (exact closed-form half-periods).
#   equianharmonic Ω = Γ(1/3)³/(2^(13/6)·π)  → 60° RHOMBIC lattice
#   lemniscatic    ω = Γ(1/4)²/(4√π)          → SQUARE lattice
# (LEMN_OM matches CEl.3's K(1/2)=K'(1/2) self-dual half-period — both are the
# lemniscatic real half-period.)
# ===========================================================================
const CPN3_EQUI_OM = 1.363034090427890310784008
const CPN3_LEMN_OM = 1.85407467730137191843385
