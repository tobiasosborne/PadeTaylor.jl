# test/_oracle_corpus_algebraic_pvi.jl
#
# Pinned mpmath / sympy oracle constants for the algebraic-PVI corpus
# (bead padetaylor-by02, epic padetaylor-25og).
#
# REGENERATION ORACLE:
#   external/probes/corpus-oracles/algebraic-pvi/capture.py  (sympy 1.14 +
#   mpmath 1.3, dps=50).  Every value below is a verbatim copy of that
#   script's labelled output.  To update: re-run `python3 capture.py`, re-pin.
#
# NON-SELF-REFERENCE (gap #9, docs/test_corpus/02_corpus_extension_plan.md
# Family F, :274-289):  two EXTERNAL algebraic PVI solutions, each gated by a
# sympy SYMBOLIC residual==0 substitution into the standard PVI ODE (NOT the
# package's own output).  The *_WPP constants are the analytic zeta-frame
# second derivative w'' = z*u' + z^2*u'' of those closed forms; the Julia test
# asserts the PACKAGE's pVI_transformed_rhs closure reproduces them.  Two
# independent typings (closed-form mpmath vs the package's FFW eq.(3) closure)
# => a real cross-check, breaking the v1 PVI self-pin.
#
#   CPvi.1  u = z^(-1/2)        PVI (a,b,g,d) = (0,   0, 1/8, -5/8)
#   CPvi.2  u = (sqrt(z)+1)/2   PVI (a,b,g,d) = (1/2, 0, 1/8, -5/8)
#
# Frame:  zeta = log z,  w(zeta) = u(z),  wp = dw/dzeta = z*u'(z),
#         WPP  = d^2 w / d zeta^2 = z*u'(z) + z^2*u''(z).
# Probes: z=4 (real) and z=2+3i (complex), sheets PLUS (+sqrt) / MINUS (-sqrt).

# ---------------------------------------------------------------------------
# CPvi.1  u = z^(-1/2),  (a,b,g,d)=(0,0,1/8,-5/8).  zeta-frame (wp, WPP) pins.
#   u(4)=+1/2 (other sheet -1/2);  u(2+3i)=0.46432545...-0.24849944...i.
# ---------------------------------------------------------------------------
const CPVI1_Z4_PLUS_U    = complex(0.5, 0.0)
const CPVI1_Z4_PLUS_WP   = complex(-0.25, 0.0)
const CPVI1_Z4_PLUS_WPP  = complex(0.125, 0.0)
const CPVI1_Z4_MINUS_U   = complex(-0.5, 0.0)
const CPVI1_Z4_MINUS_WP  = complex(0.25, 0.0)
const CPVI1_Z4_MINUS_WPP = complex(-0.125, 0.0)
const CPVI1_Z2P3I_PLUS_U    = complex(0.464325452650814958080252297505248766652260766, -0.248499440911303374762511568612851429456729049)
const CPVI1_Z2P3I_PLUS_WP   = complex(-0.232162726325407479040126148752624383326130383, 0.124249720455651687381255784306425714728364524)
const CPVI1_Z2P3I_PLUS_WPP  = complex(0.116081363162703739520063074376312191663065191, -0.0621248602278258436906278921532128573641822622)
const CPVI1_Z2P3I_MINUS_U   = complex(-0.464325452650814958080252297505248766652260766, 0.248499440911303374762511568612851429456729049)
const CPVI1_Z2P3I_MINUS_WP  = complex(0.232162726325407479040126148752624383326130383, -0.124249720455651687381255784306425714728364524)
const CPVI1_Z2P3I_MINUS_WPP = complex(-0.116081363162703739520063074376312191663065191, 0.0621248602278258436906278921532128573641822622)

# ---------------------------------------------------------------------------
# CPvi.2  u = (sqrt(z)+1)/2,  (a,b,g,d)=(1/2,0,1/8,-5/8).  zeta-frame pins.
#   u(4)=+3/2 (other sheet -1/2);  the minus sheet is NOT a negation of plus.
# ---------------------------------------------------------------------------
const CPVI2_Z4_PLUS_U    = complex(1.5, 0.0)
const CPVI2_Z4_PLUS_WP   = complex(0.5, 0.0)
const CPVI2_Z4_PLUS_WPP  = complex(0.25, 0.0)
const CPVI2_Z4_MINUS_U   = complex(-0.5, 0.0)
const CPVI2_Z4_MINUS_WP  = complex(-0.5, 0.0)
const CPVI2_Z4_MINUS_WPP = complex(-0.25, 0.0)
const CPVI2_Z2P3I_PLUS_U    = complex(1.33707461401777002022401965042452591083735434, 0.4479887380649190623578668776450217205216621)
const CPVI2_Z2P3I_PLUS_WP   = complex(0.418537307008885010112009825212262955418677169, 0.22399436903245953117893343882251086026083105)
const CPVI2_Z2P3I_PLUS_WPP  = complex(0.209268653504442505056004912606131477709338585, 0.111997184516229765589466719411255430130415525)
const CPVI2_Z2P3I_MINUS_U   = complex(-0.337074614017770020224019650424525910837354339, -0.4479887380649190623578668776450217205216621)
const CPVI2_Z2P3I_MINUS_WP  = complex(-0.418537307008885010112009825212262955418677169, -0.22399436903245953117893343882251086026083105)
const CPVI2_Z2P3I_MINUS_WPP = complex(-0.209268653504442505056004912606131477709338585, -0.111997184516229765589466719411255430130415525)
