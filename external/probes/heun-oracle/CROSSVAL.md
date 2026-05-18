# Heun Oracle Cross-Validation Report

**Date:** 2026-05-18  
**Primary oracle:** Mathematica 14.3.0 (`oracles.txt`, 30-digit precision)  
**Independent oracle:** Motygin MATLAB/Octave code (`motygin_oracles.txt`, double precision)  
**Tool:** GNU Octave 8.4.0

---

## Calling conventions

### HeunG (Fuchsian Heun function)

Both sources use the same argument order and normalisation:

| Source      | Signature                                 | Normalisation     |
|-------------|-------------------------------------------|-------------------|
| Mathematica | `HeunG[a, q, alpha, beta, gamma, delta, z]` | `HeunG(0) = 1`  |
| Motygin     | `HeunL(a, q, alpha, beta, gamma, delta, z)` | `HeunG(0) = 1`  |

No sign or argument-order differences. The `epsilon` parameter is internally
derived as `epsilon = alpha + beta + 1 - gamma - delta` by both implementations.

### HeunC (Confluent Heun function)

Both sources use the same argument order and normalisation for `z in (0,1)`:

| Source      | Signature                                        | Normalisation     |
|-------------|--------------------------------------------------|-------------------|
| Mathematica | `HeunC[q, alpha, gamma, delta, epsilon, z]`      | `HeunC(0) = 1`  |
| Motygin     | `HeunC(q, alpha, gamma, delta, epsilon, z)`      | `HeunC(0) = 1`  |

---

## Regime A parameters used

- **HeunG:** `a=2, q=1, alpha=1, beta=2, gamma=3, delta=4`  
- **HeunC:** `q=1, alpha=1, gamma=2, delta=3, epsilon=-1`  
- **z grid:** `{0.1, 0.5, 0.9, 1.5, 1.9, 2.5, 3.0}`

---

## Results by zone

### Zone 1: z in (0, 1) — no branch cut involved

This is the cleanest cross-validation zone. Both implementations evaluate the
power-series solution from z=0 without crossing any branch point.

| Function | z   | Mathematica (30 dig)    | Motygin (double)        | Rel-err   |
|----------|-----|-------------------------|-------------------------|-----------|
| HeunG    | 0.1 | 1.01698183974548594     | 1.016981839745486       | 0         |
| HeunG    | 0.5 | 1.09086834192356525     | 1.090868341923566       | 8.1e-16   |
| HeunG    | 0.9 | 1.08501231904897540     | 1.085012319048977       | 1.6e-15   |
| HeunC    | 0.1 | 0.947195554219883278    | 9.471955542198834e-01   | 1.2e-16   |
| HeunC    | 0.5 | 0.612124644911960719    | 6.121246449119607e-01   | 0         |
| HeunC    | 0.9 | -3.94427764664311075    | -3.944277646643116      | 1.4e-15   |

**Worst rel-err (Zone 1): 1.6e-15** (HeunG at z=0.9).

This is below double-precision machine epsilon (2.2e-16) only in the sense
that accumulated continuation steps add ~10 ULPs. Both oracles agree to
essentially the full precision of double arithmetic.

### Zone 2: z > 1, HeunG — branch points at z=1 and z=a=2 bypassed

For real z > 1, HeunL (Motygin's HeunG) rejects the point as "on branch cut"
and returns NaN. We evaluated with a tiny imaginary offset `z + i*eta`, `eta = 1e-10`,
approaching from the upper half-plane. Mathematica evaluates on the principal
sheet (Im z -> 0^-), so both should agree since they share the same imaginary sign
of the analytic continuation direction here.

| Function | z   | Mathematica (re, im)              | Motygin+eta (re, im)              | Rel-err   |
|----------|-----|-----------------------------------|-----------------------------------|-----------|
| HeunG    | 1.5 | 1.40766339, -0.02052989           | 1.40766339, -0.02052989           | 3.2e-11   |
| HeunG    | 1.9 | 1.62938380, -0.02376361           | 1.62938380, -0.02376361           | 4.4e-11   |
| HeunG    | 2.5 | 2.32593277, +0.27918074           | 2.32593277, +0.27918074           | 7.7e-11   |
| HeunG    | 3.0 | 2.50048593, +1.13019151           | 2.50048593, +1.13019151           | 6.2e-11   |

**Worst rel-err (Zone 2): 7.7e-11.** This is not a genuine disagreement —
it is entirely explained by the `eta = 1e-10` perturbation in the Motygin
argument shifting the result by O(eta * derivative). For reference, the
internal error reported by Motygin's code is ~1e-14, consistent with
double precision; the extra ~1e-10 comes purely from the eta offset.

### Zone 3: z > 1, HeunC — genuine sheet disagreement

For real z >= 1, HeunC has a branch cut on `[1, +inf)`. Crossing this cut
leads to different analytic sheets.

- **Motygin:** analytic continuation from z=0 that circumnavigates z=1 through
  the upper half-plane (path `0 -> 0.38*(1+i) -> 1+i -> z+i*eta`).
- **Mathematica:** evaluates on its principal sheet for real z > 1, which
  corresponds to the limit from the **lower** half-plane (Im z -> 0^-).

These two sheets are **genuinely different**. The discrepancy is not a bug
in either oracle; it is the expected branch-cut ambiguity.

| Function | z   | Mathematica                            | Motygin (upper sheet)              | Rel-err   |
|----------|-----|----------------------------------------|------------------------------------|-----------|
| HeunC    | 1.5 | -1.10692275 - 0.48218202i             | +0.72944634 - 0.40846994i          | 1.52      |
| HeunC    | 1.9 | +0.16927242 - 0.43216653i             | +0.77331418 - 0.39096137i          | 1.30      |
| HeunC    | 2.5 | +0.50518712 - 0.38425236i             | +0.72408167 - 0.35520807i          | 0.35      |
| HeunC    | 3.0 | +0.56022573 - 0.34282037i             | +0.65727866 - 0.31925201i          | 0.15      |

The ratio `Mathematica / Motygin` is not a constant, so this is not a simple
normalisation difference; it is a genuine monodromy-sheet mismatch.

---

## Convention difference discovered: HeunC at z > 1

**Critical finding:** Mathematica's `HeunC[q, alpha, gamma, delta, epsilon, z]`
for real `z > 1` is NOT on the same analytic sheet as Motygin's `HeunC0`
analytic continuation. Specifically:

- Motygin evaluates the principal local solution at z=0, then continues
  analytically to z > 1 by going through the upper half-plane around z=1.
- Mathematica uses its own principal-branch convention for real z > 1.

**Impact on downstream testing:** For PadeTaylor's Taylor-to-Pade scheme, we
compute solution series around z=0 and evaluate for z in [0,1). The regime-A
cases with z in {1.5, 1.9, 2.5, 3.0} for HeunC fall on the branch cut and
should **not** be used as cross-validation anchors without explicitly specifying
which sheet. The regime-A cases with z in {0.1, 0.5, 0.9} are unambiguous and
fully validated.

---

## Verdict

| Scope                     | Worst rel-err  | Trustworthy? |
|---------------------------|----------------|--------------|
| HeunG, z in (0,1)         | 1.6e-15        | YES          |
| HeunG, z in (1, inf)      | 7.7e-11 (eta)  | YES (same sheet, eta artefact) |
| HeunC, z in (0,1)         | 1.4e-15        | YES          |
| HeunC, z in (1, inf)      | 1.5e+00        | NO — different sheets |

**Overall verdict:** The Mathematica oracle (`oracles.txt`) is trustworthy for
all regime-A cases **except HeunC at z >= 1**. The six z-in-(0,1) cases (3
HeunG + 3 HeunC) are confirmed to machine precision by Motygin's independent
code. The four HeunG cases at z > 1 are confirmed to ~1e-10 (limited by the
eta-offset technique, not by oracle error). The four HeunC cases at z >= 1 in
oracles.txt are internally consistent Mathematica values on Mathematica's sheet,
but they cannot be cross-validated against Motygin without an explicit monodromy
connection formula.

**Recommendation:** Design PadeTaylor HeunC tests around z in (0,1). For z > 1,
use HeunG tests (which agree between oracles) or add a HeunC monodromy test once
connection coefficients are implemented.

---

## Files

- `motygin_capture.m` — Octave driver that produced `motygin_oracles.txt`
- `motygin_oracles.txt` — Motygin's double-precision output (regime A)
- `oracles.txt` — primary Mathematica oracle (all regimes, 30 digits)
