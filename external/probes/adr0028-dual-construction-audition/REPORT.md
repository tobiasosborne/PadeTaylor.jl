# ADR-0028 dual-construction audition — raw run

```
========================================================================================================
ADR-0028 DUAL-CONSTRUCTION AUDITION
  cell A = GGT diagonal (m,m) [shipped] | B_ceil = M–T ⌈m/d⌉ [ADR's cell B] | B_floor = M–T ⌊m/d⌋+shrink [square]
========================================================================================================

system           kind   d   m | err_A     err_Bceil err_Bfloor | R1_A      R1_Bf     | g_A     g_Bf   
--------------------------------------------------------------------------------------------------------
harmonic         entire  2  15 | 5.22e-09  3.33e-16  1.44e-15  | 3.26e-20  3.09e-20  | 3.3e+01 2.3e-04
exp-pair         entire  2  15 | 2.34e-09  8.87e-07  2.22e-16  | 1.66e-18  1.10e-22  | 8.3e+03 9.1e-12
weierstrass-℘    meromorphic  2  15 | 3.93e-04  1.48e-12  1.18e-11  | 7.03e-08  1.23e-13  | 5.1e+01 1.8e+02
calogero-moser   meromorphic  4  15 | 9.66e-11  2.50e-04  3.43e-13  | 2.44e-16  9.36e-25  | 3.1e+04 2.0e-14
noumi-yamada-A₄  meromorphic  5  12 | 4.66e-14  1.03e-10  5.55e-17  | 9.17e-09  1.22e-10  | 4.9e+01 1.1e+01

========================================================================================================
Q1 — B_floor (square M–T) vs A:  err_A/err_Bfloor  (>1 ⇒ B_floor wins)
========================================================================================================
  harmonic         entire       err_A=5.22e-09  err_Bfloor=1.44e-15   B_floor BETTER ×3.6e6
  exp-pair         entire       err_A=2.34e-09  err_Bfloor=2.22e-16   B_floor BETTER ×1.1e7
  weierstrass-℘    meromorphic  err_A=3.93e-04  err_Bfloor=1.18e-11   B_floor BETTER ×3.3e7
  calogero-moser   meromorphic  err_A=9.66e-11  err_Bfloor=3.43e-13   B_floor BETTER ×280.0
  noumi-yamada-A₄  meromorphic  err_A=4.66e-14  err_Bfloor=5.55e-17   B_floor BETTER ×840.0

========================================================================================================
Q2 — METRIC FIDELITY: selector(A,B_floor) vs true-error winner
    single-coeff R  vs  multi-coeff R (span 3);  Rfloor=1e-11 absolute-tie floor
========================================================================================================
  harmonic         true=square_floor  R1→diagonal     ✗   Rmulti→diagonal     ✗
  exp-pair         true=square_floor  R1→diagonal     ✗   Rmulti→diagonal     ✗
  weierstrass-℘    true=square_floor  R1→square_floor ✓   Rmulti→square_floor ✓
  calogero-moser   true=square_floor  R1→diagonal     ✗   Rmulti→diagonal     ✗
  noumi-yamada-A₄  true=square_floor  R1→square_floor ✓   Rmulti→square_floor ✓
  SCORE: single-coeff 2/5   |   multi-coeff 2/5

========================================================================================================
Q3 — CONDITIONING GAP g_A AS A 'BUILD B' GATE
========================================================================================================
  harmonic         entire        g_A=3.32e+01   A-is-worse-than-B_floor? YES (want gate→build B)
  exp-pair         entire        g_A=8.32e+03   A-is-worse-than-B_floor? YES (want gate→build B)
  weierstrass-℘    meromorphic   g_A=5.10e+01   A-is-worse-than-B_floor? YES (want gate→build B)
  calogero-moser   meromorphic   g_A=3.11e+04   A-is-worse-than-B_floor? YES (want gate→build B)
  noumi-yamada-A₄  meromorphic   g_A=4.93e+01   A-is-worse-than-B_floor? YES (want gate→build B)
  g_A where A worse: [1.5, 3.9, 1.7, 4.5, 1.7]
  g_A where A ok:    Any[]
  → if these ranges OVERLAP, g_A cannot gate the build decision.

========================================================================================================
Q5 — DETERMINISM under eps-perturbation (selector A vs B_floor, Rfloor=1e-11, multi)
========================================================================================================
  harmonic         chosen=diagonal     perturbed=diagonal     STABLE
  exp-pair         chosen=diagonal     perturbed=diagonal     STABLE
  weierstrass-℘    chosen=square_floor perturbed=square_floor STABLE
  calogero-moser   chosen=diagonal     perturbed=diagonal     STABLE
  noumi-yamada-A₄  chosen=square_floor perturbed=square_floor STABLE

========================================================================================================
Q7 — d∤m FORK:  B_ceil ⌈m/d⌉ (over-determined)  vs  B_floor ⌊m/d⌋+shrink (square)
========================================================================================================
  harmonic         d=2 m=15  err_Bceil=3.33e-16  err_Bfloor=1.44e-15   ≈ similar
  exp-pair         d=2 m=15  err_Bceil=8.87e-07  err_Bfloor=2.22e-16   DIFFER ×4.0e9 (ceil WORSE)
  weierstrass-℘    d=2 m=15  err_Bceil=1.48e-12  err_Bfloor=1.18e-11   ≈ similar
  calogero-moser   d=4 m=15  err_Bceil=2.50e-04  err_Bfloor=3.43e-13   DIFFER ×7.3e8 (ceil WORSE)
  noumi-yamada-A₄  d=5 m=12  err_Bceil=1.03e-10  err_Bfloor=5.55e-17   DIFFER ×1.9e6 (ceil WORSE)

========================================================================================================
Q6 — A/B_floor DISAGREEMENT (R_gap, pole_disagree) — does it flag ambiguous steps?
========================================================================================================
  harmonic         entire       R_gap=R1_A/R1_Bf=1.05e+00  pole_disagree=2.14e-01  (A ≥ B_floor)
  exp-pair         entire       R_gap=R1_A/R1_Bf=1.51e+04  pole_disagree=3.74e-01  (A ≥ B_floor)
  weierstrass-℘    meromorphic  R_gap=R1_A/R1_Bf=5.70e+05  pole_disagree=7.78e-01  (A ≥ B_floor)
  calogero-moser   meromorphic  R_gap=R1_A/R1_Bf=2.60e+08  pole_disagree=5.42e-01  (A ≥ B_floor)
  noumi-yamada-A₄  meromorphic  R_gap=R1_A/R1_Bf=7.54e+01  pole_disagree=9.29e-01  (A ≥ B_floor)

========================================================================================================
POLE-CROSSING PROBE — tan(z), pole at π/2 ≈ 1.5708;  A=robust_pade(m,m) vs B=robust_pade(m-1,m)
  (d=1 cell-accuracy measurement; in production d=1 short-circuits to A)
========================================================================================================
  h=1.0 (z=1.00, pre  pole) m=10 | err_A(m,m)=6.66e-16  err_B(m-1,m)=2.22e-16 | ≈ tie
  h=1.0 (z=1.00, pre  pole) m=15 | err_A(m,m)=6.66e-16  err_B(m-1,m)=2.22e-16 | ≈ tie
  h=1.3 (z=1.30, pre  pole) m=10 | err_A(m,m)=8.88e-14  err_B(m-1,m)=1.78e-15 | B wins
  h=1.3 (z=1.30, pre  pole) m=15 | err_A(m,m)=8.88e-14  err_B(m-1,m)=1.78e-15 | B wins
  h=1.8 (z=1.80, PAST pole) m=10 | err_A(m,m)=4.09e-13  err_B(m-1,m)=3.84e-13 | ≈ tie
  h=1.8 (z=1.80, PAST pole) m=15 | err_A(m,m)=5.10e-11  err_B(m-1,m)=3.84e-13 | B wins
  h=2.0 (z=2.00, PAST pole) m=10 | err_A(m,m)=1.38e-12  err_B(m-1,m)=1.07e-12 | ≈ tie
  h=2.0 (z=2.00, PAST pole) m=15 | err_A(m,m)=1.08e-10  err_B(m-1,m)=1.07e-12 | B wins
```
