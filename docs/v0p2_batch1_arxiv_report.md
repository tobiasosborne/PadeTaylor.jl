# v0.2 Literature Acquisition — Batch 1 arXiv Fetch Report

**Date:** 2026-05-19  
**Agent:** Batch 1 subagent (Claude Sonnet 4.6)  
**Scope doc:** `docs/v0p2_literature_scope.md`

## Summary

- **Fetched:** 32 / 31 attempted  
  _(Note: the scope listed 32 target IDs after excluding the pre-fetched `1502.06695`; all 32 were fetched successfully.)_
- **Failed:** 0
- **Total bytes downloaded:** 29,199,780 bytes (27.85 MB)
- **Wall-clock time:** ~229 seconds (18:10:20 – 18:14:09 CEST)
- **Pre-2007 IDs using slash-form URLs:** `math/9808003`, `math/0106208`, `math/0212117`, `math-ph/9810007`, `q-alg/9708018`, `nlin/0306020` — all fetched successfully with `https://arxiv.org/pdf/<ID>.pdf` (slash preserved).

## Per-paper results

All files verified: size > 50 KB, magic bytes = `%PDF-`, page count > 5 (via `pdfinfo`).

| Status | arXiv ID | Size (KB) | Pages | Path |
|--------|----------|-----------|-------|------|
| [OK] | 1306.6161 | 911 | 36 | `references/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41.pdf` |
| [OK] | 1101.2602 | 334 | 27 | `references/painleve_hierarchy/ClaeysGrava2011_KdV_painleve_CMP313.pdf` |
| [OK] | math/9808003 | 183 | 16 | `references/noumi_yamada/NoumiYamada1998_higher_painleve_A1l_FunkEkv41.pdf` |
| [OK] | 1612.00337 | 1759 | 29 | `references/hermite_pade/NakatsukasaSeteTrefethen2018_AAA_SISC40.pdf` |
| [OK] | 0704.2869 | 243 | 31 | `references/garnier/Sasano2007_garnier_two_variables_0704.2869.pdf` |
| [OK] | math/0106208 | 301 | 32 | `references/garnier/Mazzocco2002_geometry_garnier_classical_IMRN_math0106208.pdf` |
| [OK] | 0708.0074 | 298 | 46 | `references/noumi_yamada/Matsuda2012_rational_A4_NoumiYamada_JMP53.pdf` |
| [OK] | 0708.2960 | 486 | 75 | `references/noumi_yamada/Matsuda2007_rational_A5_NoumiYamada_0708.2960.pdf` |
| [OK] | 1811.09274 | 1928 | 26 | `references/noumi_yamada/Clarkson_etal_2020_cyclic_maya_higher_painleve_StudApplMath144.pdf` |
| [OK] | math/0212117 | 187 | 14 | `references/painleve_hierarchy/JoshiMazzocco2003_tritronquee_PII_hierarchy_math0212117.pdf` |
| [OK] | 1107.0214 | 308 | 29 | `references/painleve_hierarchy/Claeys2011_pole_free_PI_hierarchy_1107.0214.pdf` |
| [OK] | 1001.2213 | 257 | 19 | `references/painleve_hierarchy/Claeys2010_PI2_asymptotics_JPhysA43.pdf` |
| [OK] | 1807.04442 | 539 | 10 | `references/rh_numerical/KleinStoilov2018_multidomain_spectral_painleve_SIGMA14.pdf` |
| [OK] | 1111.3527 | 725 | 17 | `references/rh_numerical/ClaeysOlver2011_higher_tracy_widom_1111.3527.pdf` |
| [OK] | 0804.2543 | 541 | 43 | `references/rh_numerical/Bornemann2010_fredholm_determinants_MathComp79.pdf` |
| [OK] | math-ph/9810007 | 254 | 23 | `references/rh_numerical/KitaevKorotkin1998_schlesinger_theta_IMRN_mathph9810007.pdf` |
| [OK] | 1101.3997 | 495 | 47 | `references/rh_numerical/BertolaCafasso2012_noncommutative_PII_CMP309.pdf` |
| [OK] | 2012.05639 | 194 | 13 | `references/rh_numerical/AdlerSokolov2021_matrix_PII_TMP207.pdf` |
| [OK] | 2107.11680 | 281 | 20 | `references/rh_numerical/BobrovaSokolov2021_matrix_PIV_test_2107.11680.pdf` |
| [OK] | 1206.2446 | 4861 | 20 | `references/rh_numerical/WechslbergerBornemann2014_automatic_RH_painleve_ConstrApprox39.pdf` |
| [OK] | 1608.00958 | 635 | 45 | `references/rh_numerical/GavrylenkoLisovyy2018_isomonodromic_tau_CMP363.pdf` |
| [OK] | 1712.08546 | 481 | 26 | `references/rh_numerical/CafassoGavrylenkoLisovyy2019_widom_tau_CMP365.pdf` |
| [OK] | 1010.5725 | 161 | 18 | `references/noumi_yamada/AratynGomesZimerman2010_integrable_origins_higher_painleve_1010.5725.pdf` |
| [OK] | 2009.11668 | 431 | 38 | `references/noumi_yamada/GomezUllateGrandatiMilson2020_rational_painleve_survey_2009.11668.pdf` |
| [OK] | 1401.1408 | 1733 | 36 | `references/painleve_hierarchy/BertolaBothner2015_vorobiev_yablonski_zeros_IMRN.pdf` |
| [OK] | 1504.00440 | 5953 | 23 | `references/painleve_hierarchy/BaloghBertolaBothner2017_generalized_vorobiev_yablonski_ConstrApprox45.pdf` |
| [OK] | 1707.05222 | 1865 | 49 | `references/noumi_yamada/MasoeroRoffelsen2018_PIV_rationals_poles_SIGMA14.pdf` |
| [OK] | q-alg/9708018 | 236 | 27 | `references/noumi_yamada/NoumiYamada1999_PIV_symmetries_okamoto_NagoyaMathJ153.pdf` |
| [OK] | 1705.03295 | 575 | 35 | `references/garnier/CalligarisMazzocco2018_finite_orbits_2variable_garnier_JIntegrSyst3.pdf` |
| [OK] | 1808.09190 | 491 | 38 | `references/garnier/DiarraLoray2020_irregular_garnier_algebraic_Compositio156.pdf` |
| [OK] | nlin/0306020 | 203 | 17 | `references/garnier/KimuraMano2003_irregular_garnier_nlin0306020.pdf` |
| [OK] | 1806.08650 | 650 | 27 | `references/rh_numerical/GavrylenkoIorgovLisovyy2018_FST_system_1806.08650.pdf` |

## Notes

- Pre-2007 IDs (`math/9808003`, `math/0106208`, `math/0212117`, `math-ph/9810007`, `q-alg/9708018`, `nlin/0306020`) used slash-form URLs (`https://arxiv.org/pdf/<ID>.pdf`) and all resolved correctly without needing any fallback.
- No HTTP errors, no HTML error pages written to disk (all curl exits were 0).
- 3-second inter-request sleep maintained throughout for arXiv etiquette.
- Pre-fetched `1502.06695` (ManoTsuda2017) was correctly skipped; it is already at `references/hermite_pade/ManoTsuda2017_hermite_pade_isomonodromic_MathZ285.pdf`.
- The scope document listed 32 IDs net of the skip; all 32 are now present.
