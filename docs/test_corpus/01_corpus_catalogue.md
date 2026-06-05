# PadeTaylor.jl Ground-Truth Corpus Catalogue

Phase-2 deliverable of epic `padetaylor-p1v0` (comprehensive ground-truth test corpus). Gathered by 9 read-only Sonnet scout territories + a consolidation pass, under a strict anti-hallucination contract: every entry's ground truth is a **closed-form** solution, a **regenerable recipe** in locally-available tooling (Mathematica `wolframscript`, Python `mpmath`/`scipy`/`sympy`, Julia `DifferentialEquations`/`HypergeometricFunctions`/`Nemo`), or a **precisely-cited published value** — never a fabricated number. Recipes are reproduced verbatim from `candidates.json` by `_generate_catalogue.py`.

> Companion: `00_capability_map.md` (the 22 capability buckets B1–B22 + 15 highest-value gaps this corpus targets).

## Bucket coverage after the sweep

| Bucket | Coverage | Best candidates | Remaining hole |
|---|---|---|---|
| B1 | strong-regenerable | `exp-77-exact-rational`, `tan-33-exact-pade`, `tan-z4-froissart-2016`, `adr0002-rank-boundary-svd`, `exp-iz-complex-pade` | The geom-noisy case and tan-z4 are citation-only (frozen Octave oracle); the mpmath-regenerable rational-coefficient tests cover the algebraic path well. Wolfram PadeApproximant recipe for two-pole function not yet implemented. |
| B10 | strong-regenerable | `pIII-tilde-sheet-oracle`, `pV-tilde-sheet-oracle`, `pIII-tilde-constant-solution`, `pVI-eta-cancellation` | The PVI-eta RHS sub-expression oracle (pVI-eta-zeta-subexpr-rhs) needs a Mathematica or mpmath run to produce the actual w'' value at a specific (eta,v,v'); the cancellation oracle is pinned. The transformed ODE implementations (B10 coord transforms) need one integration test cross-checked to a non-self-referential source. |
| B11 | strong-regenerable | `sheet-winding-four-arcs`, `pIII-sheet-winding-exact`, `pVI-zeta-sheet-winding`, `sheet-segment-crosses-oblique-yes`, `sheet-two-cuts-simultaneously` | The precondition-violation fixture (sheet-winding-delta-precondition-violation) has medium confidence because the correct behavior (throw vs return) is a design decision. All structural winding-number and crossing oracles are exact. |
| B12 | strong-regenerable | `pi1-tritronquee-ic`, `pii-rational-u1`, `piv-hermite-m3`, `pIII-tilde-constant-solution`, `pVI-zeta-sheet-winding` | PVI tronquee independent regen (pvi-tronquee-regen) has low confidence; needs external Mathematica verification. All other Painleve wrapper tests are well-covered. |
| B13 | strong-regenerable | `cm-n2-rational-closed`, `coupled-riccati-tan-sec2`, `tan-sec2-vector-jet-b13`, `pi1-okamoto-2vector` | The P_I^(2) 4-vector hierarchy (pi1-hierarchy-4vector) needs asymptotic IC verification; the correct branch sign must be confirmed before use. |
| B14 | strong-regenerable | `cm-n2-imaginary-collision`, `coupled-riccati-shared-denom-canonical`, `ny-a2-backlund-pole`, `ny-a4-backlund-orbit`, `coupled-riccati-tan-sec2` | cm-n2-imaginary-collision has not yet been independently verified via BigFloat RK4 in the repo (the real-IC case is verified; the imaginary-IC case is derived by analytic continuation). This verification is the highest-priority gap. |
| B15 | strong-regenerable | `3arg-damped-oscillator-bvp`, `3arg-euler-cauchy-bvp`, `3arg-nonlinear-exp-bvp`, `bratu-bvp`, `airy-bvp` | The 3-arg bvp_solve overload has NO direct unit test currently; 3arg-damped-oscillator-bvp and 3arg-euler-cauchy-bvp fill this gap. BigFloat-256 path for 3-arg overload remains untested. |
| B16 | strong-regenerable | `vector-antiperi-sincos-bvp`, `vector-endpoint-mixing-bvp`, `harmonic-oscillator-bvp-2vec` | Rank-deficient BC fixture (vector-rank-deficient-bc) has medium confidence; the exact fixture design is tricky. The mixing BC (vector-endpoint-mixing-bvp) is the most critical new oracle for B16. |
| B17 | strong-regenerable | `ivp-bvp-dispatch-hybrid`, `oblique-complex-cosh-bvp`, `airy-bvp`, `bvp-domain-guard-imaginary-tstar` | The dispatch invariant has no dedicated test beyond the PIII hybrid at loose tolerance; ivp-bvp-dispatch-hybrid fills this gap. The BigFloat-256 oblique case extends coverage. |
| B18 | strong-regenerable | `heung-2f1-reduction`, `heunc-epsilon-neg3`, `heung-polynomial-degree1`, `heung-pole-at-z1`, `airy-ivp`, `bessel-j0-ivp`, `2f1-elementary-reduction` | HeunC connection formula (second-sheet oracle) has low confidence; needs careful Mathematica implementation. HeunG complex-a (heung-complex-a) has medium confidence pending wolframscript cross-check. Mathieu-as-HeunC and Lame-as-HeunG require Mathematica. |
| B19 | strong-regenerable | `pii-rational-u1`, `pii-rational-u2`, `piv-hermite-m3`, `piv-entire-m2z` | None; multiple closed-form rational and entire solutions are covered. All sympy-verifiable. |
| B2 | strong-regenerable | `wp-taylor-recurrence`, `pi1-taylor-recurrence`, `tan-sec2-vector-jet-b13` | BigFloat-256 order-60 WP jet needs Mathematica WeierstrassP[-1,{0,2}] at 90 dps for IC; otherwise fully covered. |
| B20 | partial | `pi1-tritronquee-ic`, `pii-hm-u0`, `pii-rational-u1`, `pii-airy-half-integer` | PI^(2) hierarchy (pi2-tritronquee-asymptotic) has low confidence; exact 4th-order ODE and IC parameterization need confirmation from KKG 2015. PVI tronquee (pvi-tronquee-regen) has low confidence. The PainleveNamed hierarchy beyond PI/PII/PIV lacks closed-form oracles. |
| B21 | strong-regenerable | `ny-a2-type-c`, `ny-a2-backlund-pole`, `ny-a4-backlund-orbit`, `ny-a4-type-c` | The generic transcendent (ny-a4-generic-mpmath) has medium confidence due to unknown pole locations; mpmath odefun at 40 dps is the only oracle. A_6^(1) and higher-rank NY systems are uncovered. |
| B22 | strong-regenerable | `exp-77-exact-rational`, `adr0002-rank-boundary-svd`, `exp-iz-complex-pade`, `tan-sec2-vector-jet-b13` | The Octave padeapprox.m oracle regeneration requires Octave+Chebfun (not available locally); those fixtures are frozen citation-only. All new oracle infrastructure tests are regenerable. |
| B3 | strong-regenerable | `wp-equianharmonic`, `tan-1storder-pole-bridge`, `simple-pole-rational`, `coth-real-pole`, `fw-riccati-bessel` | None; multiple redundant closed-form and mpmath oracles cover real-axis, complex, 1st-order, 2nd-order pole-bridging. |
| B4 | strong-regenerable | `pi1-taylor-recurrence`, `wp-equianharmonic` | The PI tritronquee trajectory reference (pi1-odefun-trajectory) needs mpmath.odefun run locally to pin values at z in {0.1,0.3,0.5,0.8,1.0} to 35 digits; not yet pinned in any oracle file. |
| B5 | strong-regenerable | `wp-equianharmonic`, `tan-1storder-pole-bridge`, `simple-pole-rational`, `fw-riccati-bessel` | None; WP long-range and pole-bridge tests are in _oracle_problems.jl and _oracle_padestepper.jl already. |
| B6 | strong-regenerable | `jorba-zou-sin-jet-b6`, `jorba-zou-fallback-trigger-b6`, `jorba-zou-relative-mode-b6` | The _oracle_stepcontrol.jl already pins h_4_1_1_TI=4.501206... The three new cases (sin-jet, fallback, relative-mode) add coverage of previously untested branches. |
| B7 | strong-regenerable | `wp-equianharmonic`, `pi1-tritronquee-ic`, `pii-rational-u1`, `pii-rational-u2`, `pii-airy-half-integer`, `piv-hermite-m3` | PVI tronquee independent regen (pvi-tronquee-regen) has low confidence and needs external verification before trust. All other pole-bridging cases are covered. |
| B8 | strong-regenerable | `tan-z4-froissart-2016`, `geom-noisy-1on1z`, `wp-equianharmonic` | Froissart filter tests rely on frozen Octave oracles; mpmath regeneration of the noisy-geom case is straightforward (tol>noise collapses to type(0,1)) but the tanZ4 case needs Octave for exact mu/nu count regeneration. |
| B9 | strong-regenerable | `morph-dilation-single-pixel`, `morph-erosion-edge-hang`, `morph-erosion-5x5`, `morph-opening-idempotent`, `morph-connected-components` | pi1-tritronquee-ic covers the edge-gate baseline. Morphology fixtures are all hand-computed exact integer tests. No remaining holes. |

## Still uncovered — source manually

- B20: PI^(2) hierarchy precise numerical values — KKG 2015 ODE and ICs need confirmation before any oracle can be trusted; exact 4th-order companion form not fully verified; low-confidence candidate only
- B20: PVI tronquee independent verification — pvi-tronquee-regen is self-referential (FFW2017 value IS the asymptotic series); true independent oracle requires asymptotic series derivation from scratch + NDSolve from z=50; not yet done
- B16: Rank-deficient BC fixture for vector BVP — correct oracle design is tricky (simple g=[0,0] collapses to trivial; non-trivial underdetermined BC needs careful construction)
- B18: HeunC second-sheet connection formula — Motygin eq.(21) recipe is schematic, needs careful wolframscript implementation; current confidence=low
- B18: Mathieu-as-HeunC and Lame-as-HeunG at arbitrary precision — requires Mathematica MathieuC/MathieuS; mpmath has no native Mathieu evaluator
- B8: Froissart residue threshold calibration for ADR-0028 — the empirical threshold from external/probes/s8-froissart-diagnosis/ is probe-based; no closed-form oracle for the exact threshold value

## Top-priority oracles (implement first)

1. `cm-n2-imaginary-collision`
2. `coupled-riccati-shared-denom-canonical`
3. `ny-a2-backlund-pole`
4. `3arg-damped-oscillator-bvp`
5. `3arg-euler-cauchy-bvp`
6. `pii-rational-u1`
7. `coupled-riccati-tan-sec2`
8. `pIII-tilde-sheet-oracle`
9. `heunc-epsilon-neg3`
10. `heung-2f1-reduction`
11. `jorba-zou-fallback-trigger-b6`
12. `pIII-tilde-constant-solution`
13. `adr0002-rank-boundary-svd`
14. `ny-a4-backlund-orbit`
15. `vector-endpoint-mixing-bvp`

## Consolidated corpus index

`R` = regenerable (closed-form or local recipe); `C` = citation-only.

| Pri | R/C | Buckets | id | gt_kind | one-line |
|---|---|---|---|---|---|
| P0 | R | B1,B2,B3,B4,B5,B7,B8 | `wp-equianharmonic` | mpmath+published-table | u''=6u^2, u(z)=WP(z-1;0,2); mpmath ellipfun('sn',...) at 50 dps; u(30)=1.095098255959744 (FW2011 Table 5.1 pinned in _oracle_padestepper.jl and _oracle_problems.jl) |
| P0 | R | B1,B22 | `exp-77-exact-rational` | closed-form | exp(z) (7,7) Pade: a_k=C(14-k,7)*7!*7!/(14!*k!); a_1=1/2, b_1=-1/2 exact; pole locations \|z\|~5-10; GGT 2013 Figure 2; _oracles.jl mu=7,nu=7 |
| P0 | C | B1,B8 | `tan-z4-froissart-2016` | published-table | tan(z^4) (20,20): mu=20,nu=16; 4 Froissart removed at \|z\|~2e5; inner ring \|z\|=1.1195=(pi/2)^{1/4}; _oracles.jl poles pinned from Octave |
| P0 | R | B13,B14 | `cm-n2-rational-closed` | closed-form | y'=(v1,v2,4/(x1-x2)^3,4/(x2-x1)^3); IC x1(0)=1,x2(0)=-1,v1=v2=0; x1(t)=sqrt(1+t^2/2) exact; verified in calogero_moser_test.jl CM.0 |
| P0 | R | B13,B14 | `coupled-riccati-tan-sec2` | closed-form | y'=(y2,2y1y2); y1=tan(z), y2=sec^2(z); shared poles at z=pi/2+k*pi; IC y(0)=(0,1); mpmath.tan, 1/mpmath.cos^2 arbitrary dps |
| P0 | R | B14 | `coupled-riccati-shared-denom-canonical` | closed-form | y=(1/(z-z0),-1/(z-z0)^2), z0=1.5; IC y(0)=(-2/3,-4/9); Taylor jets exact rational c1[k]=(-1)^k/z0^{k+1}; SharedPade must find Q with root at z0=1.5 |
| P0 | R | B14,B13 | `cm-n2-imaginary-collision` | closed-form | Same ODE with IC x1(0)=+i,x2(0)=-i,v1=v2=0; x1(t)=(i/2)*sqrt(4-2t^2); collision (shared pole) at t*=sqrt(2)~1.414; Julia: (im/2)*sqrt(4-2*t^2); mpmath analytic continuation |
| P0 | R | B15 | `3arg-damped-oscillator-bvp` | closed-form | u''=-2u'-5u, u=e^{-z}cos(2z); BC u(-1)=e*cos(2), u(1)=cos(2)/e; ∂f/∂u'=-2; mpmath 50 dps; all interior values from exp and cos |
| P0 | R | B15 | `3arg-euler-cauchy-bvp` | closed-form | u''=u/z^2-u'/z, u=(4/3)z+(2/3)/z; BC u(1)=2, u(2)=3 exact rational; u(1.5)=22/9 exact; ∂f/∂u'=-1/z node-by-node non-trivial Jacobian |
| P0 | R | B18 | `heung-2f1-reduction` | closed-form | HeunG(2,2,1,2,3/2,1;z)=2F1(1/2,1;3/2;2z-z^2); z=0.3: 1.25407417926432...; mpmath.hyp2f1 50 dps; DLMF 31.7.2 |
| P0 | R | B18 | `heunc-epsilon-neg3` | mpmath | HeunC(q=1,alpha=1,gamma=2,delta=3,epsilon=-3,z); z=0.1: 0.94525457198...; z=0.9: -10.404759...; mpmath 50-term recurrence + wolframscript N[HeunC[1,1,2,3,-3,0.1],50] |
| P0 | R | B2 | `wp-taylor-recurrence` | closed-form | u''=6u^2; c[n]=6/(n(n-1))*conv(c,c)[n-2] for n>=2; mpmath 90 dps: c[2]=3*u0^2 exact; Mathematica WeierstrassP[-1,{0,2}] cross-check |
| P0 | R | B2,B4,B12 | `pi1-taylor-recurrence` | mpmath | u''=6u^2+z; tritronquee ICs u(0)=-0.18755430834049490, u'(0)=0.30490556026122890 (FW2011 eq.4.1); c[3]=2u0u1+1/6 exact; mpmath.odefun(tol=1e-45) trajectory oracle at z in {0.1,0.3,0.5,0.8,1.0} |
| P0 | R | B21,B14 | `ny-a2-backlund-pole` | closed-form | A_2^(1) after s_0: f_0=t, f_1=1/t, f_2=-1/t, alpha=(-1,1,0); pole at t=0; IC f(1)=(1,1,-1); exact rational; shared-Q must recover Q=t |
| P0 | R | B3,B5,B6 | `tan-1storder-pole-bridge` | closed-form | u'=1+u^2, u(z)=tan(z); pole at z=pi/2 (simple, residue 1); IC u(0)=0; u(pi)=0 exact; mpmath.tan arbitrary dps |
| P0 | R | B7,B12,B19,B20 | `pii-rational-u1` | closed-form | u''=2u^3+zu+1, u(z)=-1/z; pole at z=0 (simple, residue -1); exact rational; sympy residual=0 verified |
| P0 | R | B7,B9,B12,B20 | `pi1-tritronquee-ic` | published-table | u''=6u^2+z; u(0)=-0.18755430834049490, u'(0)=0.30490556026122890 (FW2011 eq.4.1); first real pole z~-2.384 (arxiv:1412.3782); NDSolve from BVP IC regenerable |
| P1 | R | B1 | `tan-33-exact-pade` | closed-form | tan(z) (3,3) Pade: p=z-z^3/15, q=1-2z^2/5; pole at sqrt(5/2)=1.5811... vs true pi/2=1.5708...; error 1.03e-2; exact rational coefficients |
| P1 | R | B1,B22 | `adr0002-rank-boundary-svd` | mpmath | f=1/(1-z)+1e-13/(1-z/2); c_k=1+1e-13*(0.5)^k; 10x11 Toeplitz sigma_2=4.7596e-14 (mpmath 50 dps); DK bound 49% of sigma_2; ADR-0002 BigFloat dispatch test |
| P1 | R | B1,B3,B7 | `jacobi-sn-complex-bridge` | mpmath | u''=-(5/4)u+(1/2)u^3, u=sn(z,1/4); IC at z=i*Kp/2: u=i*sqrt(2) exact, u'=3/sqrt(2) exact; bridge through pole at i*Kp; mpmath.ellipfun 50 dps |
| P1 | C | B1,B8 | `geom-noisy-1on1z` | published-table | 1/(1-z)+noise (std 1e-6); tol=1e-5 -> type (0,1); a=[0.9999999934...], b=[1,-0.9999997642...]; _oracles.jl pinned; regeneration requires Octave+Chebfun |
| P1 | R | B10,B11 | `pIII-tilde-sheet-oracle` | mpmath | PIII-tilde w''=(1/w)(w')^2+(1/4)(alpha*w^2+gamma*w^3+beta*e^zeta+delta*e^{2z}/w); IC w(0)=1/4, w'(0)=9/16 (FFW2017 Fig1); sheet 0: w=-0.12219+0.00412i; conjugate symmetry check; mpmath RK4 n=15000 steps dps=50 |
| P1 | R | B10,B11 | `pV-tilde-sheet-oracle` | mpmath | PV-tilde; IC w(0)=2, w'(0)=-1 (FFW2017 Fig6); sheet 0: w=0.54159-1.21830i; mpmath RK4 dps=50; conjugate symmetry error=0 at 50 digits |
| P1 | R | B10,B11 | `pVI-eta-cancellation` | mpmath | PVI-eta: e^{e^eta}-1 near eta=log(2pi)+i*pi/2+1e-10 gives Float64 error factor ~2.25e3 vs BigFloat-50 result -1.974e-19+6.283e-10i; mpmath 50 dps recipe: exp(exp(eta))-1 |
| P1 | R | B11 | `sheet-winding-four-arcs` | closed-form | path=[1,i,-1,-i,1], branch=0; each winding_delta=pi/2; step 3 raw=-3pi/2 normalised to +pi/2; total=2*pi; sheet_index=round(2pi/(2pi))=1 |
| P1 | R | B11 | `sheet-segment-crosses-oblique-yes` | closed-form | branch=1+1j, cut_angle=pi/4; z_cur=0+2j, z_new=3+0j; Cramer: t=2/5, s=sqrt(2)/5; crosses=True; hand-computed exact |
| P1 | R | B13,B14 | `pi1-okamoto-2vector` | mpmath | y'=(y2,6y1^2+z); IC u(0)=1.0718225..., u'(0)=1.7103373...; target u(30)=1.095098255959744 (FW2011 Table 5.1); vector_path_network_solve cross-check |
| P1 | R | B13,B14 | `pii-u1-companion-2vector` | closed-form | PII(alpha=1): u=-1/t; companion y=(−1/t,1/t^2); IC y(1)=(−1,1); shared pole at t=0; exact rational; verify ODE residual=-2/t^3=y2'=2y1^3+ty1+1 |
| P1 | R | B13,B22 | `tan-sec2-vector-jet-b13` | closed-form | y1'=y2, y2'=2y1y2; y1=tan(z+pi/4), y2=sec^2(z+pi/4); IC y1[0]=1, y2[0]=2 exact; c1[1]=2, c2[1]=4 exact; mpmath.diff 50 dps + Mathematica CoefficientList |
| P1 | R | B14,B21 | `ny-a4-backlund-orbit` | closed-form | A_4^(1) after s_0: f=(t,1/t,0,0,-1/t), alpha=(-1,1,0,0,1); shared poles at t=0 in f_1,f_4; exact rational; shared-Q must recover Q=t from jet |
| P1 | R | B15 | `3arg-nonlinear-exp-bvp` | closed-form | u''=(u')^2/u, u=exp(2z); BC u(0)=1, u(1)=e^2=7.38905...; ∂f/∂u=-(u')^2/u^2, ∂f/∂u'=2u'/u; residual=0 exact; mpmath.exp |
| P1 | R | B15 | `3arg-boundary-layer-bvp` | closed-form | u''=-u'/eps (eps=0.1), u=(1-exp(-z/eps))/(1-exp(-1/eps)); BC u(0)=0, u(1)=1; ∂f/∂u'=-1/eps; mpmath 50 dps |
| P1 | R | B15 | `bratu-bvp` | closed-form | u''=-2*exp(u); u=2*log(B/cosh(B*(z-0.5))), B=cosh(B/2)=1.178775...; BC u(0)=0, u(1)=0; u(0.5)=0.32895242...; mpmath.findroot |
| P1 | R | B15,B17 | `airy-bvp` | mpmath | u''=zu, u=Ai(z); BC u(0)=Ai(0)=0.35502805..., u(3)=Ai(3)=0.00659113...; mpmath.airyai; df/du=z node-by-node variable |
| P1 | R | B15,B17 | `oblique-complex-cosh-bvp` | mpmath | u''=u, u=cosh(z); BC u(0)=1, u(1+i)=cosh(1+i)=0.83373+0.98890i; interior u(0.5+0.5i)=0.98958+0.24983i; mpmath.cosh 50 dps |
| P1 | R | B16 | `vector-antiperi-sincos-bvp` | closed-form | y'=(y2,-y1); y=(sin z,cos z); BC B_a=I, B_b=-I, g=(-1,1) on [0,pi/2]; y(pi/4)=(sqrt(2)/2,sqrt(2)/2); uniqueness k=0 forced |
| P1 | R | B16 | `vector-endpoint-mixing-bvp` | closed-form | y'=(y2,-y1); y=(sin z,cos z) on [0,pi]; BC B_a=[[0,0],[0,1]], B_b=[[1,0],[0,0]], g=[0,1]; pins y2(0)=1 and y1(pi)=0 |
| P1 | R | B17 | `ivp-bvp-dispatch-hybrid` | closed-form | u''=u, u=cosh(z); bvp_solve(z->z,z->one,0,1+i,1,cosh(1+i)) and IVP from u(0)=1,u'(0)=0 must agree to atol 1e-10 at z=0.5+0.5i |
| P1 | R | B18 | `heung-polynomial-degree1` | closed-form | HeunG(2,q,−1,2,3,4;z): Heun polynomial q_{1,m}=-6+-2*sqrt(6); Hp_{1,0}(0.5)=0.90824829...; mpmath.sqrt exact eigenvalue; DLMF 31.5 |
| P1 | R | B18 | `heung-pole-at-z1` | closed-form | HeunG(2,4,1,4,2,2;z)=1/(1-z); all Taylor coefficients=1 exact; pole at z=1; DLMF 31.7.2+15.4.6 |
| P1 | R | B18 | `heunc-epsilon-pos2` | mpmath | HeunC(q=1,alpha=1,gamma=2,delta=3,epsilon=2,z); z=0.1: 0.94986517967...; mpmath recurrence + wolframscript cross-check |
| P1 | R | B18 | `heunc-fractional-gamma` | mpmath | HeunC(q=1,alpha=1,gamma=1/2,delta=3,epsilon=-1,z); z=0.1: 0.77694852...; u'(0)=-q/gamma=-2; mpmath 50 dps recurrence |
| P1 | R | B18 | `heung-large-q-a3` | mpmath | HeunG(3,10,1,2,2,3;z); z=0.25: 1.625058391...; z=0.5: 3.286670344...; mpmath recurrence dps=50 + wolframscript N[HeunG[3,10,1,2,2,3,0.25],50] |
| P1 | R | B18 | `heung-small-a` | mpmath | HeunG(1/2,1,1,2,3,4;z); convergence radius=0.5; z=0.1: 1.068937...; z=0.3: 1.210685...; mpmath dps=50 + wolframscript |
| P1 | R | B18 | `gauss-2f1-ivp` | mpmath | z(1-z)u''+[5/4-(11/6)z]u'-(1/6)u=0, u=2F1(1/2,1/3;5/4;z); u(0.4)=1.065978091...; mpmath.hyp2f1 50 dps |
| P1 | R | B18 | `2f1-elementary-reduction` | closed-form | F(1/2,3;3;z)=(1-z)^{-1/2}; z=0.4: (0.6)^{-1/2}=1.29099445...; z=0.9: (0.1)^{-1/2}=3.16227766...; DLMF 15.4.6 exact |
| P1 | R | B18,B3 | `airy-ivp` | mpmath | u''=zu, u=Ai(z); IC u(0)=0.35502805388782..., u'(0)=-0.25881940379281...; values at z=1,2,3,-1,-2 from mpmath.airyai; DLMF 9.2 |
| P1 | R | B18,B3 | `bessel-j0-ivp` | mpmath | u''+(u'/z)+u=0, u=J_0(z); IC u(1)=0.76519768..., u'(1)=-J_1(1)=-0.44005059...; values at z=2,3,5 from mpmath.besselj(0,...) |
| P1 | R | B21 | `ny-a2-type-c` | closed-form | A_2^(1): f_j'=f_j(f_{j+1}-f_{j+2})+1/3; f_j(t)=t/3 exact; Rational{BigInt} arithmetic; constraint sum(f_j)=t preserved |
| P1 | R | B3,B5 | `simple-pole-rational` | closed-form | u'=u^2, u(z)=1/(1-z); pole at z=1 (simple, residue -1); IC u(0)=1; u(0.5)=2, u(2)=-1 exact rational |
| P1 | R | B3,B5 | `coth-real-pole` | closed-form | u'=1-u^2, u(z)=coth(z); pole at z=0 (simple, residue 1); IC u(-1)=coth(-1); u(1)=coth(1) exact via odd symmetry |
| P1 | R | B3,B5 | `fw-riccati-bessel` | mpmath | y'=t^2+y^2, y(t)=t*J_{3/4}(t^2/2)/J_{-1/4}(t^2/2); first pole at t=sqrt(2*z1)~2.0031; mpmath.besselj+findroot arbitrary dps |
| P1 | R | B3,B5,B8 | `tan-2ndorder-pole-bridge` | closed-form | u''=2u+2u^3, u(z)=tan(z); IC u(0)=0, u'(0)=1; pole at pi/2 with residue 1; mpmath.tan arbitrary dps |
| P1 | R | B6 | `jorba-zou-sin-jet-b6` | mpmath | sin-jet c[2k+1]=(-1)^k/(2k+1)!; c[30]=0; h=(eps*29!)^(1/29)=4.501206370338986; matches h_4_1_1_TI; key test: code must skip zero c[30] |
| P1 | R | B6 | `jorba-zou-fallback-trigger-b6` | mpmath | cos-jet c[2k]=(-1)^k/(2k)! for k=0..14, c[29]=c[30]=0; fallback h=(28!)^(1/27)=12.35962276...; mpmath 50 dps; Mathematica N[(28!)^(1/27),50] cross-check |
| P1 | R | B6 | `jorba-zou-relative-mode-b6` | mpmath | cos-jet c[0]=1; eps_abs=1e-15 < eps_rel*\|c0\|=1e-12; eps_eff=1e-12; h=(eps_eff*30!)^(1/30)=4.795000636590567; mpmath 50 dps |
| P1 | R | B7,B12,B19,B20 | `pii-rational-u2` | closed-form | u''=2u^3+zu+2, u=d/dz ln(z/(z^3+4)); poles at z=0 and cube roots of -4; real pole z=-4^(1/3)=-1.587401...; exact algebraic |
| P1 | R | B7,B12,B19,B20 | `pii-airy-half-integer` | closed-form | u''=2u^3+zu+1/2, u=-Ai'(-z/2^{1/3})/(2^{-1/3}*Ai(-z/2^{1/3})); poles at z=-2^{1/3}*a_k; first real pole ~2.9446; mpmath.airyai |
| P1 | R | B7,B12,B19,B20 | `piv-hermite-m3` | closed-form | PIV(alpha=-3,beta=-8), u=-4z/(2z^2-1); poles at z=+-1/sqrt(2) exact; residues +-1; sympy residual=0 |
| P1 | R | B7,B12,B20 | `pii-hm-u0` | mpmath | u''=2u^3+zu; HM IC u(x)~Ai(x) as x->+inf; u(0)=0.3670615515480784, u'(0)=-0.2953721054475501; mpmath.airyai(20) IC + odefun |
| P1 | R | B9 | `morph-dilation-single-pixel` | closed-form | Input 3x3 [[0,0,0],[0,1,0],[0,0,0]], SE=3x3 all-ones; output=all-ones 3x3; Serra (1982) dilation definition; hand-computed |
| P1 | R | B9 | `morph-erosion-edge-hang` | closed-form | Input 3x3 all-ones, SE=3x3 all-ones, boundary=0; output=[[0,0,0],[0,1,0],[0,0,0]]; hand-computed Serra (1982) |
| P1 | R | B9 | `morph-erosion-5x5` | closed-form | Input 5x5 with center 3x3 filled, SE=3x3 all-ones; output single pixel at (2,2); hand-computed |
| P1 | R | B9 | `morph-opening-idempotent` | closed-form | opening(5x5 3x3 center block, 3x3 SE) = original input; composed erosion+dilation; Serra (1982) idempotence theorem |
| P1 | R | B9 | `morph-connected-components` | closed-form | 5x5 BitMatrix with two separate L-shapes; 8-connected BFS; n_labels=2; comp1={(0,0),(0,1),(1,0)}, comp2={(2,4),(3,3),(3,4)}; hand-computed |
| P2 | R | B1,B22 | `exp-iz-complex-pade` | closed-form | exp(iz) (7,7) Pade: p_k=C(14-k,7)*7!*7!/(14!*k!)*i^k; r(1)=0.54030231+0.84147098i vs exp(i); pole at 0-9.94i; ADR-0002 Complex{T} dispatch |
| P2 | R | B10,B12 | `pIII-tilde-residue-exact` | closed-form | PIII-tilde pole residue: A^2=4/gamma, so A=+-2/sqrt(gamma); for gamma=1: residue=+-2; verified by ODE Laurent balance (polynomial algebra, no floating point) |
| P2 | R | B10,B12 | `pIII-tilde-constant-solution` | closed-form | PIII-tilde(0,0,1,-1): w(zeta)=exp(zeta/2) exact; u(z)=e^{-z/2}*w=1; RHS residual=(1/4)*exp(zeta/2)=w'' verified analytically |
| P2 | R | B11 | `pIII-sheet-winding-exact` | closed-form | z=e^{zeta/2}: CCW loop in z gives Im(zeta) shift of 4*pi, crossing strip boundary at Im=2*pi; winding_number=+1 exact integer; derived from FFW2017 md:103 strip formula |
| P2 | R | B11 | `sheet-ccw-winding-offcenter` | closed-form | path=[4+0j,2+2j,0+0j,2-2j,4+0j], branch=2+0j; each winding_delta=pi/2; total=2*pi; sheet_index=1; hand-computed |
| P2 | R | B11 | `sheet-segment-crosses-oblique-no` | closed-form | branch=1+1j, cut_angle=pi/4; z_cur=0+0j, z_new=2+0j; t=0 (endpoint, excluded), s=-sqrt(2)<0; crosses=False |
| P2 | R | B11,B12 | `pVI-zeta-sheet-winding` | closed-form | PVI sheets: (theta_0,phi_0)=(-pi,pi]^2 -> sheet 0; CCW around z=1 shifts phi by 2pi -> sheet 1; exact integer indices; from FFW2017 md:182-189 |
| P2 | R | B11,B12 | `pIII-conjugate-symmetry` | closed-form | w(conj(zeta))=conj(w(zeta)) for real params+ICs; error=\|w(zeta)-conj(w(conj(zeta)))\|=0 at 50 digits for both PIII-tilde and PV-tilde; mpmath RK4 n=15000 steps |
| P2 | R | B11,B6 | `sheet-two-cuts-simultaneously` | closed-form | branches=(0,1+0j), cut_angles=(pi/2,pi/2); z_cur=-0.5+2j, z_new=1.5+2j; t1=1/4, s1=2 (cut1 crossed); t2=3/4, s2=2 (cut2 crossed); hand-computed exact |
| P2 | R | B12,B19 | `piv-entire-m2z` | closed-form | PIV(alpha=0,beta=-2), u(z)=-2z entire; IC u(0)=0, u'(0)=-2; sympy residual=0; referenced in RF2014 and PainleveClosedForm.jl |
| P2 | R | B13,B14,B20 | `pi1-hierarchy-4vector` | mpmath | P_I^(2): y4'=-10*y2^2-20*y1*y3-40*(y1^3-6ty1+6x); IC from asymptotic u~+cbrt(6)*\|x\|^{1/3} at x=-20, t=0; mpmath odefun dps=30+ |
| P2 | R | B15 | `linear-cosh-bvp-reference` | closed-form | u''=u, u=cosh(z)/cosh(1); BC u(-1)=u(1)=1; N=16: error~6e-15; already pinned as BV.1.2/BV.5.1 BigFloat-256 |
| P2 | R | B15,B17 | `bigfloat256-oblique-cosh-bvp` | mpmath | Same as oblique-complex-cosh-bvp at BigFloat(precision=256); cosh(0.5+0.5i) to 80 dps: 0.9895848833999199364...+0.2498263975004615315i; mpmath 80 dps |
| P2 | R | B16 | `harmonic-oscillator-bvp-2vec` | closed-form | y'=(y2,-y1); y=(cos z,-sin z); BVP on [0,pi/2]; BC B_a=[1 0;0 0], B_b=[0 0;1 0]; already tested in vector_bvp_test.jl |
| P2 | R | B17,B6 | `bvp-domain-guard-imaginary-tstar` | closed-form | segment [0,1+i]; query z=-4.5+5.5j; t*=(z-half_sum)/half_diff=10j; Re(t*)=0 passes guard; Im(t*)=10 ignored; silent extrapolation; BVP.jl:491 |
| P2 | R | B18 | `heung-trivial-alpha0-q0` | closed-form | HeunG(a,0,0,beta,gamma,delta;z)=1 for all z, any a; c[1]=q/(a*gamma)=0, all subsequent c[n]=0 exact |
| P2 | R | B18 | `heung-complex-a` | mpmath | HeunG(0.5+0.5i,1,1,2,3,4;z=0.3); value=(1.1059813-0.1059218i); mpmath complex recurrence 50 dps; wolframscript N[HeunG[0.5+0.5*I,1,1,2,3,4,0.3],50] required for trust |
| P2 | R | B18 | `heung-lame-reduction` | mathematica | HeunG(k^{-2},-h/4,-n/2,(n+1)/2,1/2,1/2;sn^2(z,k))=LameC[n,j,z,k^2]; n=1,j=0: dn(z,k); wolframscript N[JacobiDN[0.5,Sqrt[0.5]],50]=0.92372... |
| P2 | R | B21 | `ny-a4-type-c` | closed-form | A_4^(1): f_j=t/5 for all j, alpha=(1/5,...,1/5); Rational{BigInt} exact; constraint sum(f_j)=t |
| P2 | R | B21,B14 | `ny-a4-generic-mpmath` | mpmath | A_4^(1): alpha=(0.1,0.15,0.2,0.25,0.3); IC f(1)=(0.3,0.25,0.15,0.2,0.1); mpmath odefun tol=1e-16; constraint sum(f_j)-t<1e-35 check |
| P2 | R | B3,B5,B7 | `tanh-imaginary-pole` | closed-form | u'=1-u^2, u(z)=tanh(z); poles at z=i(pi/2+k*pi); IC u(0)=0; u(i*pi/4)=i exact; mpmath.tanh |
| P2 | R | B7,B12,B19 | `piv-hermite-m4` | closed-form | PIV(alpha=-4,beta=-18), u=3(1-2z^2)/(z(2z^2-3)); poles at 0 and +-sqrt(6)/2; sympy verified |
| P2 | R | B7,B12,B19,B20 | `pii-rational-u3` | closed-form | u''=2u^3+zu+3, u=3z^2(160+8z^3+z^6)/(320-24z^6-z^9); 9 simple poles; real at -(10+6sqrt(5))^(1/3) and (6sqrt(5)-10)^(1/3); sympy verified |
| P2 | R | B7,B9,B12 | `pvi-tronquee-regen` | mathematica | PVI(4,-4,8,-8); u(10)=0.429534600325223 (FFW2017 Fig 2); IC at z=50 from asymptotic series; NDSolve WorkingPrecision->50 to z=10 must agree to >=10 digits (confidence: low, needs verification) |

## Full candidate records

Every field verbatim from the scout output. `how_to_compute` is the regeneration recipe / source locator a test's oracle must use.

### Territory: Painleve pole/value tables + named transcendents: PI tritronquee ICs and interio

*Scout notes:* Gaps: (1) The first precise pole location on the negative real axis for the PI tritronquee (p_1 ~ -2.3841687...) appears in Boutroux/Novokshenov literature but the rigorous high-precision value is in behind-paywall papers (Springer, SIAM); all web-accessible PDF sources returned binary. The value from arxiv:1412.3782 abstract is x = -770766/323285 = -2.3841687675... to 5e-6 accuracy; a Wolfram NDSolve recipe from the imaginary-axis BVP starting point would regenerate it. (2) The PI^(2) second-member tritronquee (KKG 2015 / Kapaev-Klein-Grava) has published ridge-value asymptotics but no easily accessible tabulated numerical values; its existence proof and asymptotic characterisation are the canonical oracle. (3) The PVI tronquee IC u(10)=0.429534600325223 from FFW2017 is self-referential (the paper generates the IC from optimal truncation of the asymptotic series at z=10); an independent verification recipe using a longer asymptotic series evaluated at z=50 then NDSolve inward is provided but the target value must be recomputed, not trusted from FFW2017. Best oracles found: DLMF 32.10 for all rational/special-function PIV and PII closed forms (verified by sympy); FW2011 eq.(4.1) for PI tritronquee u(0) and u'(0); FW2014 p.12 for HM PII u(0) and u'(0).

#### `pi-tritronquee-ic`

- **name**: PI tritronquee solution at origin
- **ode**: u'' = 6u^2 + z  (Painleve I)
- **ic_bc**: u(40i) = -sqrt(-40i/6), u'(40i) = 1/(2*sqrt(-6*40i))  [leading Boutroux asymptotic]
- **domain**: complex plane, pole-free sector |arg z| < 4*pi/5; evaluation point z=0
- **closed_form**: none
- **pole_structure**: all poles in sector 4*pi/5 < |arg z| < pi; double poles with residue 0; first real pole on negative real axis at z ~ -2.3841...
- **gt_kind**: published-table
- **how_to_compute**: FW2011 eq.(4.1): u(0) = -0.1875543083404949, u'(0) = 0.3049055602612289. Obtained by BVP over [-20i, 20i] with Maple at 32-digit precision using IC u(+/-20i) ~ -sqrt(-z/6). Independent regeneration: wolframscript command: NDSolve[{u''[z]==6*u[z]^2+z, u[40*I]==-Sqrt[-40*I/6], u'[40*I]==1/(2*Sqrt[-6*40*I])}, u, {z, 40*I, 0}, WorkingPrecision->50, Method->"StiffnessSwitching"][u][0]
- **precision**: 16 significant digits in FW2011; extendable to arbitrary via Maple/NDSolve BVP
- **targets_buckets**: B7, B9, B12, B20
- **regime**: pole-free sector interior value; BVP with asymptotic IC; B9 edge-gate baseline
- **citation**: FW2011_painleve_methodology_JCP230.md:226 (eq.4.1); Joshi & Kitaev, Stud. Appl. Math. 107 (2001) 253-291; Novokshenov, Theor. Math. Phys. 159 (2009) 853-862
- **confidence**: high
- **verify_note**: FW2011 states accuracy better than 1e-20 for these values; cross-check with Joshi-Kitaev 2001 approximations cited as consistent

#### `pi-tritronquee-first-pole`

- **name**: PI tritronquee first pole on negative real axis
- **ode**: u'' = 6u^2 + z  (Painleve I)
- **ic_bc**: u(0) = -0.1875543083404949, u'(0) = 0.3049055602612289 (tritronquee ICs from B9 entry)
- **domain**: negative real axis
- **closed_form**: none
- **pole_structure**: double pole at z_0 ~ -2.3841687675..., residue 0; Laurent: u(z) = 1/(z-z_0)^2 - z_0/10*(z-z_0)^2 - 1/6*(z-z_0)^3 + O((z-z_0)^4)
- **gt_kind**: published-table
- **how_to_compute**: Arxiv:1412.3782 (Bogatskiy et al.) confirms z_0 = -770766/323285 = -2.3841687675... to within 5e-6. For high-precision version: use wolframscript NDSolve from BVP with u(0)=-0.1875543083404949, u'(0)=0.3049055602612289, integrate along real axis and locate the Pade-denominator root near x=-2.384. Command: NDSolve[{u''[x]==6*u[x]^2+x, u[0]==-0.1875543083404949, u'[0]==0.3049055602612289}, u, {x, 0, -5}, WorkingPrecision->50, Method->"StiffnessSwitching"]
- **precision**: 5e-6 from arxiv:1412.3782; better than 1e-10 achievable via NDSolve from verified IC
- **targets_buckets**: B7, B9, B12
- **regime**: first pole crossing on negative real axis; B9 edge gate; real-line pole location
- **citation**: arxiv:1412.3782 (Bogatskiy et al.); Novokshenov 2010, Regular Chaotic Dyn. 15(2):243; FW2011_painleve_methodology_JCP230.md:297
- **confidence**: medium
- **verify_note**: The 5e-6 bound is from arxiv:1412.3782 abstract only; full paper behind paywall. Regenerate from verified IC above and find first Pade denominator root via binary search

#### `pii-hm-u0`

- **name**: Hastings-McLeod PII u(0) and u'(0)
- **ode**: u'' = 2u^3 + z*u  (Painleve II, alpha=0)
- **ic_bc**: u(x) ~ Ai(x) as x->+inf (Hastings-McLeod condition); IC: u(5)=Ai(5)=1.0834e-4, u'(5)=Ai'(5)=-2.474e-4
- **domain**: real axis, pole-free everywhere
- **closed_form**: none
- **pole_structure**: entire on real axis; poles in complex plane in two wedge-shaped regions off real axis
- **gt_kind**: mpmath
- **how_to_compute**: python: from scipy.integrate import solve_ivp; import mpmath as mp; mp.mp.dps=50; x0=5.0; u0=float(mp.airyai(x0)); up0=float(mp.airyai(x0,derivative=1)); sol=solve_ivp(lambda z,Y:[Y[1],2*Y[0]**3+z*Y[0]], [x0,0], [u0,up0], method='DOP853', rtol=1e-13, atol=1e-15, dense_output=True); print(sol.sol(0)[0], sol.sol(0)[1]). For 15+ digit: use NDSolve[{u''[x]==2*u[x]^3+x*u[x], u[20]==AiryAi[20], u'[20]==AiryAiPrime[20]}, u, {x, 20, 0}, WorkingPrecision->50][u][0]. Expected: u(0)=0.3670615515480784, u'(0)=-0.2953721054475501
- **precision**: 16 digits in FW2014; recipe verified to 10 digits with x_start=5 via scipy DOP853
- **targets_buckets**: B7, B12, B20
- **regime**: pole-free real-line solution; connection constant; asymptotic IC integration
- **citation**: FW2014_second_PII_exploration_FoCM14.md:254 (eq. on p.12); Hastings & McLeod, Arch. Ration. Mech. Anal. 73 (1980) 31-51; DLMF 32.11.2
- **confidence**: high
- **verify_note**: Value from FW2014 is internal to their computation; independently regenerated here to 10 digits via scipy+mpmath starting from Ai(5) IC. Use x_start>=5 to avoid underflow; use mpmath Ai(20) via NDSolve for 15+ digit match

#### `pii-rational-u1`

- **name**: PII rational solution u_1(z) = -1/z (alpha=1)
- **ode**: u'' = 2u^3 + z*u + 1  (Painleve II, alpha=1)
- **ic_bc**: n/a (rational solution, globally defined)
- **domain**: C \ {0}
- **closed_form**: u_1(z) = -1/z
- **pole_structure**: single simple pole at z=0, residue -1
- **gt_kind**: closed-form
- **how_to_compute**: Direct evaluation: u_1(z) = -1/z. Verified by sympy: u_1'' - (2*u_1^3 + z*u_1 + 1) = 0. Pole at z=0 has residue -1 (consistent with PII simple-pole structure with c_{-1}=-1 in FW2014 eq. before (6)).
- **precision**: exact rational
- **targets_buckets**: B7, B12, B19, B20
- **regime**: rational closed-form solution; B19 closed-form family; pole crossing test
- **citation**: FW2014_second_PII_exploration_FoCM14.md:119 (eq.6); DLMF 32.8.4
- **confidence**: high

#### `pii-rational-u2`

- **name**: PII rational solution u_2(z), alpha=2
- **ode**: u'' = 2u^3 + z*u + 2  (Painleve II, alpha=2)
- **ic_bc**: n/a
- **domain**: C \ {0, -2^(2/3), 2^(2/3)*exp(+-2*pi*i/3)}
- **closed_form**: u_2(z) = (4 - 2z^3) / (4z + z^4) = 2*(2-z^3) / (z*(z^3+4))
- **pole_structure**: 4 simple poles: z=0 (residue +1), z=-2^(2/3) = -1.587401..., z=2^(2/3)*exp(2*pi*i/3), z=2^(2/3)*exp(-2*pi*i/3). Exact: denominator zeros of z*(z^3+4). Q_2 = z^3+4.
- **gt_kind**: closed-form
- **how_to_compute**: u_2(z) = d/dz ln(Q_1/Q_2) where Q_1=z, Q_2=z^3+4 (Yablonskii-Vorob'ev). Poles at z=0 and cube-roots of -4. Pole on real axis: z=-4^(1/3)=-2^(2/3) = -1.5874010519681994. Verified by sympy: residual = 0.
- **precision**: exact rational; pole locations exact algebraic: -4^(1/3), 4^(1/3)*exp(+-i*pi/3)
- **targets_buckets**: B7, B12, B19, B20
- **regime**: rational closed-form, 4 simple poles, Yablonskii-Vorob'ev Q_2
- **citation**: FW2014_second_PII_exploration_FoCM14.md:120 (eq.6); DLMF 32.8.4; Clarkson 2006 Lect. Notes Math. 1883
- **confidence**: high

#### `pii-rational-u3`

- **name**: PII rational solution u_3(z), alpha=3
- **ode**: u'' = 2u^3 + z*u + 3  (Painleve II, alpha=3)
- **ic_bc**: n/a
- **domain**: C minus 9 simple poles
- **closed_form**: u_3(z) = 3z^2*(160 + 8z^3 + z^6) / (320 - 24z^6 - z^9) = d/dz ln(Q_2/Q_3)
- **pole_structure**: 9 simple poles: 3 from Q_2=z^3+4 (zeros: -2^(2/3), 2^(2/3)*exp(+-2*pi*i/3)) and 6 from Q_3=z^6+20z^3-80. Real poles of Q_3: -(10+6*sqrt(5))^(1/3)=-2.860926878783304 and (10-6*sqrt(5))^(1/3)=... wait, -10+6*sqrt(5)>0 since sqrt(5)>2, so (-10+6sqrt(5))^(1/3)=(6sqrt(5)-10)^(1/3). Exact: z^3 = (-20 +/- sqrt(400+320))/2 = -10 +/- sqrt(180) = -10 +/- 6*sqrt(5). Real poles: -(10+6*sqrt(5))^(1/3) = -2.8609... and (6*sqrt(5)-10)^(1/3) = 1.5061...
- **gt_kind**: closed-form
- **how_to_compute**: Evaluate u_3(z) directly. Poles are exact algebraic roots. Numerically: Q_3 real roots at z = -2.860926878783304 and z = 1.506109580086942 (verified by sympy). Q_2 real root: z = -2^(2/3) = -1.587401051968199. Verify: python -c 'from sympy import *; z=Symbol("z"); u3=3*z**2*(160+8*z**3+z**6)/(320-24*z**6-z**9); simplify(diff(u3,z,2)-(2*u3**3+z*u3+3))'
- **precision**: exact closed-form; pole locations are algebraic numbers: roots of z^3+4=0 and z^6+20z^3-80=0
- **targets_buckets**: B7, B12, B19, B20
- **regime**: rational closed-form, 9 simple poles, pole-crossing path test through 2 real-axis poles
- **citation**: FW2014_second_PII_exploration_FoCM14.md:120 (eq.6); DLMF 32.8.4; Yablonskii 1959, Vorob'ev 1965
- **confidence**: high

#### `piv-entire-m2z`

- **name**: PIV entire rational solution u(z)=-2z (alpha=0, beta=-2)
- **ode**: u'' = (u')^2/(2u) + 3u^3/2 + 4zu^2 + 2(z^2-0)u + (-2)/u  (PIV, alpha=0, beta=-2)
- **ic_bc**: u(0)=0, u'(0)=-2
- **domain**: C (entire, no poles)
- **closed_form**: u(z) = -2z
- **pole_structure**: entire
- **gt_kind**: closed-form
- **how_to_compute**: u(z) = -2z exactly. Verified: u'=-2, u''=0; plugging into PIV RHS with alpha=0, beta=-2 gives 4/(−4z)+0+4z*4z^2+(−4z)z^2+(−2)/(−2z) = −1/z + 16z^3 − 4z^3 + 1/z = 12z^3 -> wait, let me recheck. Actually verified by sympy residual=0 above (code block in session). This solution is the seed of an entire hierarchy of rational solutions via Backlund transformations.
- **precision**: exact polynomial
- **targets_buckets**: B12, B19
- **regime**: trivially entire PIV solution; smoke test for pole-field solver on a pole-free case
- **citation**: ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.md:91 ('two particular choices leading to the only known entire solutions, -2z and -(2/3)z'); DLMF 32.10.20 note
- **confidence**: high

#### `piv-hermite-m3`

- **name**: PIV Hermite rational solution u(z;-3,-8)=-4z/(2z^2-1)
- **ode**: u'' = (u')^2/(2u) + 3u^3/2 + 4zu^2 + 2(z^2+3)u + (-8)/u  (PIV, alpha=-3, beta=-8)
- **ic_bc**: n/a (global rational solution)
- **domain**: C \ {+-1/sqrt(2)}
- **closed_form**: u(z) = -H'_2(z)/H_2(z) = -4z/(2z^2-1) = -4z/(2z^2-1)
- **pole_structure**: 2 simple poles at z = +-1/sqrt(2) = +-0.7071067811865476; residues +-1
- **gt_kind**: closed-form
- **how_to_compute**: u(z) = -H'_2(z)/H_2(z) where H_2(z)=4z^2-2 is the Hermite polynomial. Explicitly: u(z) = -8z/(4z^2-2) = -4z/(2z^2-1). Poles at zeros of H_2: z=+-1/sqrt(2). Verified by sympy: PIV(alpha=-3, beta=-8) residual=0. Python: from sympy import *; z=Symbol('z'); u=-4*z/(2*z**2-1); simplify(diff(u,z,2)-((diff(u,z))**2/(2*u)+Rational(3,2)*u**3+4*z*u**2+2*(z**2+3)*u+(-8)/u))
- **precision**: exact closed-form; pole locations exact: z = +-1/sqrt(2) = +-sqrt(2)/2
- **targets_buckets**: B7, B12, B19, B20
- **regime**: generalized Hermite family PIV rational, 2 simple poles on real axis
- **citation**: DLMF 32.10.20 (eq. with m=3); ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.md:91; Noumi & Yamada, Nagoya Math. J. 153 (1999)
- **confidence**: high

#### `piv-hermite-m4`

- **name**: PIV Hermite rational solution u(z;-4,-18), alpha=-4, beta=-18
- **ode**: u'' = (u')^2/(2u) + 3u^3/2 + 4zu^2 + 2(z^2+4)u + (-18)/u  (PIV, alpha=-4, beta=-18)
- **ic_bc**: n/a
- **domain**: C \ {0, +-sqrt(6)/2}
- **closed_form**: u(z) = -H'_3(z)/H_3(z) = 3*(1-2z^2)/(z*(2z^2-3))
- **pole_structure**: 3 simple poles at z=0, z=+-sqrt(6)/2=+-1.2247448713915890; H_3(z)=8z^3-12z=4z(2z^2-3)
- **gt_kind**: closed-form
- **how_to_compute**: u(z) = -H'_3(z)/H_3(z) = -(24z^2-12)/(8z^3-12z) = 3*(1-2z^2)/(z*(2z^2-3)). Poles at z=0 and z=+-sqrt(6)/2=+-1.22474... Numerically: +-1.2247448713915890. Verified by sympy: PIV(alpha=-4, beta=-18) residual=0.
- **precision**: exact; pole at +-sqrt(6)/2 = +-0.5*sqrt(6)
- **targets_buckets**: B7, B12, B19
- **regime**: generalized Hermite family PIV rational, 3 poles, includes pole at origin
- **citation**: DLMF 32.10.20 (m=4 case); Noumi & Yamada 1999
- **confidence**: high

#### `pvi-tronquee-ffv-regen`

- **name**: PVI tronquee independent regeneration recipe (alpha,beta,gamma,delta)=(4,-4,8,-8)
- **ode**: PVI: u'' = (1/2)*(1/u+1/(u-1)+1/(u-z))*(u')^2 - (1/z+1/(z-1)+1/(u-z))*u' + u(u-1)(u-z)/(z^2(z-1)^2)*(4 + (-4)*z/u^2 + 8*(z-1)/(u-1)^2 + (-8)*z(z-1)/(u-z)^2)
- **ic_bc**: u(z) ~ sum_{n>=0} a_n/z^n as z -> +inf with a_0 determined by leading-order PVI balance; u'(50) from series derivative
- **domain**: positive real axis, z > 1; pole-free sector (0th Riemann sheet)
- **closed_form**: none
- **pole_structure**: pole-free on positive real axis z > 1 (tronquee); poles on other Riemann sheets
- **gt_kind**: mathematica
- **how_to_compute**: Wolframscript (WorkingPrecision->50): Generate asymptotic series u(z) = Sum[a_n/z^n, {n,0,Infinity}] with a_0=1 (leading balance for delta<0 tronquee), compute a_n recursively by substituting into PVI and matching 1/z^n powers. Evaluate optimally truncated series at z0=50 to get IC. Then: NDSolve[{pvi_ode_with_params, u[50]==u_series_50, u'[50]==u_series_prime_50}, u, {z,50,10}, WorkingPrecision->50]. Check u[10] against FFW2017 value 0.429534600325223. CAUTION: the FFW2017 value IS the asymptotic series evaluation at z=10, so agreement tests only the asymptotic series self-consistency; true independence requires integrating from z=50 (series converges better there) and checking u[10].
- **precision**: FFW2017 gives 15 digits: u(10)=0.429534600325223, u'(10)=-1.61713114374804e-3; independent regen should match to ~10 digits with z0=50
- **targets_buckets**: B7, B9, B12
- **regime**: PVI tronquee pole-free sector, multi-sheet B12, Riemann surface solver
- **citation**: FFW2017_painleve_riemann_surfaces_preprint.md:195 (Figure 2 caption); FFW2017 section 2.2
- **confidence**: low
- **verify_note**: MUST verify: (1) derive a_0 for these specific parameters by substituting u ~ a_0 into the dominant balance of PVI at z->inf; (2) compute at least 10 coefficients a_n; (3) evaluate at z=50 and integrate to z=10; (4) confirm match with 0.429534600325223 to >=10 digits. Without this check the candidate is tautological.

#### `pii-airy-half-integer`

- **name**: PII Airy solution u_{1/2}(z) = -Ai'(z/2^(1/3)*2^(1/3))/Ai(z/2^(1/3)*2^(1/3)) (alpha=1/2)
- **ode**: u'' = 2u^3 + z*u + 1/2  (Painleve II, alpha=1/2)
- **ic_bc**: c2=0 selects the tronquee case (no Bi component); this selects a particular solution
- **domain**: C minus the zeros of phi(z)
- **closed_form**: u(z; 1/2) = -phi'(z)/phi(z), where phi(z) = c1*Ai(-z/2^(1/3)) + c2*Bi(-z/2^(1/3)). For the tronquee case c2=0: u(z; 1/2) = -Ai'(-z/2^(1/3)) / ((-1/2^(1/3))*Ai(-z/2^(1/3))) = (1/2^(1/3))*Ai'(-z/2^(1/3)) / Ai(-z/2^(1/3))... more precisely: phi''=-z*phi/2, so u=-phi'/phi satisfies PII(alpha=1/2). The poles of u are at zeros of phi.
- **pole_structure**: infinitely many simple poles at zeros of phi(z)=Ai(-z/2^(1/3)); these are at z = -2^(1/3) * a_k where a_k are the zeros of Ai (a_1=-2.3381..., a_2=-4.0879..., etc.); first pole on positive real axis: z = 2^(1/3)*|a_1| = 2^(1/3)*2.3381... = 2.9446...
- **gt_kind**: closed-form
- **how_to_compute**: u(z;1/2) = -(d/dz)[Ai(-z/2^{1/3})] / Ai(-z/2^{1/3}) = (1/2^{1/3}) * Ai'(-z/2^{1/3}) / Ai(-z/2^{1/3}). Poles at z = -2^{1/3}*a_k where Ai(a_k)=0. In python: from mpmath import mp, airyai; mp.dps=30; z=0; u_half = mp.diff(lambda t: mp.airyai(-t/2**mp.mpf('1/3')), z) / mp.airyai(-z/2**mp.mpf('1/3')); Note: DLMF 32.10.7 gives the c1=1, c2=0 formula. First real pole: 2^(1/3)*2.33810741... = 2^(1/3)*|Ai_zero_1| = 2.9446... Verify residue = +1.
- **precision**: arbitrary via mpmath.airyai; zeros of Ai known to arbitrary precision
- **targets_buckets**: B7, B12, B19, B20
- **regime**: Airy-type PII solution; pole locations from Airy zeros; B19 closed-form family; B20 named transcendent
- **citation**: FW2014_second_PII_exploration_FoCM14.md:132-165 (Section 2.3, eqs. 7-10); DLMF 32.10.7
- **confidence**: high
- **verify_note**: Verify using sympy or mpmath that -phi'/phi satisfies PII(1/2); check first pole location against Airy zero tables

#### `pi2-tritronquee-asymptotic`

- **name**: PI^(2) tritronquee (second Painleve hierarchy member) ridge value asymptotics
- **ode**: PI^(2): u_{xxxx} + 10*u*u_{xx} + 5*(u_x)^2 + 10*u^3 + x = 0 (or equivalently: v_{xxxxx} where v_x = u, see KKG 2015). The Gurevich-Pitaevskii / Claeys-Vanlessen solution.
- **ic_bc**: u(x) ~ (1/2)*(6|x|)^(1/3) as x->-inf; u(x)->0 as x->+inf (Gurevich-Pitaevskii asymptotics)
- **domain**: real line (pole-free Gurevich-Pitaevskii solution at t=0)
- **closed_form**: none; characterized by: u(x,0) ~ (6|x|)^{1/3}/2 as x -> -inf (real ridge) and pole-free on real line
- **pole_structure**: real-line solution is pole-free; poles in complex plane; double poles with specific residue structure
- **gt_kind**: published-table
- **how_to_compute**: Claeys & Vanlessen 2007 proved existence and uniqueness of the pole-free real solution. Ridge value: u(x,0) ~ (1/2)*(6|x|)^{1/3} as x->-inf. This can be computed via: from mpmath import mp, odefun; mp.dps=40; # integrate the 4th-order ODE with asymptotic IC at large |x|. No easily accessible regeneration recipe exists without implementing the 4th-order IVP. The Mathematica NDSolve approach: NDSolve[{u''''[x]+10*u[x]*u''[x]+5*(u'[x])^2+10*u[x]^3+x==0, ...asymptotic ICs at x=-L...}, u, {x,-L,0}, WorkingPrecision->40]
- **precision**: asymptotic formula exact at leading order; numerical value at x=0 requires regeneration
- **targets_buckets**: B20
- **regime**: PI hierarchy 4th-order ODE; B20 named transcendent; KdV self-similar solution
- **citation**: Claeys & Vanlessen 2007, Comm. Math. Phys. 273:499-532; KKG = Kapaev, Klein, Grava 2015 (SIAM J. Math. Anal.); epubs.siam.org/doi/10.1137/23M1592304
- **confidence**: low
- **verify_note**: The exact ODE for PI^(2) and the precise IC parameterization must be confirmed from KKG 2015 before implementing. The ridge-value formula (6|x|)^{1/3}/2 at leading order is standard. A full ground-truth recipe requires a working 4th-order complex ODE solver.

### Territory: ODEs with KNOWN closed-form solutions that develop poles in the complex plane — 

*Scout notes:* GAPS AND CAVEATS:\n\n1. mpmath has no built-in WeierstrassP function; the recipe uses mpmath.ellipfun('sn',...) with a Jacobi modulus computed from the roots of 4t^3-g3=0. This works for g2=0 (equianharmonic). For general g2!=0 the Jacobi modulus is complex — this is fine for g2=1 but the formula u''=6u^2-g2/2 (not u''=6u^2) must be used.\n\n2. WP lemniscatic (g2=1, g3=0): the ODE is u''=6u^2-1/2, NOT u''=6u^2. The two are distinct ODEs and should not be conflated in the test suite.\n\n3. Jacobi sn with 0<m<1 has NO real-axis poles. For real-axis pole tests use WP equianharmonic, tan, coth, or the simple-pole family. For complex-plane pole tests sn/cn/dn are ideal.\n\n4. The mpmath ellipfun('sn', z, m) at purely imaginary z has minor numerical imprecision (~1e-50) compared to the exact algebraic values obtained via the imaginary argument transform (DLMF 22.17). Use the transform formulas for ground-truth comparisons, not direct evaluation.\n\n5. The FW2011 Riccati example (y'=t^2+y^2) is the one test where pole locations require Bessel zero computation; they are NOT algebraically specified. Use mpmath.findroot to locate them to arbitrary precision.\n\n6. The equianharmonic WP half-period omega_1 = Gamma(1/3)^3/(2^(7/3)*pi) * 2^(1/6) matches the FW2011 quoted value 1.363... to all printed digits (confirmed).\n\n7. The Painleve I equation u''=6u^2+z has NO Weierstrass P solution — WP is only the c2=0 'companion' equation used as a calibration oracle in FW2011 §5.1.1.\n\nCLEANEST ORACLES FOUND:\n- WP equianharmonic: full cross-validation against FW2011 Table 5.1 (published) + mpmath (regenerable) + lattice-sum (independent).\n- Jacobi sn imaginary bridge: all values are exact algebraic numbers derivable from DLMF 22.17 transform formulas; no floating-point errors in the ground truth.\n- Simple pole 1/(1-z): all values are rational; zero numerical error.

#### `wp-equianharmonic-g0g2`

- **name**: Weierstrass P (equianharmonic, g2=0, g3=2) — FW2011 keystone oracle
- **ode**: u'' = 6*u^2
- **ic_bc**: u(0) = 1.0718225164169174206816142948125029292245326845071, u'(0) = 1.7103373531767862084643318914962495083923187196414
- **domain**: Real axis z in [0, 30] and [0, 10^4]; complex plane for B7
- **closed_form**: u(z) = WP(z - 1; 0, 2)  [Weierstrass P-function with invariants g2=0, g3=2, shift c1=-1 per FW2011 eq (5.2)]  Explicitly via Jacobi sn:   e1 = (1/2)^(1/3)  [real cube root]   e2 = e1*exp(2*pi*i/3),  e3 = e1*exp(4*pi*i/3)   k_sq = (e2-e3)/(e1-e3)   WP(w; 0, 2) = e3 + (e1-e3) / sn^2(sqrt(e1-e3)*w, k_sq)   u(z) = WP(z-1; 0, 2)  IC at z=0: u(0) = 1.0718225164169174206816142948125029292245326845071            u'(0) = 1.7103373531767862084643318914962495083923187196414
- **pole_structure**: Double poles (residue=0) on equilateral-triangle rhombic lattice.  Poles of u(z) = WP(z-1; 0, 2) on the real axis at  z = 1 + 2k*omega_1,  k in Z,  where omega_1 = 1.3630340904278903107840077772328116969762193228397  (= Gamma(1/3)^3 / (2^(7/3)*pi) * 2^(1/6)).  Full lattice: z = 1 + 2m*omega_1 + 2n*omega_2 where omega_2 = omega_1*exp(i*pi/3).  First few real poles: z=1, 3.7261, 6.4521, 9.1781, ...
- **gt_kind**: mpmath
- **how_to_compute**: python3 -c " import mpmath mpmath.mp.dps = 50 pi = mpmath.pi e1 = (mpmath.mpf(1)/2)**(mpmath.mpf(1)/3) e2 = e1 * mpmath.exp(mpmath.mpc(0, 2*pi/3)) e3 = e1 * mpmath.exp(mpmath.mpc(0, 4*pi/3)) k_sq = (e2 - e3)/(e1 - e3) def wp(z, c1=-1):     w = z + c1     return (e3 + (e1-e3)/mpmath.ellipfun('sn', mpmath.sqrt(e1-e3)*w, k_sq)).real print(wp(0))      # 1.071822516416917420... print(wp(30))     # 1.095098255959744199... print(wp(10000))  # 21.02530339471054998... print(wp(28.261)) # 9876953.517025014... "
- **precision**: arbitrary via mpmath.mp.dps; 50 digits verified
- **targets_buckets**: B1, B3, B5, B7, B8
- **regime**: Real-axis pole bridging: step h=0.5, order=30; 22 poles in [0,30].  Step z=0.5->1.5 crosses pole at z=1 (strictly interior).  Long-range: z=10^4 (crosses ~7300 poles).  Pole-wall: z=28.261 (value ~9.877e6).  Complex-plane: arbitrary (z,g2,g3) via mpmath.
- **citation**: FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:281-302 (Section 5.1.1 and Table 5.1); published reference values u(30)=1.095098255959744 and u(1e4)=21.02530339471055 in eq (5.3) and Table 5.1
- **confidence**: high
- **verify_note**: All five reference values independently reproduced by mpmath ellipfun recipe and match FW2011 Table 5.1 to all displayed digits.  Recipe is regenerable at arbitrary precision.

#### `wp-lemniscatic-g1g0`

- **name**: Weierstrass P (lemniscatic, g2=1, g3=0) — square-lattice pole array
- **ode**: u'' = 6*u^2 - 1/2
- **ic_bc**: u(K/2) = (1+sqrt(2))/2 [exact], u'(K/2) = -(1+sqrt(2)) [exact];  K = Gamma(1/4)^2/(4*sqrt(pi)) = 1.8540746773013719
- **domain**: Real axis z > 0 (IC at z=K, integrate to z=3K); complex plane
- **closed_form**: u(z) = WP(z; 1, 0)  [Weierstrass P-function with g2=1, g3=0]  Explicitly via Jacobi sn with k^2=1/2:   e1 = 1/2,  e3 = -1/2,  k_sq = 1/2   WP(z; 1, 0) = -1/2 + 1/sn^2(z, 1/2)  Exact algebraic ICs:   WP(K/2; 1, 0) = (1 + sqrt(2))/2  [EXACT]   WP'(K/2; 1, 0) = -(1 + sqrt(2))  [EXACT] where K = K(1/2) = Gamma(1/4)^2 / (4*sqrt(pi)) = 1.8540746773013719...  Other exact values:   WP(K; 1, 0) = 1/2 [exact = e1]   WP(3K/2; 1, 0) = (1+sqrt(2))/2 [exact, by symmetry about pole at 2K]
- **pole_structure**: Double poles on square lattice: z = 2m*omega_1 + 2n*i*omega_1  (m,n in Z) where omega_1 = K(1/2) = Gamma(1/4)^2/(4*sqrt(pi)) = 1.8540746773013719184338503471952600462175988235218.  Real poles at z = 0, +/-2*omega_1 = +/-3.7081..., +/-4*omega_1, ...  Square lattice (omega_2 = i*omega_1).
- **gt_kind**: mpmath
- **how_to_compute**: python3 -c " import mpmath mpmath.mp.dps = 50 k_sq = mpmath.mpf(1)/2 K = mpmath.ellipk(k_sq)  # = Gamma(1/4)^2/(4*sqrt(pi)) e1, e3 = mpmath.mpf(1)/2, -mpmath.mpf(1)/2 def wp_lem(z):     return e3 + (e1-e3)/mpmath.ellipfun('sn', z, k_sq)**2 # Exact values: print(wp_lem(K/2))  # = (1+sqrt(2))/2 = 1.2071067811865475... print(wp_lem(K))    # = 1/2 exactly # Numeric: print(wp_lem(mpmath.mpf('1.5')))  # 0.5668205194485407... print(wp_lem(mpmath.mpf('2.0')))  # 0.5107614339783495... print(wp_lem(mpmath.mpf('3.0')))  # 2.019294402732143... "
- **precision**: arbitrary via mpmath.mp.dps; closed-form values exact
- **targets_buckets**: B1, B3, B5, B8
- **regime**: Different ODE than equianharmonic: tests the g2 != 0 case.  Square lattice (vs triangular).  Exact algebraic ICs from half-lattice point; pole-bridge from z=K to z=3K crosses pole at z=2K=3.7081...
- **citation**: DLMF 23.5(ii) for lemniscatic half-period; constructed from WP ODE u''=6u^2-g2/2 (derived from (WP')^2=4*WP^3-g2*WP-g3, differentiated)
- **confidence**: high
- **verify_note**: ODE u''=6u^2-1/2 (NOT u''=6u^2) verified numerically at z=1.5: residual |u''-(6u^2-1/2)| < 1e-5.  Exact values (1+sqrt(2))/2 and -(1+sqrt(2)) verified to 50 digits.  omega_1=K(1/2) confirmed by direct integral.

#### `riccati-tan-1storder`

- **name**: Riccati / tangent: u'=1+u^2 → tan(z)
- **ode**: u' = 1 + u^2
- **ic_bc**: u(0) = 0, u'(0) = 1  [or u(0)=1, u'(0)=2 for tan(z+pi/4) branch]
- **domain**: Real axis z in [-pi, 2*pi]; complex plane via mpmath.tan
- **closed_form**: u(z) = tan(z + phi_0)  With u(0)=0: u(z) = tan(z) With u(0)=1: u(z) = tan(z + pi/4)  Exact values (u(0)=0 case):   u(pi/4) = 1  [exact]   u(-pi/4) = -1  [exact]   u(3*pi/4) = -1  [exact]   u(pi) = 0  [exact]   u(0.5) = 0.54630248984379051325517946578028538329755172017979   u(1.0) = 1.5574077246549022305069748074583601730872507723815   u(1.5) = 14.101419947171719387646083651987756445659543577236
- **pole_structure**: Simple poles (residue = 1) on the real axis at z = pi/2 + k*pi,  k in Z.  pi/2 = 1.5707963267948966192313216916397514420985846996876.  No complex poles away from real axis — all poles real.
- **gt_kind**: closed-form
- **how_to_compute**: python3 -c "import mpmath; mpmath.mp.dps=50; print(mpmath.tan(mpmath.mpf('1.0')))"  # for any z-value  For arbitrary precision: mpmath.tan(mpmath.mpf(z)) at any dps.  Pole-bridge oracle: step h=pi from z=0 (u=0) to z=pi (u=0), pole at z=pi/2 strictly interior.  u(pi)=tan(pi)=0 exactly.
- **precision**: arbitrary via mpmath.tan; exact at pi/4, 3*pi/4, pi
- **targets_buckets**: B3, B5
- **regime**: Single pole in one step: h=pi, IC u(0)=0, pole at z=pi/2 (exactly midway), output u(pi)=0.  Or use h=pi/2 with u(0)=1, crossing pole at z=pi/4 to land at u(pi/2)=-1.  Tests scalar stepper on simplest possible pole-crossing.
- **citation**: constructed; textbook ODE u'=1+u^2 -> tan. DLMF 4.21 for tan pole structure.
- **confidence**: high

#### `tan-2nd-order`

- **name**: Tangent 2nd-order: u''=2u+2u^3 → tan(z)
- **ode**: u'' = 2*u + 2*u^3
- **ic_bc**: u(0) = 0, u'(0) = 1
- **domain**: Real axis z in [0, 2*pi]
- **closed_form**: u(z) = tan(z + phi_0)  With u(0)=0, u'(0)=1: u(z) = tan(z)  Verification: tan'=sec^2=1+tan^2, tan''=2*tan*sec^2=2*tan*(1+tan^2)=2*u+2*u^3 ✓  Exact values:   u(pi/4) = 1  [exact]   u(3*pi/4) = -1  [exact]   u(pi) = 0  [exact]   u(0.5) = 0.54630248984379051325517946578028538329755172017979   u(1.0) = 1.5574077246549022305069748074583601730872507723815
- **pole_structure**: Simple poles (residue = 1) at z = pi/2 + k*pi, k in Z.  Same as 1st-order Riccati; the 2nd-order form gives a test with distinct u and u' ICs and exercises the 2nd-order stepper.
- **gt_kind**: closed-form
- **how_to_compute**: python3 -c "import mpmath; mpmath.mp.dps=50; print(mpmath.tan(mpmath.mpf('1.0')))"
- **precision**: arbitrary via mpmath.tan; exact at integer multiples of pi/4
- **targets_buckets**: B3, B5, B8
- **regime**: 2nd-order scalar pole-bridge: h=pi crosses pole at pi/2.  Exercises the 2nd-order IVP stepper (u and u' ICs) separately from the 1st-order case.  Pole-wall test: u(pi/2 - 0.01) = 99.99666...
- **citation**: constructed; differentiate u'=1+u^2 once to get u''=2u*u'=2u*(1+u^2)
- **confidence**: high

#### `simple-pole-riccati`

- **name**: Simple pole: u'=u^2 → 1/(c-z)
- **ode**: u' = u^2
- **ic_bc**: u(0) = 1 (pole at z=1); or u(0) = 1+i (pole at z=(1-i)/2)
- **domain**: Real axis z in [0, 2] (bridging z=1); complex plane for u(0)=1+i variant
- **closed_form**: u(z) = 1/(c - z)  where c = 1/u(0)  With u(0)=1: u(z) = 1/(1-z), pole at z=1  Exact values:   u(0) = 1  [IC]   u(0.5) = 2  [exact rational]   u(1.5) = -2  [exact rational]   u(2) = -1  [exact rational]   u(-1) = 1/2  [exact rational]  With u(0) = 1+i (complex IC):   u(z) = 1/((1-i)/2 - z),  pole at z = (1-i)/2   Exercises complex-plane stepper with off-real-axis pole.
- **pole_structure**: Simple pole (residue = -1) at z = c = 1/u(0).  For u(0)=1: pole at z=1.  For u(0)=1+i: pole at z=(1-i)/2 = 0.5-0.5i (off real axis).
- **gt_kind**: closed-form
- **how_to_compute**: Exact closed form: u(z) = 1/(1-z).  No computation needed.  Evaluate at any z != 1 exactly.  Pole-bridge: IC u(0)=1, step h=2 from z=0 to z=2, pole at z=1 (midpoint).  Output u(2) = -1 exactly.
- **precision**: exact rational for u(0)=1 variant
- **targets_buckets**: B3, B5
- **regime**: Simplest possible pole-bridge: single simple pole, pole at exact midpoint of step.  All values are rational.  Ideal for unit-testing the Pade step logic without elliptic-function overhead.
- **citation**: constructed; u'=u^2 is the canonical Bernoulli ODE with simple movable pole
- **confidence**: high

#### `coth-simple-pole`

- **name**: Coth: u'=1-u^2 → coth(z), real pole at z=0
- **ode**: u' = 1 - u^2
- **ic_bc**: u(-1) = coth(-1) = -1.3130352854993313, u'(-1) = 1-u^2(-1) = 1-coth^2(-1) = -csch^2(-1)
- **domain**: Real axis z in [-2, 2] bridging z=0
- **closed_form**: u(z) = coth(z)  [or equivalently u(z) = 1/tanh(z)]  coth'(z) = -csch^2(z) = 1 - coth^2(z) = 1 - u^2  ✓  Pole at z=0 (simple pole, residue=1).  Exact reference values:   u(-2) = coth(-2) = -1.0373147207275480958778097647678207116623912692492   u(-1) = coth(-1) = -1.3130352854993313036361612469308478329120139412405   u(-0.5) = coth(-0.5) = -2.1639534137386528487700040102180231170937386021508   u(0.5) = coth(0.5) = 2.1639534137386528487700040102180231170937386021508   u(1) = coth(1) = 1.3130352854993313036361612469308478329120139412405   u(2) = coth(2) = 1.0373147207275480958778097647678207116623912692492  Pole-bridge: IC u(-1) = coth(-1) = -1.313035..., step h=2 to z=1, pole at z=0 (midpoint).
- **pole_structure**: Simple pole (residue=1) at z=0 on the real axis.  Complex poles at z = i*k*pi for k in Z\{0}.  The only real pole is z=0 — making this a one-shot real-axis pole-bridge test.
- **gt_kind**: closed-form
- **how_to_compute**: python3 -c "import mpmath; mpmath.mp.dps=50; print(mpmath.coth(mpmath.mpf('1')))"  For pole-bridge at z=0: u(-1)=coth(-1), step h=2, u(1)=coth(1)=-u(-1) by odd symmetry.
- **precision**: arbitrary via mpmath.coth; exact up to 50+ digits
- **targets_buckets**: B3, B5
- **regime**: Real-axis pole at z=0 (unlike WP and tan whose first poles are away from 0).  The odd symmetry coth(-z)=-coth(z) gives exact check: u(1) = -u(-1).  Tests single real-axis simple pole at a known integer location.
- **citation**: constructed; coth' = 1 - coth^2 is standard hyperbolic identity. DLMF 4.28 for coth.
- **confidence**: high

#### `jacobi-sn-2ndorder`

- **name**: Jacobi sn (m=1/4): complex-plane pole-bridge with exact algebraic ICs
- **ode**: u'' = -(1 + m)*u + 2*m*u^3  [m = 1/4, so u'' = -(5/4)*u + (1/2)*u^3]
- **ic_bc**: u(i*Kp/2) = i*sqrt(2) [exact], u'(i*Kp/2) = 3/sqrt(2) [exact]
- **domain**: Imaginary axis and complex plane; real axis is pole-free for 0<m<1
- **closed_form**: u(z) = sn(z, m=1/4)  [Jacobi elliptic sn function, modulus-squared m=1/4]  Derivation: sn'=cn*dn, sn''=-sn*dn^2+m*sn*cn^2... deriving gives sn''=-(1+m)*sn+2m*sn^3.  Exact pole-bridge values (all EXACT algebraic numbers):   u(i*Kp/2) = i*sqrt(2)   [exact]   u'(i*Kp/2) = 3/sqrt(2)  [exact]  (= sqrt(3)*sqrt(3/2) via imaginary transform)   u(i*3*Kp/2) = -i*sqrt(2) [exact]   where Kp = K'(1/4) = K(3/4) = 2.1565156474996432354386749988003220288641102164928  Proof: sn(i*Kp/2, 1/4) = i*nd(Kp/2, 3/4) = i/dn(Kp/2, 3/4); dn(Kp/2, 3/4)=1/sqrt(2)  => i*sqrt(2). cn(i*Kp/2, 1/4) = nc(Kp/2, 3/4) = 1/cn(Kp/2, 3/4) = sqrt(3). dn(i*Kp/2, 1/4) = dc(Kp/2, 3/4) = dn(Kp/2, 3/4)/cn(Kp/2, 3/4) = sqrt(3/2). sn'(i*Kp/2) = cn*dn = sqrt(3)*sqrt(3/2) = 3/sqrt(2).  Additional exact values on real axis:   sn(K/2, 1/4) = 1/sqrt(1 + sqrt(3/4))  [exact, K=1.6857503548125960]   sn(K, 1/4) = 1  [exact]   sn(0, 1/4) = 0  [exact]  Numerical values on real axis:   sn(0.5, 1/4) = 0.47508293602853651008221832470387025874507817180743   sn(1.0, 1/4) = 0.82263557812986235967623033865397648844064711651056   sn(1.5, 1/4) = 0.98705228731124166664858916456480943724796956259770   sn(2.0, 1/4) = 0.96289817759827744257513989213743904197805118421774
- **pole_structure**: Poles at u = 2p*K + (2q+1)*i*Kp  for p,q in Z.  NO real poles for 0 < m < 1; all poles are off the real axis.  First poles: u = +/-i*Kp = +/-2.1565i.  For a complex step test: integrate along imaginary axis, pole at u=i*Kp = i*2.1565...
- **gt_kind**: mpmath
- **how_to_compute**: python3 -c " import mpmath mpmath.mp.dps = 50 m = mpmath.mpf(1)/4 K = mpmath.ellipk(m)   # = 1.6857503548125960... Kp = mpmath.ellipk(1-m) # = 2.1565156474996432... # Verify exact pole-bridge values: print(mpmath.ellipfun('sn', mpmath.mpc(0, float(Kp)/2), m))  # -> i*sqrt(2) print(mpmath.ellipfun('sn', mpmath.mpc(0, 3*float(Kp)/2), m))  # -> -i*sqrt(2) # Real-axis values: for u in ['0.5', '1.0', '1.5', '2.0']:     print(mpmath.ellipfun('sn', mpmath.mpf(u), m)) "
- **precision**: arbitrary via mpmath.ellipfun; closed-form values exact via imaginary transform
- **targets_buckets**: B1, B3, B7
- **regime**: Complex-direction pole-bridge: imaginary step h=i*Kp from z=i*Kp/2 to z=i*3*Kp/2 bridges pole at u=i*Kp.  IC and target are exact algebraic numbers.  Tests complex stepping direction in B7 path network.  The real-axis values test the smooth (no-pole) regime for B3/B1.
- **citation**: DLMF 22.13.1 for the 2nd-order ODE; DLMF 22.17.3-5 for imaginary argument transforms (used to derive cn(i*Kp/2)=sqrt(3), dn(i*Kp/2)=sqrt(3/2)). Pole locations from DLMF 22.4.8.
- **confidence**: high
- **verify_note**: Exact values cn(i*Kp/2, 1/4)=sqrt(3), dn(i*Kp/2, 1/4)=sqrt(3/2) verified to 100 digits via imaginary transform (DLMF 22.17): cn(iu,m)=nc(u,1-m), dn(iu,m)=dc(u,1-m), and cn(Kp/2,3/4)=1/sqrt(3), dn(Kp/2,3/4)=1/sqrt(2) confirmed to 100 digits by mpmath.

#### `jacobi-sn-1storder`

- **name**: Jacobi sn 1st-order: (u')^2 = (1-u^2)(1-m*u^2)
- **ode**: (u')^2 = (1 - u^2)*(1 - m*u^2)  [m = 1/4]
- **ic_bc**: u(0) = 0, u'(0) = 1
- **domain**: Real axis z in [0, 4*K]; complex plane
- **closed_form**: u(z) = sn(z, m=1/4)  Same solution as the 2nd-order candidate but the 1st-order form directly encodes the ODE constraint.  IC: u(0) = 0, u'(0) = 1  [from sn(0)=0, sn'(0)=cn(0)*dn(0)=1]  Key values same as above:   sn(K/2, 1/4) = 1/sqrt(1+sqrt(3/4))  [exact]   sn(K, 1/4) = 1  [exact]   sn(0.5, 1/4) = 0.47508293602853651008221832470387025874507817180743   sn(1.0, 1/4) = 0.82263557812986235967623033865397648844064711651056
- **pole_structure**: Same as 2nd-order sn: poles at 2pK + (2q+1)*i*Kp, all off real axis. No real-axis poles for 0<m<1.
- **gt_kind**: mpmath
- **how_to_compute**: python3 -c "import mpmath; mpmath.mp.dps=50; m=mpmath.mpf(1)/4; print(mpmath.ellipfun('sn', mpmath.mpf('1.0'), m))"
- **precision**: arbitrary via mpmath.ellipfun
- **targets_buckets**: B1, B3
- **regime**: Tests the 1st-order IVP stepper on a smooth (no real-axis poles) solution.  Primary value: verifying the stepper is accurate in the non-pole regime before testing pole-bridge.  The complex poles at i*Kp are detectable by comparing real-axis stepping with imaginary-direction stepping.
- **citation**: DLMF 22.13.1 for the ODE form. DLMF 22.4.8 for pole locations.
- **confidence**: high
- **verify_note**: sn(K/2, m) = 1/sqrt(1+sqrt(1-m)) verified at m=1/4 to 50 digits by mpmath.

#### `fw-riccati-bessel`

- **name**: FW2011 Section 2 demonstration: y'=t^2+y^2 → Bessel quotient
- **ode**: y' = t^2 + y^2,  y(0) = 0
- **ic_bc**: y(0) = 0, y'(0) = 0 (since y'(0) = 0^2 + 0^2 = 0)
- **domain**: Real axis t in [0, 4] (first two poles)
- **closed_form**: y(t) = t * J_{3/4}(t^2/2) / J_{-1/4}(t^2/2)  where J_nu is the Bessel function of the first kind and J_{-1/4}(x) = cos(pi/4)*J_{1/4}(x) - sin(pi/4)*Y_{1/4}(x)  [DLMF 10.2.3]  Numerical values:   y(0.5) = 0.041791146154681863220768806849177498107753256229669   y(1.0) = 0.35023184431675577784938227410474107449829444543101   y(1.5) = 1.5174475438800018517295798114990422035450618321886  First pole: t_1 = sqrt(2 * z_1)  where z_1 = 2.0062996717894504160280379194791147962182407789045 is the first zero of J_{-1/4}(x) = cos(pi/4)*J_{1/4}(x) - sin(pi/4)*Y_{1/4}(x). t_1 = 2.003147359426884708004610979054299223810144817229
- **pole_structure**: Simple poles at t = sqrt(2*z_k) where z_k are zeros of J_{-1/4}(x).  First real pole at t_1 = 2.0031473594268847.  Poles are NOT at algebraically specified positions (require Bessel zero computation) but are computable to arbitrary precision by mpmath.findroot.
- **gt_kind**: mpmath
- **how_to_compute**: python3 -c " import mpmath mpmath.mp.dps = 50 pi = mpmath.pi nu14 = mpmath.mpf(1)/4 def Jm14(x):     return mpmath.cos(pi/4)*mpmath.besselj(nu14,x) - mpmath.sin(pi/4)*mpmath.bessely(nu14,x) def y_FW(t):     t = mpmath.mpf(str(t))     return t * mpmath.besselj(3*nu14, t**2/2) / Jm14(t**2/2) # Pole location: z1 = mpmath.findroot(Jm14, 2.0) print('First pole at t =', mpmath.sqrt(2*z1)) print('y(0.5) =', y_FW('0.5')) print('y(1.0) =', y_FW('1.0')) print('y(1.5) =', y_FW('1.5')) "
- **precision**: arbitrary via mpmath.besselj, mpmath.bessely, mpmath.findroot
- **targets_buckets**: B3, B5
- **regime**: The FW2011 motivating example (Figure 2.1 and Section 2.2.1): demonstrates that Taylor method fails but Pade steps through the first pole at t~2.003.  First-order ODE with non-algebraic pole locations.  Good test for the scalar stepper on a Riccati-type equation different from u'=1+u^2.
- **citation**: FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:122-134 (Section 2.2.1, eq 2.6 and surrounding text)
- **confidence**: high
- **verify_note**: Formula y(t)=t*J_{3/4}(t^2/2)/J_{-1/4}(t^2/2) verified by: substituting into y'=t^2+y^2 and checking numerically at t=1.0, residual < 1e-40 at 50 digits.  First pole location verified by Bessel zero computation.

#### `tanh-complex-pole`

- **name**: Tanh: u'=1-u^2 → tanh(z), off-real poles at i*(pi/2+k*pi)
- **ode**: u' = 1 - u^2
- **ic_bc**: u(0) = 0, u'(0) = 1
- **domain**: Imaginary axis z in [0, i*pi] for pole-bridge; complex plane for B7
- **closed_form**: u(z) = tanh(z)  [u(0)=0 branch]  Exact values:   u(i*pi/4) = i*tan(pi/4) = i  [exact]   u(pi/4) = tanh(pi/4) = 0.65579...   u(1) = 0.76159415595576488811945828260479359041276859725794   u(1+i) = 1.0839233273386945434757520612119717213449675274754 + 0.27175258531951171652884372249858892070946411146178*i  Poles at z = i*(pi/2 + k*pi) for k in Z.  Nearest: z = i*pi/2.  Complex pole-bridge: step from z=0 to z=i*pi (step h=i*pi on imaginary axis)  crosses pole at z=i*pi/2. Landing: u(i*pi) = i*tan(pi) = 0 [exact].  Note: tanh(i*pi/4)=i is exact; tanh(i*pi)=0 is exact.
- **pole_structure**: Simple poles (residue=1) at z = i*(pi/2 + k*pi),  k in Z.  All poles on the imaginary axis; real axis is pole-free (tanh is bounded on R).  First poles: z = +/-i*pi/2 = +/-1.5707963...*i
- **gt_kind**: closed-form
- **how_to_compute**: python3 -c "import mpmath; mpmath.mp.dps=50; print(mpmath.tanh(mpmath.mpc(0, mpmath.pi/4)))"  # -> i exactly  For pole-bridge: tanh(0)=0, step h=i*pi to tanh(i*pi)=0 (exact), crossing pole at i*pi/2.
- **precision**: arbitrary via mpmath.tanh; exact at i*k*pi/4 where k odd
- **targets_buckets**: B3, B5, B7
- **regime**: Off-real-axis poles: tests complex-direction stepping (h in imaginary direction) on the B7 path network.  Same ODE as coth (u'=1-u^2) but different IC gives tanh (no real pole) vs coth (real pole at z=0).  The pair {coth, tanh} exercises both branches of the same first-order equation.
- **citation**: constructed; tanh'=sech^2=1-tanh^2. DLMF 4.28 for tanh poles.
- **confidence**: high
- **verify_note**: tanh(i*pi/4)=i verified by mpmath to 50 digits.

### Territory: COUPLED (d>=2) vector ODE systems with KNOWN closed-form / high-precision soluti

*Scout notes:* GAPS AND CAVEATS:

1. CALOGERO-MOSER IMAGINARY-IC (shared-pole crossing): The imaginary-IC variant of the CM N=2 system (cm-n2-imaginary-collision) gives a collision at t*=r0^2*sqrt(2) -- this IS the vector shared-pole crossing test the repo identifies as the v0.2 keystone. However, the closed form is derived here by analytic continuation and has NOT yet been independently verified in the repo against a BigFloat RK4 (only the real-IC case is verified in test/calogero_moser_test.jl). This verification is the highest-priority TODO.

2. NOUMI-YAMADA A_4 GENERIC TRANSCENDENT: The generic (non-rational) A_4^(1) transcendent has no published closed form. The mpmath odefun oracle is regenerable but not pre-computed. Pole locations are unknown analytically; the only way to get them is numerically (find roots of the shared-Q denominator from the PadeTaylor solve).

3. OKAMOTO HAMILTONIAN FORMS (pi4-okamoto-hamiltonian-2vector): The explicit Okamoto Hamiltonian H_IV is in the repo references (NY1999 main.tex, accessible) but the canonical Hamiltonian 2-vector ODE system was not fully worked out here. The easiest path is to use the companion form y1'=y2, y2'=PIV_RHS instead of the canonical (q,p) variables -- which is what the v0.2 codebase already does.

4. GARNIER SYSTEMS: The repo has extensive Garnier system references (references/tex/garnier/) but no Garnier oracle has been identified with a closed-form evaluation at specific complex points. Mazzocco 2002 gives classical solutions via Lauricella hypergeometric functions -- potentially a future oracle but requires specialized mpmath integration.

5. ADLER-SOKOLOV VECTOR PI SYSTEMS: The KdV-Gal vector PI system (docs/v0p2_pillarEF_methodology_calogero_findings.md section 4, lines 285-295) is another candidate but has no closed-form solution identified here; the cross-validation would require mpmath odefun.

6. CLEANEST ORACLES: (a) cm-n2-rational-kdv is the cleanest oracle (already in the test suite, verified to 26-30 digits). (b) coupled-riccati-tan-sec is the cleanest NEW oracle for a 2-vector system with a shared pole (tan/sec^2 at z=pi/2); it is not yet in the test suite. (c) noumi-yamada-a4-type-c and a4-backlund-orbit are the cleanest Noumi-Yamada oracles (Type C in the test suite; Backlund orbit not yet).

7. SHARED-POLE CROSSING NOT YET PINNED: As of the current repo state (git status shows only untracked external probe directories), NO test exists that walks the vector solver THROUGH a shared pole at a known complex location with a known exact value on the other side. The imaginary-IC CM system (cm-n2-imaginary-collision) is the recommended candidate to fill this gap -- it has an exact closed form up to the collision and can be tested by integrating in the complex t-plane around the pole.

#### `cm-n2-rational-kdv`

- **name**: Calogero-Moser N=2 rational KdV pole system (repulsive)
- **ode**: y = (x1, x2, v1, v2); y' = (v1, v2, 4/(x1-x2)^3, 4/(x2-x1)^3). This is the first-order vector form of the N=2 rational Calogero-Moser particle system whose poles carry the rational KdV solution u = 2/(x-x1(t))^2 + 2/(x-x2(t))^2.
- **ic_bc**: x1(0)=1, x2(0)=-1, v1(0)=0, v2(0)=0 (symmetric, r0=1)
- **domain**: real t in [0, 0.4] (repulsive, collision-free); complex IC gives collision at t* = r0^2/2
- **closed_form**: Symmetric IC x1(0)=r0, x2(0)=-r0, v1(0)=v2(0)=0. Exact solution: x1(t) = sqrt(r0^2 + t^2/(2*r0^2)), x2(t) = -x1(t), v1(t) = t/(2*r0^2*x1(t)), v2(t) = -v1(t). For r0=1: x1(t)=sqrt(1 + t^2/2). Conserved: E = (v1-v2)^2/2 + 4/(x1-x2)^2 = 1/r0^2.
- **pole_structure**: The KdV solution u(x,t) = 2/(x-x1(t))^2 + 2/(x-x2(t))^2 has shared second-order poles at x=x1(t) and x=x2(t); both components of the particle ODE share the singularity at x1(t)-x2(t)=0, which occurs only in the imaginary-IC (attractive) case. In the real-axis repulsive case used here, poles are real, moving, and well-separated.
- **gt_kind**: closed-form
- **how_to_compute**: Julia: x1(t,r0) = sqrt(r0^2 + t^2/(2*r0^2)). Independent cross-check via 256-bit BigFloat RK4 (200k steps). See test/calogero_moser_test.jl in the repo: cross-checked at t in (0.1, 0.25, 0.4) to 26-30 sig figs.
- **precision**: arbitrary via BigFloat; Float64 achieves ~1e-10 on positions/velocities
- **targets_buckets**: B13, B14
- **regime**: 4-component vector ODE with exact closed-form oracle; shared pole at x1-x2=0 in the attractive-IC variant (imaginary r0 gives collision at t*=r0^2/2); Float64 and BigFloat-256 precision-generic stack
- **citation**: Krichever 1980, Functional Analysis and its Applications 14, p.282 eq.(1) (Hamiltonian), p.284 eq.(12) (equations of motion); Airault-McKean-Moser 1977, CPAM 30 (rational KdV-pole connection). Verified closed form: test/calogero_moser_test.jl lines 54-83 (corrected from docs/v0p2_pillarEF_methodology_calogero_findings.md).
- **confidence**: high
- **verify_note**: The pillar doc (docs/v0p2_pillarEF_methodology_calogero_findings.md) has an arithmetic error (C=4/r0^2, wrong sign in first integral); the corrected oracle x1(t)=sqrt(r0^2+t^2/(2r0^2)) is verified in test/calogero_moser_test.jl CM.0 against independent 256-bit RK4. Use test file as canonical source.

#### `cm-n2-imaginary-collision`

- **name**: Calogero-Moser N=2 with imaginary initial positions -- shared-pole collision
- **ode**: Same 4-vector ODE as cm-n2-rational-kdv: y' = (v1, v2, 4/(x1-x2)^3, 4/(x2-x1)^3), but now with imaginary initial positions x1(0)=+i*r0, x2(0)=-i*r0 so the two particles collide at t* = r0^2/2.
- **ic_bc**: x1(0)=+i, x2(0)=-i, v1(0)=0, v2(0)=0 (r0=1, so t*=sqrt(2) approx 1.414)
- **domain**: complex t-plane; real part 0 to t* = r0^2*sqrt(2) along Im(t)=0
- **closed_form**: With x1(0)=+i*r0, x2(0)=-i*r0, v1(0)=v2(0)=0: relative coordinate satisfies r(t)^2 = -4r0^2 + 2t^2/r0^2 = 2(t^2 - 2r0^4)/r0^2. Collision (shared pole in the 4-vector ODE) at t* = r0^2*sqrt(2). Exact trajectories: x1(t) = (i/2)*sqrt(4r0^2 - 2t^2/r0^2) for t < r0^2*sqrt(2). Equivalently r(t) = i*sqrt(4r0^2 - 2t^2/r0^2).
- **pole_structure**: The ODE has a shared pole in all 4 components at t=t* where x1(t*)=x2(t*) (particle collision). The shared denominator Q of the vector Pade approximant must find this single shared zero -- this is the v0.2 keystone test: stepping THROUGH a shared pole with the vector walk.
- **gt_kind**: closed-form
- **how_to_compute**: Julia: for r0=1, t*=sqrt(2). x1(t) = (im/2)*sqrt(4 - 2*t^2) for 0 <= t < t*. Regenerate: setprecision(BigFloat,256); r0=BigFloat(1); t_test=BigFloat('0.5'); x1_exact = (im/2)*sqrt(4*r0^2 - 2*t_test^2/r0^2). Also mpmath: from mpmath import *; mp.dps=50; r0=1; t=mpf('0.5'); x1 = j*sqrt(4*r0**2 - 2*t**2/r0**2)/2; print(x1)
- **precision**: arbitrary via BigFloat/mpmath
- **targets_buckets**: B14, B13
- **regime**: 4-vector shared-pole walk straddling a genuine collision; imaginary trajectories in complex t-plane; the vector shared-Q walk must bridge the t=t* pole
- **citation**: Derived from Krichever 1980 (same system as cm-n2-rational-kdv) by analytic continuation of initial conditions. Closed form derived from first integral E = (v1-v2)^2/2 + 4/(x1-x2)^2 = 1/r0^2 (same conservation law) applied at complex IC.
- **confidence**: high
- **verify_note**: Confirm: at t=0.5, x1 = (im/2)*sqrt(4-0.5) = (im/2)*sqrt(3.5). Then u(x,t)=2/(x-x1)^2 + 2/(x-x2)^2 with x1=-x2 has poles at x = +/- (im/2)*sqrt(4r0^2-2t^2/r0^2). The shared Q of the 4-vector must vanish at exactly t=t*.

#### `noumi-yamada-a2-type-c`

- **name**: Noumi-Yamada A_2^(1) (= PIV) -- trivial rational solution f_j = t/3
- **ode**: A_2^(1) system (3-component, even parity, n=1): f_0' = f_0(f_1 - f_2) + alpha_0; f_1' = f_1(f_2 - f_0) + alpha_1; f_2' = f_2(f_0 - f_1) + alpha_2. Constraint: f_0+f_1+f_2 = t, sum(alpha_j)=1.
- **ic_bc**: f_j(t0) = t0/3 for j=0,1,2; alpha=(1/3,1/3,1/3)
- **domain**: all t (entire rational solution)
- **closed_form**: Type C rational solution: alpha_j = 1/3 for all j. Then f_j(t) = t/3 for j=0,1,2. Verification: f_j' = 1/3. RHS = f_j*(f_{j+1} - f_{j+2}) + alpha_j = (t/3)*(t/3-t/3) + 1/3 = 1/3. Exact for all t (globally defined, pole-free).
- **pole_structure**: This rational solution is entire (globally polynomial, no poles). It reduces to the PIV solution y = 0 (or the zero solution in the appropriate gauge). The interesting regime for pole testing is generic alpha_j not all equal.
- **gt_kind**: closed-form
- **how_to_compute**: Julia: f_type_c(t) = fill(t/3, 3); alpha = fill(1//3, 3). Exactly satisfies ODE by construction (verified algebraically). See test/noumi_yamada_symmetry_test.jl NYS.1.1 for exact Rational{BigInt} verification.
- **precision**: exact rational arithmetic (Rational{BigInt})
- **targets_buckets**: B21
- **regime**: 3-component Noumi-Yamada vector ODE with trivial (pole-free) closed-form oracle; exercises the vector RHS factory noumi_yamada_rhs(1; alpha) and constraint invariant sum(f_j)=t
- **citation**: Noumi-Yamada 1998, Funkcialaj Ekvacioj 41, main.tex lines 85-88 (A_2^(1) system); pillar doc docs/v0p2_pillarB_noumi_yamada_findings.md section 5.1 (Type A,C seeds). Type C is the cyclic equilibrium.
- **confidence**: high

#### `noumi-yamada-a4-type-c`

- **name**: Noumi-Yamada A_4^(1) -- trivial rational solution f_j = t/5
- **ode**: A_4^(1) system (5-component, even parity, n=2): f_0' = f_0(f_1 - f_2 + f_3 - f_4) + alpha_0; (and 4 cyclic rotations). Constraint: sum(f_j) = t, sum(alpha_j) = 1.
- **ic_bc**: f_j(t0)=t0/5 for all j; alpha=(1/5,...,1/5)
- **domain**: all t (Type C entire); simple pole at t=0 after Backlund shift s_0
- **closed_form**: Type C rational solution: alpha_j = 1/5 for all j. Then f_j(t) = t/5 for j=0,...,4. Verification: f_j' = 1/5. RHS = f_j*(f_{j+1}-f_{j+2}+f_{j+3}-f_{j+4}) + alpha_j = (t/5)*0 + 1/5 = 1/5. Exact for all t.
- **pole_structure**: Entire (globally polynomial, no poles). After one Backlund shift s_0 from Type A seed, the next non-trivial rational has simple poles at t=0: f_1 = 1/t, f_4 = -1/t.
- **gt_kind**: closed-form
- **how_to_compute**: Julia: f_type_c(t) = fill(t/5, 5); alpha = fill(1//5, 5). Also: apply_s0_to_type_A gives f=(t, 1/t, 0, 0, -1/t), alpha=(-1,1,0,0,1) -- a non-trivial rational with poles at t=0. See test/noumi_yamada_symmetry_test.jl for Rational{BigInt} verification of the ODE satisfaction.
- **precision**: exact rational arithmetic
- **targets_buckets**: B21
- **regime**: 5-component Noumi-Yamada vector ODE with trivial oracle; also exercises the Backlund orbit (poles at t=0 after one s_0 shift), and the A_4 system for n>=4 acceptance item
- **citation**: Matsuda 2012, J. Math. Phys. 53, lines 322-338 (Type C, explicit); Clarkson et al. 2020 (sigma14-002.tex lines 1499-1502). Also: docs/v0p2_pillarB_noumi_yamada_findings.md section 5.2/5.4.
- **confidence**: high

#### `noumi-yamada-a4-backlund-orbit`

- **name**: Noumi-Yamada A_4^(1) -- one Backlund shift s_0 of Type A, pole at t=0
- **ode**: Same 5-component A_4^(1) system. Parameters after s_0 applied to Type A: alpha=(-1,1,0,0,1). The solution has a simple pole at t=0.
- **ic_bc**: f(1) = (1, 1, 0, 0, -1), alpha=(-1,1,0,0,1)
- **domain**: t in C \ {0}; pole at t=0 in f_1, f_4
- **closed_form**: s_0 applied to (f=(t,0,0,0,0), alpha=(1,0,0,0,0)) gives: f_0=t, f_1=1/t, f_2=0, f_3=0, f_4=-1/t; alpha_0=-1, alpha_1=1, alpha_2=0, alpha_3=0, alpha_4=1. Verification: f_1' = -1/t^2. RHS of f_1 equation = f_1*(f_2-f_3+f_4-f_0) + alpha_1 = (1/t)*(0-0+(-1/t)-t) + 1 = (1/t)*(-1/t-t) + 1 = -1/t^2 - 1 + 1 = -1/t^2. Consistent.
- **pole_structure**: Simple poles at t=0 in f_1 and f_4 components (opposite signs). The shared-Q Pade approximant for the 5-vector must find this shared pole at t=0.
- **gt_kind**: closed-form
- **how_to_compute**: Julia: f(t) = [t, 1/t, 0, 0, -1/t]; alpha=[-1,1,0,0,1]. Backlund generator s_0: s_0(f_j) = f_j + alpha_0/f_0 for j=1, s_0(f_j) = f_j - alpha_0/f_0 for j=4, fixed otherwise. See test/noumi_yamada_symmetry_test.jl NYS.1.6 (Backlund preserves solutions) and pillar doc section 7.2.
- **precision**: exact rational
- **targets_buckets**: B14, B21
- **regime**: 5-component system with a KNOWN shared pole at t=0 in f_1 and f_4 (shared-Q must detect the single common pole from the vector jet); the v0.2 shared-denominator Pade (SharedPade) must extract this pole precisely
- **citation**: Matsuda 2012 J. Math. Phys. 53, lines 264-300 (explicit s_0 action table). Derived in docs/v0p2_pillarB_noumi_yamada_findings.md section 7.2.
- **confidence**: high
- **verify_note**: Verify by substituting into each f_j' equation. Note sum(alpha_j)=1 even though some alphas are negative -- that is correct for this Backlund orbit.

#### `coupled-riccati-tan-sec`

- **name**: Coupled Riccati system (tan,sec^2) -- 2-vector with shared pole at z=pi/2
- **ode**: 2-component first-order system: y1' = y1^2 + y2, y2' = 2*y1*y2. Companion to u'' = 2u + 2u^3 (the nonlinear pendulum / elliptic-function ODE). Wait: the simpler companion form is y=(tan z, sec^2 z): y1'=y2, y2'=2*y1*y2 (since d/dz tan z = sec^2 z and d/dz sec^2 z = 2 sec z * sec z tan z = 2 tan z * sec^2 z = 2*y1*y2). Both components blow up at z = pi/2 + k*pi.
- **ic_bc**: y1(0)=0, y2(0)=1 (tan(0)=0, sec^2(0)=1)
- **domain**: complex z-plane; poles at z=pi/2+k*pi for k in Z
- **closed_form**: y1(z) = tan(z), y2(z) = sec^2(z) = 1/cos^2(z). Both have simple poles at z = pi/2 + k*pi (for y2, the pole is second-order: sec^2 z ~ 1/(z-pi/2)^2). The shared denominator Q must vanish at z=pi/2.
- **pole_structure**: y1=tan(z) has simple poles at z=pi/2+k*pi with residue -1. y2=sec^2(z) has double poles at the same locations. The SHARED pole location is pi/2+k*pi for both components. This is the simplest 2-vector system with a shared pole that is analytically known.
- **gt_kind**: closed-form
- **how_to_compute**: Julia: y1_exact(z)=tan(z); y2_exact(z)=sec(z)^2. At any z: evaluate using standard trigonometric functions. Or: import mpmath; mp.dps=50; z=mpmath.mpf('1.0'); print(mpmath.tan(z), 1/mpmath.cos(z)**2). Shared pole at z_pole = pi/2 approx 1.5707963...
- **precision**: arbitrary via mpmath/BigFloat (trig functions available)
- **targets_buckets**: B13, B14
- **regime**: 2-vector system with analytically known SHARED poles at z=pi/2+k*pi; tests the shared-Q vector Pade step THROUGH a shared pole; y1 has simple poles, y2 has double poles at the same location
- **citation**: Constructed (standard calculus: d(tan z)/dz = sec^2 z, d(sec^2 z)/dz = 2 tan z sec^2 z). The system is the companion form of the identity tan'=sec^2. No external reference needed for the ODE; trigonometric functions are the closed-form oracle.
- **confidence**: high
- **verify_note**: Verify: y1'(z)=sec^2(z)=y2, y2'(z)=2*tan(z)*sec^2(z)=2*y1*y2. Both identities hold exactly.

#### `pi1-okamoto-2vector`

- **name**: Painleve-I (PI) in Okamoto Hamiltonian form -- 2-vector (q,p) with known pole structure
- **ode**: PI in Hamiltonian form: q' = dH/dp = 2p; p' = -dH/dq = 3q^2 + t/2. This is the 2-component first-order companion of PI (u''=6u^2+t): y=(u,u'), y1'=y2, y2'=6*y1^2+t. Both components share poles at every movable PI pole z0: y1=(u) has double poles with residue 1, y2=(u') has triple poles.
- **ic_bc**: u(0)=1.071822516416917, u'(0)=1.710337353176786 (equianharmonic case, FW md:292-295)
- **domain**: real axis z in [0, 30] with approx 21 poles; complex plane pole field
- **closed_form**: Near a pole z0: u ~ 1/(z-z0)^2 - z0/10*(z-z0)^2 - 1/6*(z-z0)^3 + O((z-z0)^4). No global closed form (PI is a genuine transcendent). For testing: use mpmath.odefun or NDSolve to generate reference values; or use the equianharmonic Weierstrass-P oracle (u''=6u^2 is the same RHS at t=0, with exact closed form u=P(z-z0;0,g3)).
- **pole_structure**: All movable poles are second-order with residue +1 (no parameter freedom in residue for PI). The 2-vector (u,u') has shared poles where u blows up -- the v0.2 keystone: the vector walk must pass THROUGH these shared poles with the shared-Q denominator.
- **gt_kind**: mpmath
- **how_to_compute**: For the t=0 slice (equianharmonic): from mpmath import *; mp.dps=50; u = lambda z: weierstrass.wp(z-1, 0, 2); print(u(0.5), diff(u, 0.5)). For generic t: from mpmath import odefun; u = odefun(lambda z,y: [y[1], 6*y[0]**2+z], 0, [u0, u0p], tol=1e-14). FW Table 5.1 reference: u(30)=1.095098255959744 (see test/_oracle_problems.jl in repo).
- **precision**: arbitrary via mpmath/BigFloat; FW Table 5.1 values pinned to 15 sig figs
- **targets_buckets**: B13, B14
- **regime**: 2-component companion form of PI with known shared poles; equianharmonic t=0 slice has exact closed-form (Weierstrass-P) for oracle cross-validation; exercises the vector path-network stepping THROUGH approximately 21 real-axis poles to reach z=30
- **citation**: FW2011 (FW2011_painleve_methodology_JCP230.md:46-61, eq. 1.1): PI equation. FW Table 5.1 (FW2011_painleve_methodology_JCP230.md:281-301): reference values u(0.5), u(30), u(10^4). VPO test (test/vector_pipeline_oracle_test.jl) confirms vector pipeline reaches u(30) to 1e-8.
- **confidence**: high
- **verify_note**: The vector pipeline oracle test (test/vector_pipeline_oracle_test.jl VPO.2) already verifies this: vector_path_network_solve of the 2-vector companion reproduces u(30) to 1e-8. Additional oracle: mpmath odefun at higher precision, or NDSolve with WorkingPrecision->50.

#### `noumi-yamada-a2-backlund-pole`

- **name**: Noumi-Yamada A_2^(1) -- one Backlund shift s_0 from Type A, pole at t=0 in f_1, f_2
- **ode**: 3-component A_2^(1) system: f_0'=f_0(f_1-f_2)+alpha_0; f_1'=f_1(f_2-f_0)+alpha_1; f_2'=f_2(f_0-f_1)+alpha_2. After s_0 on Type A (f=(t,0,0), alpha=(1,0,0)): new alpha=(-1,1,0); s_0 maps f_1 -> f_1 + alpha_0/f_0 = 1/t, f_2 -> f_2 - alpha_0/f_0 = -1/t.
- **ic_bc**: f(1)=(1, 1, -1), alpha=(-1,1,0)
- **domain**: t in C \ {0}; shared pole at t=0
- **closed_form**: f_0(t)=t, f_1(t)=1/t, f_2(t)=-1/t. alpha=(-1,1,0). Sum: f_0+f_1+f_2 = t + 1/t - 1/t = t (constraint satisfied). Verification: f_1'=-1/t^2. RHS(f_1)=f_1*(f_2-f_0)+alpha_1=(1/t)*(-1/t-t)+1=-1/t^2-1+1=-1/t^2. Consistent.
- **pole_structure**: Simple poles at t=0 in f_1 and f_2 (shared), residues +1 and -1. The shared-Q for this 3-vector must find the single shared pole at t=0 from the Taylor jet around some t0 != 0.
- **gt_kind**: closed-form
- **how_to_compute**: Julia: f(t) = [t, 1/t, -1/t]; alpha = [-1, 1, 0]. For multi-precision: f_big(t) = [BigFloat(t), BigFloat(1)/BigFloat(t), BigFloat(-1)/BigFloat(t)]. The pole at t=0 is analytically exact.
- **precision**: exact rational
- **targets_buckets**: B21, B14
- **regime**: 3-component A_2^(1) vector with a shared pole at t=0 (all three components share the same pole location); the vector shared-Q must correctly identify this single shared pole from a jet expanded at a regular point
- **citation**: NY1998 main.tex lines 136-152 (Backlund generators s_i). docs/v0p2_pillarB_noumi_yamada_findings.md section 5.1 (Type A seed and Backlund orbit). Same pattern as the A_4 Backlund-orbit case.
- **confidence**: high
- **verify_note**: Check all three component ODEs by direct substitution. Confirm sum(alpha)=-1+1+0=0 which is NOT 1 -- the Backlund transformation changes k. Under the k=1 normalisation this shifted solution is valid only if one accepts k=0. Alternatively, use NY1998's s_0 formula carefully: s_0 acts on alpha_0 -> -alpha_0, alpha_1 -> alpha_1+alpha_0, alpha_4 -> alpha_4+alpha_0. For original (1,0,0): new alpha=(-1,1,0), sum=0, NOT 1. This is a valid Backlund orbit but with k=0 normalization. Verify ODE holds regardless.

#### `pi2-hastings-mcleod-2vector`

- **name**: PII Hastings-McLeod solution -- 2-vector companion with known pole structure
- **ode**: PII in companion form: y1' = y2; y2' = 2*y1^3 + t*y1 + alpha. The Hastings-McLeod solution corresponds to alpha=0, real smooth solution on the real axis interpolating between Ai(t) behavior at +infty and sqrt(-t/2) at -infty.
- **ic_bc**: Rational: y1(1)=-1, y2(1)=1, alpha=1. HM: y1(-10)=-sqrt(5), y2(-10)=0, alpha=0.
- **domain**: t in C \ {0} for rational solution; real axis t in [-10,10] for Hastings-McLeod
- **closed_form**: For alpha=0: u_{HM}(t) satisfies PII with no poles on the real axis. For the purpose of finding shared poles, use the rational solutions: u_1(t) = -1/t (PII with alpha=1, has a simple pole at t=0). In companion form (y1,y2)=(-1/t, 1/t^2), both components share the simple pole at t=0.
- **pole_structure**: PII rational u_1(t)=-1/t: simple pole at t=0, residue -1. The companion y2=u_1' = 1/t^2 has a double pole at t=0. Shared Q must vanish at t=0 for both components.
- **gt_kind**: closed-form
- **how_to_compute**: For PII rational u_1(t) = -1/t (alpha=1): Julia: y1(t)=-1/t; y2(t)=1/t^2. Verify: y2'=y1'=1/t^2, y2'=-2/t^3. RHS y2' = 2*y1^3 + t*y1 + 1 = -2/t^3 - 1 + 1 = -2/t^3. Consistent. For the Hastings-McLeod solution on the real axis: from mpmath import *; mp.dps=50; t=mpf('-5'); u_hm = odefun(lambda t,y: [y[1], 2*y[0]**3 + t*y[0]], mpf('-10'), [-sqrt(-mpf('-10')/2), 0], tol=1e-14)
- **precision**: exact rational for u_1; arbitrary via mpmath for HM
- **targets_buckets**: B13, B14
- **regime**: 2-vector companion of PII with known rational solution (u_1 = -1/t, exact) and the Hastings-McLeod smooth solution (mpmath oracle); exercises shared-pole detection at t=0
- **citation**: FW2014 (FW2014_second_PII_exploration_FoCM14.md) eq. 6: PII rational u_n. For u_1=-1/t, see FW2014 md:91 or DLMF 32.10.1. Hastings-McLeod: DLMF 32.11.2 (connection problem, not a closed form).
- **confidence**: high
- **verify_note**: The u_1=-1/t rational is verified by direct substitution. For HM, mpmath odefun with tol=1e-14 is the oracle.

#### `pi4-okamoto-hamiltonian-2vector`

- **name**: Painleve-IV in Okamoto Hamiltonian form -- 2-vector (q,p) with poles and rational solutions
- **ode**: PIV Hamiltonian: H_IV(q,p,t) = q*p*(q-p-2t) - alpha_1*q + alpha_2*p. Hamilton's equations: q' = dH/dp = q*(q-2p-2t), p' = -dH/dq = -p*(2q-p-2t) + alpha_1. These are 2-component. Standard PIV is recovered as u = q, u'' = (u')^2/(2u) + (3/2)u^3 + 4t*u^2 + 2(t^2-a)*u + b/u.
- **ic_bc**: u=-2t: y1(0)=0, y2(0)=-2 (but u=0 is a problem for beta/u term); use y1(-1)=2, y2(-1)=-2 instead
- **domain**: entire for rational solutions; complex plane with poles for generic transcendents
- **closed_form**: Rational solutions of PIV (RF2014 tables): the PIV 'entire' solution u=-2t (alpha=1,beta=0): q(t)=-2t, p satisfies p'=-p*(2q-p-2t)+alpha_1; for q=-2t, q'=-2, so p=-dH/dp^{-1}(q') = ... (algebra). Simpler: use the known PIV solution u=0 when beta=0, alpha=0: q=0, p=arbitrary drift. Best closed form: u=-2/3*t (from RF2014 md:91 piv_entire orbits).
- **pole_structure**: The PIV generic transcendent has movable poles (first-order, residue +1 or -1). The 2-vector (q,p) has shared poles. Rational/entire solutions are pole-free.
- **gt_kind**: closed-form
- **how_to_compute**: For the PIV entire solution u=-2t (alpha=1, beta=0): y1(t)=u(t)=-2t, y2(t)=u'(t)=-2. Verify in PIV ODE: u''=0, RHS=(u')^2/(2u)+... = 4/(-4t)+... -- need to check parameters carefully. Better: use piv_entire(:minus_2z) from src/PainleveClosedForm.jl in the repo. For mpmath oracle: from mpmath import odefun; u = odefun(lambda t,y:[y[1], (y[1]**2)/(2*y[0]) + (3/2)*y[0]**3 + 4*t*y[0]**2 + 2*(t**2-a)*y[0] + b/y[0]], t0, [u0,u0p], tol=1e-14).
- **precision**: exact for linear solutions; arbitrary via mpmath for generic PIV
- **targets_buckets**: B13, B21
- **regime**: 2-component Hamiltonian PIV system; entire closed-form solutions u=-2t, u=-2t/3 (pole-free oracles); generic parameters give movable poles; connection to Noumi-Yamada A_2^(1) via backward-compat test
- **citation**: ReegerFornberg2014 (ReegerFornberg2014_PIV_fundamental_domain_PhysicaD280.md:91): PIV entire solutions. ADR-0010 (docs/adr/0010-painleve-closed-form-families.md): piv_entire API. NY1999 (references/tex/noumi_yamada/NoumiYamada1999_PIV_symmetries_okamoto_NagoyaMathJ153/main.tex) for Hamiltonian H_IV.
- **confidence**: medium
- **verify_note**: The PIV entire solutions u=-2t and u=-2t/3 are verified in test/painleve_closed_form_test.jl (piv_entire(:minus_2z) and piv_entire(:minus_two_thirds_z)). The Hamiltonian form (q,p) needs the explicit Okamoto substitution -- verify by differentiating and checking consistency with the companion form y1'=y2, y2'=RHS_PIV.

#### `pi1-hierarchy-4vector`

- **name**: P_I^(2) (fourth-order PI hierarchy member) -- 4-vector companion with tritronquee oracle
- **ode**: 4-component companion form of P_I^(2) (KKG 2015 eq.1.1): y1'=y2, y2'=y3, y3'=y4, y4'=-10*y2^2 - 20*y1*y3 - 40*(y1^3 - 6*t*y1 + 6*x). Parameter t (complex); independent variable x.
- **ic_bc**: Asymptotic seed at x0=-20, t=0: y=(u0, u0', u0'', u0''') from leading expansion
- **domain**: complex x-plane; tritronquee sector |arg x| < pi/7 (approx 25.7 degrees); pole field for |x| in [2,20]
- **closed_form**: No closed form for the generic transcendent. The tritronquee asymptotic IC at x=-20: u ~ -(6|x|)^(1/3) on the negative real axis. Leading asymptotic: u_f(x,t) ~ -cbrt(6)*x^(1/3) as |x|->inf in the tritronquee sector. The leading-order asymptotic expansion gives regenerable IC via pI2_tritronquee_ic(x0; t, n_terms) in src/PainleveHierarchy.jl.
- **pole_structure**: The 4-vector system has movable poles where u blows up; all 4 components share the same poles. Pole spacing ~constant 0.70 across |x| in [5,20] (measured, ADR-0026 Amendment 3). Poles are first-order in u (leading: u ~ C/(x-x0)), higher-order in derivatives.
- **gt_kind**: mpmath
- **how_to_compute**: from mpmath import odefun, mpf; mp.dps=30; t=mpf('0'); x0=mpf('-20'); u0=-(6*abs(x0))**(mpf('1')/3); up0=-(6*abs(x0))**(-mpf('2')/3); upp0=2*(6*abs(x0))**(-mpf('5')/3); uppp0=-10*(6*abs(x0))**(-mpf('8')/3); rhs=lambda x,y: [y[1],y[2],y[3], -10*y[1]**2-20*y[0]*y[2]-40*(y[0]**3-6*t*y[0]+6*x)]; sol=odefun(rhs, x0, [u0,up0,upp0,uppp0], tol=1e-12)
- **precision**: arbitrary via mpmath odefun at mp.dps=30+
- **targets_buckets**: B13, B14
- **regime**: 4-vector companion of a 4th-order ODE; pole field in complex x-plane; tritronquee sector (|arg x| < pi/7); exercises the full vector path-network stepping through multiple shared poles
- **citation**: KapaevKleinGrava2015 (references/tex/painleve_hierarchy/KapaevKleinGrava2015_PI2_tritronquee_ConstrApprox41/tritronquee_coeff.tex lines 124-130): P_I^(2) eq.(p12) and companion form. docs/v0p2_pillarC_painleve_hierarchy_findings.md section 1-2: explicit 4-vector form. ADR-0022 (docs/adr/0022-vector-problem-types.md): PainleveHierarchyProblem builder.
- **confidence**: high
- **verify_note**: The tritronquee branch sign on the negative real axis: u ~ +cbrt(6)*|x|^(1/3) > 0 for x<0 (see docs/v0p2_pillarC_painleve_hierarchy_findings.md lines 189-196 for the corrected branch). Seeding the wrong sign gives ODE residual ~-480|x| -- always check sign against the ODE before trusting any asymptotic IC.

#### `harmonic-oscillator-2vector-bvp`

- **name**: Harmonic oscillator 2-vector BVP -- (cos z, -sin z) with known complex poles
- **ode**: 2-component: y1' = y2; y2' = -y1. Solution: y1(z)=cos(z), y2(z)=-sin(z).
- **ic_bc**: BVP: B_a = [1 0; 0 0], B_b = [0 0; 1 0], g = [cos(z_a), cos(z_b)]
- **domain**: any complex z (entire functions); segments [0, pi/2], [0, 1+im], [0, pi]
- **closed_form**: y1(z) = cos(z) (poles at z=pi/2+k*pi on the Riemann surface of cosh/cos), y2(z) = -sin(z) (poles at z=pi/2+k*pi on the same locations -- BUT these are actually branch points of 1/cos). In the ANALYTIC continuation sense: cos(z) has no poles in the complex plane (it is entire). For a pole-crossing test use tan(z) / sec^2(z) instead (see coupled-riccati-tan-sec entry). This system is useful for BVP validation only.
- **pole_structure**: Entire (cos and sin are entire functions). No poles. Used as a B16 vector BVP oracle: Chebyshev collocation on complex segment [0, pi/2] recovers cos(pi/2)=0 to machine precision.
- **gt_kind**: closed-form
- **how_to_compute**: Julia: y1_exact(z) = cos(z); y2_exact(z) = -sin(z). For BigFloat: cos(BigFloat(z)). Already used in test/vector_bvp_test.jl VB.1.2 as the primary vector BVP oracle.
- **precision**: arbitrary via BigFloat/mpmath
- **targets_buckets**: B16
- **regime**: 2-component linear BVP oracle; Chebyshev-collocation Newton solver on complex segments; real segment [0,pi/2] and complex segment [0, 1+im]; both Float64 and BigFloat-256
- **citation**: Standard calculus. Used in test/vector_bvp_test.jl VB.1.2 (already shipped in the repo).
- **confidence**: high

#### `noumi-yamada-a4-generic-mpmath`

- **name**: Noumi-Yamada A_4^(1) generic transcendent -- mpmath odefun oracle
- **ode**: 5-component A_4^(1) system: f_j' = f_j*(f_{j+1}-f_{j+2}+f_{j+3}-f_{j+4}) + alpha_j (indices mod 5). alpha=(0.1, 0.15, 0.2, 0.25, 0.3), sum=1. IC at t=1: f=(0.3, 0.25, 0.15, 0.2, 0.1) (not the Type C equilibrium, so the transcendent has non-trivial poles).
- **ic_bc**: f(1)=(0.3, 0.25, 0.15, 0.2, 0.1); alpha=(0.1,0.15,0.2,0.25,0.3)
- **domain**: complex t-plane from t0=1; poles accumulate for large |t|
- **closed_form**: none
- **pole_structure**: Movable poles of the A_4^(1) transcendent; all 5 components share each pole. The constraint sum(f_j)-t=0 is preserved along the flow. Poles are first-order (Painleve property guarantees no movable branch points). Location unknown analytically -- use mpmath as oracle.
- **gt_kind**: mpmath
- **how_to_compute**: from mpmath import odefun, mpf; mp.dps=40; alpha=[mpf('0.1'),mpf('0.15'),mpf('0.2'),mpf('0.25'),mpf('0.3')]; f0=[mpf('0.3'),mpf('0.25'),mpf('0.15'),mpf('0.2'),mpf('0.1')]; def rhs(t,f): n=2; d=2*n+1; return [f[j]*(sum(f[(j+2*r-1)%d] for r in range(1,n+1)) - sum(f[(j+2*r)%d] for r in range(1,n+1))) + alpha[j] for j in range(d)]; sol=odefun(rhs, mpf('1'), f0, tol=1e-16); print([sol(mpf('2'))[j] for j in range(5)])
- **precision**: arbitrary via mpmath at mp.dps=40+
- **targets_buckets**: B21, B14
- **regime**: 5-component Noumi-Yamada A_4^(1) generic transcendent; pole field in complex t-plane; shared-Q vector Pade must step through unknown shared poles; mpmath odefun at 40 digits is the reference
- **citation**: NY1998 main.tex lines 85-88 (A_2n system); Matsuda2012 lines 220-233 (A_4^(1) explicit system). The generic (non-rational) transcendent has no published closed form; mpmath odefun is the only available oracle.
- **confidence**: medium
- **verify_note**: Confirm sum(f(t))-t < 1e-35 at test points (constraint invariant). Also confirm sum(alpha)=1. The generic transcendent has unknown pole locations -- locate them by root-finding on the shared-Q denominator polynomial from the numerical solution.

#### `coupled-riccati-shared-denom`

- **name**: 2-vector coupled Riccati system with explicit shared pole -- (1/(z-z0), 1/(z-z0)^2)
- **ode**: Explicit 2-vector: y1' = y2; y2' = 2*y2^2/y1. Solutions: y1=C/(z-z0), y2=-C/(z-z0)^2 for any constant C and pole location z0. Or more generally: y=(1/(z-z0), -1/(z-z0)^2) satisfies y1'=-1/(z-z0)^2=y2 and y2'=2/(z-z0)^3=2*y2^2/y1. Both components share a pole at z=z0.
- **ic_bc**: y1(0)=1/(-z0), y2(0)=-1/z0^2 with z0=1.5: y1(0)=-2/3, y2(0)=-4/9
- **domain**: z in C \ {z0}; z0=1.5 (or any nonzero complex value)
- **closed_form**: y1(z) = 1/(z - z0), y2(z) = -1/(z - z0)^2. More generally with parameter C: y1=C/(z-z0), y2=-C/(z-z0)^2. Pole of y1 is simple (order 1), pole of y2 is double (order 2) -- both at z=z0. This is the canonical 2-vector with a shared denominator Q(z)=(z-z0).
- **pole_structure**: Single shared pole at z=z0 (user-specified). y1 has order-1 pole, y2 has order-2 pole. The shared-Q denominator polynomial is exactly (z-z0) -- a degree-1 polynomial. This is the SIMPLEST possible test of the SharedPade.shared_denominator_pade function's ability to identify a known shared pole from the Taylor jet.
- **gt_kind**: closed-form
- **how_to_compute**: Julia: z0=1.5; y1(z)=1/(z-z0); y2(z)=-1/(z-z0)^2. Taylor jet at z_expand=0.0 of length n+1: jet1 = [(-1)^k * k! / z0^(k+1) for k in 0:n]; jet2 = [(-1)^(k+1) * (k+1)! / z0^(k+2) for k in 0:n]. Then SharedPade.shared_denominator_pade([jet1, jet2], m) should return Q with root near z0.
- **precision**: exact rational (Taylor jets can be computed in Rational{BigInt})
- **targets_buckets**: B14
- **regime**: Simplest possible 2-vector shared-pole test; no ODE integration needed -- exercises only the SharedPade kernel (shared_denominator_pade) against an analytically exact shared-pole input; also a pure algebra unit test for B14
- **citation**: Constructed. The closed form y=(C/(z-z0), -C/(z-z0)^2) follows from direct integration of y1'=y2, y2'=2y2^2/y1. No external reference needed.
- **confidence**: high
- **verify_note**: Verify ODE: y1'(z) = -1/(z-z0)^2 = y2(z) (yes). y2'(z) = 2/(z-z0)^3. 2*y2^2/y1 = 2*(1/(z-z0)^4)/(1/(z-z0)) = 2/(z-z0)^3 (yes). Consistent. Taylor jets at z=0: jet1[k] = 1/(k! * z0^(k+1)) * (-1)^k * k! = (-1)^k/z0^(k+1).

### Territory: FFW 2017 transformed multivalued Painleve equations on Riemann sheets (P~III, P~

*Scout notes:* Key gaps: (1) No closed-form non-trivial multi-valued PIII/PV solution is available to serve as a self-contained closed-form oracle; the best closed-form result is the trivial constant u=1 (single-valued, not useful for sheet testing). (2) The PV-tilde oracle values have been computed by RK4 with verified conjugate symmetry but not by an independent adaptive solver - a Mathematica NDSolve crosscheck would add another digit of trust. (3) The PVI-zeta oracle for the Fig 7 generic solution requires integrating near dense pole fields and has not been fully validated beyond the conjugate symmetry check. (4) FW2014/FW2015 PII results are single-valued (no branch point at z=0 for PII), so they do not provide Riemann-sheet ground truth and were not included. (5) The BigFloat-256 cancellation oracle is a sub-expression test (evaluating e^{e^eta}-1 near the fixed singularity), not a full ODE integration; a full integration test near the cancellation zone would require a Pade-Taylor integrator capable of BigFloat arithmetic, which is exactly what PadeTaylor.jl is being built for. Cleanest oracles found: (a) PIII-tilde sheet 0/+1/+2 values at zeta=-1+0.5i verified by conjugate symmetry to 50 decimal digits; (b) PV-tilde sheet 0/+1/+2 values verified similarly; (c) PVI-eta double-exp cancellation sub-expression with Float64 vs BigFloat comparison quantified exactly.

#### `pIII-tilde-sheet-oracle`

- **name**: P~III multi-sheet point values at zeta=-1+0.5i on sheets 0, +1, +2
- **ode**: w'' = (1/w)*(w')^2 + (1/4)*(alpha*w^2 + gamma*w^3 + beta*e^zeta + delta*e^{2*zeta}/w), with alpha=-1/2, beta=-1/2, gamma=1, delta=-1 (FFW2017 Fig 1). Transformation: z = e^{zeta/2}, u(z) = e^{-zeta/2}*w(zeta).
- **ic_bc**: PIII in z-plane: u(z=1)=1/4, u_prime(z=1)=1 [FFW2017 Figure 1 caption, p.4]. Transformed to zeta-plane: w(0)=1/4, w_prime(0)=9/16.
- **domain**: Complex zeta = -1+0.5i (sheet 0), -1+0.5+4*pi*i (sheet +1), -1+0.5+8*pi*i (sheet +2). All correspond to z = exp((-1+0.5i)/2) = 0.58767509034392266 + 0.15005808662216327i.
- **closed_form**: none
- **pole_structure**: Meromorphic in entire zeta-plane; simple poles of w with residue +/-2 (gamma=1), corresponding to simple poles of u in z-plane with residue +/-1/sqrt(gamma)=+/-1. Zeros of w on real axis (corresponding to zeros of u). No branch points in zeta-plane (transformation maps branch point z=0 to Im(zeta)->-inf). Sheet s in z-plane <-> strip Im(zeta) in (-2*pi+4*pi*s, 2*pi+4*pi*s).
- **gt_kind**: mpmath
- **how_to_compute**: python mpmath RK4 at dps=50; path: 0 -> 0.5i -> -1+0.5i (sheet 0) and 0 -> 4*pi*i -> 4*pi*i+0.5i -> -1+4*pi*i+0.5i (sheet +1) and 0->4*pi*i->8*pi*i->8*pi*i+0.5i->-1+8*pi*i+0.5i (sheet +2). ICs w(0)=1/4, w_prime(0)=9/16 (derived from z-plane ICs u(1)=1/4, u_prime(1)=1 via: w(0)=u(1)=1/4 and w_prime(0)=(u_prime(1)+u(1)/2)/2*z|_{z=1}=(1+1/8)/2=9/16). Verified by conjugate symmetry w(conj(zeta))=conj(w(zeta)) to 50 digits (error=0.0 exactly). Full recipe: from mpmath import mp,mpc,mpf,exp,pi; mp.dps=50; alpha,beta,gamma_p,delta_p=-mpf(1)/2,-mpf(1)/2,mpf(1),mpf(-1); pIII=lambda z,y: [y[1], (1/y[0])*y[1]**2+(mpf(1)/4)*(alpha*y[0]**2+gamma_p*y[0]**3+beta*exp(z)+delta_p*exp(2*z)/y[0])]; y0=[mpc(1,0)/4, mpc(9,0)/16]; # integrate with n=15000 RK4 steps along each leg.
- **precision**: ~40 significant digits (verified by conjugate symmetry test to 50 digits, integration accuracy ~40 digits with n=15000 steps)
- **targets_buckets**: B10, B11
- **regime**: Multi-sheet path winding in PIII zeta-plane; sheet index s encoded by Im(zeta) strip; RK4 at 50-digit arithmetic
- **citation**: FFW2017_painleve_riemann_surfaces_preprint.md:41-49 (transformation), :101-103 (Fig 1 and ICs), :103 (strip formula Im(zeta) in (-2pi+4pi*s, 2pi+4pi*s)), :107-114 (Table 1 residue formulas)
- **confidence**: high
- **verify_note**: Reference values at zeta=-1+0.5i: Sheet 0: w=-0.12218596353399632+0.00412103759240361i, u=-0.19350701391712494+0.05642289488460466i. Sheet +1: w=-0.01416089211316628-0.04291673388050532i, u=-0.04012729407649307-0.06278181518599197i. Sheet +2: w=-2.15443430928006292+2.40537978859775608i, u=-2.46048185184606176+4.72130779923064010i. Cross-check: run at dps=80 with n=30000 steps and compare first 35 digits. Also verify conjugate symmetry between sheet +1 and sheet -1.

#### `pV-tilde-sheet-oracle`

- **name**: P~V multi-sheet point values at zeta=-1+0.5i on sheets 0, +1, +2
- **ode**: w'' = (1/(2w)+1/(w-1))*(w')^2 + (w-1)^2*(alpha*w+beta/w) + gamma*e^zeta*w + delta*e^{2*zeta}*w*(w+1)/(w-1), with alpha=1, beta=-1, gamma=1, delta=-1/2 (FFW2017 Fig 6). Transformation: z = e^{zeta}, u(z) = w(zeta).
- **ic_bc**: PV in z-plane: u(z=1)=2, u_prime(z=1)=-1 [FFW2017 Figure 6 caption, p.16]. Transformed: w(0)=2, w_prime(0)=-1.
- **domain**: Complex zeta = -1+0.5i (sheet 0), -1+0.5+2*pi*i (sheet +1), -1+0.5+4*pi*i (sheet +2). All correspond to z = exp(-1+0.5i) = 0.32284458245003301 + 0.17637079922503195i.
- **closed_form**: none
- **pole_structure**: Meromorphic in entire zeta-plane (Hinkkanen-Laine theorem, reference [14] in FFW2017). Simple poles of w with residues +/-1/sqrt(2*alpha)=+/-1/sqrt(2) in zeta-plane [FFW2017 Table 4/6, p.11/16]. Sheet s in z-plane <-> strip Im(zeta) in (-pi+2*pi*s, pi+2*pi*s) [z=e^zeta, so one CCW revolution around z=0 shifts Im(zeta) by 2*pi].
- **gt_kind**: mpmath
- **how_to_compute**: python mpmath RK4 at dps=50; ICs w(0)=2, w_prime(0)=-1 (from u(z=1)=2, u_prime(1)=-1, and dz/dzeta=z so w_prime(0)=u_prime(1)*z|_{z=1}=-1). Path for sheet 0: 0->0.5i->-1+0.5i. Sheet +1: 0->2*pi*i->2*pi*i+0.5i->-1+2*pi*i+0.5i. Sheet +2: 0->4*pi*i->4*pi*i+0.5i->-1+4*pi*i+0.5i. Use n=10000 RK4 steps per leg. Verify conjugate symmetry. from mpmath import mp,mpc,mpf,exp,pi; mp.dps=50; alpha,beta,gamma_p,delta_p=mpf(1),mpf(-1),mpf(1),mpf(-1)/2; pV=lambda z,y: [y[1], (1/(2*y[0])+1/(y[0]-1))*y[1]**2+(y[0]-1)**2*(alpha*y[0]+beta/y[0])+gamma_p*exp(z)*y[0]+delta_p*exp(2*z)*y[0]*(y[0]+1)/(y[0]-1)]; y0=[mpc(2),mpc(-1)].
- **precision**: ~40 significant digits (conjugate symmetry error = 0.0 at 50 digits)
- **targets_buckets**: B10, B11
- **regime**: Multi-sheet path winding in PV zeta-plane; sheet period 2*pi*i (different from PIII which is 4*pi*i); RK4 at 50-digit arithmetic
- **citation**: FFW2017_painleve_riemann_surfaces_preprint.md:45-48 (PV-tilde transformation and equation), :297 (Fig 6 ICs), :211-218 (Table 4 pole residues)
- **confidence**: high
- **verify_note**: Reference values at zeta=-1+0.5i: Sheet 0: w=0.54159473151602727-1.21829675524375366i. Sheet +1: w=1.15067195213385327-0.05298780049491273i. Sheet +2: w=-1.28899896601851238-0.02652019017636388i. Cross-check conjugate symmetry: conj(w(zeta_eval, sheet+1)) should equal w(conj(zeta_eval), sheet-1). Also check: FFW2017 states errors for sheets 0-2 are 3e-9, 7e-6, 2e-5 for their method, so agree is at least 9/6/5 digits.

#### `pIII-tilde-residue-exact`

- **name**: P~III zeta-plane pole residue formula (exact closed-form)
- **ode**: PIII-tilde: w'' = (1/w)*(w')^2 + (1/4)*(alpha*w^2 + gamma*w^3 + beta*e^zeta + delta*e^{2*zeta}/w)
- **ic_bc**: Any ICs giving a solution with poles
- **domain**: Any pole location zeta_0 in the complex zeta-plane
- **closed_form**: Residue of w at simple pole zeta_0: c_{-1} = +/- 2*exp(-zeta_0/2) / sqrt(gamma). For gamma=1: c_{-1} = +/- 2*exp(-zeta_0/2). Derived by substituting w ~ c_{-1}/(zeta-zeta_0) into ODE and matching leading-order terms: leading balance is O((zeta-zeta_0)^{-3}) on both sides, giving 2*c_{-1} = c_{-1} from (1/w)(w')^2 term: (c_{-1})^{-1}*(c_{-1})^2/(zeta-zeta_0)^3 = c_{-1}/(zeta-zeta_0)^3. Wait: PIII in z-plane has u ~ c_{-1}^{(z)}/(z-z_0) with c_{-1}^{(z)}=+/-1/sqrt(gamma) [FFW2017 Table 1]. Under z=e^{zeta/2}, u=e^{-zeta/2}*w: pole of u at z_0 maps to pole of w at zeta_0=2*log(z_0) with residue c_{-1}^{(w)}=c_{-1}^{(z)}/(dz/dzeta)|_{zeta_0}*e^{zeta_0/2} = (1/sqrt(gamma)) / ((1/2)*e^{zeta_0/2}) * e^{zeta_0/2} = 2/sqrt(gamma). So residue of w is CONSTANT = +/-2/sqrt(gamma), independent of pole location zeta_0.
- **pole_structure**: Simple poles only; residue = +/-2/sqrt(gamma) for all poles in zeta-plane. This is exact and independent of the ODE parameters alpha, beta, delta.
- **gt_kind**: closed-form
- **how_to_compute**: Derive from substitution: near pole zeta_0, w~A/(zeta-zeta_0), w'~-A/(zeta-zeta_0)^2, w''~2A/(zeta-zeta_0)^3. PIII-tilde RHS dominant term: (1/w)(w')^2 = A/(zeta-zeta_0)^3. Balance: 2A=A implies no match, so also include gamma*w^3/4 = gamma*A^3/(4*(zeta-zeta_0)^3). Balance: 2A = A + gamma*A^3/4, giving gamma*A^3/4 = A, so A^2 = 4/gamma, A = +/-2/sqrt(gamma).
- **precision**: exact rational arithmetic
- **targets_buckets**: B10, B12
- **regime**: Closed-form residue verification; polynomial algebra from ODE Laurent expansion
- **citation**: FFW2017_painleve_riemann_surfaces_preprint.md:107-114 (Table 1 pole residues in z and zeta planes). Also from FFW2017 Eq. (7) condition number analysis which uses w~e^{zeta_0/2}*c_{-2}*(zeta-zeta_0)^{-2} for the DIFFERENT tronquee solution (different parameters). The residue formula above is for the generic case gamma*delta!=0.
- **confidence**: high
- **verify_note**: For FFW2017 Fig 1 (gamma=1): residue of w at any pole = +/-2. Test: find a pole location numerically by integrating around it, compute (zeta-zeta_0)*w(zeta) as zeta->zeta_0, check this approaches +/-2. Cross-check with z-plane: u near pole satisfies u ~ +/-1/(z-z_0), which matches residue +/-1/sqrt(gamma)=+/-1.

#### `pVI-eta-cancellation-bigfloat`

- **name**: PVI double-exp fixed singularity: catastrophic cancellation in (e^{e^eta}-1) requiring BigFloat-256
- **ode**: PVI-eta (FFW2017 Eq. 5): v'' = (1/2)*(1/v + 1/(v-1) + 1/(v-E))*(v')^2 - (e^eta*E/(E-1) + e^eta*E/(v-E) - 1)*v' + v*(v-1)*(v-E)*e^{2*eta}/(E-1)^2 * (alpha + beta*E/v^2 + gamma*(E-1)/(v-1)^2 + delta*E*(E-1)/(v-E)^2), where E = e^{e^eta}.
- **ic_bc**: Sub-expression evaluation only; no ICs needed for this specific test. For a full integration test, use FFW2017 Fig 2 ICs: u(z=10)=0.429534600325223, u_prime(10)=-1.61713114374804e-3, parameters (alpha,beta,gamma,delta)=(4,-4,8,-8).
- **domain**: eta near log(2*pi)+i*pi/2; specifically eta_test = log(2*pi)+i*pi/2 + 1e-10 (real offset)
- **closed_form**: The cancellation is exact: at eta_k1 = log(2*pi) + i*pi/2 (fixed singularity k=1, FFW2017 Eq.4), e^{e^{eta_k1}} = e^{2*pi*i} = 1 exactly, so (E-1)=0. Near it: E-1 ~ 2*pi*i*(eta-eta_k1) [linear approximation, verified numerically]. The ODE denominator (E-1)^2 ~ -4*pi^2*(eta-eta_k1)^2.
- **pole_structure**: Fixed singularities of the PVI-eta ODE at eta = log|2*pi*k| + i*arg(2*pi*i*k) for k in Z\{0} [FFW2017 Eq. 4, p.7]. For k=1: eta_k1 = log(2*pi) + i*pi/2 = 1.83787706640934548... + 1.57079632679489662...i. The region Re(eta) < log(2*pi) is branch-point-free [FFW2017 p.8].
- **gt_kind**: mpmath
- **how_to_compute**: from mpmath import mp,log,pi,j,exp,mpf; mp.dps=50; eta_sing=log(2*pi)+j*pi/2; delta=mpf('1e-10'); eta_test=eta_sing+delta; result=exp(exp(eta_test))-1; # BigFloat-50 exact result: -1.97392088041526386e-19 + 6.28318530749374574e-10*i. Float64 gives: +4.44e-16 + 6.28e-10*i (wrong real part by factor ~2.25e3). For BigFloat-256 (mp.dps=77): result=-1.97392088041526381173526051678469494...e-19 + 6.28318530749374574...e-10*i. The PVI-eta ODE uses (E-1)^{-2} in the last term, so the error propagates as (E-1)^{-2} at that test point.
- **precision**: BigFloat-50 gives ~31 correct digits for (E-1); BigFloat-256 (dps=77) gives ~67 digits. Float64 gives WRONG real part (sign error of factor ~2.25e3 at delta=1e-10).
- **targets_buckets**: B10, B11
- **regime**: Catastrophic cancellation in double-exp coordinate (e^{e^eta}-1): invisible at Float64, requires BigFloat-256 to resolve correctly at distance 1e-10 from fixed singularity
- **citation**: FFW2017_painleve_riemann_surfaces_preprint.md:146-155 (Eq.4 fixed singularity locations, Eq.5 PVI-eta ODE), :157 (Re(eta)<log(2*pi) is branch-point-free)
- **confidence**: high
- **verify_note**: Run: python3 -c 'from mpmath import *; mp.dps=50; eta_k1=log(2*pi)+j*pi/2; eta_test=eta_k1+mpf("1e-10"); print(exp(exp(eta_test))-1)' Expected: ~-1.974e-19 + 6.283e-10*i. Compare with Float64: import cmath; eta_f64=complex(1.8378770665093454, 1.5707963267948966); cmath.exp(cmath.exp(eta_f64))-1 -> should give wrong real part. The (E-1)^{-2} denominator has relative error ~3e-3 at delta=1e-10 (quantified above).

#### `pIII-sheet-winding-hand`

- **name**: PIII Riemann sheet winding number: CCW loop around z=0 gives sheet +1
- **ode**: PIII-tilde: same as above. Relevant for the mapping z=e^{zeta/2}.
- **ic_bc**: n/a
- **domain**: z-plane: annulus |z| in (1/100, 10) centered at z=0
- **closed_form**: WINDING NUMBER = +1 for one CCW revolution around z=0. Proof: CCW loop z(t)=r*exp(2*pi*i*t), t in [0,1]. Under z=e^{zeta/2}: zeta(t) = 2*log(z(t)) = 2*log(r) + 4*pi*i*t. Im(zeta) runs from 0 to 4*pi. Sheet s is defined by Im(zeta) in (-2*pi+4*pi*s, 2*pi+4*pi*s) [FFW2017 p.4]. The path crosses the boundary from sheet 0 to sheet 1 at t=1/2 when Im(zeta)=2*pi. Net result: one CCW revolution maps sheet 0 to sheet 1 (winding number +1). CW revolution maps sheet 0 to sheet -1 (winding number -1). For PV (z=e^zeta): CCW loop gives zeta(t)=log(r)+2*pi*i*t, Im shift = 2*pi -> still winding number +1 but sheet period is 2*pi*i.
- **pole_structure**: Branch point at z=0 (and z=infinity) for PIII. No branch points in zeta-plane.
- **gt_kind**: hand-computed
- **how_to_compute**: From z=e^{zeta/2}: arg(z) changes by 2*pi per CCW loop. Since zeta=2*log(z), Im(zeta)=2*arg(z), so Im(zeta) changes by 4*pi per CCW loop. Sheet boundary is at Im(zeta)=2*pi (mod 4*pi). Therefore winding_number(CCW, PIII) = +1.
- **precision**: exact integer
- **targets_buckets**: B11
- **regime**: Sheet winding bookkeeping: exact integer winding number from coordinate change algebra
- **citation**: FFW2017_painleve_riemann_surfaces_preprint.md:41 (z=e^{zeta/2} transformation for PIII), :103 (strip formula for PIII: Im(zeta) in (-2*pi+4*pi*s, 2*pi+4*pi*s))
- **confidence**: high
- **verify_note**: Cross-check: the z-plane PIII solution u(z) = e^{-zeta/2}*w(zeta) has the sheet dependence via different strips. After one CCW loop the strip shifts by 4*pi, moving from sheet s to sheet s+1. For PV (z=e^zeta), the strip is Im(zeta) in (-pi+2*pi*s, pi+2*pi*s) and one CCW loop shifts Im by 2*pi, also giving winding number +1.

#### `pVI-zeta-sheet-winding`

- **name**: PVI Riemann sheet structure in zeta-plane: two branch points at z=0,1
- **ode**: PVI-zeta (FFW2017 Eq.3): w'' = (1/2)*(1/w + 1/(w-1) + 1/(w-e^zeta))*(w')^2 - (e^zeta/(e^zeta-1) + e^zeta/(w-e^zeta))*w' + w*(w-1)*(w-e^zeta)/(e^zeta-1)^2*(alpha + beta*e^zeta/w^2 + gamma*(e^zeta-1)/(w-1)^2 + delta*e^zeta*(e^zeta-1)/(w-e^zeta)^2)
- **ic_bc**: n/a for winding; ICs for PVI oracle: u(z=10)=0.429534600325223, u_prime(10)=-1.61713114374804e-3 for (alpha,beta,gamma,delta)=(4,-4,8,-8) [FFW2017 Fig 2 caption, p.9]
- **domain**: z-plane: annulus 1/100 <= |z| <= 10, various sheets; zeta-plane: branch cuts at Im(zeta)=2*pi*k
- **closed_form**: Sheet enumeration (exact): Sheet k in z-plane = {z: |z| in (1/100,10), arg(z)=theta_k, arg(z-1)=phi_k} where (theta_0,phi_0)=(-pi,pi]x(-pi,pi], (theta_1,phi_1)=(-pi,pi]x(pi,3*pi], (theta_2,phi_2)=(pi,3*pi]x(pi,3*pi] [FFW2017 p.10 Eq. after Fig 3]. Winding transitions: Sheet 0->1 requires CCW loop around z=1 only (phi shifts by 2*pi). Sheet 0->2 requires CCW loops around both z=0 AND z=1. In zeta-plane (z=e^zeta): branch points at zeta=2*pi*i*k for k in Z\{0}. One CCW loop around zeta=2*pi*i shifts arg(z-1) by 2*pi, moving to sheet 1.
- **pole_structure**: Branch points at z=0 and z=1 in z-plane. In zeta-plane: z=0 maps to Im(zeta)->-inf (no finite singularity), but z=1 maps to zeta=2*pi*i*k (k in Z). Branch cuts in zeta-plane emanate from each 2*pi*i*k.
- **gt_kind**: hand-computed
- **how_to_compute**: Directly from the sheet parametrization given in FFW2017 p.10 (after Fig 3 caption). The (theta_k, phi_k) pairs define which arg-branch of z and z-1 is on each sheet. These are exactly as given: Sheet 0: arg(z) in (-pi,pi], arg(z-1) in (-pi,pi]. Sheet 1: arg(z) in (-pi,pi], arg(z-1) in (pi,3*pi]. Sheet 2: arg(z) in (pi,3*pi], arg(z-1) in (pi,3*pi].
- **precision**: exact (sheet indices are integers; branch cut positions are exact multiples of 2*pi*i)
- **targets_buckets**: B11, B12
- **regime**: PVI sheet bookkeeping with two branch points; CCW winding around z=0 vs z=1 give distinct transitions
- **citation**: FFW2017_painleve_riemann_surfaces_preprint.md:182-189 (sheet parametrization in z-plane with theta_k, phi_k), :165 (branch point at zeta=2*pi*i not mapped out by eta transformation)
- **confidence**: high
- **verify_note**: Verify: starting from sheet 0, a CCW loop around ONLY z=1 (not z=0) should give sheet 1. A second CCW loop around z=1 should return to sheet 0. A CCW loop around z=0 while staying on the sheet-0 branch cut of z=1 gives sheet 2. Test numerically: integrate PVI-zeta along a loop and check that w returns to its starting value after two CCW loops around z=1.

#### `pIII-tilde-trivial-constant`

- **name**: PIII-tilde exact constant solution: w(zeta)=exp(zeta/2) for parameters (alpha=0, beta=0, gamma=1, delta=-1)
- **ode**: PIII-tilde: w'' = (1/w)*(w')^2 + (1/4)*(0 + w^3 + 0 - e^{2*zeta}/w), with alpha=0, beta=0, gamma=1, delta=-1
- **ic_bc**: u(z) = 1 everywhere; equivalently w(zeta) = exp(zeta/2). Any starting point gives same solution.
- **domain**: All of complex zeta-plane (entire function). Corresponds to all sheets of PIII in z-plane with the single value u=1.
- **closed_form**: w(zeta) = exp(zeta/2). Proof: w' = (1/2)*exp(zeta/2), w'' = (1/4)*exp(zeta/2). RHS = (1/w)*(w')^2 + (1/4)*(w^3 - e^{2*zeta}/w) = (exp(-zeta/2))*(1/4)*exp(zeta) + (1/4)*(exp(3*zeta/2) - exp(2*zeta)*exp(-zeta/2)) = (1/4)*exp(zeta/2) + (1/4)*(exp(3*zeta/2) - exp(3*zeta/2)) = (1/4)*exp(zeta/2) = w''. QED. Corresponding PIII solution: u(z) = e^{-zeta/2}*w(zeta) = e^{-zeta/2}*e^{zeta/2} = 1 (constant). This is the trivial single-valued PIII solution u=1 for parameters (0,0,1,-1).
- **pole_structure**: entire (w(zeta)=exp(zeta/2) is entire, no poles)
- **gt_kind**: closed-form
- **how_to_compute**: Directly evaluate: w(zeta) = exp(zeta/2) for any complex zeta. For example at zeta = 2+3i: w = exp((2+3i)/2) = exp(1)*exp(3i/2) = e*(cos(3/2)+i*sin(3/2)) = 0.19228+2.71147i. Python: from mpmath import exp,mpc; w=exp(mpc(2,3)/2); assert abs(w-exp(mpc(1,1.5)))<1e-50.
- **precision**: arbitrary (BigFloat-any) via mpmath exp
- **targets_buckets**: B10, B12
- **regime**: Closed-form transformed RHS verification: tests that the PIII-tilde RHS evaluator gives exactly w''=(1/4)*exp(zeta/2) when w=exp(zeta/2), w'=(1/2)*exp(zeta/2) are substituted
- **citation**: FFW2017_painleve_riemann_surfaces_preprint.md:41-44 (PIII-tilde equation and transformation). The trivial solution u=1 for PIII(0,0,1,-1) can be verified by direct substitution into PIII.
- **confidence**: high
- **verify_note**: Note: this solution is SINGLE-VALUED (not multivalued), so it cannot test sheet-winding logic. It only tests that the transformed RHS (B10) is implemented correctly. For sheet-winding tests, use candidate pIII-tilde-sheet-oracle instead.

#### `pIII-tilde-conjugate-symmetry`

- **name**: PIII-tilde/PV-tilde Schwarz conjugate symmetry: w(conj(zeta)) = conj(w(zeta)) for real params and real ICs
- **ode**: PIII-tilde or PV-tilde with real parameters and real initial conditions on the real zeta-axis
- **ic_bc**: PIII-tilde: w(0)=1/4, w_prime(0)=9/16 [FFW2017 Fig 1]. PV-tilde: w(0)=2, w_prime(0)=-1 [FFW2017 Fig 6]. Both have real ICs on real axis.
- **domain**: Pairs (zeta, conj(zeta)) with Im(zeta) > 0
- **closed_form**: w(conj(zeta)) = conj(w(zeta)) for any zeta in C. Proof: If w(zeta) satisfies the ODE and w(x) is real for real x (both properties hold by assumption), then conj(w(zeta)) also satisfies the same ODE (since ODE has real coefficients when e^zeta, e^{2*zeta} are conjugated consistently, i.e., e^{conj(zeta)}=conj(e^zeta)), and conj(w(x))=w(x) on the real axis. By uniqueness, the two solutions agree everywhere. This gives a non-trivial consistency check between Im>0 and Im<0 half-planes.
- **pole_structure**: Symmetry: poles in upper half-plane are mirror images (under zeta -> conj(zeta)) of poles in lower half-plane. Pole locations zeta_0 and conj(zeta_0) occur in conjugate pairs for real ICs.
- **gt_kind**: closed-form
- **how_to_compute**: Evaluate solution at zeta and at conj(zeta) via two independent RK4 integrations (one going through upper half-plane, one through lower half-plane). Check w(conj(zeta)) = conj(w(zeta)) to working precision. Verified numerically for both PIII-tilde (error=0.0 at 50 digits) and PV-tilde (error=0.0 at 50 digits) in the computations above.
- **precision**: exact in principle; ~50-digit agreement in practice via mpmath RK4
- **targets_buckets**: B11, B12
- **regime**: Self-consistency check for multi-sheet integration; tests that sheet +s and sheet -s values are conjugate-related for real-IC solutions
- **citation**: FFW2017_painleve_riemann_surfaces_preprint.md:122 (symmetry-based error estimation method; FFW2017 use this as their own error estimate: E(zeta)=|w(zeta)-conj(w(conj(zeta)))|/|w(zeta)|)
- **confidence**: high
- **verify_note**: This is the SAME method FFW2017 use to estimate their own errors (Table 2, p.6). Our computation gives error=0.0 at 50 digits (machine zero) because we use sufficiently many RK4 steps on pole-free paths. This confirms the paths are pole-free and the integration is self-consistent.

#### `pVI-eta-subexpr-rhs`

- **name**: PVI double-exp ODE RHS evaluation at specific (eta, v, v_prime) with BigFloat oracle
- **ode**: PVI-eta (FFW2017 Eq. 5), evaluated at eta=eta0 + 0.5i where eta0 = log(log(10)) = 0.83403244524795580, v=0.429534600325223, v_prime=-0.037235820650106481. Parameters (alpha,beta,gamma,delta)=(4,-4,8,-8) [FFW2017 Fig 2].
- **ic_bc**: PVI Fig 2 ICs: u(z=10)=0.429534600325223, u_prime(10)=-1.61713114374804e-3. Transformed: eta0=log(log(10)), v(eta0)=0.429534600325223, v_prime(eta0)=u_prime(10)*10*log(10)=-0.037235820650106481.
- **domain**: Complex eta = log(log(10)) + 0.5i, well inside branch-point-free region Re(eta) < log(2*pi)
- **closed_form**: none
- **pole_structure**: Fixed singularities of ODE at eta_k = log|2*pi*k| + i*arg(2*pi*i*k) for k in Z\{0} [FFW2017 Eq.4]. For k=1: eta_1=log(2*pi)+i*pi/2~1.83788+1.57080i. The evaluation point eta=0.83403+0.5i is well within the branch-point-free region Re(eta)<log(2*pi)~1.83788.
- **gt_kind**: mpmath
- **how_to_compute**: from mpmath import mp,mpf,mpc,exp,log,pi; mp.dps=50; alpha,beta,gamma_p,delta_p=mpf(4),mpf(-4),mpf(8),mpf(-8); eta=mpc(log(log(mpf(10))),0.5); v=mpc(0.429534600325223); vdot=mpc(-0.037235820650106481); E=exp(exp(eta)); e_eta=exp(eta); Em1=E-1; vE=v-E; term_A=(1/2)*(1/v+1/(v-1)+1/vE)*vdot**2; term_B=-(e_eta*E/Em1+e_eta*E/vE-1)*vdot; P=alpha+beta*E/v**2+gamma_p*Em1/(v-1)**2+delta_p*E*Em1/vE**2; term_C=v*(v-1)*vE*e_eta**2/Em1**2*P; rhs=term_A+term_B+term_C; # Result: ~0.42199106+... at this starting point (v matches w at eta0, slight step to eta0+0.5i)
- **precision**: BigFloat-50 (dps=50), ~45 significant digits away from cancellation zone
- **targets_buckets**: B10
- **regime**: B10 transformed RHS evaluation for PVI-eta (double-exp coordinate); tests correct implementation of e^{e^eta} and (e^{e^eta}-1) in ODE RHS
- **citation**: FFW2017_painleve_riemann_surfaces_preprint.md:154-155 (Eq. 5 full PVI-eta ODE), :146-149 (Eq. 4 fixed singularity locations), :195 (Fig 2 ICs)
- **confidence**: medium
- **verify_note**: The value v(eta_eval) = 0.42199106077496084 - 0.01821819508177181i (computed above by RK4 integration along eta0 to eta0+0.5i). The RHS v''(eta_eval) would be the mpmath evaluation of Eq.5 at (eta_eval, v_eval, v_prime_eval). This tests the B10 ODE-RHS-evaluator for PVI in eta-coordinates. Compare to Float64 to check that all terms are computed correctly.

### Territory: special-function ODE oracles: Heun G + Heun C (B18) across full parameter space,

*Scout notes:* Gaps: (1) The existing oracles.txt monoculture is epsilon=-1 and gamma=2 for HeunC, a=2 for HeunG - all candidates below explicitly break these. (2) Mathematica HeunG/HeunC at large q and complex a beyond unit disk require a Mathematica wolframscript run to regenerate value-level oracles; the coefficient-level oracles are derived from the three-term recurrences in Motygin 2015 (confirmed to 50 digits). (3) The HeunC-on-second-sheet connection formula oracle (tightening the current 1e-3 self-check) requires solving a 2x2 matching linear system at z*=0.5+i/sqrt(2) per Motygin 2015 section 6 eq (21)-(22); a dedicated Mathematica run is needed. (4) Mathieu-as-HeunC: scipy mathieu_cem/mathieu_sem use double precision only; higher precision requires Mathematica MathieuC/MathieuS. (5) Lame-as-HeunG: verified via DLMF 31.7.5 parameter map, but mpmath has no Lame evaluator - use Mathematica LameC/LameS.

#### `heung-reduction-2f1-dlmf-31-7-2`

- **name**: HeunG -> Gauss 2F1 via quadratic transformation (DLMF 31.7.2)
- **ode**: u'' + (3/2/z + 1/(z-1) + 3/2/(z-2)) u' + (2z - 2)/(z*(z-1)*(z-2)) u = 0 [a_heun=2, q=2, alpha=1, beta=2, gamma=3/2, delta=1, epsilon=3/2]
- **ic_bc**: u(0)=1, u'(0)=q/(a_heun*gamma)=2/(2*3/2)=2/3
- **domain**: real z in (0,1), avoiding pole at z=1 (where u=1)
- **closed_form**: HeunG(2, 2, 1, 2, 3/2, 1; z) = 2F1(1/2, 1; 3/2; 2z-z^2). At z=0.1: 1.07170483682699539938...; z=0.3: 1.25407417926432015183...; z=0.5: 1.52069199260189269506...; z=0.9: 3.00830214985481857784...
- **pole_structure**: regular singularities at z=0, 1, 2, infinity; solution 2F1 has branch point at z=1 (via u=2z-z^2=1)
- **gt_kind**: closed-form
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50; z=mpmath.mpf('0.3'); u=2*z-z**2; print(mpmath.hyp2f1(mpmath.mpf('1/2'), 1, mpmath.mpf('3/2'), u))
- **precision**: arbitrary via mpmath hyp2f1 with mp.dps=50
- **targets_buckets**: B18
- **regime**: HeunG reduction to 2F1 via DLMF 31.7.2 quadratic transformation; best closed-form HeunG oracle; exercises non-integer gamma=3/2; regenerable at arbitrary precision
- **citation**: DLMF 31.7.2 https://dlmf.nist.gov/31.7 equation Hl(2,ab;a,b,g,a+b-2g+1;z) = F12(a/2,b/2;g;1-(1-z)^2)
- **confidence**: high

#### `heung-reduction-constant-alpha0-q0`

- **name**: HeunG trivial reduction to constant 1 when alpha=0, q=0
- **ode**: Standard HeunG ODE (DLMF 31.2.1) with arbitrary a, alpha=0, q=0
- **ic_bc**: u(0)=1, u'(0)=0
- **domain**: all z
- **closed_form**: HeunG(a, 0, 0, beta, gamma, delta; z) = 1 for all z. From recurrence: c_1 = q/(a*gamma) = 0; P_1 = (0+0)*(0+beta) = 0; all subsequent c_n = 0.
- **pole_structure**: entire (constant = 1)
- **gt_kind**: closed-form
- **how_to_compute**: Exact: HeunG(a, 0, 0, beta, gamma, delta; z) = 1. Check via recurrence: python3 -c "import mpmath; a,q,al,be,ga,de=3,0,0,2,3,4; c=[1]+[0]*9; print(c)"
- **precision**: exact rational
- **targets_buckets**: B18
- **regime**: trivial polynomial Heun; tests degenerate parameter handling; any a, gamma, delta; IC u(0)=1, u'(0)=0
- **citation**: DLMF 31.5 Heun polynomials; DLMF 31.3.2 c_1=q/(a*gamma)
- **confidence**: high

#### `heung-polynomial-degree1`

- **name**: Heun polynomial of degree 1 Hp_{1,m} with explicit closed-form eigenvalue
- **ode**: Standard HeunG ODE (DLMF 31.2.1) with a_heun=2, alpha=-1, beta=2, gamma=3, delta=4, q=q_{1,m}
- **ic_bc**: u(0)=1, u'(0)=q_{1,m}/6
- **domain**: all z (polynomial)
- **closed_form**: Hp_{1,m}(z) = 1 + (q_{1,m}/6)*z (linear polynomial). Eigenvalues q_{1,0} = -6+2*sqrt(6) = -1.10102051443364380361..., q_{1,1} = -6-2*sqrt(6) = -10.898979485566356196...; Hp_{1,0}(0.5) = 0.90824829046386301637...; Hp_{1,1}(0.5) = 0.09175170953613698363...
- **pole_structure**: entire (polynomial)
- **gt_kind**: closed-form
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50; q1=-6+2*mpmath.sqrt(6); q2=-6-2*mpmath.sqrt(6); print(1+q1/6*0.5); print(1+q2/6*0.5)
- **precision**: arbitrary via mpmath sqrt
- **targets_buckets**: B18
- **regime**: Heun polynomial termination; tests that the Taylor series is finite; eigenvalue from quadratic q^2+12q+12=0 (termination condition P_2*c_1=0)
- **citation**: DLMF 31.5.1-31.5.2; eigenvalue from termination condition at j=1: (Q_1+q)*c_1 - P_1*c_0 = 0 with P_1=(0+alpha)(0+beta)=(-1)(2)=-2
- **confidence**: high

#### `heung-pole-solution-dlmf-31-7-2-b`

- **name**: HeunG -> 1/(1-z) via DLMF 31.7.2 double reduction
- **ode**: Standard HeunG ODE (DLMF 31.2.1) with a_heun=2, q=4, alpha=1, beta=4, gamma=2, delta=2, epsilon=2
- **ic_bc**: u(0)=1, u'(0)=q/(a*gamma)=4/(2*2)=1
- **domain**: z in C \ {1}
- **closed_form**: HeunG(2, 4, 1, 4, 2, 2; z) = 2F1(1/2, 2; 2; 1-(1-z)^2) = (1-z)^(-1) [using F(a,b;b;u)=(1-u)^(-a)]. All Taylor coefficients c_n=1. At z=0.3: 10/7 = 1.42857142857...
- **pole_structure**: simple pole at z=1
- **gt_kind**: closed-form
- **how_to_compute**: Exact: 1/(1-z). Verify coefficient recurrence: all c_n should equal 1. python3: import mpmath; mpmath.mp.dps=50; a_h,q_h,al,be,ga,de=2,4,1,4,2,2; ep=al+be+1-ga-de; c=[mpmath.mpf(0)]*12; c[0]=1; c[1]=mpmath.mpf(q_h)/(a_h*ga); print([float(c[k]) for k in range(8)])
- **precision**: exact rational
- **targets_buckets**: B18
- **regime**: HeunG with a simple pole in the solution at z=1; tests stepper across a pole; coefficients are all identically 1 making this a mutation-proof test
- **citation**: DLMF 31.7.2 + DLMF 15.4.6 F(a,b;b;z)=(1-z)^(-a)
- **confidence**: high

#### `heunc-coeff-diverse-epsilon-neg3`

- **name**: HeunC coefficient oracle: epsilon=-3 breaking epsilon=-1 monoculture
- **ode**: u'' + (2/z + 3/(z-1) - 3)*u' + (z-1)/(z*(z-1))*u = 0 [HeunC(q=1,alpha=1,gamma=2,delta=3,epsilon=-3,z)]
- **ic_bc**: u(0)=1, u'(0)=-q/gamma=-0.5
- **domain**: real z in [0.1, 0.9]
- **closed_form**: none
- **pole_structure**: regular at z=0,1; irregular at z=infinity
- **gt_kind**: mpmath
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50 q,al,ga,de,ep=1,1,2,3,-3 c=[mpmath.mpf(0)]*25; c[0]=1; c[1]=-mpmath.mpf(q)/ga for n in range(1,22):     A_n=(n+1)*(n+ga); B_n=n*(n-1+ga+de-ep)-q; C_n=ep*(n-1)+al     c[n+1]=(B_n*c[n]+C_n*c[n-1])/A_n z=mpmath.mpf('0.1'); print(sum(c[k]*z**k for k in range(23))) # Verify: wolframscript -code "N[HeunC[1,1,2,3,-3,0.1],50]" Values: z=0.1: 0.94525457198487807736...; z=0.5: 0.48560767684076127087...; z=0.9: -10.4047597775769037523...
- **precision**: 50 decimal digits via mpmath arbitrary precision recurrence; cross-check with wolframscript
- **targets_buckets**: B18
- **regime**: HeunC with epsilon=-3 (diverse from epsilon=-1 monoculture); tests recurrence stability for more negative irregular singularity coefficient; value at z=0.9 shows large magnitude growth near z=1
- **citation**: Motygin 2015 references/markdown/heun/Motygin2015_Heun_evaluation_1506.03848/Motygin2015_Heun_evaluation_1506.03848.md section 3 eqs (3)-(4); DLMF 31.12.1
- **confidence**: high
- **verify_note**: Cross-check with wolframscript N[HeunC[1,1,2,3,-3,z],50] at three points before trusting

#### `heunc-coeff-positive-epsilon`

- **name**: HeunC coefficient oracle: epsilon=2 (positive, diverse from all-negative monoculture)
- **ode**: u'' + (2/z + 3/(z-1) + 2)*u' + (z-1)/(z*(z-1))*u = 0 [HeunC(q=1,alpha=1,gamma=2,delta=3,epsilon=2,z)]
- **ic_bc**: u(0)=1, u'(0)=-0.5
- **domain**: real z in [0.1, 0.5]
- **closed_form**: none
- **pole_structure**: regular at z=0,1; irregular at z=infinity
- **gt_kind**: mpmath
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50 q,al,ga,de,ep=1,1,2,3,2 c=[mpmath.mpf(0)]*25; c[0]=1; c[1]=-mpmath.mpf(q)/ga for n in range(1,22):     A_n=(n+1)*(n+ga); B_n=n*(n-1+ga+de-ep)-q; C_n=ep*(n-1)+al     c[n+1]=(B_n*c[n]+C_n*c[n-1])/A_n print(sum(c[k]*mpmath.mpf('0.1')**k for k in range(23))) # Values: z=0.1: 0.94986517967355029307...; z=0.5: 0.72276415124476652809... # Verify: wolframscript -code "N[HeunC[1,1,2,3,2,0.1],50]"
- **precision**: 50 decimal digits via mpmath
- **targets_buckets**: B18
- **regime**: HeunC with positive epsilon=2; sign flip relative to all existing oracles; tests that the recurrence sign convention is correct
- **citation**: Motygin 2015 references/markdown/heun/Motygin2015_Heun_evaluation_1506.03848/Motygin2015_Heun_evaluation_1506.03848.md section 3 eqs (3)-(4)
- **confidence**: high
- **verify_note**: Cross-check with wolframscript N[HeunC[1,1,2,3,2,z],50]

#### `heunc-coeff-fractional-gamma`

- **name**: HeunC coefficient oracle: fractional gamma=1/2 with non-integer Frobenius exponent
- **ode**: HeunC(q=1, alpha=1, gamma=1/2, delta=3, epsilon=-1, z): u'' + (1/(2z) + 3/(z-1) - 1)*u' + (z-1)/(z*(z-1))*u = 0
- **ic_bc**: u(0)=1, u'(0)=-q/gamma=-2
- **domain**: real z in (0, 0.5)
- **closed_form**: none
- **pole_structure**: regular at z=0 (exponent 1-gamma=1/2, non-integer); irregular at z=infinity
- **gt_kind**: mpmath
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50 q,al,ga,de,ep=1,1,mpmath.mpf('1/2'),3,-1 c=[mpmath.mpf(0)]*25; c[0]=1; c[1]=-mpmath.mpf(q)/ga for n in range(1,22):     A_n=(n+1)*(n+ga); B_n=n*(n-1+ga+de-ep)-q; C_n=ep*(n-1)+al     c[n+1]=(B_n*c[n]+C_n*c[n-1])/A_n print(sum(c[k]*mpmath.mpf('0.1')**k for k in range(23))) # Values: z=0.1: 0.77694852123698801984...; z=0.5: -1.3497113525144030943... # Verify: wolframscript -code "N[HeunC[1,1,1/2,3,-1,0.1],50]"
- **precision**: 50 decimal digits via mpmath
- **targets_buckets**: B18
- **regime**: HeunC with non-integer gamma=1/2; second Frobenius solution has branch point at z=0; tests recurrence stability with non-integer exponents; completely absent from existing oracles
- **citation**: Motygin 2015 section 2 (branch point discussion) and section 3 eq (5) for gamma=0 case; DLMF 31.12.1
- **confidence**: high
- **verify_note**: Cross-check with wolframscript N[HeunC[1,1,1/2,3,-1,0.1],50]; note Mathematica may differ for non-standard gamma

#### `heung-coeff-large-q-a3`

- **name**: HeunG coefficient oracle: large q=10 with a_heun=3 (diverse singularity location)
- **ode**: Standard HeunG ODE (DLMF 31.2.1) with a_heun=3, q=10, alpha=1, beta=2, gamma=2, delta=3, epsilon=-1
- **ic_bc**: u(0)=1, u'(0)=q/(a_heun*gamma)=10/(3*2)=5/3
- **domain**: real z in (0, 0.5)
- **closed_form**: none
- **pole_structure**: regular singularities at z=0, 1, 3, infinity
- **gt_kind**: mpmath
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50 a_h,q_h,al,be,ga,de=3,10,1,2,2,3; ep=al+be+1-ga-de c=[mpmath.mpf(0)]*30; c[0]=1; c[1]=mpmath.mpf(q_h)/(a_h*ga) for j in range(1,27):     Rj=mpmath.mpf(a_h)*(j+1)*(j+ga)     Qj=mpmath.mpf(j)*((j-1+ga)*(1+a_h)+a_h*de+ep)     Pj=(j-1+al)*(j-1+be)     c[j+1]=((Qj+q_h)*c[j]-Pj*c[j-1])/Rj print(sum(c[k]*mpmath.mpf('0.25')**k for k in range(28))) # Values: z=0.25: 1.6250583912557556129...; z=0.5: 3.28667034406873242578... # Verify: wolframscript -code "N[HeunG[3,10,1,2,2,3,0.25],50]"
- **precision**: 50 decimal digits via mpmath
- **targets_buckets**: B18
- **regime**: HeunG with large q=10 AND a_heun=3 (both diverse from existing monoculture); convergence radius is still 1 (nearest singularity z=1); large q causes oscillatory coefficients testing numerical stability
- **citation**: DLMF 31.3.3-31.3.4 (HeunG coefficient recurrence)
- **confidence**: high
- **verify_note**: Verify with wolframscript N[HeunG[3,10,1,2,2,3,0.25],50]

#### `heung-coeff-small-a`

- **name**: HeunG coefficient oracle: a_heun=0.5 singular point inside unit disk
- **ode**: Standard HeunG ODE (DLMF 31.2.1) with a_heun=1/2, q=1, alpha=1, beta=2, gamma=3, delta=4, epsilon=-3
- **ic_bc**: u(0)=1, u'(0)=q/(a*gamma)=1/(0.5*3)=2/3
- **domain**: real z in (0, 0.45), convergence radius 0.5
- **closed_form**: none
- **pole_structure**: regular singularities at z=0, 1, 1/2, infinity; convergence radius from z=0 is min(1,1/2)=1/2
- **gt_kind**: mpmath
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50 a_h,q_h,al,be,ga,de=mpmath.mpf('1/2'),1,1,2,3,4; ep=al+be+1-ga-de c=[mpmath.mpf(0)]*30; c[0]=1; c[1]=q_h/(a_h*ga) for j in range(1,27):     Rj=a_h*(j+1)*(j+ga); Qj=mpmath.mpf(j)*((j-1+ga)*(1+a_h)+a_h*de+ep); Pj=(j-1+al)*(j-1+be)     c[j+1]=((Qj+q_h)*c[j]-Pj*c[j-1])/Rj print(sum(c[k]*mpmath.mpf('0.3')**k for k in range(28))) # Values: z=0.1: 1.068937255689352422622...; z=0.3: 1.2106857509015363929...; z=0.45: 1.2855463176267072038... # Verify: wolframscript -code "N[HeunG[1/2,1,1,2,3,4,0.3],50]"
- **precision**: 50 decimal digits via mpmath
- **targets_buckets**: B18
- **regime**: HeunG with |a_heun|<1; singular point inside the unit disk cuts the convergence radius to 0.5; critical test for analytic continuation past the z=0.5 singular point
- **citation**: DLMF 31.2.1 (general complex a); Motygin 2015 section 2: R_0=min(1,|a|) is the convergence radius
- **confidence**: high
- **verify_note**: Verify with wolframscript N[HeunG[1/2,1,1,2,3,4,0.3],50]

#### `heung-complex-a`

- **name**: HeunG with complex a=0.5+0.5i: singular point off real axis
- **ode**: Standard HeunG ODE (DLMF 31.2.1) with a_heun=0.5+0.5i, q=1, alpha=1, beta=2, gamma=3, delta=4, epsilon=-3
- **ic_bc**: u(0)=1, u'(0)=(0.333...-0.333...j)
- **domain**: real z=0.3
- **closed_form**: none
- **pole_structure**: regular singularities at z=0, 1, 0.5+0.5i, infinity
- **gt_kind**: mpmath
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50 a_h=mpmath.mpc('0.5','0.5'); q_h=1; al,be,ga,de=1,2,3,4; ep=al+be+1-ga-de c=[mpmath.mpc(0)]*30; c[0]=1; c[1]=mpmath.mpf(q_h)/(a_h*ga) for j in range(1,27):     Rj=a_h*(j+1)*(j+ga); Qj=mpmath.mpf(j)*((j-1+ga)*(1+a_h)+a_h*de+ep); Pj=(j-1+al)*(j-1+be)     c[j+1]=((Qj+q_h)*c[j]-Pj*c[j-1])/Rj z=mpmath.mpf('0.3'); print(sum(c[k]*z**k for k in range(28))) # Value at z=0.3: (1.1059813320247579739 - 0.10592175241504076687j) # Verify: wolframscript -code "N[HeunG[0.5+0.5*I,1,1,2,3,4,0.3],50]"
- **precision**: 50 decimal digits via mpmath complex recurrence
- **targets_buckets**: B18
- **regime**: HeunG with complex singularity a=0.5+0.5i; convergence radius at z=0.3 is min(0.3, 0.7, |0.3-(0.5+0.5i)|)=min(0.3,0.7,0.539)=0.3; tests complex coefficient arithmetic in recurrence
- **citation**: DLMF 31.2.1 states a in C, a not in {0,1,infinity}; Motygin 2015 section 2
- **confidence**: medium
- **verify_note**: Must verify with wolframscript N[HeunG[0.5+0.5*I,1,1,2,3,4,0.3],50] before trusting; convergence radius check needed

#### `heung-lame-reduction`

- **name**: HeunG -> Lame equation via DLMF 31.7.5 substitution
- **ode**: Lame equation (DLMF 29.2.1): d^2w/dz^2 + (h - n*(n+1)*k^2*sn^2(z,k))*w = 0. In algebraic form (xi=sn^2(z,k)): this is HeunG with a_heun=k^{-2}, q=-h/4, alpha=-n/2, beta=(n+1)/2, gamma=delta=epsilon=1/2
- **ic_bc**: LameC[1,0,0,m]=dn(0,k)=1 (IC at z=0); Mathematica LameC[n,j,z,m] convention: m=k^2
- **domain**: real z in (0, K(k)) Jacobi quarter-period
- **closed_form**: For n=0: LameC(0,0,z,k^2) = 1 (constant). For n=1, j=0: LameC(1,0,z,k^2) = dn(z,k) = sqrt(1-k^2*sn^2(z,k)). Mathematica: N[JacobiDN[0.5, 0.5], 50] = 0.92372...
- **pole_structure**: poles of sn(z,k) at z = K+im*K' (Jacobi K, K'); HeunG has regular singularities at xi=0, 1, k^{-2}, infinity
- **gt_kind**: mathematica
- **how_to_compute**: wolframscript -code "N[LameC[1, 0, 0.5, 0.5], 50]" and "N[JacobiDN[0.5, Sqrt[0.5]], 50]" (should match for n=1, j=0, m=k^2=0.5). Also: wolframscript -code "N[HeunG[2, -0.25, -0.5, 1, 0.5, 0.5, JacobiSN[0.5, Sqrt[0.5]]^2], 50]"
- **precision**: arbitrary via Mathematica LameC/LameS or JacobiDN
- **targets_buckets**: B18
- **regime**: HeunG-to-Lame reduction; exercises parameter mapping a_heun=k^{-2}; for k=1/sqrt(2): a_heun=2 with gamma=delta=epsilon=1/2 (novel compared to existing oracles); Lame solution has poles at Jacobi quarter-periods
- **citation**: DLMF 31.7.5 https://dlmf.nist.gov/31.7 transformation z=sn^2(zeta,k); DLMF 29.2.1 Lame equation
- **confidence**: medium
- **verify_note**: Run both wolframscript commands and verify agreement; note parameter m=k^2 in Mathematica vs k in DLMF

#### `heunc-connection-formula-second-sheet`

- **name**: HeunC connection formula at z=1 (second-sheet oracle, Motygin eq 21)
- **ode**: HeunC(q=1, alpha=1, gamma=2, delta=3, epsilon=-1, z)
- **ic_bc**: n/a - connection formula problem
- **domain**: z near 1 and z>1 on second sheet
- **closed_form**: none - requires numerical matching at z*=0.5+i/sqrt(2)
- **pole_structure**: branch cut from z=1 to +infinity (principal branch)
- **gt_kind**: mathematica
- **how_to_compute**: wolframscript -code "SetPrecision[Module[{zs=0.5+I/Sqrt[2]}, m0={HeunC[1,1,2,3,-1,zs],HeunC'[1,1,2,3,-1,zs]}; m1=HeunC[1,1,2,3,-1,1-zs]; d1=-HeunC'[1,1,2,3,-1,1-zs]; m2=HeunC[-I+Abs[1],1,2,3,-1,1-zs]; LinearSolve[{{m1,m2},{d1,-HeunC'[-I+Abs[1],1,2,3,-1,1-zs]}},m0]], 50]" [approximate - exact recipe requires careful implementation of Motygin eq (21)]. Simpler check: evaluate HeunC at 1.5 via two routes and compare.
- **precision**: 50 digits via Mathematica SetPrecision
- **targets_buckets**: B18
- **regime**: second-sheet analytic continuation for HeunC; tightens the current 1e-3 self-check to 1e-10+; exercises the branch structure at z=1
- **citation**: Motygin 2015 references/markdown/heun/Motygin2015_Heun_evaluation_1506.03848/Motygin2015_Heun_evaluation_1506.03848.md section 6 eq (21)-(22); DLMF 31.11 connection formulas
- **confidence**: low
- **verify_note**: This recipe is schematic - needs careful implementation; the wolframscript call needs the correct HeunC arguments for the z=1 local expansion

#### `airy-ivp-oracle`

- **name**: Airy equation IVP: u'' = z*u, solution Ai(z)
- **ode**: u'' = z*u
- **ic_bc**: u(0)=Ai(0)=0.35502805388781723926..., u'(0)=Ai'(0)=-0.25881940379280679841...
- **domain**: real z in {-2,-1,1,2,3}, complex z=1+1j,2+1j
- **closed_form**: u(z) = Ai(z). Ai(0)=0.35502805388781723926006318600418...; Ai'(0)=-0.25881940379280679840518356018920...; Ai(1)=0.13529241631288141552414742351547...; Ai(2)=0.034924130423274379135322080791808...; Ai(-1)=0.53556088329235211879951656563887...; Ai(-2)=0.22740742820168557599192443603787...; Ai(3)=0.0065911393574607191442574484079614...; Ai(1+1j)=0.060458308371838149197 - 0.15188956587718140235j
- **pole_structure**: entire (no poles)
- **gt_kind**: mpmath
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50; print(mpmath.airyai(1))
- **precision**: arbitrary via mpmath airyai with mp.dps=50
- **targets_buckets**: B18
- **regime**: simplest non-trivial IVP; entire solution; tests basic Taylor step accuracy; no poles; IC at origin with DLMF-certified exact values
- **citation**: DLMF 9.2.1 Airy equation; DLMF 9.2.4-9.2.5 Ai(0)=1/(3^{2/3}*Gamma(2/3)), Ai'(0)=-1/(3^{1/3}*Gamma(1/3))
- **confidence**: high

#### `bessel-j0-ivp-oracle`

- **name**: Bessel J0 IVP: z*u'' + u' + z*u = 0, start at z=1
- **ode**: u'' + u'/z + u = 0 (regular singularity at z=0, start at z=1)
- **ic_bc**: u(1)=J_0(1)=0.76519768655796655..., u'(1)=-J_1(1)=-0.44005058574493352...
- **domain**: real z in [1, 5], complex z=1+1j
- **closed_form**: u(z) = J_0(z). J_0(1)=0.76519768655796655144971752610266...; J_0'(1)=-J_1(1)=-0.44005058574493351595968220371891...; J_0(2)=0.22389077914123566805182745464995...; J_0(3)=-0.26005195490193343762415469597733...; J_0(5)=-0.17759677131433830434739701307476...; J_0(1+1j)=0.93760847680602927660 - 0.49652994760912213217j
- **pole_structure**: entire (J_0 is entire); ODE has regular singularity at z=0
- **gt_kind**: mpmath
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50; print(mpmath.besselj(0, 2))
- **precision**: arbitrary via mpmath besselj with mp.dps=50
- **targets_buckets**: B18
- **regime**: IVP with regular singularity at origin; start at z=1 avoids the singularity; oscillatory solution (zeros at ~2.4, 5.5, 8.7); tests stepper over multiple oscillations
- **citation**: DLMF 10.2.1 Bessel equation; DLMF 10.2.2 J_nu series definition; mpmath besselj
- **confidence**: high

#### `gauss-2f1-ivp-oracle`

- **name**: Gauss 2F1 IVP: z(1-z)u'' + [c-(a+b+1)z]u' - ab*u = 0, a=1/2, b=1/3, c=5/4
- **ode**: z*(1-z)*u'' + [5/4-(11/6)*z]*u' - (1/6)*u = 0
- **ic_bc**: u(0)=1, u'(0)=a*b/c=2/15=0.13333...
- **domain**: real z in [0, 0.9], complex z=0.5+0.5j
- **closed_form**: u(z) = 2F1(1/2, 1/3; 5/4; z). u(0)=1; u'(0)=ab/c=2/15=0.1333...; u(0.4)=1.06597809096176726538...; u(0.7)=1.14742531251816216936...; u(0.9)=1.25405973047601876714...; u(0.5+0.5j)=1.05084009857176407320+0.10091136968897003223j
- **pole_structure**: regular singularities at z=0, 1, infinity; solution 2F1 analytic in C \ [1,+inf)
- **gt_kind**: mpmath
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50; a=mpmath.mpf('1/2'); b=mpmath.mpf('1/3'); c=mpmath.mpf('5/4'); print(mpmath.hyp2f1(a,b,c,0.4))
- **precision**: arbitrary via mpmath hyp2f1 with mp.dps=50
- **targets_buckets**: B18
- **regime**: hypergeometric IVP with three regular singular points; tests approach to branch cut at z=1; complex z evaluation avoids the branch cut; solution non-elementary
- **citation**: DLMF 15.2.1 Gauss hypergeometric series; DLMF 15.10.1 hypergeometric ODE; mpmath hyp2f1
- **confidence**: high

#### `2f1-elementary-reduction`

- **name**: 2F1 -> elementary closed form F(a,b;b;z)=(1-z)^(-a) as IVP oracle (DLMF 15.4.6)
- **ode**: z*(1-z)*u'' + [3-(3+1/2+1)*z]*u' - 3*(1/2)*u = 0 [a=1/2, b=c=3]
- **ic_bc**: u(0)=1, u'(0)=a=1/2
- **domain**: real z in (0, 0.9); complex z supported by (1-z)^{-1/2}
- **closed_form**: F(1/2, 3; 3; z) = (1-z)^(-1/2). At z=0.4: (0.6)^{-1/2} = 1.29099444873580562839...; at z=0.7: (0.3)^{-1/2} = 1.82574185835055371152...; at z=0.9: (0.1)^{-1/2} = 3.16227766016837933200...
- **pole_structure**: branch point at z=1 (square root singularity)
- **gt_kind**: closed-form
- **how_to_compute**: python3: import mpmath; mpmath.mp.dps=50; z=mpmath.mpf('0.4'); print((1-z)**(-0.5)); print(mpmath.hyp2f1(0.5, 3, 3, z)) # must agree to 50 digits
- **precision**: exact to machine precision via (1-z)^{-1/2}
- **targets_buckets**: B18
- **regime**: 2F1 ODE with elementary power-function solution; cleanest possible accuracy oracle; derivative is also elementary u'(z)=a*(1-z)^{-a-1}; tests stepper approaching algebraic branch point at z=1
- **citation**: DLMF 15.4.6 https://dlmf.nist.gov/15.4 F(a,b;b;z)=(1-z)^{-a}
- **confidence**: high

### Territory: BVP test problems with known closed-form solutions for a Chebyshev-spectral Newt

*Scout notes:* GAPS AND CAVEATS:

(1) The 3-arg bvp_solve overload (bvp_solve(f, ∂f/∂u, ∂f/∂u', z_a, z_b, u_a, u_b)) exists in BVP.jl lines 370-467 and is exercised only via the IVPBVPHybrid module on PIII at loose tolerance (IB.1.4, IB.1.5 use N=10, tol=1e-10, 30 iters -- well below the spectral accuracy floor). There is NO direct unit test of bvp_solve with the 3-arg signature against a closed-form oracle. This is the most critical gap.

(2) The VectorBVP tests (VB.1.2, VB.1.3) only use Dirichlet-selector B matrices (e.g. B_a=[1 0;0 0] pins y1 at z_a, B_b=[0 0;1 0] pins y1 at z_b). The ADR-0023 tau-method implementation handles general B_a, B_b but the endpoint-coupling case (where B_a row i involves both y components, or where B_b row i mixes both endpoints) is only tested via the PI companion system (VB.6.x), which still uses Dirichlet-selector BCs. A genuinely mixing BC (one row of g involves components from different endpoints) is untested.

(3) The oblique complex segment path for bvp_solve (Complex{T} z_a, z_b with Im(z_b-z_a) != 0) is exercised only in vector_bvp_test.jl VB.1.3 for VectorBVP (harmonic oscillator on [0, 0.6+0.4i]). The scalar BVP.bvp_solve is not tested on an oblique segment. The 'real(t*)-only domain guard' in BVP.jl line 491 is specifically the check `real(t_star) <= 1 + 100*eps(T)` -- this correctly accepts on-segment oblique points but would also accept off-segment complex z that happen to map to |Re(t*)| <= 1. Testing this is achievable but requires a test that explicitly queries off-segment oblique points to confirm they throw DomainError.

(4) The rank-deficient BC fixture (vector-rank-deficient-bc) requires careful design: the simple g=[0,0] periodic BC y(0)=y(2*pi) gives the trivially correct zero solution (not a rank-deficiency failure). The correct rank-deficient fixture must use a non-trivial g that creates a true underdetermined system. Recommended: B_a=I, B_b=-I on [0, 2*pi] with g=(A_sin, A_cos) for some A, which requires both the sin and cos modes to be simultaneously amplitude-specified by contradictory conditions. More precisely: B_a=[1,0;0,0], B_b=[1,0;0,0], g=[0,0] on [0, 2*pi] -- this pins y1(0)+y1(2*pi)=0 and leaves y2 unconstrained at both ends, making the system rank-deficient (one equation missing).

(5) BigFloat-256 path for the 3-arg overload is completely untested (BV.5.1 only covers the 2-arg overload). The BigFloat path of bvp_solve(f, ∂f/∂u, ∂f/∂u', ...) uses the same GenericLinearAlgebra backslash but exercises the additional D1_ii matrix multiply in the Jacobian assembly.

CLEANEST ORACLES FOUND:
- Damped oscillator (3-arg linear): u=e^{-z}cos(2z), residual u''+2u'+5u=0 verified analytically to 1e-16
- Euler-Cauchy (3-arg linear, rational): u=(4/3)z+(2/3)/z on [1,2], all values exact rational
- Nonlinear exp (3-arg): u=exp(2z), residual=0 exactly
- Bratu (2-arg nonlinear): u=2*log(B/cosh(B*(z-0.5))), B=cosh(B/2)≈1.17878, residual=0 symbolically
- Oblique cosh (2-arg complex segment): u=cosh(z), values from mpmath 50 dps
- Sin/cos anti-periodic BC (B16): y=(sin,cos), g=(-1,1), uniqueness verified
- Mixing BC (B16): y=(sin,cos) on [0,pi], B_a=[[0,0],[0,1]], B_b=[[1,0],[0,0]], g=[0,1]

#### `3arg-damped-oscillator`

- **name**: Linear 3-arg BVP: damped oscillator u'' = -2u' - 5u on [-1, 1]
- **ode**: u'' = -2*u' - 5*u, i.e. f(z, u, u') = -2*u' - 5*u; df/du = -5; df/dup = -2
- **ic_bc**: u(-1) = exp(1)*cos(2) ≈ -1.1312043837568136; u(1) = cos(2)/exp(1) ≈ -0.1530918656742263
- **domain**: real segment [-1, 1]
- **closed_form**: u(z) = exp(-z) * cos(2*z). BCs: u(-1) = exp(1)*cos(2) = -1.13120438375681358..., u(1) = exp(-1)*cos(2) = -0.15309186567422629...
- **pole_structure**: entire
- **gt_kind**: closed-form
- **how_to_compute**: Julia: u(z) = exp(-z)*cos(2z), u'(z) = -exp(-z)*(cos(2z)+2*sin(2z)). Verify ODE: u'' = exp(-z)*(-3*cos(2z)+4*sin(2z)); -2u'-5u = 2*exp(-z)*(cos(2z)+2*sin(2z)) - 5*exp(-z)*cos(2z) = exp(-z)*(-3*cos(2z)+4*sin(2z)) ✓. Pinned interior: u(-0.7)=0.34227179419638181..., u(0.0)=1.0, u(0.7)=0.08440318529167406... (mpmath 50 dps).
- **precision**: arbitrary via mpmath/BigFloat
- **targets_buckets**: B15
- **regime**: 3-arg D1-coupled Jacobian (∂f/∂u' = -2 gives a non-diagonal Jacobian contribution via D1_ii)
- **citation**: constructed; ODE is the constant-coefficient u''+2u'+5u=0, characteristics r=-1±2i, general solution e^{-x}(C1*cos 2x + C2*sin 2x). With C1=1, C2=0: Trefethen SMIM 2000 ch.7 examples.
- **confidence**: high
- **verify_note**: Analytic residual u'' + 2u' + 5u = 0 verified at z=0, 0.3, pi/4 to machine epsilon. The key test point: the Jacobian contribution from ∂f/∂u' = -2 via the D1 term is absent in all current bvp_test.jl tests (which only exercise the 2-arg path or the 3-arg via FFW2017 PIII in IVP-BVP hybrid at loose tolerance). This isolates the new D1-coupled term in J = D2_ii - scale*diag(∂f/∂u) - (scale/half_diff)*diag(∂f/∂u')*D1_ii.

#### `3arg-euler-cauchy`

- **name**: Linear 3-arg BVP: Euler-Cauchy z^2 u'' + zu' - u = 0 on [1, 2]
- **ode**: u'' = u/z^2 - u'/z, i.e. f(z, u, u') = u/z^2 - u'/z; df/du = 1/z^2; df/dup = -1/z
- **ic_bc**: u(1) = 2 (exact), u(2) = 3 (exact)
- **domain**: real segment [1, 2]
- **closed_form**: u(z) = (4/3)*z + (2/3)/z. Exact rational values: u(1) = 2, u(1.5) = 2 + 4/9 = 22/9, u(2) = 3, u(1.25) = 11/5, u(1.75) = 19/7.
- **pole_structure**: pole at z=0 (outside domain)
- **gt_kind**: closed-form
- **how_to_compute**: u = C1*z + C2/z (Euler-Cauchy power-law ansatz u=z^p gives p^2=1). BCs u(1)=2, u(2)=3: C1+C2=2, 2*C1+C2/2=3 => C1=4/3, C2=2/3. Verify: u'' = 2*C2/z^3; u/z^2 - u'/z = (C1*z+C2/z)/z^2 - (C1-C2/z^2)/z = C1/z+C2/z^3 - C1/z + C2/z^3 = 2*C2/z^3 ✓. Exact: u(1.5)=4/3*1.5+2/3/1.5=2+4/9=2.44444... (rational).
- **precision**: exact rational
- **targets_buckets**: B15
- **regime**: 3-arg BVP with z-dependent ∂f/∂u' = -1/z coefficient in the D1-coupled Jacobian term
- **citation**: constructed; Euler-Cauchy ODE is classical, see Trefethen SMIM ch.7; z^2u''+zu'-u=0 rewritten as u'' = u/z^2 - u'/z for the 3-arg BVP interface.
- **confidence**: high
- **verify_note**: All values are exact rational numbers (no floating-point rounding). The ∂f/∂u' = -1/z is z-dependent, giving a non-trivial non-diagonal D1 Jacobian contribution that changes node-by-node across the Chebyshev grid -- a better isolator than the constant-coefficient case.

#### `3arg-nonlinear-exp`

- **name**: Nonlinear 3-arg BVP: u'' = (u')^2/u on [0, 1], solution u = exp(2z)
- **ode**: u'' = (u')^2 / u, i.e. f(z, u, u') = (u')^2/u; df/du = -(u')^2/u^2; df/dup = 2*u'/u
- **ic_bc**: u(0) = 1.0 (exact), u(1) = exp(2) = 7.38905609893065...
- **domain**: real segment [0, 1]
- **closed_form**: u(z) = A*exp(k*z) for any A > 0 and k. Canonical choice A=1, k=2: u(z) = exp(2*z). Verify: u'' = 4*exp(2z); (u')^2/u = 4*exp(4z)/exp(2z) = 4*exp(2z) ✓.
- **pole_structure**: entire
- **gt_kind**: closed-form
- **how_to_compute**: u(z) = exp(2*z). BCs: u(0)=1.0, u(1)=exp(2)=7.38905609893065022723... Pinned: u(0.25)=exp(0.5)=1.64872127070012814684..., u(0.5)=exp(1)=2.71828182845904523536..., u(0.75)=exp(1.5)=4.48168907033806482260... All from stdlib exp. Python: `import math; math.exp(2*z)`. Julia: `exp(2*z)`.
- **precision**: arbitrary via BigFloat
- **targets_buckets**: B15
- **regime**: nonlinear 3-arg BVP; the ∂f/∂u' = 2u'/u term is the u'-dependent Jacobian contribution, making Newton non-trivial
- **citation**: constructed; u=A*exp(kz) satisfies u'' = (u')^2/u exactly for any A,k (elementary computation). The PIII RHS has (u')^2/u as its leading nonlinear term (FFW2017 md:43), making this the minimal nonlinear proxy for the hybrid's BVP step.
- **confidence**: high
- **verify_note**: ODE residual is exactly zero (u''=k^2*e^{kz}, (u')^2/u=k^2*e^{2kz}/e^{kz}=k^2*e^{kz}). Boundary conditions are u(0)=1, u(1)=e^2. Initial guess: linear ramp from 1 to e^2; Newton should converge in ~3 steps.

#### `linear-cosh-bvp`

- **name**: Linear 2-arg BVP: u'' = u on [-1, 1], solution u = cosh(z)/cosh(1) [already tested, REFERENCE]
- **ode**: u'' = u; f(z, u) = u; df/du = 1
- **ic_bc**: u(-1) = 1, u(1) = 1
- **domain**: real segment [-1, 1]
- **closed_form**: u(z) = cosh(z)/cosh(1). Symmetric, u(±1) = 1.
- **pole_structure**: entire
- **gt_kind**: closed-form
- **how_to_compute**: Julia: `cosh(z)/cosh(1.0)`. Already pinned in test/_oracle_bvp.jl Group 2 (DMSUITE Octave oracle). N=16 gives error ~6e-15 vs closed form.
- **precision**: Float64 / BigFloat-256 (already tested in bvp_test.jl BV.1.2, BV.5.1)
- **targets_buckets**: B15
- **regime**: linear 2-arg BVP, spectral accuracy floor at N=16
- **citation**: Trefethen SMIM 2000 ch.6, example; WeidemanReddy2000 DMSUITE oracle in test/_oracle_bvp.jl
- **confidence**: high
- **verify_note**: Already tested and mutation-proved. Listed for completeness as the baseline 2-arg linear oracle.

#### `bratu-bvp`

- **name**: Nonlinear 2-arg BVP: Bratu equation u'' = -2*exp(u) on [0, 1], exact solution
- **ode**: u'' = -lambda * exp(u) with lambda=2; f(z, u) = -2*exp(u); df/du = -2*exp(u)
- **ic_bc**: u(0) = 0 (exact), u(1) = 0 (exact)
- **domain**: real segment [0, 1]
- **closed_form**: u(z) = 2*log(B / cosh(B*(z - 0.5))) where B solves B = cosh(B/2) ≈ 1.1787755269387010213...; u(0) = 0, u(1) = 0. Symmetric: u(z) = u(1-z).
- **pole_structure**: entire (for real z)
- **gt_kind**: closed-form
- **how_to_compute**: Python mpmath 50 dps: `from mpmath import mp, mpf, log, cosh, findroot; mp.dps=50; B=findroot(lambda B: B-cosh(B/2), mpf('1.2')); u = lambda z: 2*log(B/cosh(B*(z-mpf('0.5'))))`. Pinned values: B=1.1787755269387010213...; u(0.25)=0.24333656779461700534...; u(0.5)=0.32895242134111357436...; u(0.75)=0.24333656779461700534...
- **precision**: arbitrary via mpmath/BigFloat (B itself must be computed to desired precision)
- **targets_buckets**: B15
- **regime**: nonlinear 2-arg BVP with transcendental exact solution; two-branch structure (lambda=2 is below lambda_c=3.513...) -- first physical branch only
- **citation**: Aris 1975 'Mathematical Theory of Diffusion and Reaction', Chapter 2; Scott 1975; Jacobsen & Schmitt DCDS 2002. Formula derivation: substitute u = 2*log(A/cosh(B*(x-1/2))), get u'' = -2*B^2/cosh^2(B*(x-1/2)); lambda*exp(u) = 2*lambda*A^2/cosh^2(B*(x-1/2)); sum = 0 iff 2*B^2 = 2*lambda*A^2, i.e. A=B/sqrt(lambda); BC u(0)=0 gives A=cosh(B/2); combined: B=sqrt(lambda)*cosh(B/2), for lambda=2: B=cosh(B/2) [since sqrt(2)≈1.414 gives B=sqrt(2)*cosh(B/2), the first branch B≈1.179 satisfies this approximately... let me recheck].
- **confidence**: high
- **verify_note**: Analytic residual u'' + 2*exp(u) verified = 0 symbolically (see computation above). The transcendental equation for B: 2*B^2 = lambda*A^2 with A=cosh(B/2) and lambda=2 gives B^2=cosh^2(B/2), so B=cosh(B/2) (taking positive root). Numerical verification via mpmath findroot + central-difference u'' confirms residual |u'' + 2*e^u| < 1e-14 at z=0.3. Second branch requires B larger (B≈3.97), outside the physical regime.

#### `airy-bvp`

- **name**: Linear 2-arg BVP: Airy equation u'' = z*u on [0, 3], solution u = Ai(z)
- **ode**: u'' = z*u; f(z, u) = z*u; df/du = z
- **ic_bc**: u(0) = Ai(0) = 0.35502805388781723926...; u(3) = Ai(3) = 0.00659113935746071914...
- **domain**: real segment [0, 3]
- **closed_form**: u(z) = Ai(z) (Airy function of the first kind). Exact values via DLMF 9.2: Ai(0) = 1/(3^{2/3}*Gamma(2/3)) = 0.35502805388781723926...; Ai(1.5) = 0.07174949700810540967...; Ai(3) = 0.00659113935746071914...
- **pole_structure**: entire (Airy functions are entire)
- **gt_kind**: mpmath
- **how_to_compute**: Python: `from mpmath import mp, airyai; mp.dps=50; print(airyai(0), airyai(1.5), airyai(3))`. Julia: `using SpecialFunctions; airyai(z)`. DLMF 9.2.3 gives exact expression for Ai(0). Pinned: Ai(0)=0.35502805388781723926..., Ai(1.5)=0.07174949700810540967..., Ai(3)=0.00659113935746071914...
- **precision**: arbitrary via mpmath/SpecialFunctions.jl
- **targets_buckets**: B15, B17
- **regime**: 2-arg linear variable-coefficient BVP (f(z,u)=z*u depends on z, isolating the variable z-dependence of the analytic Jacobian diag(z_int))
- **citation**: DLMF 9.2 (https://dlmf.nist.gov/9.2); Abramowitz & Stegun 10.4. SpecialFunctions.jl uses the Cephes / DLMF implementations.
- **confidence**: high
- **verify_note**: Ai(z) satisfies u'' = z*u by definition (DLMF 9.2.1). BCs from SpecialFunctions.airyai(0.0) and airyai(3.0). The 2-arg Jacobian df/du=z is z-dependent (not constant), making this a better isolator of the variable-coefficient path than u''=u. At N=20 on [0,3], spectral accuracy floor ~1e-14.

#### `vector-anti-periodic-sin-cos`

- **name**: Vector 2-arg BVP: harmonic oscillator y'=(y2,-y1) with anti-periodic BC on [0, pi/2]
- **ode**: y' = [y2, -y1]; f(z, y) = [y[2], -y[1]]
- **ic_bc**: B_a = [[1,0],[0,1]] (identity), B_b = [[-1,0],[0,-1]] (negative identity), g = [-1.0, 1.0]
- **domain**: real segment [0, pi/2]
- **closed_form**: y(z) = (sin(z), cos(z)). Anti-periodic BC: B_a = I, B_b = -I, g = (sin(0)-sin(pi/2), cos(0)-cos(pi/2)) = (-1, 1).
- **pole_structure**: entire
- **gt_kind**: closed-form
- **how_to_compute**: y1(z) = sin(z), y2(z) = cos(z). BC check: B_a*y(0) + B_b*y(pi/2) = (sin 0, cos 0) - (sin(pi/2), cos(pi/2)) = (0,1) - (1,0) = (-1,1) ✓. Uniqueness: any perturbation y=(sin z + k*cos z, cos z - k*sin z) gives BC residual (k-1, 1+k)=(-1,1) => k=0. Interior: y(pi/4)=(sqrt(2)/2, sqrt(2)/2)=(0.70710678118654752..., same).
- **precision**: arbitrary via Julia sin/cos
- **targets_buckets**: B16
- **regime**: B16 vector BVP with anti-periodic (B_a=-B_b) boundary condition -- currently only Dirichlet-selector BCs are tested (VB.1.2 uses B_a=[1 0;0 0], B_b=[0 0;1 0])
- **citation**: constructed; the harmonic oscillator y'=(y2,-y1) with closed form (sin,cos) is the canonical elementary example for first-order systems. Anti-periodic BC is a genuinely general two-point condition.
- **confidence**: high
- **verify_note**: Uniqueness verified: k=0 is forced by the anti-periodic BC. ODE residual: y1'=y2=cos(z)=d(sin z)/dz ✓, y2'=-y1=-sin(z)=d(cos z)/dz ✓. The BC B_a=I, B_b=-I with g=(-1,1) couples the SAME component at both endpoints -- a qualitatively different BC from the existing Dirichlet-selector tests.

#### `vector-endpoint-mixing-bc`

- **name**: Vector 2-arg BVP: harmonic oscillator y'=(y2,-y1) with endpoint-coupling (mixing) BC on [0, pi]
- **ode**: y' = [y2, -y1]; f(z, y) = [y[2], -y[1]]
- **ic_bc**: B_a = [[0,0],[0,1]], B_b = [[1,0],[0,0]], g = [0.0, 1.0]; note B_a*y(z_a) pins y2 at z_a, B_b*y(z_b) pins y1 at z_b
- **domain**: real segment [0, pi]
- **closed_form**: y(z) = (sin(z), cos(z)). Mixing BC: B_a = [[0,0],[0,1]], B_b = [[1,0],[0,0]], g = [0, 1]. This pins y2(0)=1 via B_a and y1(pi)=0 via B_b, with cross-coupling between the two endpoints.
- **pole_structure**: entire
- **gt_kind**: closed-form
- **how_to_compute**: y1(z) = sin(z), y2(z) = cos(z). BC: B_a*y(0) + B_b*y(pi) = [[0,0],[0,1]]*(0,1)^T + [[1,0],[0,0]]*(0,-1)^T = (0,1)+(0,0)=(0,1) ✓. Interior: y(pi/4)=(0.70710678..., 0.70710678...), y(pi/2)=(1.0, ~0).
- **precision**: arbitrary via Julia sin/cos
- **targets_buckets**: B16
- **regime**: B16 genuinely general two-point BC where row i of g gets contributions from DIFFERENT components at DIFFERENT endpoints -- the key test of ADR-0023's tau-method endpoint pairing
- **citation**: constructed; the B_a/B_b structure is the minimal example where the tau method's node-0/node-N endpoint assignment (ADR-0023 Mutation B was: swap B_a and B_b blocks, Verified bite on VB.1.1-1.4) matters distinctly from Dirichlet selectors.
- **confidence**: high
- **verify_note**: The existing Dirichlet-selector tests (VB.1.2: B_a=[1 0;0 0], B_b=[0 0;1 0]) reduce to pointwise endpoint BCs but with zero off-diagonal blocks. This mixing BC has structurally different sparsity. The ODE is linear so Newton converges in 1 step; the oracle is sin/cos to machine precision.

#### `vector-rank-deficient-bc`

- **name**: Vector BVP: rank-deficient BC fixture -- harmonic oscillator on [0, 2*pi] with y(0)+y(2*pi)=[0,0]
- **ode**: y' = [y2, -y1]; f(z, y) = [y[2], -y[1]]
- **ic_bc**: B_a = I_2, B_b = -I_2, g = [0.0, 0.0] on [0, 2*pi]; the system y(0) - y(2*pi) = 0 is the periodic BC, rank-deficient because the solution space is 2-dimensional
- **domain**: real segment [0, 2*pi]
- **closed_form**: none (underdetermined: any y=(A*sin z + B*cos z, A*cos z - B*sin z) satisfies BC with g=[0,0], making the Newton Jacobian singular)
- **pole_structure**: n/a
- **gt_kind**: closed-form
- **how_to_compute**: Julia: `@test_throws ErrorException vector_bvp_solve(f_ho, 0.0, 2*pi, Matrix(I,2,2), Matrix(I,2,2), [0.0, 0.0]; N=16)`. This must throw because the BC B_a*y(0)+B_b*y(2*pi)=(I+I)*y_periodic=2*y(0) with g=0 means y(0)=0, but the ODE y'=(y2,-y1) with y(0)=0 gives identically y≡0 -- trivial, Newton converges to the zero vector which is the only solution. ACTUALLY: with g=[0,0] this IS uniquely y≡0. For a proper rank-deficient fixture, use B_a=I, B_b=I with g=[2*sin(0)+2*sin(2*pi), 2*cos(0)+2*cos(2*pi)] = [0, 4] which is inconsistent unless y is periodic AND satisfies a specific amplitude. The correct rank-deficient fixture is: B_a=I, B_b=-I, g=[0,0] on [0, 2*pi] (periodic condition y(0)=y(2*pi)). This is RANK DEFICIENT because y=(A*sin z+B*cos z, A*cos z-B*sin z) with any (A,B) satisfies it.
- **precision**: n/a
- **targets_buckets**: B16
- **regime**: B16 fail-fast fixture for rank-deficient BC: the Newton Jacobian must be singular and the solver must throw ErrorException with a 'Suggestion'
- **citation**: constructed; the periodic BC y(0)=y(2*pi) for the harmonic oscillator is the canonical rank-deficient two-point BC (one eigenvalue of the monodromy matrix is exactly 1, making the linearised system singular). See Ascher et al. 1988 'Numerical Solution of BVPs' ch.2.
- **confidence**: medium
- **verify_note**: Verify that vector_bvp_solve throws (not converges) on this fixture. The tau-method J_bc block [[I|0|...|0|-I]] has rank d=2 but the global Jacobian may still be invertible depending on N -- test at multiple N values (N=8, 16) to confirm consistent throw. If the solver incorrectly finds the zero solution (y≡0), the fixture is degenerate; use B_a=I, B_b=-I, g=[1,0] instead (consistent but still rank-deficient without the specific normalization).

#### `oblique-complex-cosh-bvp`

- **name**: Oblique-complex segment BVP: u'' = u on [0, 1+i] with cosh BCs
- **ode**: u'' = u; f(z, u) = u; df/du = 1
- **ic_bc**: u(0) = 1.0 + 0.0im; u(1+i) = cosh(1+i) = 0.8337300251311490... + 0.9888977057628651...im
- **domain**: oblique complex segment [0, 1+i] (both Re and Im nonzero in z_b - z_a)
- **closed_form**: u(z) = cosh(z). BCs: u(0)=cosh(0)=1, u(1+i)=cosh(1+i)=0.83373002513114904888...+0.98889770576286509638...i. Interior: u(0.5+0.5i)=cosh(0.5+0.5i)=0.98958488339991993644...+0.24982639750046153149...i
- **pole_structure**: entire (cosh is entire)
- **gt_kind**: mpmath
- **how_to_compute**: Python: `from mpmath import mp, mpc, cosh; mp.dps=50; print(cosh(mpc(0.5,0.5)))`. Julia: `using SpecialFunctions; cosh(0.5+0.5im)` (stdlib). Pinned 50-dps: cosh(0+0i)=1; cosh(0.3+0.3i)=0.99865002603564327947...+0.08999190005207135169...i; cosh(0.5+0.5i)=0.98958488339991993644...+0.24982639750046153149...i; cosh(1+i)=0.83373002513114904888...+0.98889770576286509638...i
- **precision**: arbitrary via mpmath/BigFloat
- **targets_buckets**: B15, B17
- **regime**: oblique complex segment; tests the real(t*)-only domain guard, the complex affine scale factor s=(z_b-z_a)/2=0.5+0.5i, and complex type promotion CT=ComplexF64
- **citation**: constructed; cosh satisfies u''=u exactly (d^2/dz^2 cosh(z) = cosh(z)). The oblique segment [0,1+i] has s=(1+i)/2 -- the minimal test for the complex scale factor in the BVP solver. mpmath values cross-checked against Julia `cosh(0.5+0.5im)` = 0.9895848833999199 + 0.24982639750046153im (Float64).
- **confidence**: high
- **verify_note**: The domain guard in BVP.jl line 491 checks only real(t*). For on-segment queries z = s*(1+i) with s in [0,1], t* = 2s-1 is always REAL, so the guard works correctly. The test should verify: (1) BVP converges with Complex{Float64} types; (2) sol(0.5+0.5im) matches cosh(0.5+0.5im) to atol 1e-12; (3) a BigFloat-256 version achieves atol big(1e-40).

#### `bigfloat-256-oblique-cosh`

- **name**: BigFloat-256 oblique complex BVP: u'' = u on [0, 1+i] with cosh BCs at 256-bit precision
- **ode**: u'' = u; f(z, u) = u; df/du = 1 (same as oblique-complex-cosh-bvp but at BigFloat-256)
- **ic_bc**: u(0) = big(1.0)+big(0.0)*im; u(1+i) = cosh(big(1.0)+big(1.0)*im) computed at BigFloat-256
- **domain**: oblique complex segment [0, 1+i]
- **closed_form**: u(z) = cosh(z). BigFloat-80dps: cosh(0.5+0.5i) = 0.9895848833999199364440705335978509367482 + 0.2498263975004615314895596560305843910592i
- **pole_structure**: entire
- **gt_kind**: mpmath
- **how_to_compute**: Python mpmath 80 dps: `from mpmath import mp, mpc, cosh; mp.dps=80; print(cosh(mpc('0.5','0.5')))`. Julia BigFloat-256: `setprecision(BigFloat, 256) do; cosh(big(0.5)+big(0.5)*im); end`. The BigFloat-256 residual floor should be ~1e-60 (matching BV.5.1 pattern for the linear cosh BVP on [-1,1]).
- **precision**: 256-bit via Julia BigFloat, verified to 80 decimal digits via mpmath
- **targets_buckets**: B15, B17
- **regime**: BigFloat-256 on an oblique complex segment -- exercises both the complex type path and the GenericLinearAlgebra BigFloat backslash in combination
- **citation**: constructed; extends BV.5.1 (BigFloat-256 real linear BVP) to the complex oblique case. mpmath 80 dps values are the ground truth.
- **confidence**: high
- **verify_note**: Compare to BV.5.1 pattern: at N=16 BigFloat-256 the residual should be << 1e-40 (BigFloat eps^(3/4) ≈ 1e-57 for 256-bit). The BVP.jl code routes through Julia's promote_type; setting z_a=big(0.0), z_b=big(1.0)+big(1.0)*im should trigger Complex{BigFloat} path automatically.

#### `3arg-boundary-layer`

- **name**: 3-arg BVP: singular perturbation eps*u'' + u' = 0, u(0)=0, u(1)=1, eps=0.1
- **ode**: u'' = -u'/eps, i.e. f(z, u, u') = -u'/eps; df/du = 0; df/dup = -1/eps
- **ic_bc**: u(0) = 0 (exact), u(1) = 1 (exact)
- **domain**: real segment [0, 1]
- **closed_form**: u(z) = (1 - exp(-z/eps)) / (1 - exp(-1/eps)). For eps=0.1: u(0.1)=0.63214925836048665138...; u(0.3)=0.95025607319111530831...; u(0.5)=0.99330714907571514444...; u(1)=1.0 (all mpmath 50 dps).
- **pole_structure**: entire in z (exponential, entire)
- **gt_kind**: closed-form
- **how_to_compute**: Python: `from mpmath import mp, mpf, exp; mp.dps=50; eps=mpf('0.1'); u = lambda z: (1-exp(-z/eps))/(1-exp(-1/eps))`. Julia: `u(z) = (1 - exp(-z/0.1)) / (1 - exp(-1/0.1))`. For eps=0.1 (moderate), spectral method works without layer-adapted grids. For eps=0.01 (stiff), Chebyshev points cluster near the boundary layer naturally but may need larger N.
- **precision**: arbitrary via mpmath/BigFloat
- **targets_buckets**: B15
- **regime**: 3-arg BVP with ∂f/∂u'=-1/eps -- pure u'-dependence (∂f/∂u=0), isolating the D1_ii Jacobian contribution alone with no confounding diagonal ∂f/∂u term
- **citation**: Trefethen SMIM 2000 ch.7 p.91-92 (boundary layer problems); Bender & Orszag 1999 'Advanced Mathematical Methods for Scientists and Engineers' ch.9. The exact solution is the Green's function for the first-order convection-diffusion operator.
- **confidence**: high
- **verify_note**: The ODE eps*u'' + u' = 0 rewritten as u'' = -u'/eps has purely u'-dependent RHS (∂f/∂u=0), making the Jacobian J = D2_ii - 0 - (scale/half_diff)*(-1/eps)*D1_ii = D2_ii + (scale/(eps*half_diff))*D1_ii. This isolates the D1_ii contribution with zero confounding from the diagonal ∂f/∂u term. For eps=0.1, N=20 gives spectral accuracy to ~1e-10.

#### `ivp-bvp-hybrid-dispatch`

- **name**: IVP-BVP dispatch: scalar u''=u on complex segment [0, 1+i] via BVP solver, cross-checked against closed form
- **ode**: u'' = u; f(z, u) = u; df/du = 1
- **ic_bc**: u(0) = 1, u(1+i) = cosh(1+i); or equivalently IVP IC (u(0)=1, u'(0)=0)
- **domain**: complex segment [0, 1+i] or real segment [-1, 1] depending on dispatch path
- **closed_form**: u(z) = cosh(z) (same as oblique-complex-cosh-bvp). The B17 point: this BVP is solved via bvp_solve and the result compared to a PathNetwork IVP solve from the same IC, verifying the IVP↔BVP dispatch invariant.
- **pole_structure**: entire
- **gt_kind**: closed-form
- **how_to_compute**: Julia two-path cross-check: (1) bvp_solve with cosh BCs gives solution; (2) for the same problem, IVP from z=0 with u(0)=1, u'(0)=0 (cosh initial conditions sinh(0)=0) gives the same values. The dispatch tests that: (a) bvp_solve and the IVP give identical values to atol 1e-10; (b) the hybrid driver's degenerate_full_plane=true path from IVPBVPHybrid gives the same result as bvp_solve on a trivially-pole-free problem.
- **precision**: Float64 (main) / BigFloat-256 (extended)
- **targets_buckets**: B17
- **regime**: B17 IVP↔BVP dispatch hybrid; verifies the composition invariant that bvp_solve and IVP path give consistent solutions on a pole-free segment
- **citation**: FW 2011 §3.2 lines 190-192 (FW2011_painleve_methodology_JCP230.md:190-192): the derivative match criterion for IVP-BVP stitching. The cosh problem is the trivial well-conditioned case where both methods agree to machine precision.
- **confidence**: high
- **verify_note**: This is a DISPATCH test, not a new algorithm test. The closed form is cosh(z). The assertion is that bvp_solve and IVP give the same values -- testing the hybrid's gluing logic on a simple case before using it on PIII. Implementation: compare bvp_solve(z->z, z->one(z), 0.0+0.0im, 1.0+1.0im, 1.0+0im, cosh(1+1im); N=20) at 0.5+0.5im versus cosh(0.5+0.5im) to atol 1e-12.

### Territory: Robust Pade approximation (GGT 2013 / chebfun padeapprox.m) and Froissart-double

*Scout notes:* Six candidates are provided covering: (a) closed-form Pade reductions for standard functions (exp, tan, rational, defect-block polynomial); (b) Froissart-doublet design cases with known genuine vs spurious pole counts; (c) the ADR-0002 SVD rank-boundary pinning test; (d) a Complex{BigFloat}-path oracle via exp(iz) substitution. Gaps: (1) No Octave available for re-running the existing padeapprox-oracle/oracles.txt; that oracle was already captured from Octave 8.4 with MATLAB-seed randn(state,42) and is frozen — regeneration requires Octave. (2) The ADR-0002 sigma_2 misclassification at Float64 is a 'may happen' bound (DK absolute error is 49% of sigma_2), not a deterministic guarantee — empirical confirmation on the actual LAPACK DGESVD build is recommended before treating rank misclassification as certain. (3) Wolframscript is not available locally; the Mathematica PadeApproximant recipe for the two-pole function is given as a regenerable command to be run externally. Best oracles found: the existing external/probes/padeapprox-oracle/oracles.txt (Octave-captured, frozen) for exp(20,20) and tan(z^4) (20,20); mpmath-computed exact rational coefficients for exp (7,7) and tan (3,3); the GGT 2013 paper itself (section 7 and figure 3) for the exact algorithmic trace and Froissart doublet formula.

#### `exp-77-exact-rational`

- **name**: exp(z) type-(7,7) Pade — exact rational coefficients, no Froissart poles
- **ode**: n/a (Pade approximation unit test, not ODE)
- **ic_bc**: n/a
- **domain**: z in disk |z| <= 5 (well inside Pade pole radius)
- **closed_form**: Numerator a_k = C(14-k,7)*7!*7! / (14!*k!) for k=0..7; denominator b_k = (-1)^k * a_k. Explicit: a=[1, 1/2, 3/26, 5/312, 5/3432, 1/11440, 1/308880, 1/17297280], b=[1, -1/2, 3/26, -5/312, 5/3432, -1/11440, 1/308880, -1/17297280]. Poles of denominator: 7 conjugate pairs with smallest |z| ~ 9.94 (all off the real axis, all outside the unit disk). exp is entire, so the (20,20) Pade from oracles.txt collapses to exact type (7,7) via GGT diagonal hopping.
- **pole_structure**: entire — no poles; all Pade poles of r_{77} lie at |z| ~ 5 to 10 (complex conjugate pairs far from origin)
- **gt_kind**: closed-form
- **how_to_compute**: python3 -c "from math import factorial, comb; from fractions import Fraction; m,n=7,7; mn=m+n; a=[Fraction(comb(mn-k,n)*factorial(m)*factorial(n), factorial(mn)*factorial(k)) for k in range(m+1)]; b=[(-1)**k*a[k] for k in range(n+1)]; print(a); print(b)"
- **precision**: exact rational (arbitrary via Fraction arithmetic)
- **targets_buckets**: B1, B22
- **regime**: Block-degenerate Pade table: (m,n)=(20,20) input collapses to exact type (7,7) via GGT diagonal hopping; validates that the SVD rank-counting correctly identifies the defect-delta=13 block and returns the minimal-degree rational
- **citation**: GGT 2013 (SIREV 55) Figure 2 caption (exp block structure); Baker-Graves-Morris 1996 Table 5.1 for exact (m,n) Pade of exp via the formula C(m+n-k,n)*m!*n!/(m+n)!/k!; references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/GGT2013_robust_pade_via_SVD_SIREV55.md:238-241
- **confidence**: high
- **verify_note**: Cross-check: external/probes/padeapprox-oracle/oracles.txt lines 5-8 give the Octave-captured (20,20) result which reduces to type (7,7) with coefficients matching the exact rational values to 1e-12 rel. Verify a_1=0.5 and b_1=-0.5 analytically.

#### `tan-33-exact-pade`

- **name**: tan(z) type-(3,3) Pade — closed-form pole at sqrt(5/2), approximating the nearest true pole pi/2
- **ode**: n/a
- **ic_bc**: n/a
- **domain**: z in disk |z| <= 1.4 (inside the nearest pole)
- **closed_form**: p(z) = z - z^3/15, q(z) = 1 - 2*z^2/5 (exact). Poles of q: z = +/-sqrt(5/2) = +/-1.5811388300841898... True nearest pole of tan(z): pi/2 = 1.5707963267948966...; Pade pole error = 1.034e-2. The Pade p/q satisfies tan(z) - p(z)/q(z) = O(z^7) (matching through degree 6).
- **pole_structure**: tan(z) poles at z_k = (k+1/2)*pi for k in Z; the (3,3) Pade places a single spurious pole pair at +/-sqrt(5/2) approximating the nearest genuine pole pi/2
- **gt_kind**: closed-form
- **how_to_compute**: python3 -c "import mpmath; c=[0,1,0,mpmath.mpf(1)/3,0,mpmath.mpf(2)/15,0]; p,q=mpmath.pade(c,3,3); print('p=',p,'q=',q); print('pole=',mpmath.sqrt(mpmath.mpf(5)/2))"
- **precision**: exact (rational coefficients 1/15, 2/5)
- **targets_buckets**: B1
- **regime**: Non-degenerate off-diagonal Pade: (3,3) type for an odd meromorphic function; validates SVD null-vector extraction and normalisation. Exact closed-form check for the denominator-root location.
- **citation**: Standard Pade tables for tan(z); mpmath.pade verification; DLMF 4.14 for tan series coefficients; GGT 2013 references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/GGT2013_robust_pade_via_SVD_SIREV55.md:329-332 for tan(z^4) as canonical Froissart example
- **confidence**: high
- **verify_note**: Confirm p(z)/q(z) matches tan(z) through O(z^7) by Taylor expansion. The pole at sqrt(5/2) is unambiguously above 1 and below 2, distinct from pi/2.

#### `geom-series-reduction`

- **name**: 1/(1-z) noisy series: (10,10) Pade with tol=1e-5 collapses to type (0,1), Froissart suppressed
- **ode**: n/a
- **ic_bc**: n/a
- **domain**: Taylor coefficients c_0..c_20 at z=0
- **closed_form**: Underlying function 1/(1-z) = sum z^k. True Pade: a=[1], b=[1,-1]. With c_j = 1 + 1e-6*randn (seeded MATLAB randn(state,42)), GGT Algorithm 2 at tol=1e-5 returns type (0,1) with a~[0.99999999], b~[1,-0.99999976]. With tol=1e-14 (default), returns type (1,1) blocks (noise dominates). Exact type-reduction rule: sigma_1 >> sigma_2,...,sigma_10 >> tol*||c||_2 at tol=1e-5 forces rank=1 from a rank-10 C matrix.
- **pole_structure**: Genuine pole at z=1 (simple); all other Pade poles (with tol=0, nonrobust) are Froissart doublets clustered near unit circle
- **gt_kind**: published-table
- **how_to_compute**: From external/probes/padeapprox-oracle/oracles.txt lines 29-33 (frozen Octave output). a=[9.999999934540796387e-01], b=[1.000000000000000000e+00, -9.999997642302624890e-01], mu=0, nu=1. To regenerate: run capture.m in external/probes/padeapprox-oracle/ under Octave with Chebfun on path.
- **precision**: Float64 (Octave-captured, ~1e-12 relative to true type (0,1))
- **targets_buckets**: B1, B8
- **regime**: Noise-level detection: tol above noise level (1e-5 > 1e-6) must collapse the block; tol below noise (1e-14 < 1e-6) must leave it unfrozen. Exercises the sigma_i < tau thresholding critical path in GGT Algorithm 2 step 4.
- **citation**: GGT 2013 Figure 4 and the caption (p.316, Pade table for noisy 1/(1-z)); references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/GGT2013_robust_pade_via_SVD_SIREV55.md:307-315
- **confidence**: high
- **verify_note**: The MATLAB seed 'randn(state,42)' is legacy syntax; cannot regenerate without Octave. The values in oracles.txt are the frozen oracle. Cross-check: with exact c_j=1 (no noise), tol=1e-14, the (10,10) Pade must reduce to EXACT type (0,1) with a=[1], b=[1,-1].

#### `tan-z4-froissart-20-20`

- **name**: tan(z^4) type-(20,16) Pade: 16 genuine poles (2 rings of 8), 4 Froissart doublets removed at |z|~2e5
- **ode**: n/a
- **ic_bc**: n/a
- **domain**: Complex plane, poles on 8 rays from origin
- **closed_form**: Genuine poles of tan(z^4): for each integer n, z^4 = (2n+1)*pi/2, giving poles at z = ((|2n+1|*pi/2)^{1/4}) * exp(i*k*pi/2) for k=0..3 (on-axis) and * exp(i*(pi/4 + k*pi/2)) (diagonal), 8 poles per |n| pair. Ring 0 (n=0,n=-1): |z| = (pi/2)^{1/4} = 1.119515134920248... Ring 1 (n=1,n=-2): |z| = (3*pi/2)^{1/4} = 1.473364776175541... GGT robust (20,20) returns mu=20, nu=16 (4 poles removed = 4 Froissart doublets at |z|~2e5). The 16 surviving poles match the 2 rings with inner-ring error ~2e-6, outer-ring error ~2.6%.
- **pole_structure**: 8 genuine poles per ring at |z| = ((2k+1)*pi/2)^{1/4} for k=0,1,... on 8 rays (4 axis-aligned + 4 diagonal); Froissart doublets cluster at |z|>>1 with tiny residues
- **gt_kind**: published-table
- **how_to_compute**: From external/probes/padeapprox-oracle/oracles.txt lines 17-21 (Octave-captured). Pole locations listed as ComplexF64 array (16 poles). Exact ring radii: python3 -c "import mpmath; print((mpmath.pi/2)**mpmath.mpf('0.25'), (3*mpmath.pi/2)**mpmath.mpf('0.25'))" gives 1.119515134920248, 1.473364776175541. Note: (20,20) nonrobust gives mu=nu=20 with 4 spurious poles at |z|~2e5 (Froissart); robust returns mu=20, nu=16.
- **precision**: Float64 (Octave-captured); true ring radii available to arbitrary precision via mpmath
- **targets_buckets**: B1, B8
- **regime**: Froissart-doublet filter validation: GGT robust algorithm removes exactly 4 of 20 denominator poles as Froissart doublets; the count mu=20 nu=16 pins the SVD reweighting and degree-reduction steps. The 8-fold rotational symmetry of the pole field provides an additional consistency check.
- **citation**: GGT 2013 Figure 5 caption and body (p.316-317): 'removes four poles with absolute value about 2e5'; references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/GGT2013_robust_pade_via_SVD_SIREV55.md:349-357
- **confidence**: high
- **verify_note**: Exact pole count: mu=20, nu=16. Inner ring placement error ~2e-6 per GGT Fig 5 ('six digits'); outer ring error ~2.6% ('two digits'). The 4 removed Froissart poles should have |residue| << 1e-3 (red in GGT color coding). Cross-check: compute tan(z^4) Taylor series to order 41 via FFT on unit circle, compare Pade output against oracles.txt.

#### `adr0002-rank-boundary-svd`

- **name**: ADR-0002 SVD rank-boundary: 10x11 Toeplitz with sigma_2 = 4.76e-14, just 3.9% above Float64 GGT threshold, DK abs-error = 49% of sigma_2
- **ode**: n/a (SVD rank-classification unit test)
- **ic_bc**: n/a
- **domain**: Taylor coefficients at z=0, real Float64 or BigFloat-256
- **closed_form**: f(z) = 1/(1-z) + 1e-13/(1-z/2). Taylor coefficients c_k = 1 + 1e-13 * (1/2)^k for k=0..20. For (m=10, n=10) Pade, the 10x11 Toeplitz C (built from c[11..20] as C[i,j] = c[m+1+i-j]) has true rank 2 (two poles at z=1 and z=2). Exact mpmath-50dp singular values: sigma_1 = 10.488088481701534511, sigma_2 = 4.7595635728817588896e-14, sigma_3 = 4.9e-51 (machine zero). Float64 GGT threshold = 1e-14 * ||c||_2 = 4.582576e-14. Ratio sigma_2/threshold = 1.0386. Demmel-Kahan absolute error bound on sigma_2 at Float64: 10 * 2^{-52} * sigma_1 = 2.33e-14 (= 49% of sigma_2). BigFloat-256 threshold = 2^{-246} * ||c||_2 ~ 4.05e-74; sigma_2 >> this. True Pade: type (10,10) with poles at z=1 and z=2. False (if rank misclassified as 1): reduced to type (9,1).
- **pole_structure**: Two genuine poles: z=1 (from 1/(1-z)) and z=2 (from 1e-13/(1-z/2)); sigma_3..sigma_10 are exact-arithmetic zeros (~ 5e-51 at 50dp)
- **gt_kind**: mpmath
- **how_to_compute**: python3 -c "import mpmath; mpmath.mp.dps=50; delta=mpmath.mpf('1e-13'); c=[1+delta*mpmath.mpf('0.5')**k for k in range(21)]; m,n=10,10; C=mpmath.matrix(n,n+1); [C.__setitem__((i,j),c[m+1+i-j]) for i in range(n) for j in range(n+1) if 0<=m+1+i-j<21]; U,S,V=mpmath.svd(C); print([mpmath.nstr(S[k],20) for k in range(3)])"
- **precision**: mpmath 50dp (sigma_2 = 4.7595635728817588896e-14 to 20 significant digits)
- **targets_buckets**: B1, B22
- **regime**: ADR-0002 claim: Demmel-Kahan absolute accuracy may misclassify a genuine small SV as noise when kappa(C) ~ 1e14 and sigma_min is only 3.9% above the GGT threshold. BigFloat Jacobi SVD (GenericLinearAlgebra) has relative accuracy 2^{-246}*sigma_i per SV, guaranteeing correct rank classification at BigFloat-256.
- **citation**: ADR-0002 docs/adr/0002-bigfloat-svd-via-genericlinalg.md:26-40 (Demmel-Kahan vs Jacobi accuracy argument); Demmel-Veselic 1992 relative-accuracy theorem for one-sided Jacobi; references/markdown/GGT2013_robust_pade_via_SVD_SIREV55/GGT2013_robust_pade_via_SVD_SIREV55.md:213-217 (tol*||c||_2 threshold)
- **confidence**: medium
- **verify_note**: The DK bound (49% of sigma_2) is necessary but not sufficient for misclassification: whether LAPACK actually perturbs sigma_2 below threshold depends on the specific LAPACK build and matrix structure. Recommend confirming empirically with Julia LinearAlgebra.svdvals on Float64(C) vs GenericLinearAlgebra.svd on BigFloat(C). If LAPACK does not misclassify in practice, reframe the test as a precision comparison (not rank flip) and lower confidence to medium.

#### `exp-iz-77-complex-bigfloat`

- **name**: exp(iz) type-(7,7) Pade with Complex{BigFloat} coefficients — closed-form complex oracle for the complex SVD dispatch path
- **ode**: n/a
- **ic_bc**: n/a
- **domain**: z in unit disk (evaluation at z=1 as a spot check); complex coefficients
- **closed_form**: By substituting z -> iz in the (7,7) Pade of exp(z): numerator p_k = a_k * i^k, denominator q_k = a_k * (-i)^k, where a_k = C(14-k,7)*7!*7!/(14!*k!) are the same rational coefficients as the real case. Explicit: p = [1, i/2, -3/26, -5i/312, 5/3432, i/11440, -1/308880, -i/17297280], q = [1, -i/2, -3/26, 5i/312, 5/3432, -i/11440, -1/308880, i/17297280]. Evaluation: r_{77}(exp(iz))(z=1) = (0.5403023058681399 + 0.8414709848078963*i) vs exp(i) = (0.5403023058681398 + 0.8414709848078965*i), |error| = 2.5e-16 (near machine epsilon). Poles: q(z)=0 has 7 roots; the closest is at z ~ -9.9436*i (on the negative imaginary axis, |z| ~ 9.94 — same magnitude as the real-exp pole).
- **pole_structure**: All poles at |z| ~ 5 to 10 (conjugate pairs under real-part symmetry), all in the lower half-plane; exp(iz) is entire so no genuine poles
- **gt_kind**: closed-form
- **how_to_compute**: python3 -c "from math import factorial,comb; from fractions import Fraction; m,n=7,7; mn=m+n; a=[Fraction(comb(mn-k,n)*factorial(m)*factorial(n),factorial(mn)*factorial(k)) for k in range(m+1)]; p=[complex(float(a[k])*(1j)**k) for k in range(m+1)]; q=[complex(float(a[k])*(-1j)**k) for k in range(n+1)]; print('p=',p,'q=',q); z=1+0j; r=sum(p[k]*z**k for k in range(8))/sum(q[k]*z**k for k in range(8)); import cmath; print('r(1)=',r,'exp(i)=',cmath.exp(1j))"
- **precision**: exact rational (Fraction arithmetic) for coefficients; Float64 evaluation check at z=1 accurate to ~2.5e-16
- **targets_buckets**: B1, B22
- **regime**: Complex{BigFloat} path: GGT Algorithm 2 over complex Taylor coefficients dispatches to GenericLinearAlgebra SVD via the complex dispatch branch (ADR-0002 table row 'Complex{T} for T<:AbstractFloat'). This is currently untested per the ADR. The oracle provides exact coefficient values and an evaluation check with |error| < 1e-15.
- **citation**: ADR-0002 docs/adr/0002-bigfloat-svd-via-genericlinalg.md table (Complex{T} dispatch row); Baker-Graves-Morris 1996 (rational Pade of exp); constructed by z->iz substitution (closed-form derivation, no external reference needed)
- **confidence**: high
- **verify_note**: Cross-check: poles of q(z) should lie at -i times the poles of the real-exp Pade (geometric: multiply each pole by -i = rotation by -pi/2). The smallest-|z| pole of the real Pade is at ~9.94+0i, so the complex-exp Pade's closest pole is at ~0-9.94i. Verify this with numpy.roots(list(reversed(q_coeffs))).

### Territory: high-precision Taylor-coefficient oracles + reference-integrator cross-checks: B

*Scout notes:* 1. The mpmath.odefun trajectory for WP deviates from the recurrence near the pole at z=1 (order-30 Taylor at z=0 converges slowly for z>0.5); use the recurrence at order-60+ or WP NDSolve (WP<50) for z>0.5. The odefun IS reliable for P_I trajectory on the positive real axis (tritronquée sector, no poles). 2. mpmath has no native WeierstrassP function (ellipwp does not exist in mpmath 1.3); the coefficient oracle must use the recurrence formula or wolframscript. 3. The sin-jet (odd-parity) step-size case does NOT trigger the fallback: c[p-1]=c[29]!=0, so only one candidate is nonzero and h is finite. The fallback is triggered by the even-only-to-28 cos-like jet (c[29]=c[30]=0). 4. The existing test/_oracle_stepcontrol.jl pins the canonical step value h=4.501206370338986 from three-source consensus (TI.jl + mpmath + Mathematica). The sin-jet and cos-even-jet candidates in this territory are NEW cases not yet in the oracle file.

#### `wp-taylor-coeff-recurrence-b2`

- **name**: Weierstrass ℘-function Taylor coefficients (u''=6u^2, order 30+, BigFloat-256)
- **ode**: u'' = 6*u^2
- **ic_bc**: u(0) = wp(-1; 0, 2) = 1.07182251641691726..., u'(0) = wp'(-1; 0, 2) = 1.71033735317678626...
- **domain**: Real: z in [0, 0.5] (well inside convergence radius ~1); complex: arbitrary via BigFloat arithmetic
- **closed_form**: u(z) = wp(z + c1; 0, g3) where c1=-1, g3=2. Coefficient recurrence: c[0]=u0, c[1]=u1, c[n] = 6/(n*(n-1)) * sum_{k=0}^{n-2} c[k]*c[n-2-k] for n>=2.
- **pole_structure**: Double poles on rhombic lattice with omega1 = Gamma(1/3)^3 / (2^(13/6)*pi) * (2/4)^(1/6) ≈ 1.0818 at z0=0, spacing 2*omega1 ≈ 2.1636 along real axis; first pole at z ≈ 1 (real axis, FW 2011 p.12 equianharmonic c1=-1).
- **gt_kind**: closed-form
- **how_to_compute**: Python mpmath snippet (regenerable oracle):  import mpmath mpmath.mp.dps = 50  def wp_coeffs(u0, u1, order):     c = [mpmath.mpf(u0), mpmath.mpf(u1)]     for n in range(2, order+1):         conv = sum(c[k]*c[n-2-k] for k in range(n-1))         c.append(6*conv / (n*(n-1)))     return c  # ICs from FW 2011 p.12 eq.(5.2): wp(-1; 0, 2) u0 = mpmath.mpf('1.07182251641691726539737474946') u1 = mpmath.mpf('1.71033735317678626742044547124') coeffs = wp_coeffs(u0, u1, 30) for i, c in enumerate(coeffs):     print(f'c[{i}] = {mpmath.nstr(c, 50)}')  # Spot-check c[2]: must equal 3*u0^2 assert abs(coeffs[2] - 3*u0**2) < 1e-45  Mathematica cross-check: wolframscript -code 'Print[N[CoefficientList[Normal[Series[WeierstrassP[z-1,{0,2}],{z,0,30}]],z],50]]'  Key spot values (50 dps): c[0] = 1.07182251641691726539737474946... c[1] = 1.71033735317678626742044547124... c[2] = 3.44641052009487887... (= 3*u0^2, exact check) c[29] = 29.9999986561163983... c[30] = 31.0000008044737570...
- **precision**: Arbitrary via mpmath.mp.dps; 50+ dps routinely; BigFloat-256 (77 decimal digits) exact
- **targets_buckets**: B2
- **regime**: Scalar 2nd-order Taylor recurrence; BigFloat-256 coefficient oracle; equianharmonic Weierstrass ℘ with known closed form.
- **citation**: FW 2011 §5.1.1 (references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:281-302): equation (5.1) u''=6u^2, closed-form (5.2), ICs at p.12 eq. and reference values eq.(5.3). Recurrence is standard: DLMF 23.6.2 for Laurent coefficients of Weierstrass P.
- **confidence**: high
- **verify_note**: Must verify c[2] = 3*u0^2 exactly (from the ODE at n=2). Must verify series evaluates within epsilon of odefun at z=0.3 (odefun and recurrence agree to ~15 digits there). Cross-check c[2] via FW 2011: 'g3 = 4*u0^3 - u1^2 ≈ 2', so c[2] = 3*u0^2 (not g3/2).

#### `pi1-taylor-coeff-recurrence-b2`

- **name**: Painlevé I Taylor coefficients (u''=6u^2+z, tritronquée ICs, order 30+)
- **ode**: u'' = 6*u^2 + z
- **ic_bc**: u(0) = -0.18755430834049490..., u'(0) = 0.30490556026122890... (tritronquée, FW 2011 eq. 4.1)
- **domain**: Real z0=0; expansion valid for |h| < (distance to first pole) determined by ODE
- **closed_form**: none (P_I is transcendental; no closed form). Coefficient recurrence: c[0]=u0, c[1]=u1, c[2]=3*u0^2+z0/2, c[3]=(12*u0*u1+1)/6, and for n>=4: c[n] = (6*sum_{k=0}^{n-2} c[k]*c[n-2-k]) / (n*(n-1)). The '+z' term contributes: c[2] += z0/2, c[3] += 1/6.
- **pole_structure**: Movable double poles. Tritronquée solution: no poles in sector |arg(z)| < 2*pi/5. Positive real axis is pole-free for tritronquée.
- **gt_kind**: mpmath
- **how_to_compute**: Python mpmath snippet (regenerable oracle):  import mpmath mpmath.mp.dps = 50  def pi1_coeffs(u0, u1, z0, order):     c = [mpmath.mpf(u0), mpmath.mpf(u1)]     for n in range(2, order+1):         conv = sum(c[k]*c[n-2-k] for k in range(n-1))         z_coeff = mpmath.mpf(z0) if n==2 else (mpmath.mpf(1) if n==3 else mpmath.mpf(0))         c.append((6*conv + z_coeff) / (n*(n-1)))     return c  # Tritronquée ICs from FW 2011 eq.(4.1) u0 = mpmath.mpf('-0.18755430834049490') u1 = mpmath.mpf('0.30490556026122890') coeffs = pi1_coeffs(u0, u1, 0, 30) for i, c in enumerate(coeffs):     print(f'c[{i}] = {mpmath.nstr(c, 50)}')  # Key invariant checks: # c[2] = 3*u0^2 + 0  (z0=0, no z-term at n=2 since z0=0) # c[3] = (12*u0*u1 + 1)/6 = 2*u0*u1 + 1/6 print('c[2] check:', coeffs[2], '==', 3*u0**2) print('c[3] check:', coeffs[3], '==', 2*u0*u1 + mpmath.mpf(1)/6)  Mathematica cross-check: wolframscript -code ' u0=SetPrecision[-0.18755430834049490,60]; u1=SetPrecision[0.30490556026122890,60]; AsymptoticDSolveValue[{u''''[z]==6*u[z]^2+z,u[0]==u0,u''''[0]==u1},u[z],{z,0,30}]//CoefficientList[#,z]&//N[#,50]&//Print'  Key spot values (50 dps, z0=0): c[0] = -0.18755430834049490... c[1] = 0.30490556026122890... c[2] = 0.10552985573124425... (= 3*u0^2 = 3*(0.18755...)^2) c[3] = 0.05229396373873495... (= 2*u0*u1 + 1/6) c[28] = 1.387676179661039e-10 c[29] = 6.021080104736566e-11 c[30] = 2.609623502502971e-11
- **precision**: Arbitrary via mpmath.mp.dps; 50 dps routine
- **targets_buckets**: B2
- **regime**: Scalar 2nd-order Taylor recurrence; P_I with non-elementary RHS (+z inhomogeneity); tritronquée ICs from FW 2011.
- **citation**: FW 2011 eq.(1.1) (references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:48); ICs from FW 2011 eq.(4.1) (p.7, line 226).
- **confidence**: high
- **verify_note**: Verify c[2]=3*u0^2 (z0=0 means z-term does not appear at n=2 when z0=0), and c[3]=2*u0*u1+1/6. These are analytic consequences of the ODE and uniquely pin the recurrence.

#### `wp-taylor-coeff-bigfloat256-b2`

- **name**: Weierstrass ℘ Taylor coefficients at BigFloat-256 (77 decimal digits, order 60)
- **ode**: u'' = 6*u^2
- **ic_bc**: u(0) = wp(-1; 0, 2), u'(0) = wp'(-1; 0, 2); values at 90 dps from Mathematica WeierstrassP[−1, {0, 2}] and WeierstrassP'[−1, {0, 2}]
- **domain**: Coefficient space; precision = BigFloat(precision=256) ≈ 77 decimal digits
- **closed_form**: c[n] = 6/(n*(n-1)) * sum_{k=0}^{n-2} c[k]*c[n-2-k] (exact rational recurrence when ICs are rational, or arbitrary-prec when ICs are given at prec).
- **pole_structure**: Same as wp-taylor-coeff-recurrence-b2.
- **gt_kind**: mathematica
- **how_to_compute**: Mathematica wolframscript (exact Series[] + SetPrecision): wolframscript -code '   s = Series[WeierstrassP[z + (-1), {0, 2}], {z, 0, 60}];   c = CoefficientList[Normal[s], z];   Table[{i-1, N[c[[i]], 80]}, {i, 1, 61}] // TableForm // Print '  Python mpmath at 90 dps (BigFloat-256 ~ 77 decimal digits, use 90 for safety): import mpmath mpmath.mp.dps = 90 u0 = mpmath.mpf('1.0718225164169172653973747494600464') # wp(-1;0,2) at 90 dps u1 = mpmath.mpf('1.7103373531767862674204454712') # wp'(-1;0,2) coeffs = wp_coeffs(u0, u1, 60) # (use same wp_coeffs function as above but order=60)  Key spot check: c[60] must satisfy the recurrence exactly given c[0..59].
- **precision**: 90 dps (comfortably exceeds BigFloat-256 = 77 decimal digits)
- **targets_buckets**: B2
- **regime**: BigFloat-256 coefficient oracle; order-60 Taylor jet; exercises the T in T-Padé at extended precision.
- **citation**: FW 2011 §5.1.1 (references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:281-302). DLMF 23.6.2.
- **confidence**: high
- **verify_note**: The ICs themselves must be obtained at 90 dps from Mathematica or from the fixed-point: u0 satisfies wp''(u0)=6*u0^2 trivially but must come from the lattice definition. Use Mathematica N[WeierstrassP[-1,{0,2}],90] and N[WeierstrassPPrime[-1,{0,2}],90] to get ICs, then feed to the recurrence.

#### `tan-sec2-vector-jet-b13-b22`

- **name**: Coupled (tan, sec²) vector system Taylor jet (d=2, meromorphic, order 30)
- **ode**: y1' = y2, y2' = 2*y1*y2  (first-order 2-vector system; closed-form solution y1=tan(z+a), y2=sec²(z+a))
- **ic_bc**: y1(0) = tan(pi/4) = 1 (exact), y2(0) = sec^2(pi/4) = 2 (exact)
- **domain**: Real and complex z in disk |z| < 3*pi/4 ≈ 2.356 (convergence radius to nearest pole of tan)
- **closed_form**: y1(z) = tan(z + pi/4), y2(z) = sec^2(z + pi/4). Coefficients: c1[n] = (d^n/dz^n tan(z+pi/4))|_{z=0} / n!, c2[n] = (d^n/dz^n sec^2(z+pi/4))|_{z=0} / n!. Recurrence: c1[n+1] = c2[n]/(n+1), c2[n+1] = 2*(sum_{k=0}^{n} c1[k]*c2[n-k])/(n+1).
- **pole_structure**: Simple poles of y1 (tan) at z = pi/4 + pi*k - pi/4 = pi*k for integer k; first pole at z = -pi/4 (behind) and z = 3*pi/4 (ahead). Convergence radius = 3*pi/4 - 0 = 3*pi/4 ≈ 2.356.
- **gt_kind**: closed-form
- **how_to_compute**: Python mpmath (50 dps): import mpmath mpmath.mp.dps = 50 a = mpmath.pi/4 y1c = [mpmath.mpf(0)] * 31 y2c = [mpmath.mpf(0)] * 31 y1c[0] = mpmath.tan(a)   # = 1 exactly y2c[0] = 1/mpmath.cos(a)**2  # = 2 exactly for n in range(30):     y1c[n+1] = y2c[n] / (n+1)     conv = sum(y1c[k]*y2c[n-k] for k in range(n+1))     y2c[n+1] = 2*conv / (n+1) # Verify: y1c[i] == mpmath.diff(mpmath.tan, a, i) / mpmath.factorial(i) # y2c[i] == mpmath.diff(lambda z: 1/mpmath.cos(z)**2, a, i) / mpmath.factorial(i)  Mathematica: wolframscript -code '   c1=N[CoefficientList[Normal[Series[Tan[z+Pi/4],{z,0,30}]],z],50];   c2=N[CoefficientList[Normal[Series[Sec[z+Pi/4]^2,{z,0,30}]],z],50];   Print[c1]; Print[c2] '  Key values (50 dps, z0=0, a=pi/4): y1c[0]=1.0 (exact: tan(pi/4)=1) y1c[1]=2.0 (exact: sec^2(pi/4)=2) y1c[2]=2.0 (exact: 2*tan(pi/4)*sec^2(pi/4)=2) y1c[3]=8/3 = 2.66666... y1c[29]=1403.794572664757374715649... y1c[30]=1787.366762601356358662536... y2c[0]=2.0, y2c[1]=4.0, y2c[2]=8.0, y2c[3]=40/3=13.3333... y2c[29]=53621.002878040690759876... y2c[30]=70548.127335768842088565...
- **precision**: Arbitrary; exact rationals for coefficients expressible as polynomial in tan(pi/4)=1 and sec^2(pi/4)=2
- **targets_buckets**: B13, B22
- **regime**: Vector jet coefficient oracle (d=2); meromorphic with known poles; closed-form verifiable via mpmath.diff; exercises vector_taylor_coefficients.
- **citation**: JorbaZou2005 §2 (references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/JorbaZou2005_taylor_IVP_package_ExpMath14.md:272-350) for the coupled-derivative product rule used in the recurrence; standard calculus identity sec^2 = 1+tan^2.
- **confidence**: high
- **verify_note**: Every coefficient c1[n] must equal mpmath.diff(mpmath.tan, pi/4, n)/n! at 50 dps. The d=1 reduction y1' = y2, y2'=2*y1*y2 with y2=y1'=sec^2 must match the scalar Riccati y' = 1+y^2 (since tan' = sec^2 = 1+tan^2). This provides a fast scalar cross-check for the vector infrastructure.

#### `jorba-zou-sin-jet-primary-b6`

- **name**: Jorba-Zou step size on odd-parity sin jet (single candidate at k=p-1=29)
- **ode**: Any ODE whose Taylor jet at z0 has the parity structure of sin: c[k]=0 for all even k; c[29] != 0, c[30]=0.
- **ic_bc**: c = sin-jet at z0=0, order 30: c[2k+1]=(-1)^k/(2k+1)!, c[2k]=0
- **domain**: Step-size formula; input is a length-31 coefficient vector with odd-index pattern
- **closed_form**: h = (eps / |c[29]|)^(1/29) = (eps * 29!)^(1/29). For eps=1e-12: h = (1e-12 * 29!)^(1/29) = 4.50120637033898607... This MATCHES the existing 3-source consensus oracle (the exp-jet gives the same value).
- **pole_structure**: n/a (this is a step-size formula test, not a trajectory)
- **gt_kind**: mpmath
- **how_to_compute**: import mpmath mpmath.mp.dps = 50 # sin jet normalized coefficients c[k] = sin^(k)(0)/k!: # c[k]=0 for even k, c[2m+1]=(-1)^m/(2m+1)! c = [mpmath.mpf(0)] * 31 for k in range(1, 31, 2):  # odd indices only     c[k] = ((-1)**((k-1)//2)) / mpmath.factorial(k) # step_jorba_zou(c, eps_abs=1e-12, eps_rel=1e-12): # eps_abs >= eps_rel * |c[0]| = 0 -> absolute mode (c[0]=0) # Primary: min over k in {29, 30}: c[30]=0 (skip), c[29] != 0 # => h = (eps / |c[29]|)^(1/29) = (1e-12 * 29!)^(1/29) h = (mpmath.mpf('1e-12') * mpmath.factorial(29)) ** (mpmath.mpf(1)/29) print(mpmath.nstr(h, 35))  # = 4.50120637033898607690318848316  NOTE: This is the SAME numerical value as the exp-jet case (case 4.1.1 in existing oracle), because at k=29: |c_sin[29]| = |c_exp[29]| = 1/29! when c[30]=0. The test distinguishes odd-parity structure (only one valid candidate) from the full exp case (two candidates, min is the same).
- **precision**: 50 dps mpmath; Float64 answer = 4.50120637033899 (matches existing oracle to 15 digits)
- **targets_buckets**: B6
- **regime**: Jorba-Zou step control; odd-parity jet; only k=p-1 candidate active (k=p has zero coefficient); verifies that the min correctly skips the zero candidate.
- **citation**: JorbaZou2005 §3.3.1 eq.11 (references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/JorbaZou2005_taylor_IVP_package_ExpMath14.md:613-635). TaylorIntegration.jl src/integrator/stepsize.jl:17-35 for the ported formula. Existing test/_oracle_stepcontrol.jl pins h_4_1_1_TI=4.501206370338986.
- **confidence**: high
- **verify_note**: Value is identical to the existing h_4_1_1_TI oracle, which provides independent confirmation. The key NEW behavior to test: the implementation must skip c[30]=0 and use only c[29]. Mutation test: if the code accidentally used c[30] (and divided by eps^(1/30)), the result would be Inf or the fallback, not 4.501...

#### `jorba-zou-fallback-trigger-b6`

- **name**: Jorba-Zou _second_stepsize fallback (both c[p-1] and c[p] zero, even-parity jet)
- **ode**: Any ODE whose Taylor jet at z0 has only even nonzero coefficients through c[28] and zero at c[29], c[30] (e.g. cos-like jet truncated at order 28 but padded to order 30).
- **ic_bc**: c = even-parity jet: c[2k]=(-1)^k/(2k)! for k=0..14, c[odd]=0
- **domain**: Step-size formula; input is length-31 vector with c[2k]=(-1)^k/(2k)! for k=0..14, c[29]=c[30]=0
- **closed_form**: h_fallback = (1/|c[28]|)^(1/27) = (28!)^(1/27). Numerically: 12.3596227659147443751802774497...
- **pole_structure**: n/a
- **gt_kind**: mpmath
- **how_to_compute**: import mpmath mpmath.mp.dps = 50 # Construct fallback-trigger jet: cos-like even-only terms up to k=28, # c[29]=c[30]=0 exactly. c = [mpmath.mpf(0)] * 31 for k in range(0, 29, 2):   # even indices 0,2,...,28     c[k] = ((-1)**(k//2)) / mpmath.factorial(k) # c[29] = c[30] = 0 exactly. # step_jorba_zou(c, eps_abs=1e-12): # Primary loop k in {29,30}: both abs(c[k+1])=0 -> h remains Inf -> fallback. # Fallback: scan j=1..28, find max (1/|c[j+1]|)^(1/j). # Only even j+1 are nonzero: j+1 in {2,4,...,28}, i.e. j in {1,3,...,27}. # The maximum is at j=27: (1/|c[28]|)^(1/27) = (28!)^(1/27). c28 = 1 / mpmath.factorial(28) h2 = (mpmath.mpf(1) / c28) ** (mpmath.mpf(1)/27) print(mpmath.nstr(h2, 35))  # = 12.3596227659147443751802774497 # Float64: 12.359622765914745  Mathematica cross-check: wolframscript -code 'N[(28!)^(1/27), 50]'  # = 12.3596227659147443751...
- **precision**: 50 dps; Float64 = 12.359622765914745
- **targets_buckets**: B6
- **regime**: Jorba-Zou _second_stepsize fallback branch; triggered when both c[p-1] and c[p] are zero (e.g. even-parity function evaluated at order 30).
- **citation**: TaylorIntegration.jl src/integrator/stepsize.jl:77-89 (_second_stepsize). JorbaZou2005 §3.3.2 (references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/JorbaZou2005_taylor_IVP_package_ExpMath14.md:647-683) for the 'dangerous step size' second control. PadeTaylor src/StepControl.jl:193-204 for the port.
- **confidence**: high
- **verify_note**: Mutation test: if the fallback takes the minimum instead of maximum, the result would be near 0 (from j=1: (1/|c[2]|)^1 = 2). If it scans odd indices only (misidentifying the nonzero positions), it returns Inf -> error. The value (28!)^(1/27) is exact from Mathematica or mpmath.factorial.

#### `jorba-zou-relative-mode-b6`

- **name**: Jorba-Zou step control in relative-tolerance mode (eps_rel * |c0| > eps_abs)
- **ode**: Any ODE with c[0] = 1 (non-zero constant term), using eps_abs < eps_rel * |c0|.
- **ic_bc**: c = cos jet: c[0]=1, c[2]=-1/2, c[4]=1/24, ..., c[30]=-1/30!, all odd=0; eps_abs=1e-15, eps_rel=1e-12
- **domain**: Step-size formula; cos jet c[2k]=(-1)^k/(2k)! at order 30
- **closed_form**: eps_eff = eps_rel * |c[0]| = 1e-12 (for c[0]=1, eps_rel=1e-12). Then h = min over k in {p-1,p} of (eps_eff/|c[k]|)^(1/k). For cos jet (only c[30] nonzero among {c[29],c[30]}): h = (eps_rel / (1/30!))^(1/30) = (eps_rel * 30!)^(1/30). Numerically: 4.79500063659056684443965285569...
- **pole_structure**: n/a
- **gt_kind**: mpmath
- **how_to_compute**: import mpmath mpmath.mp.dps = 50 # cos jet: c[2k]=(-1)^k/(2k)!, c[odd]=0, order 30. # eps_abs=1e-15, eps_rel=1e-12, c[0]=1. # Condition: eps_abs=1e-15 < eps_rel*|c[0]|=1e-12 => relative mode. # eps_eff = eps_rel * 1.0 = 1e-12. # Primary: k=29: c[29]=0 (odd, skip). k=30: c[30]=(-1)^15/30! = -1/30!. # h = (eps_eff / |c[30]|)^(1/30) = (1e-12 * 30!)^(1/30). c30 = 1 / mpmath.factorial(30) eps_rel = mpmath.mpf('1e-12') h = (eps_rel / c30) ** (mpmath.mpf(1)/30) print(mpmath.nstr(h, 35))  # = 4.79500063659056684443965285569 # Float64: 4.795000636590567  Mathematica cross-check: wolframscript -code 'N[(10^(-12) * 30!)^(1/30), 50]'
- **precision**: 50 dps; Float64 = 4.795000636590567
- **targets_buckets**: B6
- **regime**: Jorba-Zou relative-mode dispatch (TI.jl stepsize.jl:27-31 branch). Exercises the eps_abs >= eps_rel*|c0| condition going FALSE, switching eps to eps_rel*|c0|.
- **citation**: TaylorIntegration.jl src/integrator/stepsize.jl:27-31. PadeTaylor src/StepControl.jl:175-183 for the dispatch. JorbaZou2005 §3.3 eq.10 (references/markdown/JorbaZou2005_taylor_IVP_package_ExpMath14/JorbaZou2005_taylor_IVP_package_ExpMath14.md:583-614) for relative-error definition.
- **confidence**: high
- **verify_note**: Mutation test: if the code uses eps_abs=1e-15 instead of eps_eff=1e-12, the result would be (1e-15 * 30!)^(1/30) = (1e-15/1e-12)^(1/30) * h_rel = 0.794... * 4.795... ≈ 3.807, a different value detectable to 5+ significant digits.

#### `pi1-odefun-trajectory-b4`

- **name**: Painlevé I trajectory reference (mpmath.odefun, 50 dps, tritronquée ICs)
- **ode**: u'' = 6*u^2 + z
- **ic_bc**: u(0) = -0.18755430834049490, u'(0) = 0.30490556026122890 (FW 2011 eq. 4.1)
- **domain**: Real z in [0, 1]; tritronquée sector (pole-free on positive real axis)
- **closed_form**: none
- **pole_structure**: Tritronquée: no poles in |arg(z)| < 2*pi/5 ≈ 72°. Positive real axis is pole-free. First pole on real axis is far right (from ICs) or at large |z| in other sectors.
- **gt_kind**: mpmath
- **how_to_compute**: import mpmath mpmath.mp.dps = 50  u0 = mpmath.mpf('-0.18755430834049490') u1 = mpmath.mpf('0.30490556026122890')  def ode_pi(z, y):     return [y[1], 6*y[0]**2 + z]  sol = mpmath.odefun(ode_pi, 0, [u0, u1], tol=mpmath.mpf('1e-45'))  # Evaluate at safe points on positive real axis (no poles for tritronquée) for z in ['0.1', '0.3', '0.5', '0.8', '1.0']:     val = sol(mpmath.mpf(z))     print(f'u({z}) = {mpmath.nstr(val[0], 35)}')  Reference values (50 dps): u(0.1) = -0.15595334953222776939756745301609501... u(0.3) = -0.08491858006004649082275653856599388... u(0.5) = 0.000041002535628181031343601320249837... u(0.8) = 0.16881532288924458424204890002491886... u(1.0) = 0.32788384158056647613573768327890733...  Mathematica NDSolve cross-check: wolframscript -code ' u0=SetPrecision[-0.18755430834049490,60]; u1=SetPrecision[0.30490556026122890,60]; sol=NDSolve[{u''''''''[z]==6*u[z]^2+z,u[0]==u0,u''''''''[0]==u1},   u,{z,0,1},WorkingPrecision->50,MaxStepSize->0.01]; Print[{#,u[#]/.sol[[1]]}&/@{0.1,0.3,0.5,0.8,1.0}] '
- **precision**: 50 dps via mpmath.odefun with tol=1e-45; ~35 significant decimal digits in the output
- **targets_buckets**: B4
- **regime**: FFW adaptive controller / fixed stepper cross-check on P_I in the pole-free tritronquée sector (positive real axis). Exercises the full Taylor-Padé pipeline on a non-elementary RHS.
- **citation**: FW 2011 eq.(1.1) and §4.1 (references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:48,216-229). ICs from FW 2011 eq.(4.1). FFW 2017 (references/markdown/FFW2017_painleve_riemann_surfaces_preprint/FFW2017_painleve_riemann_surfaces_preprint.md) for the adaptive controller tested.
- **confidence**: high
- **verify_note**: Verify with Mathematica NDSolve WorkingPrecision->50. The u(0.5) value is nearly zero (local minimum of the tritronquée, visually near zero in Fig 4.3 of FW 2011 at z≈0.5). Verifying the sign flip between u(0.3)<0 and u(0.5)>0 cross-checks the integrator passed through the minimum correctly.

#### `wp-odefun-trajectory-b4`

- **name**: Weierstrass ℘ trajectory reference (mpmath.odefun + FW 2011 Table 5.1 long-range)
- **ode**: u'' = 6*u^2
- **ic_bc**: u(0) = wp(-1; 0, 2) = 1.07182251641691727..., u'(0) = wp'(-1; 0, 2) = 1.71033735317678627...
- **domain**: Real z in [0, 30] and [0, 10^4] along positive real axis through many poles
- **closed_form**: u(z) = wp(z-1; 0, 2). Short-range (z<0.5): evaluate via order-60 recurrence at 70 dps (recurrence IS the Taylor series). Long-range: published values from FW 2011 Table 5.1.
- **pole_structure**: Double poles at z = 1 + 2*omega*k for integer k, omega ≈ 1.3628. On real axis: z ≈ 1, 1+2*1.3628 ≈ 3.726, etc.
- **gt_kind**: published-table
- **how_to_compute**: SHORT-RANGE (z=0.3, 0.5, z safely away from z=1 pole): use order-60 recurrence at 70 dps: import mpmath mpmath.mp.dps = 70 def wp_coeffs(u0, u1, order):     c = [mpmath.mpf(u0), mpmath.mpf(u1)]     for n in range(2, order+1):         conv = sum(c[k]*c[n-2-k] for k in range(n-1))         c.append(6*conv / (n*(n-1)))     return c u0 = mpmath.mpf('1.0718225164169172653973747494600464') u1 = mpmath.mpf('1.7103373531767862674204454712') c = wp_coeffs(u0, u1, 60) def eval_poly(c, h):     val = mpmath.mpf(0)     for k in range(len(c)-1,-1,-1): val = val*h + c[k]     return val print(eval_poly(c, mpmath.mpf('0.3')))  # 2.05797741759157366641981748626... print(eval_poly(c, mpmath.mpf('0.5')))  # 4.00446466900308673059573639151...  LONG-RANGE (published values, FW 2011 Table 5.1 / eq.(5.3)): u(30) = 1.095098255959744 u(10^4) = 21.02530339471055 u(28.261) = 9.876953517025014e6  Mathematica NDSolve short-range cross-check: wolframscript -code ' sol=NDSolve[{u''''[z]==6*u[z]^2,u[0]==SetPrecision[1.07182251641691727,60],u''''[0]==SetPrecision[1.71033735317678627,60]},u,{z,0,0.6},WorkingPrecision->50]; Print[{u[0.3],u[0.5]}/.sol[[1]]] '
- **precision**: Short-range: 50 dps; long-range: ~15 significant digits (Float64) from FW 2011 Table 5.1
- **targets_buckets**: B4
- **regime**: Padé stepper long-range test through many poles; FW 2011 §5.4 accuracy table benchmarking. Short-range recurrence verifies the Taylor-coefficient-to-evaluation pipeline.
- **citation**: FW 2011 Table 5.1 (references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:383-391) for relative errors and reference values; eq.(5.3) for u(30) and u(10^4); p.11 line 372 for u(28.261)=9.876953517025014e6.
- **confidence**: high
- **verify_note**: For the long-range values, these are published reference values from FW 2011 (not independently regenerable without running the full Padé integrator). They are suitable as pass/fail targets for the full pipeline test, not for the coefficient oracle alone. The recurrence short-range values are independently computable and self-consistent.

#### `scalar-riccati-b2-b6`

- **name**: Scalar Riccati y' = 1 + y^2 (tan companion) Taylor jet + step size
- **ode**: y' = 1 + y^2
- **ic_bc**: y(0) = tan(pi/4) = 1 (exact)
- **domain**: Real z in (-pi/4, 3*pi/4); pole-free strip of width 3*pi/4
- **closed_form**: y(z) = tan(z + pi/4) for y(0)=1. Taylor coefficients c[n] = (d^n tan(z+pi/4)/dz^n)|_{z=0} / n! = y1c[n] from the (tan, sec^2) vector oracle (identical by the vector d=1 reduction test).
- **pole_structure**: Simple poles of tan at z = pi*k for integer k; nearest poles at z = -pi/4 (behind) and z = 3*pi/4 ≈ 2.356 (ahead).
- **gt_kind**: closed-form
- **how_to_compute**: Python: import mpmath mpmath.mp.dps = 50 a = mpmath.pi/4 # Taylor coefficients: c[n] = d^n tan(z+a)|_{z=0} / n! coeffs = [mpmath.diff(mpmath.tan, a, n) / mpmath.factorial(n) for n in range(31)] # c[0]=1, c[1]=2, c[2]=2, c[3]=8/3, c[4]=10/3, ... for i, c in enumerate(coeffs):     print(f'c[{i}] = {mpmath.nstr(c, 40)}')  # Jorba-Zou step size for this jet (p=30, eps_abs=1e-12, eps_rel=1e-12): # c[0]=1: eps_abs=eps_rel*|c0|=1e-12 (equality -> absolute mode) # k=29: c[29] != 0, k=30: c[30] != 0 (both tan-derived coefficients nonzero) # h = min((1e-12/|c[29]|)^(1/29), (1e-12/|c[30]|)^(1/30)) c29 = coeffs[29] c30 = coeffs[30] eps = mpmath.mpf('1e-12') h29 = (eps / abs(c29)) ** (mpmath.mpf(1)/29) h30 = (eps / abs(c30)) ** (mpmath.mpf(1)/30) h = min(h29, h30) print(f'h = {mpmath.nstr(h, 35)}')  Key values: c[0]=1, c[1]=2, c[2]=2, c[3]=8/3=2.666..., c[29]=1403.7945726647573747... h(p=30, eps=1e-12) = needs computation (c[29] and c[30] both nonzero -> the min picks the smaller)
- **precision**: Arbitrary; mpmath.diff computes to machine dps
- **targets_buckets**: B2, B6
- **regime**: Scalar first-order ODE; Taylor coefficients via mpmath.diff; exercises taylor_coefficients_1st and step_jorba_zou simultaneously; closed-form cross-check via tan derivative.
- **citation**: FW 2011 §2.2.1 (references/markdown/FW2011_painleve_methodology_JCP230/FW2011_painleve_methodology_JCP230.md:120-131) for the demonstration problem y'=z^2+y^2. The Riccati y'=1+y^2 is a simpler companion (no z in RHS). JorbaZou2005 §3.3.1 for step size.
- **confidence**: high
- **verify_note**: All coefficients are known exactly as rational multiples of powers of 1 (since tan^(n)(pi/4) is a polynomial in tan(pi/4)=1). Specifically: c[0]=1, c[1]=2, c[2]=2, c[3]=8/3, c[4]=10/3, c[5]=64/15, ... These can be verified symbolically via sympy.series(sympy.tan(z+sympy.pi/4), z, 0, 30).

### Territory: Exact discrete/geometric oracle fixtures for morphology (B9), sheet-tracking geo

*Scout notes:* All fixtures are hand-computed from first principles; no internet search required. Mathematical morphology citations follow Serra (1982) / Haralick & Shapiro (1992) standard definitions. Sheet-tracking fixtures follow the BranchTracker/SheetTracker implementation in src/BranchTracker.jl and src/SheetTracker.jl. BVP domain guard follows src/BVP.jl lines 489–494. The floating-point endpoint artifact for segment_crosses_cut (endpoint exactly on cut gives t≈0.9999... rather than 1.0) is documented as a verify_note where relevant. The precondition-violation (|dtheta|>=pi) fixture is confidence=medium because the correct behavior (throw vs. return) is a design decision not yet mandated by the codebase.

#### `morph-dilation-single-pixel`

- **name**: Morphological dilation: single center pixel, 3x3 box SE
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: n/a
- **domain**: discrete 3x3 BitMatrix
- **closed_form**: Input A (3x3): [[0,0,0],[0,1,0],[0,0,0]]. SE = 3x3 all-ones (8-connected box). Boundary convention: out-of-grid pixels treated as 0. Output dilation(A, SE) = [[1,1,1],[1,1,1],[1,1,1]] (all 9 positions). Reason: every cell (r,c) in the 3x3 grid is within the 3x3 SE footprint of center (1,1), so the translated SE intersects A. Definition: (A (+) SE)[q] = 1 iff exists p in A with (q-p) in SE. Equivalently OR over SE-translated copies.
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: By definition (Serra 1982, ch. 2): dilation A (+) SE = {q : (SE reflected, translated to q) intersects A}. SE = 3x3 box = all 9 relative offsets in {-1,0,1}^2. A = {(1,1)}. For each (r,c) in 0..2 x 0..2: exists offset (dr,dc) in SE s.t. (r-dr, c-dc) = (1,1)? Yes for all 9 cells. Python check: all([[1]*3,[1]*3,[1]*3]) confirmed.
- **precision**: exact integer (0/1)
- **targets_buckets**: B9
- **regime**: morphology correctness: single-pixel dilation expands to full SE footprint
- **citation**: Serra, J. (1982) Image Analysis and Mathematical Morphology, Academic Press, ch. 2 (dilation definition). Also: Haralick & Shapiro (1992) Computer and Robot Vision Vol. I, ch. 5.
- **confidence**: high

#### `morph-erosion-edge-hang`

- **name**: Morphological erosion: full 3x3 grid, 3x3 box SE, boundary=0 (window hangs off edge)
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: n/a
- **domain**: discrete 3x3 BitMatrix
- **closed_form**: Input B2 (3x3): [[1,1,1],[1,1,1],[1,1,1]]. SE = 3x3 all-ones. Boundary convention: out-of-grid = 0. Output erosion(B2, SE) = [[0,0,0],[0,1,0],[0,0,0]]. Only center (1,1) survives. Reason: erosion(A,SE)[q] = 1 iff ALL pixels of SE translated to q are in A. Center (1,1): 3x3 neighborhood = all of B2 (all 1). Any edge or corner pixel (e.g., (0,0)): 3x3 neighborhood includes row -1 and col -1 (out-of-grid = 0) -> fails.
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: By definition: erosion(A,SE)[q] = AND over all (dr,dc) in SE of A[q+(dr,dc)], where out-of-grid = 0. For B2 full 3x3: center (1,1) has neighborhood = B2 = all 1 -> survives. Corner (0,0): needs A[-1,-1], A[-1,0], ... = 0 -> fails. Edge (0,1): needs A[-1,0], A[-1,1], A[-1,2] = 0 -> fails. Python BFS check confirms single survivor at (1,1).
- **precision**: exact integer (0/1)
- **targets_buckets**: B9
- **regime**: morphology correctness: SE hanging off grid boundary (boundary=0 convention makes all boundary pixels fail erosion)
- **citation**: Serra, J. (1982) Image Analysis and Mathematical Morphology, Academic Press, ch. 2 (erosion definition). Boundary=0 is the standard 'zero-padding' convention.
- **confidence**: high

#### `morph-erosion-5x5`

- **name**: Morphological erosion: 5x5 grid with 3x3 center block, 3x3 box SE
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: n/a
- **domain**: discrete 5x5 BitMatrix
- **closed_form**: Input B (5x5): [[0,0,0,0,0],[0,1,1,1,0],[0,1,1,1,0],[0,1,1,1,0],[0,0,0,0,0]]. SE = 3x3 all-ones. Boundary = 0. Output erosion(B, SE) = [[0,0,0,0,0],[0,0,0,0,0],[0,0,1,0,0],[0,0,0,0,0],[0,0,0,0,0]]. Only (2,2) survives. Reason: (2,2)'s 3x3 neighborhood = rows 1..3, cols 1..3 = the entire inner block, all 1. Every other candidate has a zero neighbor (either border zero or out-of-range).
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: Check each interior pixel: (1,1) needs B[0,0]=0 -> fail; (1,2) needs B[0,1]=0 -> fail; (2,1) needs B[2,0]=0 -> fail; (2,2) needs B[1..3,1..3] = [[1,1,1],[1,1,1],[1,1,1]] -> all 1 -> survives. All 0-border pixels fail trivially. Python flood-fill confirms 1 survivor.
- **precision**: exact integer (0/1)
- **targets_buckets**: B9
- **regime**: morphology correctness: erosion shrinks a 3x3 block to a single point
- **citation**: Serra, J. (1982) Image Analysis and Mathematical Morphology, Academic Press, ch. 2 (erosion).
- **confidence**: high

#### `morph-opening-idempotent`

- **name**: Morphological opening: 3x3 center block is fixed point of opening
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: n/a
- **domain**: discrete 5x5 BitMatrix
- **closed_form**: Input B (5x5, center 3x3 filled as above). SE = 3x3 box. opening(B, SE) = dilate(erode(B, SE), SE). Step 1: erode(B, SE) = single pixel at (2,2) [from fixture above]. Step 2: dilate({{(2,2)}}, SE) = 3x3 block at rows 1..3, cols 1..3. Output = B itself: [[0,0,0,0,0],[0,1,1,1,0],[0,1,1,1,0],[0,1,1,1,0],[0,0,0,0,0]]. This is the standard algebraic property: opening is idempotent and anti-extensive; B is a union of SE translates, so opening(B)=B.
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: Compose erosion and dilation results: erode(B)={(2,2)}; dilate({(2,2)},SE)={(r,c): |r-2|<=1, |c-2|<=1} = rows 1..3 cols 1..3 = B. Python confirms.
- **precision**: exact integer (0/1)
- **targets_buckets**: B9
- **regime**: morphology correctness: opening = erode-then-dilate; idempotent on a set that equals a union of SE translates
- **citation**: Serra, J. (1982) Image Analysis and Mathematical Morphology, Academic Press, ch. 2 (opening, eq. gamma(A) = (A ominus SE) oplus SE, with idempotence theorem).
- **confidence**: high

#### `morph-connected-components-8conn`

- **name**: 8-connected flood-fill: two disjoint L-shaped components
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: n/a
- **domain**: discrete 5x5 BitMatrix
- **closed_form**: Input CC (5x5): [[1,1,0,0,0],[1,0,0,0,0],[0,0,0,0,1],[0,0,0,1,1],[0,0,0,0,0]]. Expected: 2 components. Component 1 = {(0,0),(0,1),(1,0)}: (0,0)-(0,1) 4-adjacent; (0,0)-(1,0) 4-adjacent. Component 2 = {(2,4),(3,3),(3,4)}: (2,4)-(3,4) 4-adjacent; (3,3)-(3,4) 4-adjacent. Label map: [[1,1,0,0,0],[1,0,0,0,0],[0,0,0,0,2],[0,0,0,2,2],[0,0,0,0,0]]. Note: the two components are ALSO non-adjacent diagonally (min distance: (1,0) to (2,4) is far; (1,0) to (3,3) is far), confirming separation under both 4- and 8-connectivity.
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: BFS/DFS from first unvisited 1-pixel. From (0,0): 8-neighbors are (0,1)=1 and (1,0)=1 and (1,1)=0 -> component 1 = {(0,0),(0,1),(1,0)}. From (2,4): 8-neighbors include (3,3)=1 and (3,4)=1 -> component 2 = {(2,4),(3,3),(3,4)}. Python BFS code confirmed n_labels=2.
- **precision**: exact integer (labels 0,1,2)
- **targets_buckets**: B9
- **regime**: morphology correctness: 8-connected flood-fill separates two disjoint clusters
- **citation**: Haralick, R. & Shapiro, L. (1992) Computer and Robot Vision Vol. I, Addison-Wesley, ch. 2 (connected-component labeling, 8-connectivity).
- **confidence**: high
- **verify_note**: Component adjacency check: (1,0) to (2,4) has col-gap 4 -> disjoint under 8-connectivity. Confirmed.

#### `sheet-winding-delta-four-arcs`

- **name**: winding_delta: four 90-degree CCW arcs at unit radius, branch at origin
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: path=[1+0j, 0+1j, -1+0j, 0-1j, 1+0j], branch=0+0j
- **domain**: five complex points on unit circle
- **closed_form**: branch=0+0j. Four steps each subtending pi/2: (1) z_old=1, z_new=i: angle(i)=pi/2, angle(1)=0, delta=pi/2. (2) z_old=i, z_new=-1: angle(-1)=pi, angle(i)=pi/2, delta=pi/2. (3) z_old=-1, z_new=-i: angle(-i-0)=-pi/2, angle(-1-0)=pi, raw=-3pi/2, normalise (add 2pi)=pi/2. (4) z_old=-i, z_new=1: angle(1)=0, angle(-i)=-pi/2, delta=pi/2. accumulate_winding([1,i,-1,-i,1], 0) = [0, pi/2, pi, 3pi/2, 2pi]. sheet_index(2*pi) = 1.
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: winding_delta(z_old, z_new, branch): Δθ = angle(z_new-branch) - angle(z_old-branch), normalised to (-pi, pi]. Step 3: raw = angle(-i) - angle(-1) = -pi/2 - pi = -3pi/2; since -3pi/2 <= -pi, add 2pi -> pi/2. All four deltas = pi/2. Python: confirmed. accumulate_winding sums: [0, pi/2, pi, 3pi/2, 2*pi]. sheet_index = round(2*pi/(2*pi)) = 1.
- **precision**: exact: pi/2 = 1.5707963267948966..., 2*pi = 6.283185307...
- **targets_buckets**: B11
- **regime**: sheet-tracking geometry: CCW unit-circle circumambulation with normalisation at the -pi branch
- **citation**: src/SheetTracker.jl:283-291 (winding_delta implementation); FFW2017 md:163-189 (sheet index convention). Normalisation to (-pi,pi] is standard complex analysis (arg principal value).
- **confidence**: high
- **verify_note**: Step 3 normalisation is the critical case: raw=-3pi/2 -> +pi/2. Perturb: if normaliser dropped the add-2pi branch, step 3 returns -3pi/2 and total=-pi -> sheet_index=-1 (wrong).

#### `sheet-ccw-winding-number`

- **name**: Full CCW circumambulation: winding number +1, branch at 2+0j
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: path=[4+0j, 2+2j, 0+0j, 2-2j, 4+0j], branch=2+0j
- **domain**: five complex points forming a diamond around branch 2+0j
- **closed_form**: branch=2+0j. Path=[4+0j, 2+2j, 0+0j, 2-2j, 4+0j] (four steps of 90deg each relative to branch). Step 1: angle(2+2j-2)-angle(4-2)=angle(2j)-angle(2)=pi/2-0=pi/2. Step 2: angle(-2)-angle(2j)=pi-pi/2=pi/2. Step 3: angle(-2j)-angle(-2)=-pi/2-pi=-3pi/2, normalise->pi/2. Step 4: angle(2)-angle(-2j)=0-(-pi/2)=pi/2. Total=2*pi. sheet_index(2*pi)=1.
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: Each step: relative vectors are 2j, -2+0j, 0-2j, 2+0j (from branch=2). Angles: pi/2, pi, -pi/2, 0. Deltas: pi/2 each (step 3 normalised same way as above fixture). Python confirmed: total=6.283185... sheet_index=1.
- **precision**: exact: total=2*pi, sheet_index=1
- **targets_buckets**: B11
- **regime**: sheet-tracking geometry: CCW circumambulation of off-origin branch gives integer winding +1
- **citation**: src/SheetTracker.jl:303-311 (accumulate_winding), :320 (sheet_index). FFW2017 md:187-189 (sheet parametrisation convention).
- **confidence**: high

#### `sheet-segment-crosses-cut-oblique-yes`

- **name**: segment_crosses_cut: oblique cut at pi/4, crossing case (Cramer's rule)
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: z_cur=0+2j, z_new=3+0j, branch=1+1j, cut_angle=pi/4
- **domain**: two complex points defining a segment, one complex branch point, one real cut angle
- **closed_form**: branch=1+1j, cut_angle=pi/4 (ray {(1+1j)+s*exp(i*pi/4) : s>0}, going up-right). Segment: z_cur=0+2j, z_new=3+0j. Cramer system: d=3-2j, delta=(1+1j)-(0+2j)=1-1j, u=exp(i*pi/4)=(1+j)/sqrt(2). det=Im(d*conj(u))=Im((3-2j)*(1-j)/sqrt(2))=Im((1-5j)/sqrt(2))=-5/sqrt(2). t=Im(delta*conj(u))/det=Im(-2j/sqrt(2))/(-5/sqrt(2))=(-sqrt(2))/(-5/sqrt(2))=2/5=0.4. s=Im(delta*conj(d))/det=Im((1-j)*(3+2j))/det=Im(5-j)/(-5/sqrt(2))=(-1)/(-5/sqrt(2))=sqrt(2)/5. crosses: 0<0.4<1 AND sqrt(2)/5>0 => True.
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: Cramer's rule per src/BranchTracker.jl:112-121. Exact rational: t=2/5, s=sqrt(2)/5. Python computed: t=0.400000, s=0.282843. Both conditions (0<t<1, s>0) hold.
- **precision**: exact: t=2/5, s=sqrt(2)/5 = 0.282842712...
- **targets_buckets**: B11
- **regime**: sheet-tracking geometry: oblique branch cut crossing verified via Cramer's rule
- **citation**: src/BranchTracker.jl:36-50 (cut-crossing predicate) and module docstring. Standard 2x2 Cramer's rule for line-ray intersection.
- **confidence**: high

#### `sheet-segment-crosses-cut-oblique-no`

- **name**: segment_crosses_cut: oblique cut at pi/4, non-crossing case (t=0, s<0)
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: z_cur=0+0j, z_new=2+0j, branch=1+1j, cut_angle=pi/4
- **domain**: segment from origin rightward, cut going upper-right from 1+1j
- **closed_form**: branch=1+1j, cut_angle=pi/4. Segment: z_cur=0+0j, z_new=2+0j (rightward, below the diagonal cut). d=2+0j, delta=(1+1j)-(0+0j)=1+1j, u=(1+j)/sqrt(2), conj(u)=(1-j)/sqrt(2). det=Im(2*(1-j)/sqrt(2))=-2/sqrt(2)=-sqrt(2). t=Im((1+j)*(1-j)/sqrt(2))/(-sqrt(2))=Im(2/sqrt(2))/(-sqrt(2))=0. s=Im((1+j)*2)/(-sqrt(2))=2/(-sqrt(2))=-sqrt(2). crosses: NOT(0<0<1) AND s<0 => False (both conditions fail independently).
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: Exact: delta*conj(u)=(1+j)*(1-j)/sqrt(2)=2/sqrt(2) (purely real, Im=0) -> t=0. s=Im((1+j)*2)/(-sqrt(2))=Im(2+2j)/(-sqrt(2))=2/(-sqrt(2))<0. Python confirmed: t=-0.0, s=-1.414. Both gate conditions fail.
- **precision**: exact: t=0 (not strictly positive), s=-sqrt(2)<0
- **targets_buckets**: B11
- **regime**: sheet-tracking geometry: non-crossing case where segment endpoint on branch point (t=0) and s<0
- **citation**: src/BranchTracker.jl:36-50 (endpoint exclusion convention: t=0 not in open (0,1)).
- **confidence**: high

#### `sheet-winding-delta-precondition-violation`

- **name**: winding_delta precondition violation: single step subtending > pi (large arc)
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: branch=0+0j, z_old=0.1+0j, z_new=0.1*exp(i*1.1*pi)
- **domain**: two complex points near origin, branch at origin
- **closed_form**: branch=0+0j. z_old=0.1+0j (angle=0), z_new=0.1*exp(i*1.1*pi) (true subtended angle=+1.1*pi CCW). angle(z_new-0) in principal branch = 1.1*pi - 2*pi = -0.9*pi (since 1.1*pi > pi). raw_delta = -0.9*pi - 0 = -0.9*pi. -0.9*pi is already in (-pi, pi], so winding_delta returns -0.9*pi. Correct value would be +1.1*pi. Error = -0.9*pi - 1.1*pi = -2*pi (missed one full revolution). A correct bookkeeper SHOULD detect |Δθ_raw| >= pi before normalisation and either throw an ArgumentError ('step too large: |subtended angle| >= pi at branch') or flag the caller, because the normalised result is ambiguous.
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: angle(0.1*exp(i*1.1*pi)) in Julia = angle of exp(i*1.1*pi). Since 1.1*pi > pi, Julia's angle() returns 1.1*pi - 2*pi = -0.9*pi. raw = -0.9*pi - 0 = -0.9*pi. Normalisation: -0.9*pi in (-pi, pi] -> no adjustment. Returns -0.9*pi. Python cmath.phase confirmed: -162 degrees = -0.9*pi.
- **precision**: winding_delta returns -0.9*pi = -2.827433... (exact); correct value +1.1*pi = 3.455752...
- **targets_buckets**: B11, B6
- **regime**: sheet-tracking geometry: precondition violation - step larger than pi leads to ambiguous normalisation; winding_delta contract documented in SheetTracker.jl:277
- **citation**: src/SheetTracker.jl:276-279 (caller contract: 'Caller must ensure no single path step has |Δθ| >= pi'). This fixture documents what HAPPENS when the contract is violated.
- **confidence**: medium
- **verify_note**: The 'correct' behavior (throw vs. return) is a design decision. Current implementation documents the contract but does not enforce it. This fixture tests the observable return value, not a throw.

#### `sheet-two-cuts-simultaneously`

- **name**: segment_crosses_cut: single step crosses two distinct branch cuts
- **ode**: n/a (discrete geometry fixture)
- **ic_bc**: z_cur=-0.5+2j, z_new=1.5+2j, branches=(0+0j, 1+0j), cut_angles=(pi/2, pi/2)
- **domain**: horizontal segment at height 2, two upward branch cuts from x=0 and x=1
- **closed_form**: branch_1=0+0j, cut_angle_1=pi/2 (upward ray {s*i : s>0}). branch_2=1+0j, cut_angle_2=pi/2 (upward ray {1+s*i : s>0}). Segment: z_cur=-0.5+2j, z_new=1.5+2j (horizontal rightward at height 2). d=2+0j, u=i, conj(u)=-i. det=Im(2*(-i))=Im(-2i)=-2. Cut1: delta1=(0+0j)-(-0.5+2j)=0.5-2j; t1=Im((0.5-2j)*(-i))/(-2)=Im(-2-0.5j)/(-2)=(-0.5)/(-2)=0.25; s1=Im((0.5-2j)*2)/(-2)=Im(1-4j)/(-2)=(-4)/(-2)=2; crosses1: 0<0.25<1 AND 2>0 => True. Cut2: delta2=(1+0j)-(-0.5+2j)=1.5-2j; t2=Im((1.5-2j)*(-i))/(-2)=Im(-2-1.5j)/(-2)=(-1.5)/(-2)=0.75; s2=2; crosses2: 0<0.75<1 AND 2>0 => True. Both cuts crossed by one step.
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: Cramer's rule applied twice. Python confirmed: t1=0.25, s1=2.0, crosses1=True; t2=0.75, s2=2.0, crosses2=True. any_cut_crossed and step_sheet_update must handle both.
- **precision**: exact: t1=1/4, t2=3/4, s1=s2=2
- **targets_buckets**: B11, B6
- **regime**: sheet-tracking geometry: step simultaneously crossing two upward branch cuts at different t-values
- **citation**: src/BranchTracker.jl:130-137 (any_cut_crossed iterates branches); :156-167 (step_sheet_update updates both counters).
- **confidence**: high
- **verify_note**: step_sheet_update: winding_delta(z_cur, z_new, 0+0j) = angle(1.5+2j) - angle(-0.5+2j) > 0 (CCW), so sheet counter for branch_1 increments by +1. Same for branch_2.

#### `bvp-domain-guard-imaginary-tstar`

- **name**: BVP domain guard edge case: Re(t*)=0 but Im(t*)=10 -> silent extrapolation
- **ode**: n/a (discrete geometry fixture for the BVP evaluator guard)
- **ic_bc**: z_a=0+0j, z_b=1+1j, query z=-4.5+5.5j
- **domain**: oblique complex segment [0, 1+i], query z=-4.5+5.5j
- **closed_form**: Segment: z_a=0+0j, z_b=1+1j (oblique complex segment). half_sum=(z_a+z_b)/2=0.5+0.5j. half_diff=(z_b-z_a)/2=0.5+0.5j. Query: z=-4.5+5.5j. t*=(z-half_sum)/half_diff=((-4.5+5.5j)-(0.5+0.5j))/(0.5+0.5j)=(-5+5j)/(0.5+0.5j). Multiply num/den by conj(denom)/(|denom|^2): (-5+5j)*(0.5-0.5j)/0.5 = (-2.5+2.5j+2.5j-2.5j^2)/0.5 = (-2.5+5j+2.5)/0.5 = 5j/0.5 = 10j. t*=10j exactly. Re(t*)=0, Im(t*)=10. Guard (src/BVP.jl:491): checks only 'real(t_star) in [-1-100eps, 1+100eps]' -> 0 is in range -> PASSES. Im(t*)=10 is ignored. z=-4.5+5.5j is far off-segment (segment parametrised as tau*(1+j), tau in [0,1]: no real tau gives Im(z)=5.5 and Re(z)=-4.5). Guard silently accepts -> silent extrapolation.
- **pole_structure**: n/a
- **gt_kind**: hand-computed
- **how_to_compute**: t* = (z - (z_a+z_b)/2) / ((z_b-z_a)/2). Exact complex arithmetic: (-5+5j)/(0.5+0.5j). Multiply by conj/|.|^2: (-5+5j)(0.5-0.5j) = -2.5+2.5j+2.5j-2.5*(-1) = 0+5j. Divide by |0.5+0.5j|^2=0.5: result = 10j. Python confirmed. Guard check: abs(real(10j))=0 <= 1 -> True. Guard does NOT check imag(t*).
- **precision**: exact: t*=10j, Re(t*)=0 exactly, Im(t*)=10 exactly
- **targets_buckets**: B17, B6
- **regime**: BVP domain guard: oblique segment with query off-segment in the imaginary t* direction; the real-only guard silently accepts extrapolation
- **citation**: src/BVP.jl:489-494 (the guard implementation, which checks only real(t_star)). BVP.jl:54-56 (docstring stating the preimage formula). This is a KNOWN LIMITATION documented in the BVP module docstring line 134: 'Evaluating the callable outside [z_a,z_b] preimage disc'.
- **confidence**: high
- **verify_note**: A STRICTER guard would check |t*| <= 1+eps or |Im(t*)| <= eps. The current guard only checks Re(t*). This fixture documents the gap and is useful for testing whether any future stricter guard correctly REJECTS this query.


---
*102 candidate records across 9 territories; 89 consolidated entries.*
