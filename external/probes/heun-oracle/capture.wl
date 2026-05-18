(* capture.wl — Heun oracle capture for PadeTaylor.jl

   Run:  wolframscript -file capture.wl 2>/dev/null | grep -v "^Prefetching" > oracles.txt
   Or simply:  ./run-capture.sh

   Emits one record per line:
     FuncName|key=val|...|z=val|value=re±imi

   All numeric inputs are 50-digit arbitrary-precision Mathematica numbers.
   Values are emitted at 30 significant digits.
*)

(* ── Formatter ────────────────────────────────────────────────────────────── *)

formatReal[x_] := ToString[NumberForm[N[x, 30], 30, ExponentFunction -> (Null &)]]

formatComplex[z_] := Module[{re, im, sign},
  re   = Re[z]; im = Im[z];
  sign = If[im >= 0, "+", "-"];
  formatReal[re] <> sign <> formatReal[Abs[im]] <> "i"
]

(* ── Single-record emitter (errors → "ERROR:<msg>" record, no truncation) ── *)

emit[func_String, paramStr_String, zStr_String, val_] :=
  Module[{vstr},
    vstr = If[NumericQ[val] || (Head[val] === Complex && NumericQ[Re[val]] && NumericQ[Im[val]]),
              formatComplex[val],
              "ERROR:nonnumeric:" <> ToString[InputForm[val]]];
    Print[func <> "|" <> paramStr <> "|z=" <> zStr <> "|value=" <> vstr]
  ]

(* ── Header ─────────────────────────────────────────────────────────────── *)

Print["# Heun oracle outputs — PadeTaylor.jl probe"];
Print["# Mathematica version : " <> $Version];
Print["# Date                : " <> DateString["ISODate"]];
Print["# Working precision   : 50 digits (Mathematica arbitrary precision)"];
Print["# Output precision    : 30 significant digits"];
Print["# HeunG signature     : HeunG[a, q, alpha, beta, gamma, delta, z]"];
Print["#   Fuchsian constraint: epsilon = alpha+beta+1-gamma-delta (internal)"];
Print["#   Normalisation      : HeunG[a,q,al,be,ga,de, 0] = 1"];
Print["# HeunC signature     : HeunC[q, alpha, gamma, delta, epsilon, z]"];
Print["#   Normalisation      : HeunC[q,al,ga,de,ep, 0] = 1"];
Print["# Record format       : Func|key=val|...|z=val|value=re±imi"];
Print["#"];

(* Helper: stringify a rational as a decimal for the record's z field. *)
fmtZ[zv_] := formatReal[zv]

(* ── Regime A: easy reference points ─────────────────────────────────────── *)
Print["#"];
Print["# === Regime A: easy HeunG, HeunC across z grid ==="];

aG = SetPrecision[2, 50]; qG = SetPrecision[1, 50];
alG = SetPrecision[1, 50]; beG = SetPrecision[2, 50];
gaG = SetPrecision[3, 50]; deG = SetPrecision[4, 50];
psG = "func=HeunG|a=2|q=1|alpha=1|beta=2|gamma=3|delta=4";

zsA = {1/10, 1/2, 9/10, 3/2, 19/10, 5/2, 3};
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunG[aG, qG, alG, beG, gaG, deG, zPrec];
    emit["HeunG", psG, fmtZ[zv], val]
  ],
  {zv, zsA}
];

qC = SetPrecision[1, 50]; alC = SetPrecision[1, 50];
gaC = SetPrecision[2, 50]; deC = SetPrecision[3, 50];
epC = SetPrecision[-1, 50];
psC = "func=HeunC|q=1|alpha=1|gamma=2|delta=3|epsilon=-1";

Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunC[qC, alC, gaC, deC, epC, zPrec];
    emit["HeunC", psC, fmtZ[zv], val]
  ],
  {zv, zsA}
];

(* ── Regime B: HeunG special-case identity checks ────────────────────────── *)
Print["#"];
Print["# === Regime B: special-case identities (DLMF Ch 31 has no num tables) ==="];

(* B1: q=alpha=0 → HeunG ≡ 1 along positive real axis on |z|<1 sheet. *)
psB1 = "func=HeunG|a=2|q=0|alpha=0|beta=2|gamma=3|delta=4";
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunG[SetPrecision[2,50], SetPrecision[0,50], SetPrecision[0,50],
                SetPrecision[2,50], SetPrecision[3,50], SetPrecision[4,50], zPrec];
    emit["HeunG", psB1, fmtZ[zv], val]
  ],
  {zv, {1/4, 1/2, 3/4}}
];

(* B2: epsilon=0 from (alpha=1,beta=0,gamma=1,delta=1): epsilon = 1+0+1-1-1 = 0. *)
psB2 = "func=HeunG|a=2|q=1|alpha=1|beta=0|gamma=1|delta=1";
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunG[SetPrecision[2,50], SetPrecision[1,50], SetPrecision[1,50],
                SetPrecision[0,50], SetPrecision[1,50], SetPrecision[1,50], zPrec];
    emit["HeunG", psB2, fmtZ[zv], val]
  ],
  {zv, {1/4, 1/2, 3/4}}
];

(* B3: balanced parameters, a=4. epsilon = 1+2+1-2-2 = 0. *)
psB3 = "func=HeunG|a=4|q=2|alpha=1|beta=2|gamma=2|delta=2";
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunG[SetPrecision[4,50], SetPrecision[2,50], SetPrecision[1,50],
                SetPrecision[2,50], SetPrecision[2,50], SetPrecision[2,50], zPrec];
    emit["HeunG", psB3, fmtZ[zv], val]
  ],
  {zv, {1/4, 1/2, 3/4}}
];

(* ── Regime C: hard parameter regimes ──────────────────────────────────── *)
Print["#"];
Print["# === Regime C: hard cases (a near 1, large |z|, large q) ==="];

(* C1: a = 1.001 (singularities at z=1 and z=a nearly coincident) *)
psC1 = "func=HeunG|a=1.001|q=0.5|alpha=1|beta=2|gamma=3|delta=4";
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunG[SetPrecision[1001/1000,50], SetPrecision[1/2,50], SetPrecision[1,50],
                SetPrecision[2,50], SetPrecision[3,50], SetPrecision[4,50], zPrec];
    emit["HeunG", psC1, fmtZ[zv], val]
  ],
  {zv, {1/4, 1/2, 3/4}}
];

(* C2: q = 3/2 (heuristic "near apparent singularity") *)
psC2 = "func=HeunG|a=2|q=1.5|alpha=1|beta=2|gamma=3|delta=4";
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunG[SetPrecision[2,50], SetPrecision[3/2,50], SetPrecision[1,50],
                SetPrecision[2,50], SetPrecision[3,50], SetPrecision[4,50], zPrec];
    emit["HeunG", psC2, fmtZ[zv], val]
  ],
  {zv, {1/4, 1/2, 3/4}}
];

(* C3: large |z| from Regime A baseline *)
psC3 = "func=HeunG|a=2|q=1|alpha=1|beta=2|gamma=3|delta=4";
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunG[aG, qG, alG, beG, gaG, deG, zPrec];
    emit["HeunG", psC3, fmtZ[zv], val]
  ],
  {zv, {5, 10}}
];

(* C4: HeunC with large |q|=20 *)
psC4 = "func=HeunC|q=20|alpha=1|gamma=2|delta=3|epsilon=-1";
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunC[SetPrecision[20,50], SetPrecision[1,50], SetPrecision[2,50],
                SetPrecision[3,50], SetPrecision[-1,50], zPrec];
    emit["HeunC", psC4, fmtZ[zv], val]
  ],
  {zv, {1/4, 1/2, 3/4}}
];

(* ── Regime D: Teukolsky-Schwarzschild scalar l=2 n=0 QNM HeunC ─────────── *)
(* Parameters from docs/teukolsky_heun_mapping.md (Fiziev 2009 Eq II.6, 2M=1 units).
   For l=2, s=0, n=0: omega ≈ (0.967284, -0.193532) in 2M=1 units (Leaver 1985 Table I).
   alpha = +2 I omega
   gamma = +2 I omega
   delta = +2 omega^2
   eta   = -l(l+1) + s^2 + 2 omega^2     -- Mathematica's HeunC uses q (DLMF q ≡ Fiziev's η? confirm convention)
   We expose these as "raw Fiziev" parameters; Mathematica's HeunC[q,alpha,gamma,delta,epsilon,z]
   needs the convention bridge — for now we emit values at the Fiziev parameters and let
   the Julia side translate. *)
Print["#"];
Print["# === Regime D: Teukolsky Schwarzschild scalar QNM l=2 n=0 (Fiziev params) ==="];
Print["# Convention: Fiziev 2009 Eq II.6, 2M=1 units; omega from Leaver 1985 Table I."];

omegaQNM = SetPrecision[967284/1000000, 50] - I*SetPrecision[193532/1000000, 50];
alphaT = 2*I*omegaQNM;
gammaT = 2*I*omegaQNM;
deltaT = 2*omegaQNM^2;
etaT   = SetPrecision[-6, 50] + 2*omegaQNM^2;   (* -l(l+1) + s^2 + 2 omega^2, l=2,s=0 *)

(* Mathematica HeunC convention: HeunC[q, alpha, gamma, delta, epsilon, z].
   Fiziev's eta corresponds to Mathematica's q (per DLMF 31.12).
   Fiziev's delta corresponds to Mathematica's epsilon (per same).
   For now we sample with the variable substitution z_M ≡ r (physical radial coord). *)
psD = "func=HeunC|fiziev_q=eta|fiziev_alpha=2Iom|fiziev_gamma=2Iom|fiziev_delta=2om2|omega=0.967284-0.193532i|l=2|s=0|n=0";
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunC[etaT, alphaT, gammaT, deltaT, etaT, zPrec];   (* placeholder: needs translation per worklog 051 *)
    emit["HeunC", psD, fmtZ[zv], val]
  ],
  {zv, {3/2, 2, 5/2}}    (* sample at r=1.5, 2, 2.5 in 2M=1 units; r=1 is horizon *)
];

(* ── Regime E: off-real complex z (branch-structure probe) ──────────────── *)
Print["#"];
Print["# === Regime E: complex z (branch-structure probe) ==="];

psEG = "func=HeunG|a=2|q=1|alpha=1|beta=2|gamma=3|delta=4";
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunG[aG, qG, alG, beG, gaG, deG, zPrec];
    emit["HeunG", psEG, formatComplex[N[zv,30]], val]
  ],
  {zv, {1/2 + I/2, 3/2 + I/4, 5/2 - I/2}}
];

psEC = "func=HeunC|q=1|alpha=1|gamma=2|delta=3|epsilon=-1";
Do[
  Module[{zPrec = SetPrecision[zv, 50], val},
    val = HeunC[qC, alC, gaC, deC, epC, zPrec];
    emit["HeunC", psEC, formatComplex[N[zv,30]], val]
  ],
  {zv, {1/2 + I/2, 3/2 + I/4}}
];

Print["#"];
Print["# === capture complete ==="];
