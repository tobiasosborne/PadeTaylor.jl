(* capture.wl -- Phase 5 PadeStepper oracle (Mathematica primary).
   Pin closed-form Weierstrass-P values and high-precision NDSolve
   values for Painleve I, both starting from the FW 2011 line-292/295
   pinned ICs.

   Usage:  wolframscript -file capture.wl
   Outputs: oracles_wolfram.txt (Julia-readable literals).

   Variable-name discipline: NO underscores in Mathematica-side names
   (Mathematica parses `u_05` as Pattern[u, Blank[]]).  Output keys are
   underscored Julia-friendly identifiers.
*)

ClearAll["Global`*"];

cOne = -1;
cTwo = 2;
uIC  = SetPrecision[1071822516416917/10^15,  50];
upIC = SetPrecision[1710337353176786/10^15,  50];

strm = OpenWrite["oracles_wolfram.txt"];

WriteString[strm,
  "# Phase 5 PadeStepper oracle (wolframscript primary).\n",
  "# Captured by: wolframscript -file external/probes/padestepper-oracle/capture.wl\n",
  "# Format: Julia-readable literals.\n\n",
  "# FW 2011 ICs (verbatim string-match from FW2011_*.md:292-295):\n",
  "u_0_FW   = " <> ToString[CForm[N[uIC,  17]]] <> "\n",
  "up_0_FW  = " <> ToString[CForm[N[upIC, 17]]] <> "\n",
  "c_1_FW   = -1\n",
  "c_2_FW   = 2\n\n"
];

(* Closed-form: weier[zarg] = WeierstrassP(zarg + cOne; 0, cTwo). *)
weierU[zarg_]  := WeierstrassP[zarg + cOne, {0, cTwo}];
weierUp[zarg_] := WeierstrassPPrime[zarg + cOne, {0, cTwo}];

(* Sanity at z=0. *)
sU0  = N[weierU[0]  - uIC,  17];
sUp0 = N[weierUp[0] - upIC, 17];
WriteString[strm,
  "# Sanity check (should both be < 1e-15):\n",
  "u0_residual_at_z0  = " <> ToString[CForm[sU0]]  <> "\n",
  "up0_residual_at_z0 = " <> ToString[CForm[sUp0]] <> "\n\n"
];

(* ---- 5.1.1: weier at z=0.5 (closed-form + NDSolve cross-check) ---- *)
uHalf  = N[weierU[1/2],  17];
upHalf = N[weierUp[1/2], 17];

(* Independent check via NDSolve on u'' = 6 u^2 (same equianharmonic ODE). *)
weierSolNd = First[NDSolve[
  {y''[t] == 6 y[t]^2, y[0] == uIC, y'[0] == upIC},
  y, {t, 0, 1/2},
  WorkingPrecision -> 50
]];
ifnW            = y /. weierSolNd;
uHalfFromNd     = N[ifnW[1/2],                   17];
upHalfFromNd    = N[Derivative[1][ifnW][1/2],    17];
weierMethodDiff = N[Abs[uHalf - uHalfFromNd],    17];

WriteString[strm,
  "# Case 5.1.1 -- one step h=0.5 on u''=6u^2 from FW ICs.\n",
  "# Source A (primary): closed-form WeierstrassP[z + c1, {0, c2}].\n",
  "# Source B (cross-check): NDSolve at WorkingPrecision=50.\n",
  "u_5_1_1_at_05            = " <> ToString[CForm[uHalf]]            <> "\n",
  "up_5_1_1_at_05           = " <> ToString[CForm[upHalf]]           <> "\n",
  "u_5_1_1_at_05_via_NDSolve  = " <> ToString[CForm[uHalfFromNd]]    <> "\n",
  "up_5_1_1_at_05_via_NDSolve = " <> ToString[CForm[upHalfFromNd]]   <> "\n",
  "u_5_1_1_method_diff       = " <> ToString[CForm[weierMethodDiff]] <> "\n\n"
];

(* ---- Sanity at z=1: pole? ---- *)
invUatOne = N[1/weierU[1], 17];
WriteString[strm,
  "# Sanity at z=1: 1/u(1) (small ~ on a pole, since c1=-1 and z+c1=0).\n",
  "inv_u_at_z1 = " <> ToString[CForm[invUatOne]] <> "\n\n"
];

(* ---- 5.1.2: continuation to z=0.9 (pole at z=1 avoided) ---- *)
uNine  = N[weierU[9/10],  17];
upNine = N[weierUp[9/10], 17];
WriteString[strm,
  "# Case 5.1.2 -- continuation step h=0.4 from z=0.5 to z=0.9.\n",
  "z_target_5_1_2 = 0.9\n",
  "u_5_1_2_at_09  = " <> ToString[CForm[uNine]]  <> "\n",
  "up_5_1_2_at_09 = " <> ToString[CForm[upNine]] <> "\n\n"
];

(* ---- 5.1.3: PI step from FW ICs to z=0.5 via NDSolve ---- *)
piSol = First[NDSolve[
  {y''[t] == 6 y[t]^2 + t, y[0] == uIC, y'[0] == upIC},
  y, {t, 0, 1/2},
  WorkingPrecision -> 50
]];
uPiHalf  = N[(y[1/2]   /. piSol),                17];
upPiHalf = N[(y'[1/2]  /. piSol /. (y -> y)),    17];
(* The above ReplaceAll trick: y'[1/2]/.piSol applies the
   InterpolatingFunction, then we do the no-op replace to materialise. *)
(* Cleaner alternative: extract the InterpolatingFunction directly. *)
ifn = y /. piSol;
uPiHalfClean  = N[ifn[1/2],            17];
upPiHalfClean = N[Derivative[1][ifn][1/2], 17];

WriteString[strm,
  "# Case 5.1.3 -- PI step (u''=6u^2+z) from FW ICs to z=0.5; NDSolve WorkingPrecision=50.\n",
  "u_5_1_3_PI_at_05  = " <> ToString[CForm[uPiHalfClean]]  <> "\n",
  "up_5_1_3_PI_at_05 = " <> ToString[CForm[upPiHalfClean]] <> "\n",
  "# Diff from Weier at z=0.5 (the +z RHS contribution):\n",
  "diff_PI_minus_weier_at_05  = " <> ToString[CForm[N[uPiHalfClean  - uHalf,  17]]] <> "\n",
  "diff_PI_minus_weierP_at_05 = " <> ToString[CForm[N[upPiHalfClean - upHalf, 17]]] <> "\n\n"
];

(* ---- 5.1.4: near-pole step from z=0.9 (h=0.05 -> z=0.95). ---- *)
WriteString[strm,
  "# Case 5.1.4 -- near-pole step from z=0.9 to z=0.95 (pole at z=1, h=0.05).\n",
  "z_step_start_5_1_4 = 0.9\n",
  "h_step_5_1_4       = 0.05\n",
  "z_target_5_1_4     = 0.95\n",
  "u_state_at_09      = " <> ToString[CForm[uNine]]  <> "\n",
  "up_state_at_09     = " <> ToString[CForm[upNine]] <> "\n"
];

uNineFive  = N[weierU[95/100],  17];
upNineFive = N[weierUp[95/100], 17];
WriteString[strm,
  "u_5_1_4_at_095     = " <> ToString[CForm[uNineFive]]  <> "\n",
  "up_5_1_4_at_095    = " <> ToString[CForm[upNineFive]] <> "\n",
  "u_5_1_4_magnitude  = " <> ToString[CForm[N[Abs[uNineFive], 17]]] <> "\n"
];

Close[strm];
Print["Wrote oracles_wolfram.txt"];
