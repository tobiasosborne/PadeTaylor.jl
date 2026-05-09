(* capture.wl -- Phase 6 Problems oracle (Mathematica primary).
   Pin closed-form WeierstrassP values for two distinct purposes:

   (A) The pole-bridge demo points z ∈ {0.5, 0.95, 1.05, 1.4}, used by
       the revised Phase-6 v1 acceptance (worklog 004).  These pin the
       Padé-extends-past-Taylor-radius demonstration.  z=0.5 and z=0.95
       are BEFORE the lattice pole at z=1; z=1.05 and z=1.4 are AFTER.

   (B) The deferred FW 2011 Table 5.1 reference values at z ∈ {30, 10⁴,
       28.261, 7.123} (lines 301, 372 of the markdown).  Kept intact for
       the v2 long-range integration work.

   Usage:  wolframscript -file capture.wl
   Outputs: oracles_wolfram.txt (Julia-readable literals).

   Variable-name discipline: NO underscores in Mathematica-side names
   (Mathematica parses `u_30` as Pattern[u, Blank[]]; friction recorded
   in worklog 003 §F2).  Only Julia output keys are underscored.

   The c1 = -1, c2 = g3 = 2 lattice setup is FW 2011 line 292.
*)

ClearAll["Global`*"];

cOne = -1;
cTwo = 2;
uIC  = SetPrecision[1071822516416917/10^15,  50];
upIC = SetPrecision[1710337353176786/10^15,  50];

(* WeierstrassP returns Complex even at real z; cast to Re because
   the lattice we use (g_2 = 0, g_3 = 2) makes ℘ real on the real axis
   between poles. *)
weierU[zarg_]  := Re[WeierstrassP[zarg + cOne, {0, cTwo}]];
weierUp[zarg_] := Re[WeierstrassPPrime[zarg + cOne, {0, cTwo}]];

strm = OpenWrite["oracles_wolfram.txt"];
WriteString[strm,
  "# Phase 6 Problems oracle (wolframscript primary).\n",
  "# Captured by: wolframscript -file external/probes/problems-oracle/capture.wl\n",
  "# Format: Julia-readable literals.\n\n",
  "# (A) Pole-bridge demo (worklog 004; new Phase-6 v1 acceptance):\n",
  "#       z ∈ {0.5, 0.95, 1.05, 1.4}; pole at z=1.\n",
  "# (B) FW 2011 Table 5.1 reference values (deferred to v2):\n",
  "#       u(30)     = 1.095098255959744                  (line 301)\n",
  "#       u(10^4)   = 21.02530339471055                   (line 301)\n",
  "#       u(28.261) = 9.876953517025014e6                 (line 372)\n\n"
];

(* ===================================================================
   (A) Pole-bridge demo points -- new Phase-6 v1 acceptance
   =================================================================== *)

(* ---- z = 0.5 (BEFORE pole): closed-form + NDSolve cross-check ---- *)
uHalf  = N[weierU[1/2],  17];
upHalf = N[weierUp[1/2], 17];
solHalf = First[NDSolve[
  {y''[t] == 6 y[t]^2, y[0] == uIC, y'[0] == upIC},
  y, {t, 0, 1/2},
  WorkingPrecision -> 50
]];
ifnHalf       = y /. solHalf;
uHalfFromNd   = N[ifnHalf[1/2],                  17];
upHalfFromNd  = N[Derivative[1][ifnHalf][1/2],   17];

WriteString[strm,
  "# z = 0.5 (BEFORE pole at z=1) -- closed-form + NDSolve cross-check.\n",
  "u_at_0_5             = " <> ToString[CForm[uHalf]]        <> "\n",
  "up_at_0_5            = " <> ToString[CForm[upHalf]]       <> "\n",
  "u_at_0_5_via_NDSolve  = " <> ToString[CForm[uHalfFromNd]]  <> "\n",
  "up_at_0_5_via_NDSolve = " <> ToString[CForm[upHalfFromNd]] <> "\n\n"
];

(* ---- z = 0.95 (NEAR pole, BEFORE): closed-form + NDSolve cross-check ---- *)
uNineFive  = N[weierU[95/100],  17];
upNineFive = N[weierUp[95/100], 17];
solNineFive = First[NDSolve[
  {y''[t] == 6 y[t]^2, y[0] == uIC, y'[0] == upIC},
  y, {t, 0, 95/100},
  WorkingPrecision -> 50
]];
ifnNineFive       = y /. solNineFive;
uNineFiveFromNd   = N[ifnNineFive[95/100],                  17];
upNineFiveFromNd  = N[Derivative[1][ifnNineFive][95/100],   17];

WriteString[strm,
  "# z = 0.95 (NEAR pole, just BEFORE z=1) -- closed-form + NDSolve.\n",
  "u_at_0_95             = " <> ToString[CForm[uNineFive]]        <> "\n",
  "up_at_0_95            = " <> ToString[CForm[upNineFive]]       <> "\n",
  "u_at_0_95_via_NDSolve  = " <> ToString[CForm[uNineFiveFromNd]]  <> "\n",
  "up_at_0_95_via_NDSolve = " <> ToString[CForm[upNineFiveFromNd]] <> "\n\n"
];

(* ---- z = 1.05 (PAST pole): closed-form only as primary, plus 80-dps ---- *)
uOnePointZeroFive  = N[weierU[105/100],  17];
upOnePointZeroFive = N[weierUp[105/100], 17];
WriteString[strm,
  "# z = 1.05 (PAST pole at z=1) -- closed-form sole primary source;\n",
  "# mpmath.odefun cannot integrate through the pole.  Padé bridges it.\n",
  "u_at_1_05  = " <> ToString[CForm[uOnePointZeroFive]]  <> "\n",
  "up_at_1_05 = " <> ToString[CForm[upOnePointZeroFive]] <> "\n"
];

(* 80-dps strings at z=1.05 for the BigFloat-256 round-trip test. *)
uOnePointZeroFiveHp  = N[weierU[105/100],  80];
upOnePointZeroFiveHp = N[weierUp[105/100], 80];
WriteString[strm,
  "# z = 1.05 at 80-digit precision (for BigFloat-256 round-trip).\n",
  "u_at_1_05_80dps_str  = \"" <> ToString[NumberForm[uOnePointZeroFiveHp,  {78, 77}, ExponentFunction -> (Null &)]] <> "\"\n",
  "up_at_1_05_80dps_str = \"" <> ToString[NumberForm[upOnePointZeroFiveHp, {78, 77}, ExponentFunction -> (Null &)]] <> "\"\n\n"
];

(* ---- z = 1.4 (FURTHER PAST pole): closed-form + NDSolve restart from 1.05 ---- *)
uOnePointFour  = N[weierU[14/10],  17];
upOnePointFour = N[weierUp[14/10], 17];
WriteString[strm,
  "# z = 1.4 (further past pole) -- closed-form sole primary source.\n",
  "u_at_1_4  = " <> ToString[CForm[uOnePointFour]]  <> "\n",
  "up_at_1_4 = " <> ToString[CForm[upOnePointFour]] <> "\n"
];

(* NDSolve restart attempt: start at z=1.05 from closed-form IC, push to 1.4. *)
uIcRestart  = SetPrecision[uOnePointZeroFiveHp,  50];
upIcRestart = SetPrecision[upOnePointZeroFiveHp, 50];
restartResult = Quiet[Check[
  Module[{sol, ifn},
    sol = First[NDSolve[
      {y''[t] == 6 y[t]^2, y[105/100] == uIcRestart, y'[105/100] == upIcRestart},
      y, {t, 105/100, 14/10},
      WorkingPrecision -> 50
    ]];
    ifn = y /. sol;
    {N[ifn[14/10], 17], N[Derivative[1][ifn][14/10], 17]}
  ],
  $Failed
]];

If[restartResult === $Failed || Head[restartResult] =!= List,
  WriteString[strm,
    "# NDSolve restart at z=1.05 did not converge; closed-form is sole source for z > 1.\n\n"
  ],
  Module[{uRest, upRest},
    {uRest, upRest} = restartResult;
    WriteString[strm,
      "# NDSolve restart from closed-form IC at z=1.05, integrated to z=1.4.\n",
      "u_at_1_4_via_NDSolve_restart  = " <> ToString[CForm[uRest]]  <> "\n",
      "up_at_1_4_via_NDSolve_restart = " <> ToString[CForm[upRest]] <> "\n\n"
    ]
  ]
];

(* ===================================================================
   (B) FW 2011 Table 5.1 reference values -- deferred to v2 (kept intact)
   =================================================================== *)

(* ---- Test 6.1.1: u(30) ---- *)
u30  = N[weierU[30],  17];
up30 = N[weierUp[30], 17];
WriteString[strm,
  "# Test 6.1.1 -- u(30) and u'(30); FW reference 1.095098255959744.\n",
  "u_at_30        = " <> ToString[CForm[u30]]  <> "\n",
  "up_at_30       = " <> ToString[CForm[up30]] <> "\n",
  "u_at_30_FW_ref = 1.095098255959744\n\n"
];

(* ---- Test 6.1.2: u(10^4) ---- *)
u10000  = N[weierU[10000],  17];
up10000 = N[weierUp[10000], 17];
WriteString[strm,
  "# Test 6.1.2 -- u(10^4) and u'(10^4); FW reference 21.02530339471055.\n",
  "u_at_10000        = " <> ToString[CForm[u10000]]  <> "\n",
  "up_at_10000       = " <> ToString[CForm[up10000]] <> "\n",
  "u_at_10000_FW_ref = 21.02530339471055\n\n"
];

(* ---- Test 6.1.3: u(28.261) -- pole wall ---- *)
u28261  = N[weierU[28261/1000],  17];
up28261 = N[weierUp[28261/1000], 17];
WriteString[strm,
  "# Test 6.1.3 -- u(28.261) (high on pole wall); FW ref 9.876953517025014e6.\n",
  "u_at_28261        = " <> ToString[CForm[u28261]]  <> "\n",
  "up_at_28261       = " <> ToString[CForm[up28261]] <> "\n",
  "u_at_28261_FW_ref = 9.876953517025014e6\n\n"
];

(* ---- Test 6.1.4: u(30) at high precision (for BigFloat-256 test) ---- *)
WriteString[strm,
  "# Test 6.1.4 -- u(30) at 80-digit precision (for BigFloat-256 BigFloat 1e-65 tol).\n",
  "# Source: N[WeierstrassP[30 + c1, {0, c2}], 80].\n"
];
u30hp  = N[weierU[30],  80];
up30hp = N[weierUp[30], 80];
WriteString[strm,
  "u_at_30_80dps_str  = \"" <> ToString[NumberForm[u30hp,  {78, 77}, ExponentFunction -> (Null &)]] <> "\"\n",
  "up_at_30_80dps_str = \"" <> ToString[NumberForm[up30hp, {78, 77}, ExponentFunction -> (Null &)]] <> "\"\n\n"
];

(* ---- Test 6.1.5: dense interpolation at z = 7.123 ---- *)
u7123  = N[weierU[7123/1000],  17];
up7123 = N[weierUp[7123/1000], 17];
WriteString[strm,
  "# Test 6.1.5 -- dense interp target u(7.123) for problem 6.1.1.\n",
  "u_at_7_123  = " <> ToString[CForm[u7123]]  <> "\n",
  "up_at_7_123 = " <> ToString[CForm[up7123]] <> "\n\n"
];

(* ---- FW 2011 ICs (re-pinned for test convenience) ---- *)
WriteString[strm,
  "# FW 2011 ICs (FW2011_*.md:292-295):\n",
  "u_0_FW   = 1.071822516416917\n",
  "up_0_FW  = 1.710337353176786\n",
  "c_1_FW   = -1\n",
  "c_2_FW   = 2\n"
];

Close[strm];
Print["Wrote oracles_wolfram.txt"];
