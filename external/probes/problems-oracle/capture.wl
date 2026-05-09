(* capture.wl -- Phase 6 Problems oracle (Mathematica primary).
   Pin closed-form WeierstrassP values matching FW 2011 Table 5.1
   reference values (lines 301, 372 of the markdown).

   Usage:  wolframscript -file capture.wl
   Outputs: oracles_wolfram.txt (Julia-readable literals).

   Variable-name discipline: NO underscores in Mathematica-side names
   (Mathematica parses `u_30` as Pattern[u, Blank[]]).

   The c1 = -1, c2 = g3 = 2 lattice setup is FW 2011 line 292.  We pin
   four values at full Float64 precision plus one at ≈80 decimal digits
   for the BigFloat-256 test case.
*)

ClearAll["Global`*"];

cOne = -1;
cTwo = 2;

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
  "# FW 2011 Table 5.1 cross-validation targets, c1 = -1, c2 = g3 = 2.\n",
  "# Reference values pinned in FW 2011:\n",
  "#   u(30)     = 1.095098255959744                  (line 301)\n",
  "#   u(10^4)   = 21.02530339471055                   (line 301)\n",
  "#   u(28.261) = 9.876953517025014e6                 (line 372)\n\n"
];

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
