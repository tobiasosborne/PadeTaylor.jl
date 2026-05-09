(* capture.wl -- Phase 3 Coefficients oracle (Mathematica primary).
   Pin Taylor coefficient values for the 4 numerical test cases.

   Usage:  wolframscript -file capture.wl
   Outputs: oracles.txt (Julia-readable literals; copied to
            test/_oracle_coefficients.jl).

   ASCII-only output to avoid wolframscript UTF-8 glitches.

   Cases pinned:
     3.1.1  exp(z) Taylor coeffs, order 10 -- Float64.
            Source: Series[Exp[z], {z, 0, 10}]; matches the well-known
            c_k = 1/k! (any standard analysis text, e.g. Rudin 8.6).

     3.1.2  Taylor coeffs of t*J_{3/4}(t^2/2)/J_{-1/4}(t^2/2), order 14.
            FW 2011 sec. 2.2.1 closed-form for dy/dz = z^2 + y^2, y(0) = 0.
            Source citation: references/markdown/FW2011_*/FW2011_*.md
            describes this closed form in section 2.2.

     3.1.3  Taylor coeffs of u(z) where u'' = 6 u^2 with FW 2011-pinned
            equianharmonic-Weierstrass ICs, order 30 -- Float64.
            ICs from FW 2011 lines 292-295 (verbatim string-match):
              u(0)  ~= 1.071822516416917
              u'(0) ~= 1.710337353176786
              c_1 = -1, c_2 = g_3 = 2
            Source: AsymptoticDSolveValue on the ODE with these ICs;
            cross-checked against Series[WeierstrassP[z + c_1, {0, c_2}]].

     3.1.4  exp(z) Taylor coeffs, order 60 -- scientific-notation strings
            for parse(BigFloat, _) at precision 256. 1/60! ~= 1.2e-82, so
            decimal-only formatting truncates; we emit mantissa+exponent.
            Source: N[1/k!, 80] for k = 0..60.

   Precision: 35 digits for Float64-targeted cases (well above Float64
   resolution of ~16 digits), 80 digits for the BigFloat-256 case
   (above the ~77-digit BigFloat-256 resolution).
*)

ClearAll["Global`*"];

(* --------- Output helpers --------- *)

(* Format a Float64 value as a Julia-readable literal with 17 significant
   digits (Julia's BigFloat-clean Float64 round-trip). Always include a
   decimal so Julia parses as Float64 not Int. *)
juliaFloat[x_] := Module[{val, s},
  val = N[x, 17];
  s = ToString[CForm[N[val * 1.0, 17]]];
  (* Ensure a decimal point or exponent appears so Julia infers Float64. *)
  If[StringFreeQ[s, "." | "e" | "E"], s = s <> ".0"];
  s
];

(* Format a high-precision real as a quoted scientific-notation string for
   parse(BigFloat, _, base=10). Uses CForm to force mantissa-exponent form
   so values < 10^-78 (e.g. 1/60!) round-trip correctly. *)
juliaBigFloatString[x_, digits_] := Module[{n, s},
  n = N[x, digits];
  s = If[n == 0,
    "0.0",
    ToString[CForm[n]]
  ];
  "\"" <> s <> "\""
];

writeFloat64Vec[strm_, name_, v_] := Module[{},
  WriteString[strm, name <> " = ["];
  Do[
    If[k > 1, WriteString[strm, ", "]];
    WriteString[strm, juliaFloat[v[[k]]]],
    {k, Length[v]}];
  WriteString[strm, "]\n"]
];

writeBigFloatStrVec[strm_, name_, v_, digits_] := Module[{},
  WriteString[strm, name <> " = ["];
  Do[
    If[k > 1, WriteString[strm, ",\n  "]];
    WriteString[strm, juliaBigFloatString[v[[k]], digits]],
    {k, Length[v]}];
  WriteString[strm, "]\n"]
];

(* --------- Open output --------- *)

strm = OpenWrite["oracles.txt"];
WriteString[strm,
  "# Phase 3 Coefficients oracle --wolframscript output.\n",
  "# Captured by: wolframscript -file external/probes/coefficients-oracle/capture.wl\n",
  "# Format: Julia-readable literals (top-level assignments).\n",
  "# Cross-validated against capture.py (Python sympy/mpmath).\n\n"
];

(* --------- 3.1.1 --exp(z) Taylor, order 10, Float64 --------- *)

WriteString[strm, "# Case 3.1.1 --exp(z) Taylor coeffs, order 10 (Float64).\n"];
WriteString[strm, "# Source: Series[Exp[z], {z, 0, 10}]; matches c_k = 1/k!.\n"];
expCoefs = CoefficientList[Normal[Series[Exp[z], {z, 0, 10}]], z];
expCoefs = PadRight[expCoefs, 11];
writeFloat64Vec[strm, "c_3_1_1_exp_o10", expCoefs];
WriteString[strm, "\n"];

(* --------- 3.1.2 --Bessel ratio Taylor, order 14, Float64 --------- *)

WriteString[strm, "# Case 3.1.2 --t*J_{3/4}(t^2/2)/J_{-1/4}(t^2/2), order 14.\n"];
WriteString[strm, "# Closed-form solution to dy/dz = z^2 + y^2, y(0)=0 per FW 2011 sec. 2.2.1.\n"];
besselExpr = t BesselJ[3/4, t^2/2] / BesselJ[-1/4, t^2/2];
(* The ratio is analytic at t=0 (singularities cancel); leading coeff is t^3/3. *)
besselSeries = Series[besselExpr, {t, 0, 14}];
besselCoefs = CoefficientList[Normal[besselSeries], t];
besselCoefs = PadRight[besselCoefs, 15];
writeFloat64Vec[strm, "c_3_1_2_bessel_o14", besselCoefs];
WriteString[strm, "\n"];

(* --------- 3.1.3 --u'' = 6u² with FW 2011 pinned ICs, order 30 --------- *)

WriteString[strm, "# Case 3.1.3 --u(z) with u''=6u^2, u(0)=1.071822516416917, u'(0)=1.710337353176786.\n"];
WriteString[strm, "# ICs string-matched from FW 2011 markdown lines 292-295.\n"];
WriteString[strm, "# Source: AsymptoticDSolveValue (Mathematica's ODE-Taylor solver).\n"];
u0Val = SetPrecision[1.071822516416917, 30];
up0Val = SetPrecision[1.710337353176786, 30];
weierSol = AsymptoticDSolveValue[
  {u''[z] == 6 u[z]^2, u[0] == u0Val, u'[0] == up0Val},
  u[z], {z, 0, 30}];
weierCoefs = CoefficientList[weierSol, z];
weierCoefs = PadRight[weierCoefs, 31];
writeFloat64Vec[strm, "c_3_1_3_weierstrass_o30", weierCoefs];
WriteString[strm, "\n"];

(* Cross-check: g_3 = 4 u_0^3 - u'_0^2 should be ≈ 2 (FW 2011 line 292). *)
g3val = 4 u0Val^3 - up0Val^2;
WriteString[strm, "# Cross-check: g_3 = 4*u_0^3 - u'_0^2 (FW 2011 line 292: c_2 = g_3 = 2)\n"];
WriteString[strm, "c_3_1_3_g3_check = " <> juliaFloat[g3val] <> "\n"];

(* Cross-check 2: independent path via WeierstrassP[z + c_1, {0, c_2}]. *)
weierSeries2 = Series[WeierstrassP[z + (-1), {0, 2}], {z, 0, 30}];
weierCoefs2 = CoefficientList[Normal[weierSeries2], z];
weierCoefs2 = PadRight[weierCoefs2, 31];
WriteString[strm, "# Cross-check: same coeffs via Series[WeierstrassP[z-1, {0,2}], ...].\n"];
writeFloat64Vec[strm, "c_3_1_3_weierstrass_o30_via_P", weierCoefs2];
WriteString[strm, "\n"];

(* --------- 3.1.4 --1/k! at 80 digits for BigFloat-256 test, order 60 --------- *)

WriteString[strm, "# Case 3.1.4 --1/k! at 80 decimal digits, k=0..60.\n"];
WriteString[strm, "# For BigFloat(precision=256) exp(z) Taylor coefficient test.\n"];
expBigCoefs = Table[1/Factorial[k], {k, 0, 60}];
writeBigFloatStrVec[strm, "c_3_1_4_exp_o60_str", expBigCoefs, 78];
WriteString[strm, "\n"];

Close[strm];
Print["Wrote oracles.txt"];
