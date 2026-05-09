(* capture.wl -- Phase 4 StepControl oracle, wolframscript triangulation.

   Third-source cross-check for the consensus formula
       h = min over k in {p-1, p} of (eps / |x[k]|)^(1/k)
   on the test inputs of cases 4.1.1, 4.1.3, 4.1.4.

   Usage:  wolframscript -file capture.wl
   Outputs: oracles_wolfram.txt (Julia-readable literals).
*)

ClearAll["Global`*"];

strm = OpenWrite["oracles_wolfram.txt"];
WriteString[strm,
  "# Phase 4 StepControl oracle (wolframscript triangulation).\n",
  "# Cross-validates capture.jl (TI.jl) and capture.py (mpmath/sympy)\n\n"
];

(* ---- 4.1.1 -- step_jorba_zou expected ---- *)
WriteString[strm, "# Case 4.1.1 -- (eps * k!)^(1/k) for k in {29, 30}, eps = 1e-12.\n"];
WriteString[strm, "# Source: Mathematica N[..., 50] from the consensus formula.\n"];

eps = 10^-12;
p = 30;
candidate[k_] := N[(eps * Factorial[k])^(1/k), 50];
h29 = candidate[29];
h30 = candidate[30];
hmin = Min[h29, h30];

WriteString[strm, "h_4_1_1_at_29_wl = " <> ToString[CForm[N[h29, 17]]] <> "\n"];
WriteString[strm, "h_4_1_1_at_30_wl = " <> ToString[CForm[N[h30, 17]]] <> "\n"];
WriteString[strm, "h_4_1_1_min_wl   = " <> ToString[CForm[N[hmin, 17]]] <> "\n"];
WriteString[strm, "h_4_1_1_min_wl_50dps = \"" <> ToString[NumberForm[hmin, {48, 47}, ExponentFunction -> (Null &)]] <> "\"\n\n"];

(* ---- 4.1.3 -- roots of 1 - z/2 ---- *)
WriteString[strm, "# Case 4.1.3 -- roots of Q(z) = 1 - z/2.\n"];
roots413 = z /. NSolve[1 - z/2 == 0, z];
WriteString[strm, "# roots: " <> ToString[roots413] <> "\n"];
WriteString[strm, "step_4_1_3_expected_wl = 2.0\n\n"];

(* ---- 4.1.4 -- roots of 1 - 0.6 z + 0.1 z^2 ---- *)
WriteString[strm, "# Case 4.1.4 -- roots of Q(z) = 1 - 0.6 z + 0.1 z^2.\n"];
roots414 = z /. NSolve[1 - 6/10 z + 1/10 z^2 == 0, z];
WriteString[strm, "# roots: " <> ToString[roots414] <> "\n"];
(* Real-axis projection of pole at 3+i: Re((3+i)*Conjugate[1])/|1| = 3 *)
WriteString[strm, "step_4_1_4_expected_wl = 3.0\n"];

Close[strm];
Print["Wrote oracles_wolfram.txt"];
