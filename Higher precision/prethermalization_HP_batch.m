(* ================================================================ *)
(*  HIGH-PRECISION PRETHERMALIZATION WORKFLOW                         *)
(*                                                                    *)
(*  Single initial condition, arbitrary precision.                    *)
(*  Parallelization: over time-point chunks (not ICs).               *)
(*  Precision report after every numerical step.                     *)
(*                                                                    *)
(*  Modules:                                                          *)
(*    1. Setup + global precision                                     *)
(*    2. Helpers (in .wl)                                             *)
(*    3. Parameters                                                   *)
(*    4. Main loop (L-first, SRN+LRN per L)                          *)
(*    5. Plotting                                                     *)
(* ================================================================ *)


(* ================================================================ *)
(* 1. SETUP                                                           *)
(* ================================================================ *)

(* >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> *)
(*  GLOBAL WORKING PRECISION — set this once at the top.             *)
(*  All downstream computations inherit this value.                  *)
(*  Rule of thumb: to reach time t_max with |E_k|~O(L),             *)
(*  you need wp > Log10(|E_max| * t_max) + desired_digits.           *)
(*  E.g., t_max=10^12, |E|~10, desired=8 digits => wp > 21.         *)
$WorkingPrecision = 30;
(* <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< *)

$baseDir = If[$InputFileName =!= "" && $InputFileName =!= None,
  DirectoryName[$InputFileName],
  NotebookDirectory[]
];
SetDirectory[$baseDir];

Get[FileNameJoin[{$baseDir, "QMB_prethermalization_HP.wl"}]];

$resultsBase = FileNameJoin[{$baseDir, "results_HP"}];
createDirs[regime_, LL_] := Module[{base},
  base = FileNameJoin[{$resultsBase, regime, "L" <> ToString[LL]}];
  If[!DirectoryQ[base],
    CreateDirectory[base, CreateIntermediateDirectories -> True]];
  base
];
$plotDir = FileNameJoin[{$resultsBase, "plots"}];
If[!DirectoryQ[$plotDir],
  CreateDirectory[$plotDir, CreateIntermediateDirectories -> True]];

(* Parallel kernels *)
nKernels = Min[$ProcessorCount - 1, 14];
If[Length[Kernels[]] < nKernels,
  LaunchKernels[nKernels - Length[Kernels[]]]];
ParallelEvaluate[
  Get[FileNameJoin[{$baseDir, "QMB_prethermalization_HP.wl"}]]];

Print["SETUP | wp=", $WorkingPrecision,
  ", Kernels=", Length[Kernels[]],
  ", Dir=", $baseDir];


(* ================================================================ *)
(* 2. HELPERS                                                         *)
(* ================================================================ *)

(* --- Hamiltonian setup with precision tracking --- *)
setupForHG$HP[hval_, gval_, LL_, wp_] :=
  Module[{Ha, Heven, Hodd,
    evenEigvals, evenEigvecs, oddEigvals, oddEigvecs,
    evenMap, oddMap, Me, Mo, p},

    Print["  Building H..."];
    Ha = IsingHamiltonianExact[
      SetPrecision[hval, wp], SetPrecision[gval, wp],
      SetPrecision[1, wp], LL,
      BoundaryConditions -> "Open"];
    ReportPrecision["Hamiltonian (exact)", Ha];

    Print["  Sector projection..."];
    {Heven, Hodd} = SectorHamiltoniansHP[Ha, LL, wp];
    ReportPrecision["Heven", Heven];
    ReportPrecision["Hodd", Hodd];

    Print["  Diagonalizing even sector (dim=",
      Length[Heven], ")..."];
    {evenEigvals, evenEigvecs} = DiagonalizeSectorHP[Heven, wp];
    p = ReportPrecision["Even eigenvalues", evenEigvals];
    ReportPrecision["Even eigenvectors", evenEigvecs];

    Print["  Diagonalizing odd sector (dim=",
      Length[Hodd], ")..."];
    {oddEigvals, oddEigvecs} = DiagonalizeSectorHP[Hodd, wp];
    p = ReportPrecision["Odd eigenvalues", oddEigvals];
    ReportPrecision["Odd eigenvectors", oddEigvecs];

    (* Orthonormality check *)
    Module[{GE, errE},
      GE = evenEigvecs . ConjugateTranspose[evenEigvecs];
      errE = Max[Abs[GE - IdentityMatrix[Length[evenEigvals]]]];
      Print["  Orthonormality error (even): ",
        ScientificForm[errE, 3],
        If[errE < 10^(-wp + 5), " [OK]", " [WARNING]"]];
    ];

    {evenMap, oddMap} = buildMaps[LL];

    <|"eVals" -> evenEigvals, "eVecs" -> evenEigvecs,
      "oVals" -> oddEigvals,  "oVecs" -> oddEigvecs,
      "evenMap" -> evenMap, "oddMap" -> oddMap|>
  ];


(* --- Single-IC entropy computation with parallelized time chunks --- *)
computeEntropyHP[setup_Association, ttlist_, LL_, wp_, seed_: 1234] :=
  Module[{evenMap, oddMap, Me, Mo,
    psi0, psiE, psiO, ckE, ckO,
    dA = 2^(LL/2), dB = 2^(LL/2),
    Nt, cSize, nChunks, chunkRanges, results},

    evenMap = setup["evenMap"]; oddMap = setup["oddMap"];
    Me = makeExpansionMatrixExact[evenMap, LL, "even"];
    Mo = makeExpansionMatrixExact[oddMap, LL, "odd"];

    (* Generate initial state at controlled precision *)
    BlockRandom[SeedRandom[seed]];
    psi0 = Flatten[KroneckerProduct[
      HaarStateHP[dA, wp], HaarStateHP[dB, wp]], 1];
    ReportPrecision["Initial state", psi0];

    (* Verify S(t=0) = 0 *)
    Module[{s0},
      s0 = EntanglementEntropySVDHP[psi0, dA, dB, wp];
      Print["  S(t=0) = ", ScientificForm[s0, 3],
        If[Abs[s0] < 10^(-wp + 5), " [OK]", " [WARNING]"]];
    ];

    (* Project into sectors *)
    psiE = ConjugateTranspose[Me] . psi0;
    psiO = ConjugateTranspose[Mo] . psi0;
    ReportPrecision["Projected even", psiE];
    ReportPrecision["Projected odd", psiO];

    (* Norm check *)
    Module[{nrm},
      nrm = Re[Conjugate[psiE] . psiE + Conjugate[psiO] . psiO];
      Print["  Norm (even+odd) = ", NumberForm[nrm, {wp, wp - 2}],
        If[Abs[nrm - 1] < 10^(-wp + 3), " [OK]", " [WARNING]"]];
    ];

    (* Overlap coefficients (computed once) *)
    ckE = ComputeOverlapsHP[setup["eVecs"], psiE];
    ckO = ComputeOverlapsHP[setup["oVecs"], psiO];
    ReportPrecision["ckE", ckE];
    ReportPrecision["ckO", ckO];

    (* --- Parallel time evolution --- *)
    (* Chunk the time list for parallel dispatch *)
    Nt = Length[ttlist];
    cSize = Max[1, Ceiling[Nt / Length[Kernels[]]]];
    chunkRanges = Partition[Range[Nt], UpTo[cSize]];
    nChunks = Length[chunkRanges];
    Print["  Evolving ", Nt, " time points in ", nChunks,
      " chunks (cSize=", cSize, ")..."];

    (* Inject all data via With for zero-overhead dispatch *)
    results = With[
      {wEVecs = setup["eVecs"], wEVals = setup["eVals"],
       wOVecs = setup["oVecs"], wOVals = setup["oVals"],
       wCkE = ckE, wCkO = ckO,
       wMe = Me, wMo = Mo,
       wTlist = ttlist, wDA = dA, wDB = dB, wWP = wp,
       wChunks = chunkRanges},

      ParallelTable[
        Module[{indices, entropies},
          indices = wChunks[[ci]];
          entropies = Table[
            Module[{t, evolE, evolO, fullState},
              t = wTlist[[idx]];
              evolE = EvolveSectorHP[wEVecs, wEVals, wCkE, t];
              evolO = EvolveSectorHP[wOVecs, wOVals, wCkO, t];
              fullState = wMe . evolE + wMo . evolO;
              EntanglementEntropySVDHP[fullState, wDA, wDB, wWP]
            ],
            {idx, indices}
          ];
          entropies
        ],
        {ci, nChunks}
      ]
    ];

    (* Flatten results back to ordered list *)
    Flatten[results]
  ];


(* --- Smoothing --- *)
smoothEntropy[s_, w_] := Module[{n = Length[s]},
  If[w <= 1 || n <= w, Return[s]];
  Table[Mean[s[[Max[1, i - Floor[(w-1)/2]] ;;
    Min[n, i + Floor[(w-1)/2]]]]], {i, n}]
];

(* --- Color function --- *)
cf = Blend[{{0, Blue}, {0.5, Purple}, {1, Red}}, #] &;
cfBar[logMin_, logMax_] :=
  Function[x, cf[Rescale[x, {logMin, logMax}, {0, 1}]]];

Print["HELPERS | Defined."];


(* ================================================================ *)
(* 3. PARAMETERS                                                      *)
(* ================================================================ *)

wp = $WorkingPrecision;

Llist = {8, 10};
Nic   = 1;  (* single IC *)
seed  = 1234;

(* Time grid at working precision *)
tlist = Table[SetPrecision[10^\[Tau], wp],
  {\[Tau], -1, 12, 1/200}];
Nt = Length[tlist];

(* SRN: fix g, scan h *)
gSRNlist = SetPrecision[{1}, wp];
hSRN = Table[SetPrecision[10^x, wp], {x, -6, -1, 1/4}];

(* LRN: fix h, scan g *)
hLRNlist = SetPrecision[{1}, wp];
gLRN = Table[SetPrecision[10^x, wp], {x, -6, -1, 1/4}];

maWindow = 100;

Print["PARAMETERS | wp=", wp, ", Llist=", Llist, ", Nt=", Nt];
Print["  SRN: ", Length[gSRNlist], " g x ", Length[hSRN], " h"];
Print["  LRN: ", Length[hLRNlist], " h x ", Length[gLRN], " g"];


(* ================================================================ *)
(* 4. MAIN LOOP                                                       *)
(* ================================================================ *)

Do[
  Print["\n###################################################"];
  Print["###  L = ", LL, "  (wp = ", wp, ")  ###"];
  Print["###################################################"];

  Spage = PageEntropy[LL/2, LL/2];
  $SRNspage[LL] = Spage;
  $LRNspage[LL] = Spage;

  (* --- 4a: SRN for this L --- *)
  Do[
    Print["\n========== SRN | L=", LL, ", g=", gFixed, " =========="];
    dataDir = createDirs["SRN", LL];
    resultsSRN = Association[];

    Do[
      hval = hSRN[[ih]];
      Print["\n  h=", hval, " (", ih, "/", Length[hSRN], ")"];
      t0 = AbsoluteTime[];

      Print["  --- Hamiltonian & diagonalization ---"];
      setup = setupForHG$HP[hval, gFixed, LL, wp];

      Print["  --- Time evolution ---"];
      entropyList = computeEntropyHP[setup, tlist, LL, wp, seed];
      ReportPrecision["Entropy list", entropyList];

      (* Spot-check: precision of last entropy value *)
      Print["  S(t_max) = ", entropyList[[-1]],
        ", prec = ", Precision[entropyList[[-1]]]];

      resultsSRN[hval] = <|"meanEntropy" -> N[entropyList]|>;

      Export[FileNameJoin[{dataDir,
        "meanEntropy_h" <> ToString[AccountingForm[N[hval], 8]] <>
        ".m"}], N[entropyList]];

      Print["  Total: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
      , {ih, Length[hSRN]}
    ];

    Export[FileNameJoin[{dataDir,
      "allResults_g" <> ToString[AccountingForm[N[gFixed], 8]] <>
      ".m"}], resultsSRN];
    $SRNdata[{LL, gFixed}] = resultsSRN;

    , {gFixed, gSRNlist}
  ];
  Print["\n--- SRN complete for L=", LL, " ---"];

  (* --- 4b: LRN for this L --- *)
  Do[
    Print["\n========== LRN | L=", LL, ", h=", hFixed, " =========="];
    dataDir = createDirs["LRN", LL];
    resultsLRN = Association[];

    Do[
      gval = gLRN[[ig]];
      Print["\n  g=", gval, " (", ig, "/", Length[gLRN], ")"];
      t0 = AbsoluteTime[];

      Print["  --- Hamiltonian & diagonalization ---"];
      setup = setupForHG$HP[hFixed, gval, LL, wp];

      Print["  --- Time evolution ---"];
      entropyList = computeEntropyHP[setup, tlist, LL, wp, seed];
      ReportPrecision["Entropy list", entropyList];

      resultsLRN[gval] = <|"meanEntropy" -> N[entropyList]|>;

      Export[FileNameJoin[{dataDir,
        "meanEntropy_g" <> ToString[AccountingForm[N[gval], 8]] <>
        ".m"}], N[entropyList]];

      Print["  Total: ", Round[AbsoluteTime[] - t0, 0.1], " s"];
      , {ig, Length[gLRN]}
    ];

    Export[FileNameJoin[{dataDir,
      "allResults_h" <> ToString[AccountingForm[N[hFixed], 8]] <>
      ".m"}], resultsLRN];
    $LRNdata[{LL, hFixed}] = resultsLRN;

    , {hFixed, hLRNlist}
  ];
  Print["\n--- LRN complete for L=", LL, " ---"];
  Print["=== L=", LL, " FULLY DONE ==="];

  , {LL, Llist}
];


(* ================================================================ *)
(* 5. PLOTTING                                                        *)
(* ================================================================ *)

(* --- 5a: SRN entropy plots --- *)
Do[
  Module[{res, paramList, nH, hMin, hMax, logMin, logMax,
    curves, colors, Sp, plotTitle, plt},

    res = $SRNdata[{LL, gFixed}];
    If[!AssociationQ[res], Continue[]];

    paramList = Sort[Keys[res]];
    nH = Length[paramList];
    If[nH == 0, Continue[]];

    hMin = Min[N[paramList]]; hMax = Max[N[paramList]];
    logMin = Log10[hMin]; logMax = Log10[hMax];
    Sp = $SRNspage[LL];

    curves = Table[
      Transpose[{N[tlist],
        smoothEntropy[
          N[res[paramList[[i]]]["meanEntropy"]], maWindow]}],
      {i, nH}];

    colors = Table[
      cf[Rescale[Log10[N[paramList[[i]]]], {logMin, logMax}, {0, 1}]],
      {i, nH}];

    plotTitle = "HP (wp=" <> ToString[wp] <> "), L=" <>
      ToString[LL] <> ", g=" <>
      ToString[NumberForm[N[gFixed], {6, 5}]] <> " (SRN)";

    plt = Legended[
      ListLogLinearPlot[curves,
        PlotRange -> {{N[tlist[[1]]], N[tlist[[-1]]]},
          {0, Sp + 0.3}},
        PlotStyle -> (Directive[#, AbsoluteThickness[1.2]] & /@
          colors),
        Frame -> True,
        FrameLabel -> {Style["Time (t)", 14],
          Style["Entropy S(t)", 14]},
        FrameTicksStyle -> 12,
        PlotLabel -> Style[plotTitle, 12, Bold],
        ImageSize -> 700, AspectRatio -> 0.55,
        Epilog -> {
          {Black, Dashed, AbsoluteThickness[2],
            InfiniteLine[{0, Sp}, {1, 0}]},
          {Darker[Green], Dashed, AbsoluteThickness[2],
            InfiniteLine[{0, Log[2]}, {1, 0}]}}
      ],
      BarLegend[{cfBar[logMin, logMax], {logMin, logMax}},
        LegendLabel -> Style["Field (h)", 12],
        LabelStyle -> 11,
        "Ticks" -> Table[
          {v, Superscript["10", ToString[Round[v]]]},
          {v, Ceiling[logMin], Floor[logMax]}]]
    ];

    Export[FileNameJoin[{$plotDir,
      "HP_Entropy_L" <> ToString[LL] <>
      "_g" <> ToString[NumberForm[N[gFixed], {6, 5}]] <>
      "_SRN.png"}], plt, ImageResolution -> 300];
  ];
  , {LL, Llist}, {gFixed, gSRNlist}
];

(* --- 5b: LRN entropy plots --- *)
Do[
  Module[{res, paramList, nG, gMin, gMax, logMin, logMax,
    curves, colors, Sp, plotTitle, plt},

    res = $LRNdata[{LL, hFixed}];
    If[!AssociationQ[res], Continue[]];

    paramList = Sort[Keys[res]];
    nG = Length[paramList];
    If[nG == 0, Continue[]];

    gMin = Min[N[paramList]]; gMax = Max[N[paramList]];
    logMin = Log10[gMin]; logMax = Log10[gMax];
    Sp = $LRNspage[LL];

    curves = Table[
      Transpose[{N[tlist],
        smoothEntropy[
          N[res[paramList[[i]]]["meanEntropy"]], maWindow]}],
      {i, nG}];

    colors = Table[
      cf[Rescale[Log10[N[paramList[[i]]]], {logMin, logMax}, {0, 1}]],
      {i, nG}];

    plotTitle = "HP (wp=" <> ToString[wp] <> "), L=" <>
      ToString[LL] <> ", h=" <>
      ToString[NumberForm[N[hFixed], {6, 5}]] <> " (LRN)";

    plt = Legended[
      ListLogLinearPlot[curves,
        PlotRange -> {{N[tlist[[1]]], N[tlist[[-1]]]},
          {0, Sp + 0.3}},
        PlotStyle -> (Directive[#, AbsoluteThickness[1.2]] & /@
          colors),
        Frame -> True,
        FrameLabel -> {Style["Time (t)", 14],
          Style["Entropy S(t)", 14]},
        FrameTicksStyle -> 12,
        PlotLabel -> Style[plotTitle, 12, Bold],
        ImageSize -> 700, AspectRatio -> 0.55,
        Epilog -> {
          {Black, Dashed, AbsoluteThickness[2],
            InfiniteLine[{0, Sp}, {1, 0}]}}
      ],
      BarLegend[{cfBar[logMin, logMax], {logMin, logMax}},
        LegendLabel -> Style["Field (g)", 12],
        LabelStyle -> 11,
        "Ticks" -> Table[
          {v, Superscript["10", ToString[Round[v]]]},
          {v, Ceiling[logMin], Floor[logMax]}]]
    ];

    Export[FileNameJoin[{$plotDir,
      "HP_Entropy_L" <> ToString[LL] <>
      "_h" <> ToString[NumberForm[N[hFixed], {6, 5}]] <>
      "_LRN.png"}], plt, ImageResolution -> 300];
  ];
  , {LL, Llist}, {hFixed, hLRNlist}
];

Print["\n=== HP WORKFLOW COMPLETE ==="];
Print["Results: ", $resultsBase];
Print["Plots:   ", $plotDir];
