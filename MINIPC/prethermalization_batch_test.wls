#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ================================================================ *)
(*  PRETHERMALIZATION: PRODUCTION WORKFLOW                            *)
(*                                                                    *)
(*  Structure:                                                        *)
(*    1. Setup                                                        *)
(*    2. Helpers                                                      *)
(*    3. Parameters                                                   *)
(*    4. Main loop (SRN + LRN)                                        *)
(*    5. Plotting (entropy evolution)                                  *)
(*    6. Melting time (extraction + power-law plots)                   *)
(* ================================================================ *)


(* ================================================================ *)
(* 1. SETUP                                                           *)
(* ================================================================ *)

(* Machine-independent base directory:                               *)
(*   - In wolframscript: uses $InputFileName                         *)
(*   - In notebook: falls back to NotebookDirectory[]                *)
$baseDir = If[$InputFileName =!= "" && $InputFileName =!= None,
  DirectoryName[$InputFileName],
  NotebookDirectory[]
];
SetDirectory[$baseDir];

Get[FileNameJoin[{$baseDir, "QMB_prethermalization.wl"}]];

(* --- Create directory tree --- *)
$resultsBase = FileNameJoin[{$baseDir, "results"}];
createDirs[regime_, LL_] := Module[{base},
  base = FileNameJoin[{$resultsBase, regime, "L" <> ToString[LL]}];
  If[!DirectoryQ[base], CreateDirectory[base, CreateIntermediateDirectories -> True]];
  base
];
$plotDir = FileNameJoin[{$resultsBase, "plots"}];
If[!DirectoryQ[$plotDir], CreateDirectory[$plotDir, CreateIntermediateDirectories -> True]];

Print["SETUP | Kernels: ", Length[Kernels[]], "  | Dir: ", $baseDir];


(* --- Parallel kernels --- *)
LaunchKernels[14];
ParallelEvaluate[Get[FileNameJoin[{$baseDir, "QMB_prethermalization.wl"}]]];


(* ================================================================ *)
(* 2. HELPERS                                                         *)
(* ================================================================ *)

(* --- Hamiltonian setup for given (h, g, L) --- *)
setupForHG[hval_, gval_, LL_] := Module[
  {Ha, Heven, Hodd,
   evenEigvals, evenEigvecs, oddEigvals, oddEigvecs,
   evenMap, oddMap, Me, Mo,
   specFull, specBlock, specErr,
   JJ = 1., ddim = 2^LL, ePrec = 25},

  Ha = N[Normal[IsingHamiltonian[hval, gval, JJ, LL,
    BoundaryConditions -> "Open"]]];

  {evenMap, oddMap} = buildMaps[LL];
  Me = makeExpansionMatrix[evenMap, LL, "even"];
  Mo = makeExpansionMatrix[oddMap, LL, "odd"];

  {Heven, Hodd} = SectorHamiltonians[Ha, LL];
  {evenEigvals, evenEigvecs} = DiagonalizeSector[Heven];
  {oddEigvals, oddEigvecs}   = DiagonalizeSector[Hodd];

  specFull  = Sort[Eigenvalues[Ha]];
  specBlock = Sort[Join[evenEigvals, oddEigvals]];
  specErr   = Max[Abs[specFull - specBlock]];
  Print["  Spectrum error: ", ScientificForm[specErr, 3],
    If[specErr < 1.*^-9, " [OK]", " [WARNING]"]];

  <|"eVals" -> evenEigvals, "eVecs" -> evenEigvecs,
    "oVals" -> oddEigvals,  "oVecs" -> oddEigvecs,
    "eValsHP" -> SetPrecision[evenEigvals, ePrec],
    "oValsHP" -> SetPrecision[oddEigvals, ePrec],
    "Me" -> Me, "Mo" -> Mo|>
];

(* --- Batch entropy over Nic initial conditions --- *)
batchEntropy[setup_Association, NNic_, ttlist_, LL_] :=
  Module[
    {eVals, eVecs, oVals, oVecs, eValsHP, oValsHP,
     mMe, mMo, mMeT, mMoT, eVecsT, oVecsT,
     ddA = 2^(LL/2), ddB = 2^(LL/2), NNt = Length[ttlist],
     cSize = 1000, tThresh = 1.*^5, ePrec = 25},

    eVals = setup["eVals"]; eVecs = setup["eVecs"];
    oVals = setup["oVals"]; oVecs = setup["oVecs"];
    eValsHP = setup["eValsHP"]; oValsHP = setup["oValsHP"];
    mMe = setup["Me"]; mMo = setup["Mo"];
    mMeT = Transpose[mMe]; mMoT = Transpose[mMo];
    eVecsT = Transpose[eVecs]; oVecsT = Transpose[oVecs];

    (* Use With to inject values directly into ParallelTable.       *)
    (* This avoids DistributeDefinitions overhead entirely \[LongDash]         *)
    (* With substitutes at parse time before dispatch to kernels.    *)

    With[
      {wEVals = eVals, wEVecs = eVecs, wOVals = oVals, wOVecs = oVecs,
       wEValsHP = eValsHP, wOValsHP = oValsHP,
       wMe = mMe, wMo = mMo, wMeT = mMeT, wMoT = mMoT,
       wEVecsT = eVecsT, wOVecsT = oVecsT,
       wTlist = ttlist, wCSize = cSize, wTThresh = tThresh,
       wEPrec = ePrec, wLL = LL, wDA = ddA, wDB = ddB, wNt = NNt},

      ParallelTable[
        Module[
          {psi0, psiE, psiO, ckE, ckO, entropies,
           iEnd, tChunk, nChunk, phaseE, phaseO, tHP,
           evolE, evolO, fullChunk},

          psi0 = N[RandomChainProductState[wLL]];
          psiE = wMeT . psi0;
          psiO = wMoT . psi0;
          ckE = Conjugate[wEVecs] . psiE;
          ckO = Conjugate[wOVecs] . psiO;
          entropies = ConstantArray[0., wNt];

          Do[
            iEnd   = Min[iStart + wCSize - 1, wNt];
            tChunk = wTlist[[iStart ;; iEnd]];
            nChunk = iEnd - iStart + 1;
            If[Max[Abs[tChunk]] < wTThresh,
              phaseE = ckE * Exp[-I Outer[Times, wEVals, tChunk]];
              phaseO = ckO * Exp[-I Outer[Times, wOVals, tChunk]];,
              tHP = SetPrecision[tChunk, wEPrec];
              phaseE = N[ckE * Exp[-I Outer[Times, wEValsHP, tHP]]];
              phaseO = N[ckO * Exp[-I Outer[Times, wOValsHP, tHP]]];
            ];
            evolE = wEVecsT . phaseE;
            evolO = wOVecsT . phaseO;
            fullChunk = Transpose[wMe . evolE + wMo . evolO];
            Do[entropies[[iStart + j - 1]] =
                EntanglementEntropySVD[fullChunk[[j]], wDA, wDB];,
              {j, nChunk}];,
            {iStart, 1, wNt, wCSize}];
          entropies],
        {ic, NNic},
        Method -> "CoarsestGrained"]
    ]  (* end With *)
  ];

(* --- Moving average (preserves length) --- *)
smoothEntropy[s_, w_] := Module[{n = Length[s]},
  If[w <= 1 || n <= w, Return[s]];
  Table[Mean[s[[Max[1, i - Floor[(w-1)/2]] ;; Min[n, i + Floor[(w-1)/2]]]]], {i, n}]
];

(* --- Melting time extraction --- *)
extractMeltingTime[meanS_, tgrid_, Sthreshold_] :=
  Module[{idx},
    idx = SelectFirst[Range[Length[tgrid]], meanS[[#]] > Sthreshold &];
    If[MissingQ[idx], Missing["NotReached"], tgrid[[idx]]]
  ];

(* --- Color function: blue -> purple -> red --- *)
cf = Blend[{{0, Blue}, {0.5, Purple}, {1, Red}}, #] &;

(* --- Rescaled color function for BarLegend --- *)
(* BarLegend passes raw tick values (logMin..logMax) to the function.
   We must rescale to [0,1] before calling cf. *)
cfBar[logMin_, logMax_] :=
  Function[x, cf[Rescale[x, {logMin, logMax}, {0, 1}]]];

(* --- Power-law fit on a log-log subrange --- *)
fitPowerLaw[data_, {pMin_, pMax_}] :=
  Module[{subset, logData, lm},
    subset = Select[data, pMin <= #[[1]] <= pMax && NumericQ[#[[2]]] &];
    If[Length[subset] < 3,
      Return[<|"slope" -> Missing[], "intercept" -> Missing[],
               "model" -> None, "R2" -> Missing[]|>]];
    logData = {Log10[#[[1]]], Log10[#[[2]]]} & /@ subset;
    lm = LinearModelFit[logData, x, x];
    <|"slope" -> lm["BestFitParameters"][[2]],
      "intercept" -> lm["BestFitParameters"][[1]],
      "model" -> lm,
      "R2" -> lm["RSquared"]|>
  ];

(* --- Find crossover in log-log data --- *)
findCrossover[data_] :=
  Module[{logData, n, bestIdx = 0, bestResid = Infinity, resid, lm1, lm2},
    logData = {Log10[#[[1]]], Log10[#[[2]]]} & /@ data;
    n = Length[logData];
    Do[
      If[i < 4 || n - i < 4, Continue[]];
      lm1 = LinearModelFit[logData[[1 ;; i]], x, x];
      lm2 = LinearModelFit[logData[[i ;; n]], x, x];
      resid = Total[lm1["FitResiduals"]^2] + Total[lm2["FitResiduals"]^2];
      If[resid < bestResid, bestResid = resid; bestIdx = i];,
      {i, 4, n - 3}];
    If[bestIdx > 0, 10^logData[[bestIdx, 1]], Missing["NoCrossover"]]
  ];

Print["HELPERS | All functions defined."];


(* ================================================================ *)
(* 3. PARAMETERS                                                      *)
(* ================================================================ *)

Llist = {4,6,8,10,12,14};
J     = 1.;
Nic   = 400;

(* --- Time grid: 10^{-1} to 10^{12} --- *)
tlist = Developer`ToPackedArray[N[Table[10^i, {i, -0.5, 10, 0.0005}]]];
Nt = Length[tlist];

(* --- SRN regime: fix g, scan h --- *)
(* gSRN will be varied over gSRNlist; for each gSRN, scan hSRN *)
gSRNlist = N[Table[10^i, {i, -2, 0.5, 0.1}]];
gSRNlist=AppendTo[gSRNlist,2.];
hSRN     = Sort[DeleteDuplicates[N[Table[10^i, {i, -6, -1, 0.05}]]]];

(* --- LRN regime: fix h, scan g --- *)
(* hLRN will be varied over hLRNlist; for each hLRN, scan gLRN *)
hLRNlist = N[Table[10^i, {i, -4, 0.5, 0.1}]];
gLRN     = Sort[DeleteDuplicates[N[Table[10^i, {i, -6, -1, 0.05}]]]];

(* --- Smoothing window --- *)
maWindow = 200;

Print["PARAMETERS | Llist=", Llist, ", Nt=", Nt];
Print["  SRN: ", Length[gSRNlist], " g-values x ", Length[hSRN], " h-values"];
Print["  LRN: ", Length[hLRNlist], " h-values x ", Length[gLRN], " g-values"];


(* ================================================================ *)
(* 4. MAIN LOOP                                                       *)
(*   For each system size L, run ALL SRN then ALL LRN before          *)
(*   moving to the next L.                                            *)
(* ================================================================ *)

Do[
  Print["\n\n###################################################"];
  Print["###  L = ", LL, "  ###"];
  Print["###################################################"];

  Spage = PageEntropy[LL/2, LL/2];
  $SRNspage[LL] = Spage;
  $LRNspage[LL] = Spage;

  (* --- 4a: SRN for this L --- *)

  Do[
    Print["\n========== SRN | L=", LL,
      ", g=", gFixed, " =========="];

    dataDir = createDirs["SRN", LL];
    resultsSRN = Association[];

    Do[
      hval = hSRN[[ih]];
      Print["  h=", hval, " (", ih, "/", Length[hSRN], ")"];
      t0 = AbsoluteTime[];

      setup = setupForHG[hval, gFixed, LL];
      entropyMatrix = batchEntropy[setup, Nic, tlist, LL];
      meanS = Mean[entropyMatrix];
      Clear[entropyMatrix];  (* free RAM immediately *)

      resultsSRN[hval] = <|"meanEntropy" -> meanS|>;

      Export[FileNameJoin[{dataDir,
        "meanEntropy_h" <> ToString[AccountingForm[hval, 8]] <> ".m"}],
        meanS];

      Print["    ", Round[AbsoluteTime[] - t0, 0.1], " s"];
      , {ih, Length[hSRN]}
    ];

    Export[FileNameJoin[{dataDir,
      "allResults_g" <> ToString[AccountingForm[gFixed, 8]] <> ".m"}],
      resultsSRN];

    $SRNdata[{LL, gFixed}] = resultsSRN;

    , {gFixed, gSRNlist}
  ];

  Print["\n--- SRN complete for L=", LL, " ---"];

  (* --- 4b: LRN for this L --- *)

  Do[
    Print["\n========== LRN | L=", LL,
      ", h=", hFixed, " =========="];

    dataDir = createDirs["LRN", LL];
    resultsLRN = Association[];

    Do[
      gval = gLRN[[ig]];
      Print["  g=", gval, " (", ig, "/", Length[gLRN], ")"];
      t0 = AbsoluteTime[];

      setup = setupForHG[hFixed, gval, LL];
      entropyMatrix = batchEntropy[setup, Nic, tlist, LL];
      meanS = Mean[entropyMatrix];
      Clear[entropyMatrix];  (* free RAM immediately *)

      resultsLRN[gval] = <|"meanEntropy" -> meanS|>;

      Export[FileNameJoin[{dataDir,
        "meanEntropy_g" <> ToString[AccountingForm[gval, 8]] <> ".m"}],
        meanS];

      Print["    ", Round[AbsoluteTime[] - t0, 0.1], " s"];
      , {ig, Length[gLRN]}
    ];

    Export[FileNameJoin[{dataDir,
      "allResults_h" <> ToString[AccountingForm[hFixed, 8]] <> ".m"}],
      resultsLRN];

    $LRNdata[{LL, hFixed}] = resultsLRN;

    , {hFixed, hLRNlist}
  ];

  Print["\n--- LRN complete for L=", LL, " ---"];
  Print["=== L=", LL, " FULLY DONE ==="];

  , {LL, Llist}
];


(* ================================================================ *)
(* 5. PLOTTING (entropy evolution)                                    *)
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

    hMin = Min[paramList]; hMax = Max[paramList];
    logMin = Log10[hMin]; logMax = Log10[hMax];
    Sp = $SRNspage[LL];

    curves = Table[
      Transpose[{tlist,
        smoothEntropy[res[paramList[[i]]]["meanEntropy"], maWindow]}],
      {i, nH}];

    colors = Table[
      cf[Rescale[Log10[paramList[[i]]], {logMin, logMax}, {0, 1}]],
      {i, nH}];

    plotTitle = "L=" <> ToString[LL] <>
      ", g=" <> ToString[NumberForm[gFixed, {6, 5}]] <> " (SRN)";

    plt = Legended[
      ListLogLinearPlot[curves,
        PlotRange -> {{tlist[[1]], tlist[[-1]]}, {0, Sp + 0.3}},
        PlotStyle -> (Directive[#, AbsoluteThickness[1]] & /@ colors),
        Frame -> True,
        FrameLabel -> {Style["Time (t)", 14],
          Style["Entropy S(t)", 14]},
        FrameTicksStyle -> 12,
        PlotLabel -> Style[plotTitle, 14, Bold],
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
      "Entropy_L" <> ToString[LL] <>
      "_g" <> ToString[NumberForm[gFixed, {6, 5}]] <> "_SRN.png"}],
      plt, ImageResolution -> 300];
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

    gMin = Min[paramList]; gMax = Max[paramList];
    logMin = Log10[gMin]; logMax = Log10[gMax];
    Sp = $LRNspage[LL];

    curves = Table[
      Transpose[{tlist,
        smoothEntropy[res[paramList[[i]]]["meanEntropy"], maWindow]}],
      {i, nG}];

    colors = Table[
      cf[Rescale[Log10[paramList[[i]]], {logMin, logMax}, {0, 1}]],
      {i, nG}];

    plotTitle = "L=" <> ToString[LL] <>
      ", h=" <> ToString[NumberForm[hFixed, {6, 5}]] <> " (LRN)";

    plt = Legended[
      ListLogLinearPlot[curves,
        PlotRange -> {{tlist[[1]], tlist[[-1]]}, {0, Sp + 0.3}},
        PlotStyle -> (Directive[#, AbsoluteThickness[1]] & /@ colors),
        Frame -> True,
        FrameLabel -> {Style["Time (t)", 14],
          Style["Entropy S(t)", 14]},
        FrameTicksStyle -> 12,
        PlotLabel -> Style[plotTitle, 14, Bold],
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
      "Entropy_L" <> ToString[LL] <>
      "_h" <> ToString[NumberForm[hFixed, {6, 5}]] <> "_LRN.png"}],
      plt, ImageResolution -> 300];
  ];
  , {LL, Llist}, {hFixed, hLRNlist}
];

Print["PLOTTING | All entropy evolution plots saved to ", $plotDir];


(* ================================================================ *)
(* 6. MELTING TIME                                                    *)
(* ================================================================ *)

(* --- 6a: SRN melting times (threshold = Log[2] always) --- *)

SthreshSRN = Log[2];

Do[
  Module[{res, meltData, meltClean},
    res = $SRNdata[{LL, gFixed}];
    If[!AssociationQ[res], Continue[]];

    meltData = Table[
      {hSRN[[i]],
       extractMeltingTime[res[hSRN[[i]]]["meanEntropy"], tlist, SthreshSRN]},
      {i, Length[hSRN]}];
    meltClean = Select[meltData, NumericQ[Last[#]] &];

    $SRNmelt[{LL, gFixed}] = meltClean;
    Export[FileNameJoin[{createDirs["SRN", LL],
      "meltingData_g" <> ToString[AccountingForm[gFixed, 8]] <> ".m"}],
      meltClean];
    Print["SRN melt | L=", LL, ", g=", gFixed,
      ": ", Length[meltClean], "/", Length[hSRN], " points"];
  ];
  , {LL, Llist}, {gFixed, gSRNlist}
];


(* --- 6b: LRN melting times --- *)
(* >>> SET SthreshLRN AFTER INSPECTING THE LRN ENTROPY PLOTS <<< *)
(* One value per system size in Llist = {6, 8, 10, 12}.          *)
(* The threshold increases with L because S_page increases.       *)
SthreshLRN = <|6 -> 1.35, 8 -> 1.8, 10 -> 1.3, 12 -> 1.3|>;  (* <-- EDIT THESE *)

Do[
  Module[{res, meltData, meltClean, thresh},
    res = $LRNdata[{LL, hFixed}];
    If[!AssociationQ[res], Continue[]];

    thresh = SthreshLRN[LL];
    meltData = Table[
      {gLRN[[i]],
       extractMeltingTime[res[gLRN[[i]]]["meanEntropy"], tlist, thresh]},
      {i, Length[gLRN]}];
    meltClean = Select[meltData, NumericQ[Last[#]] &];

    $LRNmelt[{LL, hFixed}] = meltClean;
    Export[FileNameJoin[{createDirs["LRN", LL],
      "meltingData_h" <> ToString[AccountingForm[hFixed, 8]] <> ".m"}],
      meltClean];
    Print["LRN melt | L=", LL, ", h=", hFixed,
      ", thresh=", thresh,
      ": ", Length[meltClean], "/", Length[gLRN], " points"];
  ];
  , {LL, Llist}, {hFixed, hLRNlist}
];


(* --- 6c: SRN melting time power-law plots --- *)
(* Red vertical line at gFixed (the g value used for this SRN run) *)

Do[
  Module[{meltClean, crossover, fit1, fit2,
    xMin, xMax, fitLine1, fitLine2, epilogElems,
    insetText, plt},

    meltClean = $SRNmelt[{LL, gFixed}];
    If[!ListQ[meltClean] || Length[meltClean] < 5, Continue[]];

    crossover = findCrossover[meltClean];
    xMin = Min[meltClean[[All, 1]]];
    xMax = Max[meltClean[[All, 1]]];

    If[NumericQ[crossover],
      fit1 = fitPowerLaw[meltClean, {xMin, crossover}];
      fit2 = fitPowerLaw[meltClean, {crossover, xMax}];,
      fit1 = fitPowerLaw[meltClean, {xMin, xMax}];
      fit2 = <|"slope" -> Missing[]|>;
    ];

    (* Red vertical line at gFixed *)
    epilogElems = {
      {Red, AbsoluteThickness[1.5],
       InfiniteLine[{Log[gFixed], 0}, {0, 1}]}};

    fitLine1 = If[NumericQ[fit1["slope"]],
      LogLogPlot[10^(fit1["intercept"] + fit1["slope"] Log10[x]),
        {x, xMin, If[NumericQ[crossover], crossover, xMax]},
        PlotStyle -> {Blue, Dashed, AbsoluteThickness[1.5]}],
      Graphics[]];
    fitLine2 = If[NumericQ[fit2["slope"]],
      LogLogPlot[10^(fit2["intercept"] + fit2["slope"] Log10[x]),
        {x, If[NumericQ[crossover], crossover, xMin], xMax},
        PlotStyle -> {Darker[Green], Dashed, AbsoluteThickness[1.5]}],
      Graphics[]];

    insetText = {};
    If[NumericQ[fit1["slope"]],
      AppendTo[insetText,
        Style[Row[{Subscript["\[Alpha]", "1"], " = ",
          NumberForm[fit1["slope"], {4, 2}],
          " (", Superscript["R", "2"], "=",
          NumberForm[fit1["R2"], {3, 2}], ")"}],
          12, Blue]]];
    If[NumericQ[fit2["slope"]],
      AppendTo[insetText,
        Style[Row[{Subscript["\[Alpha]", "2"], " = ",
          NumberForm[fit2["slope"], {4, 2}],
          " (", Superscript["R", "2"], "=",
          NumberForm[fit2["R2"], {3, 2}], ")"}],
          12, Darker[Green]]]];

    plt = Show[
      ListLogLogPlot[meltClean,
        PlotStyle -> {Black, PointSize[0.015]},
        PlotRange -> All, Frame -> True,
        FrameLabel -> {Style["h", 14, Italic],
          Style["Melting Time (\!\(\*SubscriptBox[\(T\), \(M\)]\))", 14]},
        FrameTicksStyle -> 12,
        PlotLabel -> Style[
          "L = " <> ToString[LL] <> ", J = 1, g = " <>
          ToString[gFixed], 14],
        ImageSize -> 600, AspectRatio -> 0.65,
        Epilog -> epilogElems],
      fitLine1, fitLine2,
      If[Length[insetText] > 0,
        Graphics[Inset[Column[insetText, Spacings -> 0.3],
          Scaled[{0.75, 0.85}]]],
        Graphics[]]
    ];

    Export[FileNameJoin[{$plotDir,
      "meltingtime_L" <> ToString[LL] <>
      "_g" <> ToString[NumberForm[gFixed, {6, 5}]] <> "_SRN.png"}],
      plt, ImageResolution -> 300];
  ];
  , {LL, Llist}, {gFixed, gSRNlist}
];


(* --- 6d: LRN melting time power-law plots --- *)
(* Single power-law fit only (no crossover for LRN).               *)
(* Red vertical line at hFixed (the h value used for this LRN run) *)

Do[
  Module[{meltClean, fit1,
    xMin, xMax, fitLine1, epilogElems,
    insetText, plt},

    meltClean = $LRNmelt[{LL, hFixed}];
    If[!ListQ[meltClean] || Length[meltClean] < 3, Continue[]];

    xMin = Min[meltClean[[All, 1]]];
    xMax = Max[meltClean[[All, 1]]];

    (* Single power-law over full range *)
    fit1 = fitPowerLaw[meltClean, {xMin, xMax}];

    (* Red vertical line at hFixed *)
    epilogElems = {
      {Red, AbsoluteThickness[1.5],
       InfiniteLine[{Log[hFixed], 0}, {0, 1}]}};

    fitLine1 = If[NumericQ[fit1["slope"]],
      LogLogPlot[10^(fit1["intercept"] + fit1["slope"] Log10[x]),
        {x, xMin, xMax},
        PlotStyle -> {Blue, Dashed, AbsoluteThickness[1.5]}],
      Graphics[]];

    insetText = {};
    If[NumericQ[fit1["slope"]],
      AppendTo[insetText,
        Style[Row[{"\[Alpha]", " = ",
          NumberForm[fit1["slope"], {4, 2}],
          " (", Superscript["R", "2"], "=",
          NumberForm[fit1["R2"], {3, 2}], ")"}],
          12, Blue]]];

    plt = Show[
      ListLogLogPlot[meltClean,
        PlotStyle -> {Black, PointSize[0.015]},
        PlotRange -> All, Frame -> True,
        FrameLabel -> {Style["g", 14, Italic],
          Style["Melting Time (\!\(\*SubscriptBox[\(T\), \(M\)]\))", 14]},
        FrameTicksStyle -> 12,
        PlotLabel -> Style[
          "L = " <> ToString[LL] <> ", J = 1, h = " <>
          ToString[hFixed], 14],
        ImageSize -> 600, AspectRatio -> 0.65,
        Epilog -> epilogElems],
      fitLine1,
      If[Length[insetText] > 0,
        Graphics[Inset[Column[insetText, Spacings -> 0.3],
          Scaled[{0.75, 0.85}]]],
        Graphics[]]
    ];

    Export[FileNameJoin[{$plotDir,
      "meltingtime_L" <> ToString[LL] <>
      "_h" <> ToString[NumberForm[hFixed, {6, 5}]] <> "_LRN.png"}],
      plt, ImageResolution -> 300];
  ];
  , {LL, Llist}, {hFixed, hLRNlist}
];


Print["\n=== WORKFLOW COMPLETE ==="];
Print["Results: ", $resultsBase];
Print["Plots:   ", $plotDir];
