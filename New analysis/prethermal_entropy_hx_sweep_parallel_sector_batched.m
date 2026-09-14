(* ::Package:: *)

(* ::Title:: *)
(*Setup*)


qmbInitPath = "C:\\Users\\Miguel\\Github\\libs\\QMB\\Kernel\\init.m";
If[!FileExistsQ[qmbInitPath],
    Print["ERROR: QMB init.m was not found: ", qmbInitPath]; Abort[];
];
Get[qmbInitPath];
SetDirectory[NotebookDirectory[]];


(* ::Title::Closed:: *)
(*Definitions and compiled kernels*)


ClearAll[reflectionSectorBases];
reflectionSectorBases[L_Integer] := Module[
    {dim, perm, plusRules = {}, minusRules = {}, np = 0, nm = 0, j, k},
    dim = 2^L;
    perm = Table[
        1 + FromDigits[Reverse[IntegerDigits[j - 1, 2, L]], 2],
        {j, dim}
    ];
    Do[
        k = perm[[j]];
        Which[
            j < k,
            np++;
            AppendTo[plusRules, {j, np} -> 1/Sqrt[2]];
            AppendTo[plusRules, {k, np} -> 1/Sqrt[2]];
            nm++;
            AppendTo[minusRules, {j, nm} -> 1/Sqrt[2]];
            AppendTo[minusRules, {k, nm} -> -1/Sqrt[2]],
            j == k,
            np++;
            AppendTo[plusRules, {j, np} -> 1]
        ],
        {j, dim}
    ];
    <|
        "Even" -> SparseArray[plusRules, {dim, np}],
        "Odd"  -> SparseArray[minusRules, {dim, nm}]
    |>
];


ClearAll[sortedEigensystem];
sortedEigensystem[M_] := Module[{evals, evecs, ord},
    {evals, evecs} = Eigensystem[Normal[M]];
    ord = Ordering[evals];
    {
        Developer`ToPackedArray[N[evals[[ord]]]],
        Developer`ToPackedArray[N[evecs[[ord]]]]
    }
];


(* One uniform-grid phase step. *)
ClearAll[phaseStepC];
phaseStepC = Compile[
    {
        {aR, _Real, 1}, {aI, _Real, 1},
        {cosStep, _Real, 1}, {sinStep, _Real, 1}
    },
    Module[{newR, newI},
        newR = aR*cosStep + aI*sinStep;
        newI = aI*cosStep - aR*sinStep;
        {newR, newI}
    ],
    CompilationTarget -> "C",
    RuntimeOptions -> "Speed"
];


(* Generate a batch of amplitudes on a uniform linear time grid. *)
ClearAll[phaseBatchC];
phaseBatchC = Compile[
    {
        {aR0, _Real, 1}, {aI0, _Real, 1},
        {cosStep, _Real, 1}, {sinStep, _Real, 1},
        {nBatch, _Integer}
    },
    Module[
        {d = Length[aR0], r, im, curR = aR0, curI = aI0,
         nextR, nextI, i, j},
        r = Table[0., {Length[aR0]}, {nBatch}];
        im = Table[0., {Length[aR0]}, {nBatch}];
        For[j = 1, j <= nBatch, j++,
            For[i = 1, i <= d, i++,
                r[[i, j]] = curR[[i]];
                im[[i, j]] = curI[[i]];
            ];
            If[j < nBatch,
                nextR = curR*cosStep + curI*sinStep;
                nextI = curI*cosStep - curR*sinStep;
                curR = nextR;
                curI = nextI;
            ];
        ];
        {r, im}
    ],
    CompilationTarget -> "C",
    RuntimeOptions -> "Speed"
];


(* Independent phases at arbitrary times: used for logarithmic grids. *)
ClearAll[phaseAtTimesBatchC];
phaseAtTimesBatchC = Compile[
    {
        {cR, _Real, 1}, {cI, _Real, 1},
        {energies, _Real, 1}, {times, _Real, 1}
    },
    Module[{d = Length[energies], b = Length[times], r, im, ct, st, i, j},
        r = Table[0., {Length[energies]}, {Length[times]}];
        im = Table[0., {Length[energies]}, {Length[times]}];
        For[j = 1, j <= b, j++,
            For[i = 1, i <= d, i++,
                ct = Cos[energies[[i]]*times[[j]]];
                st = Sin[energies[[i]]*times[[j]]];
                r[[i, j]] = cR[[i]]*ct + cI[[i]]*st;
                im[[i, j]] = cI[[i]]*ct - cR[[i]]*st;
            ];
        ];
        {r, im}
    ],
    CompilationTarget -> "C",
    RuntimeOptions -> "Speed"
];


(* Independent phase at one arbitrary time. *)
ClearAll[reseedC];
reseedC = Compile[
    {
        {cR, _Real, 1}, {cI, _Real, 1},
        {energies, _Real, 1}, {t, _Real}
    },
    Module[{ct, st},
        ct = Cos[energies*t];
        st = Sin[energies*t];
        {cR*ct + cI*st, cI*ct - cR*st}
    ],
    CompilationTarget -> "C",
    RuntimeOptions -> "Speed"
];


ClearAll[entropyFromSVC];
entropyFromSVC = Compile[
    {{sv, _Real, 1}},
    Module[{sum = 0., p = 0., i, n},
        n = Length[sv];
        For[i = 1, i <= n, i++,
            p = sv[[i]]*sv[[i]];
            If[p > 1.*^-14, sum -= p*Log[p]];
        ];
        sum
    ],
    CompilationTarget -> "C",
    RuntimeOptions -> "Speed"
];


ClearAll[entropyFast];
entropyFast[psiR_, psiI_, dA_Integer, dB_Integer] := Module[{mat, sv},
    mat = ArrayReshape[psiR + I psiI, {dA, dB}];
    sv = SingularValueList[mat];
    entropyFromSVC[sv]
];


(* A trajectory must be a finite machine-real vector of the expected length. *)
ClearAll[validEntropyVectorQ];
validEntropyVectorQ[x_, n_Integer] :=
    TrueQ[VectorQ[x, MachineNumberQ[#] && TrueQ[Im[#] == 0.] &]] &&
    Length[x] == n;

(* Dense reconstruction remains inside the smaller reflection sectors.
   The sparse B+ and B- maps lift the states only after the fast
   matrix-matrix products. *)
ClearAll[entropyBatchFromSectors];
entropyBatchFromSectors[
    WP_, WM_, BpN_, BmN_, batchP_, batchM_, dA_Integer, dB_Integer
] := Module[{uPR, uPI, uMR, uMI, psiR, psiI, b, entropies, normErrors},
    uPR = WP . batchP[[1]];
    uPI = WP . batchP[[2]];
    uMR = WM . batchM[[1]];
    uMI = WM . batchM[[2]];
    psiR = BpN . uPR + BmN . uMR;
    psiI = BpN . uPI + BmN . uMI;
    b = Dimensions[psiR][[2]];
    entropies = Developer`ToPackedArray[
        Table[
            entropyFast[psiR[[All, j]], psiI[[All, j]], dA, dB],
            {j, b}
        ]
    ];
    normErrors = Abs[Total[psiR*psiR + psiI*psiI] - 1.];
    <|"entropy" -> entropies, "maxNormError" -> Max[normErrors]|>
];


(* Batched uniform-grid evolution with periodic absolute reseeding. *)
ClearAll[linearEntropyEvolutionSector];
linearEntropyEvolutionSector[
    WP_, WM_, BpN_, BmN_, evalsP_, evalsM_, coeffP_, coeffM_,
    t0_Real, dt_Real, nSteps_Integer, batchSize_Integer,
    resetEvery_Integer, dA_Integer, dB_Integer
] := Module[
    {nSamples = nSteps + 1, entropy, cPR, cPI, cMR, cMI,
     cosP, sinP, cosM, sinM, aPR, aPI, aMR, aMI,
     batchP, batchM, values, pos = 1, step, len,
     stepsUntilReset, currentTime, batchResult, maxNormError = 0.},

    entropy = Developer`ToPackedArray[ConstantArray[0., nSamples]];
    cPR = Developer`ToPackedArray[Re[coeffP]];
    cPI = Developer`ToPackedArray[Im[coeffP]];
    cMR = Developer`ToPackedArray[Re[coeffM]];
    cMI = Developer`ToPackedArray[Im[coeffM]];
    cosP = Developer`ToPackedArray[Cos[evalsP*dt]];
    sinP = Developer`ToPackedArray[Sin[evalsP*dt]];
    cosM = Developer`ToPackedArray[Cos[evalsM*dt]];
    sinM = Developer`ToPackedArray[Sin[evalsM*dt]];

    {aPR, aPI} = reseedC[cPR, cPI, evalsP, t0];
    {aMR, aMI} = reseedC[cMR, cMI, evalsM, t0];

    While[pos <= nSamples,
        step = pos - 1;
        If[step > 0 && resetEvery > 0 && Mod[step, resetEvery] == 0,
            currentTime = t0 + step*dt;
            {aPR, aPI} = reseedC[cPR, cPI, evalsP, N[currentTime]];
            {aMR, aMI} = reseedC[cMR, cMI, evalsM, N[currentTime]];
        ];

        stepsUntilReset = If[
            resetEvery > 0,
            resetEvery - Mod[step, resetEvery],
            nSamples
        ];

        len = Min[batchSize, nSamples - pos + 1, stepsUntilReset];
        batchP = phaseBatchC[aPR, aPI, cosP, sinP, len];
        batchM = phaseBatchC[aMR, aMI, cosM, sinM, len];
        batchResult = entropyBatchFromSectors[
            WP, WM, BpN, BmN, batchP, batchM, dA, dB
        ];
        values = batchResult["entropy"];
        If[!validEntropyVectorQ[values, len],
            Return[Failure["InvalidEntropyBatch", <|"startIndex" -> pos|>]]
        ];
        maxNormError = Max[maxNormError, batchResult["maxNormError"]];
        entropy[[pos ;; pos + len - 1]] = values;

        If[pos + len <= nSamples && Mod[step + len, resetEvery] != 0,
            {aPR, aPI} = phaseStepC[
                batchP[[1, All, -1]], batchP[[2, All, -1]], cosP, sinP
            ];
            {aMR, aMI} = phaseStepC[
                batchM[[1, All, -1]], batchM[[2, All, -1]], cosM, sinM
            ];
        ];
        pos += len;
    ];
    <|"entropy" -> entropy, "maxNormError" -> maxNormError|>
];


(* Batched arbitrary/logarithmic grid: every sample is generated directly
   from the original spectral coefficients, so recurrence drift is absent. *)
ClearAll[arbitraryEntropyEvolutionSector];
arbitraryEntropyEvolutionSector[
    WP_, WM_, BpN_, BmN_, evalsP_, evalsM_, coeffP_, coeffM_,
    times_, batchSize_Integer, dA_Integer, dB_Integer
] := Module[
    {nSamples = Length[times], entropy, cPR, cPI, cMR, cMI,
     pos = 1, len, timesChunk, batchP, batchM, values,
     batchResult, maxNormError = 0.},

    entropy = Developer`ToPackedArray[ConstantArray[0., nSamples]];
    cPR = Developer`ToPackedArray[Re[coeffP]];
    cPI = Developer`ToPackedArray[Im[coeffP]];
    cMR = Developer`ToPackedArray[Re[coeffM]];
    cMI = Developer`ToPackedArray[Im[coeffM]];

    While[pos <= nSamples,
        len = Min[batchSize, nSamples - pos + 1];
        timesChunk = Developer`ToPackedArray[
            N[times[[pos ;; pos + len - 1]]]
        ];
        batchP = phaseAtTimesBatchC[cPR, cPI, evalsP, timesChunk];
        batchM = phaseAtTimesBatchC[cMR, cMI, evalsM, timesChunk];
        batchResult = entropyBatchFromSectors[
            WP, WM, BpN, BmN, batchP, batchM, dA, dB
        ];
        values = batchResult["entropy"];
        If[!validEntropyVectorQ[values, len],
            Return[Failure["InvalidEntropyBatch", <|"startIndex" -> pos|>]]
        ];
        maxNormError = Max[maxNormError, batchResult["maxNormError"]];
        entropy[[pos ;; pos + len - 1]] = values;
        pos += len;
    ];
    <|"entropy" -> entropy, "maxNormError" -> maxNormError|>
];


(* Direct one-time check using the same eigensystem but no recurrence. *)
ClearAll[directEntropyAndNormSector];
directEntropyAndNormSector[
    WP_, WM_, BpN_, BmN_, evalsP_, evalsM_, coeffP_, coeffM_,
    t_Real, dA_Integer, dB_Integer
] := Module[{cPR, cPI, cMR, cMI, aP, aM, psiR, psiI, s, normError},
    cPR = Developer`ToPackedArray[Re[coeffP]];
    cPI = Developer`ToPackedArray[Im[coeffP]];
    cMR = Developer`ToPackedArray[Re[coeffM]];
    cMI = Developer`ToPackedArray[Im[coeffM]];
    aP = reseedC[cPR, cPI, evalsP, t];
    aM = reseedC[cMR, cMI, evalsM, t];
    psiR = BpN . (WP . aP[[1]]) + BmN . (WM . aM[[1]]);
    psiI = BpN . (WP . aP[[2]]) + BmN . (WM . aM[[2]]);
    s = entropyFast[psiR, psiI, dA, dB];
    normError = Abs[Total[psiR*psiR + psiI*psiI] - 1.];
    {s, normError}
];


ClearAll[sampledEigensystemResidual];
sampledEigensystemResidual[M_, evals_, evecs_] := Module[{idx, scale},
    idx = DeleteDuplicates[
        Round[Subdivide[1, Length[evals], Min[4, Max[1, Length[evals] - 1]]]]
    ];
    scale = Max[1., Norm[M, "Frobenius"]];
    Max[
        Table[
            Norm[M . evecs[[j]] - evals[[j]]*evecs[[j]]]/scale,
            {j, idx}
        ]
    ]
];


(* Formatting helpers: force numerical labels and explicit 10^x notation. *)
ClearAll[numericString];
numericString[x_?NumericQ] := ToString[InputForm[N[x]]];

ClearAll[exponentString];
exponentString[x_?NumericQ] := Module[{xr = Round[x]},
    If[
        Abs[x - xr] < 10^-10,
        ToString[xr],
        ToString[NumberForm[N[x], {8, 3}, NumberPadding -> {"", ""}]]
    ]
];

ClearAll[powerOfTenLabel];
powerOfTenLabel[x_?NumericQ, size_: 22] := Style[
    Superscript["10", exponentString[x]], size, Black
];


(* ::Title::Closed:: *)
(*Parameters*)


L = 8;
J = 1.;
hz = 1.;


(* hx = 10^alpha. *)
hxExponents = N[Range[-2, 0, 2/10]];
hxList = Developer`ToPackedArray[N[10^hxExponents]];


(* Time-grid mode:
   "Linear" -> equally spaced in t, arbitrary real dt.
   "Log"    -> equally spaced in log10(t). *)
timeGridType = "Linear";


(* LINEAR GRID PARAMETERS.
   For the long production target use, for example:
       tmax = 10^7;
       dt   = 1.;
   dt need not be an integer. *)
tmin = 0.;
tmax = 10^6;
dt = 1.;


(* LOG GRID PARAMETERS.
   logTMin=-1, logTMax=10, pointsPerDecade=1000 gives
   11001 positive samples from 10^-1 to 10^10. *)
logTMin = -1.;
logTMax = 10.;
pointsPerDecade = 2000;
includeZeroInLogGrid = True;


(* Linear recurrence reseeding interval in number of samples. *)
resetEvery = 100000;


(* Matrix-matrix batch size.  These conservative defaults limit memory
   while giving BLAS substantially more work per call. *)
batchSize = Which[
    L <= 10, 64,
    L <= 12, 32,
    True, 16
];


(* Independent hx values are distributed over eight local kernels. *)
requestedKernels = 8;


(* Accuracy/validation tolerances. *)
hermiticityTolerance = 10^-12;
reflectionTolerance = 10^-12;
reconstructionTolerance = 10^-10;
eigensystemResidualTolerance = 10^-10;
validationEntropyTolerance = 10^-8;
validationNormTolerance = 10^-10;


(* Fixed RPS across the entire hx sweep. *)
seed = 12345;


If[!EvenQ[L] || L <= 0,
    Print["ERROR: L must be a positive even integer for the half-chain cut."];
    Abort[];
];

If[!AllTrue[{requestedKernels, batchSize, resetEvery}, IntegerQ[#] && # > 0 &],
    Print["ERROR: requestedKernels, batchSize and resetEvery must be positive."];
    Abort[];
];

LA = L/2;
dA = 2^LA;
dB = 2^(L - LA);


(* Construct only the time data needed by the selected grid. *)
Switch[
    timeGridType,

    "Linear",
    If[dt <= 0. || tmax < tmin,
        Print["ERROR: Linear grid requires dt > 0 and tmax >= tmin."];
        Abort[];
    ];
    nSteps = Round[(tmax - tmin)/dt];
    If[
        Abs[tmin + nSteps*dt - tmax] > 10^-10 Max[1., Abs[tmax]],
        Print["ERROR: tmin, tmax and dt do not define an integer number of linear steps."];
        Abort[];
    ];
    nSamples = nSteps + 1;
    timeGrid = None;
    maximumRequestedTime = tmax,

    "Log",
    If[logTMax <= logTMin || !IntegerQ[pointsPerDecade] || pointsPerDecade < 1,
        Print["ERROR: Log grid requires logTMax > logTMin and pointsPerDecade >= 1."];
        Abort[];
    ];
    nLogIntervals = Max[1, Round[(logTMax - logTMin)*pointsPerDecade]];
    positiveLogTimeGrid = Developer`ToPackedArray[
        N[10.^Subdivide[logTMin, logTMax, nLogIntervals]]
    ];
    timeGrid = Developer`ToPackedArray[
        If[
            TrueQ[includeZeroInLogGrid],
            Join[{0.}, positiveLogTimeGrid],
            positiveLogTimeGrid
        ]
    ];
    nSamples = Length[timeGrid];
    nSteps = nSamples - 1;
    maximumRequestedTime = Max[timeGrid],

    _,
    Print["ERROR: timeGridType must be \"Linear\" or \"Log\"."];
    Abort[];
];


If[maximumRequestedTime >= 10^9,
    Print[
        "WARNING: requested times reach ", maximumRequestedTime,
        ". Direct spectral phases remove recurrence drift, but machine-precision ",
        "eigenvalue errors still accumulate approximately as t deltaE. ",
        "Representative long-time points near 10^10 should be independently ",
        "validated at higher precision."
    ];
];


(* ============================================================ *)
(* Parameter-resolved output directory                           *)
(* ============================================================ *)

ClearAll[tagNumber];
tagNumber[x_?NumericQ] := Module[{s},
    s = ToString[InputForm[N[Round[x, 10^-12]]]];
    s = StringReplace[s, RegularExpression["\\.$"] -> ""];
    StringReplace[s, {"-" -> "m", "." -> "p", " " -> ""}]
];

hxExponentMin = First[hxExponents];
hxExponentMax = Last[hxExponents];
hxExponentStep = If[
    Length[hxExponents] > 1,
    hxExponents[[2]] - hxExponents[[1]],
    0.
];


timeTag = Switch[
    timeGridType,
    "Linear",
    StringRiffle[
        {
            "grid_linear",
            "tmin_" <> tagNumber[tmin],
            "tmax_" <> tagNumber[tmax],
            "dt_" <> tagNumber[dt]
        },
        "_"
    ],
    "Log",
    StringRiffle[
        {
            "grid_log",
            "logt_" <> tagNumber[logTMin] <> "_to_" <> tagNumber[logTMax],
            "ppd_" <> ToString[pointsPerDecade],
            "zero_" <> ToString[includeZeroInLogGrid]
        },
        "_"
    ]
];


configTag = StringRiffle[
    {
        "L_" <> ToString[L],
        "J_" <> tagNumber[J],
        "hz_" <> tagNumber[hz],
        "hxexp_" <> tagNumber[hxExponentMin] <>
            "_to_" <> tagNumber[hxExponentMax] <>
            "_step_" <> tagNumber[hxExponentStep],
        timeTag,
        "seed_" <> ToString[seed]
    },
    "_"
];


configDir = FileNameJoin[
    {NotebookDirectory[], "prethermal_entropy_hx_sweeps", configTag}
];

If[!DirectoryQ[configDir],
    CreateDirectory[configDir, CreateIntermediateDirectories -> True]
];

existingRunDirs = Select[FileNames["run_*", configDir], DirectoryQ];
existingRunNumbers = Cases[
    FileNameTake /@ existingRunDirs,
    s_String :> Quiet@Check[
        ToExpression[StringReplace[s, StartOfString ~~ "run_" -> ""]],
        Nothing
    ]
];

runIndex = If[existingRunNumbers === {}, 1, Max[existingRunNumbers] + 1];
outputDir = FileNameJoin[
    {configDir, "run_" <> IntegerString[runIndex, 10, 3]}
];
CreateDirectory[outputDir];

Print["Output configuration: ", configTag];
Print["Output run directory: ", outputDir];
Print["Time-grid type: ", timeGridType, "   samples per hx = ", nSamples];
Print["Batch size: ", batchSize];


(* ::Title::Closed:: *)
(*Fixed reflection bases and fixed RPS*)


bases = reflectionSectorBases[L];
Bp = bases["Even"];
Bm = bases["Odd"];

(* Numerical sparse bases for production linear algebra. *)
BpN = N[Bp];
BmN = N[Bm];


(* One-time basis validation. *)
basisOrthogonalityError = Max[
    Norm[
        ConjugateTranspose[BpN] . BpN -
        IdentityMatrix[Dimensions[BpN][[2]]],
        "Frobenius"
    ],
    Norm[
        ConjugateTranspose[BmN] . BmN -
        IdentityMatrix[Dimensions[BmN][[2]]],
        "Frobenius"
    ],
    Norm[
        ConjugateTranspose[BpN] . BmN,
        "Frobenius"
    ]
];

Print["Reflection-basis orthogonality error = ", basisOrthogonalityError];
If[!TrueQ[basisOrthogonalityError <= reconstructionTolerance], Abort[]];


SeedRandom[seed];
psi0 = Developer`ToPackedArray[N[RandomChainProductState[L]]];

If[Length[psi0] =!= 2^L,
    Print["ERROR: RandomChainProductState returned the wrong dimension."];
    Abort[];
];

psi0NormError = Abs[Norm[psi0]^2 - 1.];
If[psi0NormError > 10^-12,
    Print[
        "ERROR: initial RPS is not normalized. Norm-squared error = ",
        psi0NormError
    ];
    Abort[];
];


psi0P = Developer`ToPackedArray[
    N[ConjugateTranspose[BpN] . psi0]
];
psi0M = Developer`ToPackedArray[
    N[ConjugateTranspose[BmN] . psi0]
];

sectorNormError = Abs[Norm[psi0P]^2 + Norm[psi0M]^2 - 1.];
Print["Fixed RPS sector norm error = ", sectorNormError];
If[!TrueQ[sectorNormError <= reconstructionTolerance], Abort[]];


manifest = <|
    "L" -> L,
    "J" -> N[J],
    "hz" -> N[hz],
    "seed" -> seed,
    "psi0" -> psi0,
    "hxExponents" -> hxExponents,
    "hxList" -> hxList,
    "timeGridType" -> timeGridType,
    "tmin" -> N[tmin],
    "tmax" -> N[tmax],
    "dt" -> N[dt],
    "logTMin" -> N[logTMin],
    "logTMax" -> N[logTMax],
    "pointsPerDecade" -> pointsPerDecade,
    "includeZeroInLogGrid" -> includeZeroInLogGrid,
    "nSamples" -> nSamples,
    "nSteps" -> nSteps,
    "resetEvery" -> resetEvery,
    "batchSize" -> batchSize,
    "requestedKernels" -> requestedKernels,
    "configTag" -> configTag,
    "runIndex" -> runIndex,
    "outputDir" -> outputDir,
    "entropyLogBase" -> E,
    "entropyProbabilityCutoff" -> 1.*^-14,
    "basisOrthogonalityError" -> basisOrthogonalityError,
    "psi0NormError" -> psi0NormError,
    "sectorNormError" -> sectorNormError,
    "sectorDimensions" -> {
        Dimensions[BpN][[2]],
        Dimensions[BmN][[2]]
    },
    "WolframVersion" -> $Version,
    "SystemID" -> $SystemID
|>;

If[timeGridType === "Log",
    manifest = Append[manifest, "timeGrid" -> timeGrid]
];

manifestPath = FileNameJoin[{outputDir, "manifest.wxf"}];
manifestExportResult = Export[manifestPath, manifest, "WXF"];
If[!StringQ[manifestExportResult] || !FileExistsQ[manifestPath],
    Print["ERROR: manifest export failed."];
    Abort[];
];


(* ::Title:: *)
(*Automatic hx sweep*)


(* One complete hx trajectory.  All temporary symbols are local so
   different hx values can be evaluated safely on different kernels. *)
ClearAll[runOneHx];
runOneHx[k_Integer] := Module[
    {
        hxExponent, hx, runStart, H, hScale, hermiticityError,
        crossBlock, reflectionBlockError, Hp, Hm, diagSeconds,
        eigensystems, evalsP, evecsP, evalsM, evecsM,
        eigResidualP, eigResidualM, coeffP, coeffM,
        maxImagV, WP, WM, reconstructedPsi0, reconstructionError,
        weightsP, weightsM, meanEnergy, energyVariance,
        evolutionSeconds, evolutionResult, entropy, evolutionNormError,
        validationIndices, validationTimes,
        validationResults, validationEntropyDifferences,
        validationNormErrors, maxValidationEntropyDifference,
        maxValidationNormError, runFile, resultAssociation,
        computeSeconds, exportSeconds, exportResult, totalSeconds
    },

    hxExponent = hxExponents[[k]];
    hx = hxList[[k]];
    runStart = AbsoluteTime[];


    (* --------------------------------------------------------- *)
    (* Hamiltonian checks before any Chop or block diagonalization *)
    (* --------------------------------------------------------- *)

    H = N[IsingHamiltonian[hx, hz, J, L]];
    hScale = Max[1., Norm[H, "Frobenius"]];

    hermiticityError =
        Norm[H - ConjugateTranspose[H], "Frobenius"]/hScale;

    If[hermiticityError > hermiticityTolerance,
        Return[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Hermiticity check failed.",
            "hermiticityError" -> hermiticityError
        |>]
    ];

    crossBlock = ConjugateTranspose[BpN] . H . BmN;
    reflectionBlockError = Norm[crossBlock, "Frobenius"]/hScale;

    If[reflectionBlockError > reflectionTolerance,
        Return[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Reflection block check failed.",
            "reflectionBlockError" -> reflectionBlockError
        |>]
    ];

    Hp = ConjugateTranspose[BpN] . H . BpN;
    Hm = ConjugateTranspose[BmN] . H . BmN;


    (* --------------------------------------------------------- *)
    (* Separate parity-sector diagonalizations                   *)
    (* --------------------------------------------------------- *)

    diagSeconds = First[AbsoluteTiming[
        eigensystems = {
            sortedEigensystem[Hp],
            sortedEigensystem[Hm]
        };
    ]];

    {{evalsP, evecsP}, {evalsM, evecsM}} = eigensystems;

    eigResidualP = sampledEigensystemResidual[Hp, evalsP, evecsP];
    eigResidualM = sampledEigensystemResidual[Hm, evalsM, evecsM];

    If[Max[eigResidualP, eigResidualM] > eigensystemResidualTolerance,
        Return[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Sampled eigensystem residual check failed.",
            "eigResidualP" -> eigResidualP,
            "eigResidualM" -> eigResidualM
        |>]
    ];


    (* --------------------------------------------------------- *)
    (* Fixed RPS expressed in the two new eigenbases             *)
    (* --------------------------------------------------------- *)

    coeffP = Conjugate[evecsP] . psi0P;
    coeffM = Conjugate[evecsM] . psi0M;

    maxImagV = Max[
        Max[Abs[Im[evecsP]]],
        Max[Abs[Im[evecsM]]]
    ];

    If[maxImagV > 10^-10,
        Return[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "The sector eigenbasis is not numerically real.",
            "maxImagV" -> maxImagV
        |>]
    ];


    (* --------------------------------------------------------- *)
    (* Sector-sized reconstruction matrices: no full D x D Vfull *)
    (* --------------------------------------------------------- *)

    WP = Developer`ToPackedArray[Re[Transpose[evecsP]]];
    WM = Developer`ToPackedArray[Re[Transpose[evecsM]]];

    reconstructedPsi0 =
        BpN . (WP . coeffP) +
        BmN . (WM . coeffM);

    reconstructionError = Norm[reconstructedPsi0 - psi0];

    If[reconstructionError > reconstructionTolerance,
        Return[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Initial-state reconstruction check failed.",
            "reconstructionError" -> reconstructionError
        |>]
    ];


    (* --------------------------------------------------------- *)
    (* Save spectral weights and energy distribution of the RPS  *)
    (* --------------------------------------------------------- *)

    weightsP = Developer`ToPackedArray[Abs[coeffP]^2];
    weightsM = Developer`ToPackedArray[Abs[coeffM]^2];

    meanEnergy = Total[weightsP*evalsP] + Total[weightsM*evalsM];
    energyVariance =
        Total[weightsP*(evalsP - meanEnergy)^2] +
        Total[weightsM*(evalsM - meanEnergy)^2];


    (* --------------------------------------------------------- *)
    (* Long-time entropy evolution                               *)
    (* --------------------------------------------------------- *)

    evolutionSeconds = First[AbsoluteTiming[
        evolutionResult = Switch[
            timeGridType,

            "Linear",
            linearEntropyEvolutionSector[
                WP, WM, BpN, BmN,
                evalsP, evalsM, coeffP, coeffM,
                N[tmin], N[dt], nSteps,
                batchSize, resetEvery, dA, dB
            ],

            "Log",
            arbitraryEntropyEvolutionSector[
                WP, WM, BpN, BmN,
                evalsP, evalsM, coeffP, coeffM,
                timeGrid, batchSize, dA, dB
            ],
            _, Failure["InvalidTimeGridType", <|"value" -> timeGridType|>]
        ];
    ]];

    If[!AssociationQ[evolutionResult],
        Return[<|"success" -> False, "index" -> k, "hx" -> hx,
            "message" -> "Evolution did not return a result association. Reevaluate all definitions and worker setup.",
            "resultHead" -> Head[evolutionResult]|>]
    ];
    entropy = Lookup[evolutionResult, "entropy", Missing["NotAvailable"]];
    If[!validEntropyVectorQ[entropy, nSamples],
        Return[<|"success" -> False, "index" -> k, "hx" -> hx,
            "message" -> "Evolution did not return the expected finite real entropy vector.",
            "resultHead" -> Head[entropy], "expectedLength" -> nSamples|>]
    ];
    evolutionNormError = Lookup[evolutionResult, "maxNormError", Infinity];
    If[!TrueQ[evolutionNormError <= validationNormTolerance],
        Return[<|"success" -> False, "index" -> k, "hx" -> hx,
            "message" -> "Norm drift in the evolved batches exceeded tolerance.",
            "maxEvolutionNormError" -> evolutionNormError|>]
    ];
    Clear[evolutionResult];


    (* --------------------------------------------------------- *)
    (* Lightweight direct-phase checks at selected samples       *)
    (* --------------------------------------------------------- *)

    validationIndices = DeleteDuplicates[
        Clip[
            Switch[
                timeGridType,
                "Linear",
                {
                    1,
                    1 + Min[nSteps, Max[0, resetEvery - 1]],
                    1 + Min[nSteps, resetEvery],
                    1 + Floor[nSteps/2],
                    nSamples
                },
                "Log",
                {
                    1,
                    1 + Floor[(nSamples - 1)/3],
                    1 + Floor[2 (nSamples - 1)/3],
                    nSamples
                }
            ],
            {1, nSamples}
        ]
    ];

    validationTimes = Switch[
        timeGridType,
        "Linear",
        N[tmin + (validationIndices - 1)*dt],
        "Log",
        N[timeGrid[[validationIndices]]]
    ];

    validationResults = Table[
        directEntropyAndNormSector[
            WP, WM, BpN, BmN,
            evalsP, evalsM, coeffP, coeffM,
            N[validationTimes[[j]]], dA, dB
        ],
        {j, Length[validationIndices]}
    ];

    validationEntropyDifferences = Abs[
        entropy[[validationIndices]] - validationResults[[All, 1]]
    ];
    validationNormErrors = validationResults[[All, 2]];
    maxValidationEntropyDifference = Max[validationEntropyDifferences];
    maxValidationNormError = Max[validationNormErrors];

    If[
        maxValidationEntropyDifference > validationEntropyTolerance ||
        maxValidationNormError > validationNormTolerance,
        Return[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Long-time validation check failed.",
            "maxValidationEntropyDifference" -> maxValidationEntropyDifference,
            "maxValidationNormError" -> maxValidationNormError
        |>]
    ];


    (* --------------------------------------------------------- *)
    (* Independent result file for this hx                       *)
    (* --------------------------------------------------------- *)

    runFile = FileNameJoin[
        {outputDir, "run_" <> IntegerString[k, 10, 2] <> ".wxf"}
    ];

    computeSeconds = N[AbsoluteTime[] - runStart];

    resultAssociation = <|
        "index" -> k,
        "L" -> L,
        "J" -> N[J],
        "hz" -> N[hz],
        "hxExponent" -> N[hxExponent],
        "hx" -> N[hx],
        "timeGridType" -> timeGridType,
        "tmin" -> N[tmin],
        "tmax" -> N[tmax],
        "dt" -> N[dt],
        "logTMin" -> N[logTMin],
        "logTMax" -> N[logTMax],
        "pointsPerDecade" -> pointsPerDecade,
        "includeZeroInLogGrid" -> includeZeroInLogGrid,
        "nSamples" -> nSamples,
    "nSteps" -> nSteps,
        "entropy" -> entropy,
        "evalsP" -> evalsP,
        "evalsM" -> evalsM,
        "weightsP" -> weightsP,
        "weightsM" -> weightsM,
        "meanEnergy" -> meanEnergy,
        "energyVariance" -> energyVariance,
        "hermiticityError" -> hermiticityError,
        "reflectionBlockError" -> reflectionBlockError,
        "eigResidualP" -> eigResidualP,
        "eigResidualM" -> eigResidualM,
        "reconstructionError" -> reconstructionError,
        "maxImagV" -> maxImagV,
        "validationIndices" -> validationIndices,
        "validationTimes" -> validationTimes,
        "validationEntropyDifferences" -> validationEntropyDifferences,
        "validationNormErrors" -> validationNormErrors,
        "maxValidationEntropyDifference" -> maxValidationEntropyDifference,
        "maxValidationNormError" -> maxValidationNormError,
        "maxEvolutionNormError" -> evolutionNormError,
        "diagonalizationSeconds" -> diagSeconds,
        "evolutionSeconds" -> evolutionSeconds,
        "computeSecondsBeforeExport" -> computeSeconds
    |>;

    If[timeGridType === "Log",
        resultAssociation = Append[resultAssociation, "timeGrid" -> timeGrid]
    ];

    exportSeconds = First[AbsoluteTiming[
        exportResult = Export[runFile, resultAssociation, "WXF"];
    ]];
    totalSeconds = computeSeconds + exportSeconds;

    If[!StringQ[exportResult] || !FileExistsQ[runFile],
        Return[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Export failed.",
            "runFile" -> runFile
        |>]
    ];

    <|
        "success" -> True,
        "index" -> k,
        "hxExponent" -> hxExponent,
        "hx" -> hx,
        "diagonalizationSeconds" -> diagSeconds,
        "evolutionSeconds" -> evolutionSeconds,
        "exportSeconds" -> exportSeconds,
        "totalSeconds" -> totalSeconds,
        "hermiticityError" -> hermiticityError,
        "reflectionBlockError" -> reflectionBlockError,
        "eigResidualP" -> eigResidualP,
        "eigResidualM" -> eigResidualM,
        "reconstructionError" -> reconstructionError,
        "maxValidationEntropyDifference" -> maxValidationEntropyDifference,
        "maxValidationNormError" -> maxValidationNormError,
        "maxEvolutionNormError" -> evolutionNormError,
        "runFile" -> runFile
    |>
];


(* ------------------------------------------------------------ *)
(* Launch and initialize eight local worker kernels             *)
(* ------------------------------------------------------------ *)

CloseKernels[];
LaunchKernels[requestedKernels];
activeKernels = Length[Kernels[]];

Print["Parallel kernels active = ", activeKernels];
If[activeKernels < requestedKernels,
    Print["WARNING: fewer kernels were launched than requested."]
];

(* The private QMB package is explicitly loaded on every worker. *)
If[activeKernels == 0,
    Print["ERROR: no parallel worker kernels are available."]; Abort[];
];
(* With inserts the literal path before ParallelEvaluate sends the expression. *)
workerPackageStatus = With[{initPath = qmbInitPath},
    ParallelEvaluate[Quiet[Check[Get[initPath]; True, False]]]
];
If[!AllTrue[workerPackageStatus, TrueQ],
    Print["ERROR: QMB package initialization failed on a worker."];
    Abort[];
];

DistributeDefinitions[
    runOneHx, validEntropyVectorQ,
    sortedEigensystem,
    phaseStepC,
    phaseBatchC,
    phaseAtTimesBatchC,
    reseedC,
    entropyFromSVC,
    entropyFast,
    entropyBatchFromSectors,
    linearEntropyEvolutionSector,
    arbitraryEntropyEvolutionSector,
    directEntropyAndNormSector,
    sampledEigensystemResidual,
    BpN, BmN, psi0, psi0P, psi0M,
    hxExponents, hxList,
    L, J, hz, dA, dB,
    timeGridType, tmin, tmax, dt, nSteps, nSamples, timeGrid,
    logTMin, logTMax, pointsPerDecade, includeZeroInLogGrid,
    resetEvery, batchSize, outputDir,
    hermiticityTolerance, reflectionTolerance,
    reconstructionTolerance, eigensystemResidualTolerance,
    validationEntropyTolerance, validationNormTolerance
];


(* ------------------------------------------------------------ *)
(* Coarse-grained parallel hx sweep                             *)
(* ------------------------------------------------------------ *)

(* Small worker smoke test: catches missing/stale definitions before the sweep. *)
workerChecks = ParallelEvaluate[
    Quiet[Check[
        Module[{id = IdentityMatrix[2], lin, arb},
            (* A single two-dimensional sector with zero odd-sector weight. *)
            lin = linearEntropyEvolutionSector[
                id, id, N[id], N[id], {0., 1.}, {0., 1.},
                {1. + 0. I, 0. + 0. I}, {0. + 0. I, 0. + 0. I},
                0., 0.25, 4, 3, 2, 1, 2
            ];
            arb = arbitraryEntropyEvolutionSector[
                id, id, N[id], N[id], {0., 1.}, {0., 1.},
                {1. + 0. I, 0. + 0. I}, {0. + 0. I, 0. + 0. I},
                {0., 0.1, 1.}, 2, 1, 2
            ];
            AssociationQ[lin] && AssociationQ[arb] &&
            validEntropyVectorQ[lin["entropy"], 5] &&
            validEntropyVectorQ[arb["entropy"], 3] &&
            Max[Abs[lin["entropy"]]] < 10^-12 &&
            Max[Abs[arb["entropy"]]] < 10^-12
        ], False]]
];
If[!AllTrue[workerChecks, TrueQ],
    Print["ERROR: numerical worker startup check failed: ", workerChecks]; Abort[];
];

(* Store results inside the timed expression: a final semicolon is safe. *)
sweepSeconds = First[AbsoluteTiming[
    runSummary = ParallelMap[
        runOneHx,
        Range[Length[hxList]],
        Method -> "CoarsestGrained"
    ];
]];

If[!ListQ[runSummary] || Length[runSummary] != Length[hxList],
    Print["ERROR: the sweep did not return one result per hx."]; Abort[];
];
failedRuns = Select[runSummary,
    !AssociationQ[#] || !TrueQ[Lookup[#, "success", False]] &
];

If[failedRuns =!= {},
    Print["ERROR: one or more hx runs failed validation."];
    Print[failedRuns];
    Export[FileNameJoin[{outputDir, "failed_runs.wxf"}], failedRuns, "WXF"];
    Abort[];
];


Do[
    Print[""];
    Print["============================================================"];
    Print[
        "Run ", r["index"], "/", Length[hxList],
        "   exponent = ", r["hxExponent"],
        "   hx = ", r["hx"]
    ];
    Print["Diagonalization time [s] = ", r["diagonalizationSeconds"]];
    Print["Evolution time [s] = ", r["evolutionSeconds"]];
    Print["Export time [s] = ", r["exportSeconds"]];
    Print["Total worker run time [s] = ", r["totalSeconds"]];
    Print["Reconstruction error = ", r["reconstructionError"]];
    Print[
        "Max direct-phase entropy check = ",
        r["maxValidationEntropyDifference"]
    ];
    Print[
        "Max norm error at validation samples = ",
        r["maxValidationNormError"]
    ];
    Print["Saved: ", r["runFile"]],
    {r, runSummary}
];


(* ::Title::Closed:: *)
(*Sweep summary*)


summaryHeader = {
    "index",
    "hxExponent",
    "hx",
    "diagonalizationSeconds",
    "evolutionSeconds",
    "exportSeconds",
    "totalSeconds",
    "hermiticityError",
    "reflectionBlockError",
    "eigResidualP",
    "eigResidualM",
    "reconstructionError",
    "maxValidationEntropyDifference",
    "maxValidationNormError",
    "file"
};

summaryRows = (
    {
        #["index"],
        #["hxExponent"],
        #["hx"],
        #["diagonalizationSeconds"],
        #["evolutionSeconds"],
        #["exportSeconds"],
        #["totalSeconds"],
        #["hermiticityError"],
        #["reflectionBlockError"],
        #["eigResidualP"],
        #["eigResidualM"],
        #["reconstructionError"],
        #["maxValidationEntropyDifference"],
        #["maxValidationNormError"],
        #["runFile"]
    } &
) /@ runSummary;

Export[
    FileNameJoin[{outputDir, "run_summary.csv"}],
    Prepend[summaryRows, summaryHeader]
];

Export[
    FileNameJoin[{outputDir, "run_summary.wxf"}],
    runSummary,
    "WXF"
];

Print[""];
Print["============================================================"];
Print["Sweep complete."];
Print["Configuration: ", configTag];
Print["Run index: ", runIndex];
Print["Output directory: ", outputDir];
Print["Number of hx values: ", Length[hxList]];
Print["Requested parallel kernels: ", requestedKernels];
Print["Active parallel kernels: ", activeKernels];
Print["Whole parallel sweep wall time [s] = ", sweepSeconds];
Print["============================================================"];


(* ::Title::Closed:: *)
(*Plotting and export*)


(* Plotting remains export-only.  No thinning, interpolation, or joining
   is applied to the raw data.  A guard prevents accidental rendering of
   hundreds of millions of points unless explicitly requested. *)

exportRawLinearPlot = True;
exportRawLogPlot = True;
exportMovingAverageLogPlot = True;

movingAverageWindow = 200;

maxAutomaticPlotPoints = 5*10^6;
forceHugePlotExport = False;

totalRequestedPlotPoints = Length[hxList]*nSamples;

If[
    totalRequestedPlotPoints > maxAutomaticPlotPoints &&
    !TrueQ[forceHugePlotExport],

    Print[
        "Plot export skipped automatically because the combined dataset contains ",
        totalRequestedPlotPoints,
        " points. No thinning is performed. Set forceHugePlotExport = True ",
        "if you explicitly want Mathematica to render every point."
    ],


    (* ======================================================== *)
    (* Load all completed hx runs                              *)
    (* ======================================================== *)

    dataDir = outputDir;

    files = Sort[
        FileNames[
            "run_" ~~ DigitCharacter ~~ DigitCharacter ~~ ".wxf",
            dataDir
        ]
    ];

    runs = Import[#, "WXF"] & /@ files;


    (* ======================================================== *)
    (* Plot settings                                           *)
    (* ======================================================== *)

    plotDir = FileNameJoin[{dataDir, "plots"}];
    If[!DirectoryQ[plotDir], CreateDirectory[plotDir]];

    cf = Blend[
        {
            {0., Blue},
            {0.5, Purple},
            {1., Red}
        },
        #
    ] &;

    plotHxExponentMin = Min[runs[[All, "hxExponent"]]];
    plotHxExponentMax = Max[runs[[All, "hxExponent"]]];

    colors = If[plotHxExponentMin == plotHxExponentMax,
        ConstantArray[cf[0.5], Length[runs]],
        cf /@ Rescale[runs[[All, "hxExponent"]],
            {plotHxExponentMin, plotHxExponentMax}]
    ];


    (* Numerical 10^x color-bar labels; never rational fractions. *)
    barTicks = Table[
        {x, powerOfTenLabel[N[x], 24]},
        {x, Ceiling[plotHxExponentMin], Floor[plotHxExponentMax], 1}
    ];

    bar = BarLegend[
        {cf, {plotHxExponentMin, plotHxExponentMax}},
        LegendLabel -> Style[Subscript["h", "x"], 30, Black],
        LabelStyle -> {Black, 20},
        ColorFunctionScaling -> True,
        LegendMarkerSize -> {12, 360},
        Ticks -> barTicks
    ];

    plotTitle = Style[
        Row[
            {
                "RPS, L = ", N[L],
                ", J = ", N[J],
                ", ", Subscript["h", "z"], " = ", N[hz]
            }
        ],
        34,
        Black
    ];


    (* ======================================================== *)
    (* Raw data coordinates                                    *)
    (* ======================================================== *)

    rawData = Table[
        runTimes = If[
            runs[[k]]["timeGridType"] === "Linear",
            Developer`ToPackedArray[
                N[
                    runs[[k]]["tmin"] +
                    Range[0, runs[[k]]["nSamples"] - 1]*runs[[k]]["dt"]
                ]
            ],
            runs[[k]]["timeGrid"]
        ];

        Transpose[{runTimes, runs[[k]]["entropy"]}],
        {k, Length[runs]}
    ];


    (* ======================================================== *)
    (* RAW: linear-time points only                            *)
    (* ======================================================== *)

    If[TrueQ[exportRawLinearPlot],
        rawLinearPlot = Legended[
            ListPlot[
                rawData,
                Joined -> False,
                PlotStyle -> (
                    Directive[#, Opacity[0.75], PointSize[0.002]] &
                ) /@ colors,
                PlotRange -> {All, {0., All}},
                Frame -> True,
                Axes -> False,
                ImageSize -> 1000,
                PlotTheme -> "Detailed",
                FrameStyle -> Directive[Black, 26],
                PlotLabel -> plotTitle,
                FrameLabel -> {
                    Style["t", 30, Black],
                    Style["S(t)", 30, Black]
                },
                GridLines -> None
            ],
            Placed[bar, Right]
        ];

        Export[
            FileNameJoin[{plotDir, "entropy_hx_raw_linear.png"}],
            rawLinearPlot,
            ImageResolution -> 200
        ];
    ];


    (* ======================================================== *)
    (* Positive times and explicit decade ticks                *)
    (* ======================================================== *)

    rawDataPositiveTime = (
        Select[#, First[#] > 0. &] &
    ) /@ rawData;

    If[AllTrue[rawDataPositiveTime, Length[#] == 0 &],
        exportRawLogPlot = False;
        exportMovingAverageLogPlot = False;
        rawDataPositiveTime = {{{1., 0.}}}; (* ticks only; no log plot exported *)
    ];
    minPositivePlotTime = Min[
        Flatten[rawDataPositiveTime[[All, All, 1]]]
    ];

    maxPlotTime = Max[
        Flatten[rawDataPositiveTime[[All, All, 1]]]
    ];

    logTimeTicks = Table[
        {10.^x, powerOfTenLabel[N[x], 22]},
        {
            x,
            Ceiling[Log10[minPositivePlotTime]],
            Floor[Log10[maxPlotTime]],
            1
        }
    ];


    (* ======================================================== *)
    (* RAW: logarithmic-time points only                       *)
    (* ======================================================== *)

    If[TrueQ[exportRawLogPlot],
        rawLogPlot = Legended[
            ListLogLinearPlot[
                rawDataPositiveTime,
                Joined -> False,
                PlotStyle -> (
                    Directive[#, Opacity[0.75], PointSize[0.002]] &
                ) /@ colors,
                PlotRange -> {All, {0., All}},
                Frame -> True,
                Axes -> False,
                ImageSize -> 1000,
                PlotTheme -> "Detailed",
                FrameStyle -> Directive[Black, 26],
                PlotLabel -> plotTitle,
                FrameLabel -> {
                    Style["t", 30, Black],
                    Style["S(t)", 30, Black]
                },
                FrameTicks -> {
                    {Automatic, None},
                    {logTimeTicks, None}
                },
                GridLines -> None
            ],
            Placed[bar, Right]
        ];

        Export[
            FileNameJoin[{plotDir, "entropy_hx_raw_logtime.png"}],
            rawLogPlot,
            ImageResolution -> 200
        ];
    ];


    (* ======================================================== *)
    (* Full-window moving average, as in the working workflow.    *)
    (* Both time and entropy are averaged; raw data are unchanged. *)
    (* ======================================================== *)

    If[TrueQ[exportMovingAverageLogPlot] && nSamples >= movingAverageWindow,
        movingAverageDataPositiveTime = (
            Select[#, First[#] > 0. &] &
        ) /@ movingAverageData;
        
        movingAverageData = Map[
    Function[data,
        Module[{times, values, sums, n, first},
            times = data[[All, 1]];
            values = data[[All, 2]];
            sums = Prepend[Accumulate[values], 0.];
            n = Length[values];

            Transpose[{
                times,
                Table[
                    first = Max[1, i - movingAverageWindow + 1];
                    (sums[[i + 1]] - sums[[first]])/
                        (i - first + 1),
                    {i, n}
                ]
            }]
        ]
    ],
    rawData
];

        movingAverageTitle = Style[
            Row[
                {
                    "RPS, L = ", N[L],
                    ", J = ", N[J],
                    ", ", Subscript["h", "z"], " = ", N[hz],
                    ", moving average = ", N[movingAverageWindow]
                }
            ],
            34,
            Black
        ];

        movingAverageLogPlot = Legended[
            ListLogLinearPlot[
                movingAverageDataPositiveTime,
                Joined -> False,
                PlotStyle -> (
                    Directive[#, Opacity[0.90], PointSize[0.0025]] &
                ) /@ colors,
                PlotRange -> {All, {0., All}},
                Frame -> True,
                Axes -> False,
                ImageSize -> 1000,
                PlotTheme -> "Detailed",
                FrameStyle -> Directive[Black, 26],
                PlotLabel -> movingAverageTitle,
                FrameLabel -> {
                    Style["t", 30, Black],
                    Style["S(t)", 30, Black]
                },
                FrameTicks -> {
                    {Automatic, None},
                    {logTimeTicks, None}
                },
                GridLines -> None
            ],
            Placed[bar, Right]
        ];

        Export[
            FileNameJoin[
                {
                    plotDir,
                    "entropy_hx_moving_average_" <>
                    ToString[movingAverageWindow] <>
                    "_logtime.png"
                }
            ],
            movingAverageLogPlot,
            ImageResolution -> 200
        ];
    ];


    (* No graphics returned to the notebook. *)
    Clear[rawLinearPlot, rawLogPlot, movingAverageLogPlot];
];


(* Optional cleanup after all work and exports. *)
CloseKernels[];
