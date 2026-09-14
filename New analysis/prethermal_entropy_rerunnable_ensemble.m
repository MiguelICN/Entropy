(* ::Package:: *)

(* ::Title:: *)
(* Setup*)


qmbInitPath = "C:\\Users\\Miguel\\Github\\libs\\QMB\\Kernel\\init.m";
If[!FileExistsQ[qmbInitPath], Print["ERROR: QMB init.m not found: ", qmbInitPath]; Abort[]];
Get[qmbInitPath];
SetDirectory[NotebookDirectory[]];



(* ::Title::Closed:: *)
(* Definitions and compiled kernels *)


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



(* Filesystem helpers: only letters, digits and underscores in parameter tags. *)
ClearAll[tagNumber, makeRunDirectory, checkedExport];
tagNumber[x_?NumericQ] := Module[{s},
    s = ToString[InputForm[N[x]]];
    s = StringReplace[s, {"*^" -> "e", "-" -> "m", "." -> "p", "+" -> ""}];
    s = StringReplace[s, RegularExpression["[^A-Za-z0-9_]" ] -> ""];
    s = StringReplace[s, {"pe" -> "e", RegularExpression["p$"] -> ""}];
    StringTake[s, UpTo[20]]
];
checkedExport[path_String, data_, format_String] := Module[{r},
    r = Quiet[Check[Export[path, data, format], $Failed]];
    If[!StringQ[r] || !FileExistsQ[path],
        Failure["ExportFailed", <|"path" -> path|>], r]
];
makeRunDirectory[root_String, size_Integer, coupling_, field_, exponents_List] := Module[
    {tag, config, directory, probe, attempt = 0, result},
    tag = "L" <> ToString[size] <> "_J" <> tagNumber[coupling] <>
        "_hz" <> tagNumber[field] <> "_hx" <> tagNumber[Min[exponents]] <>
        "to" <> tagNumber[Max[exponents]];
    config = FileNameJoin[{root, "entropy_sweeps", tag}];
    If[!DirectoryQ[config],
        Quiet[Check[CreateDirectory[config, CreateIntermediateDirectories -> True], Null]]
    ];
    If[!DirectoryQ[config], Return[Failure["CreateDirectoryFailed", <|"path" -> config|>]]];
    While[attempt < 10,
        attempt++;
        directory = FileNameJoin[{config, "run_" <> StringTake[CreateUUID[], 12]}];
        result = Quiet[Check[CreateDirectory[directory], $Failed]];
        If[StringQ[result] && DirectoryQ[directory], Break[]]
    ];
    If[!StringQ[result] || !DirectoryQ[directory],
        Return[Failure["CreateRunFailed", <|"path" -> config|>]]];
    probe = FileNameJoin[{directory, "write_probe.tmp"}];
    result = checkedExport[probe, "write test", "Text"];
    If[FailureQ[result], Return[result]];
    DeleteFile[probe];
    directory
];

ClearAll[runOneHx];
runOneHx[cfg_Association, k_Integer] := Catch[Module[
    {
        L = cfg["L"], J = cfg["J"], hz = cfg["hz"],
        hxExponents = cfg["hxExponents"], hxList = cfg["hxList"],
        dA = cfg["dA"], dB = cfg["dB"], BpN = cfg["BpN"], BmN = cfg["BmN"],
        timeGridType = cfg["timeGridType"], tmin = cfg["tmin"], tmax = cfg["tmax"],
        dt = cfg["dt"], nSteps = cfg["nSteps"], nSamples = cfg["nSamples"],
        timeGrid = cfg["timeGrid"], resetEvery = cfg["resetEvery"],
        batchSize = cfg["batchSize"], outputDir = cfg["outputDir"],
        logTMin = cfg["logTMin"], logTMax = cfg["logTMax"],
        pointsPerDecade = cfg["pointsPerDecade"], includeZeroInLogGrid = cfg["includeZeroInLogGrid"],
        hermiticityTolerance = cfg["hermiticityTolerance"], reflectionTolerance = cfg["reflectionTolerance"],
        reconstructionTolerance = cfg["reconstructionTolerance"],
        eigensystemResidualTolerance = cfg["eigensystemResidualTolerance"],
        validationEntropyTolerance = cfg["validationEntropyTolerance"],
        validationNormTolerance = cfg["validationNormTolerance"],
        stateCount = cfg["stateCount"], member, psi0, psi0P, psi0M,
        meanEntropy, worstReconstruction = 0., worstEntropyDifference = 0.,
        worstDirectNorm = 0., worstEvolutionNorm = 0., ensembleSeconds,
        hxExponent, hx, runStart, H, hScale, hermiticityError,
        crossBlock, reflectionBlockError, Hp, Hm, diagSeconds,
        eigensystems, evalsP, evecsP, evalsM, evecsM,
        eigResidualP, eigResidualM, coeffP, coeffM,
        maxImagV, WP, WM, reconstructedPsi0, reconstructionError,
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
    If[!DirectoryQ[outputDir], Throw[<|"success" -> False, "index" -> k,
        "message" -> "Output directory is not visible on this kernel.", "path" -> outputDir|>, "runHxFailure"]];


    (* --------------------------------------------------------- *)
    (* Hamiltonian checks before any Chop or block diagonalization *)
    (* --------------------------------------------------------- *)

    H = N[IsingHamiltonian[hx, hz, J, L]];
    hScale = Max[1., Norm[H, "Frobenius"]];

    hermiticityError =
        Norm[H - ConjugateTranspose[H], "Frobenius"]/hScale;

    If[hermiticityError > hermiticityTolerance,
        Throw[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Hermiticity check failed.",
            "hermiticityError" -> hermiticityError
        |>, "runHxFailure"]
    ];

    crossBlock = ConjugateTranspose[BpN] . H . BmN;
    reflectionBlockError = Norm[crossBlock, "Frobenius"]/hScale;

    If[reflectionBlockError > reflectionTolerance,
        Throw[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Reflection block check failed.",
            "reflectionBlockError" -> reflectionBlockError
        |>, "runHxFailure"]
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
        Throw[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Sampled eigensystem residual check failed.",
            "eigResidualP" -> eigResidualP,
            "eigResidualM" -> eigResidualM
        |>, "runHxFailure"]
    ];


    (* --------------------------------------------------------- *)
    (* Fixed RPS expressed in the two new eigenbases             *)
    (* --------------------------------------------------------- *)

    WP = Developer`ToPackedArray[Re[Transpose[evecsP]]];
    WM = Developer`ToPackedArray[Re[Transpose[evecsM]]];
    meanEntropy = Developer`ToPackedArray[ConstantArray[0., nSamples]];
    ensembleSeconds = 0.;
    Do[
        psi0 = cfg["initialStates"][[member]];
        psi0P = cfg["initialStatesP"][[member]];
        psi0M = cfg["initialStatesM"][[member]];
    coeffP = Conjugate[evecsP] . psi0P;
    coeffM = Conjugate[evecsM] . psi0M;

    maxImagV = Max[
        Max[Abs[Im[evecsP]]],
        Max[Abs[Im[evecsM]]]
    ];

    If[maxImagV > 10^-10,
        Throw[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "The sector eigenbasis is not numerically real.",
            "maxImagV" -> maxImagV
        |>, "runHxFailure"]
    ];


    (* --------------------------------------------------------- *)
    (* Sector-sized reconstruction matrices: no full D x D Vfull *)
    (* --------------------------------------------------------- *)


    reconstructedPsi0 =
        BpN . (WP . coeffP) +
        BmN . (WM . coeffM);

    reconstructionError = Norm[reconstructedPsi0 - psi0];

    If[reconstructionError > reconstructionTolerance,
        Throw[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Initial-state reconstruction check failed.",
            "reconstructionError" -> reconstructionError
        |>, "runHxFailure"]
    ];


    (* --------------------------------------------------------- *)
    (* Evolve this member before updating the running mean       *)
    (* --------------------------------------------------------- *)

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
        Throw[<|"success" -> False, "index" -> k, "hx" -> hx,
            "message" -> "Evolution did not return a result association. Reevaluate all definitions and worker setup.",
            "resultHead" -> Head[evolutionResult]|>, "runHxFailure"]
    ];
    entropy = Lookup[evolutionResult, "entropy", Missing["NotAvailable"]];
    If[!validEntropyVectorQ[entropy, nSamples],
        Throw[<|"success" -> False, "index" -> k, "hx" -> hx,
            "message" -> "Evolution did not return the expected finite real entropy vector.",
            "resultHead" -> Head[entropy], "expectedLength" -> nSamples|>, "runHxFailure"]
    ];
    evolutionNormError = Lookup[evolutionResult, "maxNormError", Infinity];
    If[!TrueQ[evolutionNormError <= validationNormTolerance],
        Throw[<|"success" -> False, "index" -> k, "hx" -> hx,
            "message" -> "Norm drift in the evolved batches exceeded tolerance.",
            "maxEvolutionNormError" -> evolutionNormError|>, "runHxFailure"]
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
        Throw[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Long-time validation check failed.",
            "maxValidationEntropyDifference" -> maxValidationEntropyDifference,
            "maxValidationNormError" -> maxValidationNormError
        |>, "runHxFailure"]
    ];



        (* Online arithmetic mean of ENTROPIES; no member trajectory is exported. *)
        meanEntropy += (entropy - meanEntropy)/member;
        ensembleSeconds += evolutionSeconds;
        worstReconstruction = Max[worstReconstruction, reconstructionError];
        worstEntropyDifference = Max[worstEntropyDifference, maxValidationEntropyDifference];
        worstDirectNorm = Max[worstDirectNorm, maxValidationNormError];
        worstEvolutionNorm = Max[worstEvolutionNorm, evolutionNormError];
        Clear[entropy];
        , {member, stateCount}
    ];
    entropy = Developer`ToPackedArray[meanEntropy];
    Clear[meanEntropy];
    evolutionSeconds = ensembleSeconds;
    reconstructionError = worstReconstruction;
    maxValidationEntropyDifference = worstEntropyDifference;
    maxValidationNormError = worstDirectNorm;
    evolutionNormError = worstEvolutionNorm;

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
        "stateCount" -> stateCount,
        "stateSeeds" -> cfg["stateSeeds"],
        "calculationMode" -> cfg["calculationMode"],
        "observable" -> "Arithmetic mean of pure-state entanglement entropies (nats)",
        "evalsP" -> evalsP,
        "evalsM" -> evalsM,
        "hermiticityError" -> hermiticityError,
        "reflectionBlockError" -> reflectionBlockError,
        "eigResidualP" -> eigResidualP,
        "eigResidualM" -> eigResidualM,
        "reconstructionError" -> reconstructionError,
        "maxImagV" -> maxImagV,
        "validationIndices" -> validationIndices,
        "validationTimes" -> validationTimes,
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
        exportResult = checkedExport[runFile, resultAssociation, "WXF"];
    ]];
    totalSeconds = computeSeconds + exportSeconds;

    If[!StringQ[exportResult] || !FileExistsQ[runFile],
        Throw[<|
            "success" -> False,
            "index" -> k,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "message" -> "Export failed.",
            "runFile" -> runFile
        |>, "runHxFailure"]
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
], "runHxFailure"];



(* Configuration is rebuilt from the current parameter values on EVERY call. *)
ClearAll[requireSweep, buildSweepConfig, initializeSweepWorkers, runEntropySweep];
requireSweep[condition_, message_String] := If[!TrueQ[condition],
    Throw[Failure["InvalidSweep", <|"message" -> message|>], "sweep"]
];
buildSweepConfig[mode_String] := Catch[Module[
    {cfg, count, seeds, states, bases, bp, bm, ts = None, steps, samples,
     dir, ortho, nlog, values},
    requireSweep[MemberQ[{"Single", "Ensemble"}, mode], "Unknown calculation mode."];
    requireSweep[IntegerQ[L] && EvenQ[L] && L > 0, "L must be a positive even integer."];
    requireSweep[VectorQ[hxExponents, NumericQ] && Length[hxExponents] > 0,
        "hxExponents must be a nonempty numerical list."];
    requireSweep[AllTrue[{batchSize, resetEvery, requestedKernels}, IntegerQ[#] && # > 0 &],
        "Batch size, reset interval and requested kernels must be positive integers."];
    count = If[mode === "Single", 1, ensembleSize];
    requireSweep[IntegerQ[count] && count >= 1, "ensembleSize must be a positive integer."];
    seeds = If[mode === "Single", {seed}, ensembleSeed + Range[0, count - 1]];
    requireSweep[VectorQ[seeds, IntegerQ], "State seeds must be integers."];
    Switch[timeGridType,
        "Linear",
        requireSweep[NumericQ[tmin] && NumericQ[tmax] && NumericQ[dt] && dt > 0 && tmax >= tmin,
            "Linear times require tmax >= tmin and dt > 0."];
        steps = Round[(tmax - tmin)/dt];
        requireSweep[Abs[(tmax - tmin)/dt - steps] <= 10^-7,
            "The linear endpoints must define an integer number of steps."];
        samples = steps + 1,
        "Log",
        requireSweep[logTMax > logTMin && IntegerQ[pointsPerDecade] && pointsPerDecade >= 1,
            "Invalid logarithmic grid parameters."];
        nlog = Max[1, Round[(logTMax - logTMin)*pointsPerDecade]];
        ts = N[10.^Subdivide[logTMin, logTMax, nlog]];
        If[TrueQ[includeZeroInLogGrid], ts = Prepend[ts, 0.]];
        ts = Developer`ToPackedArray[ts];
        requireSweep[VectorQ[ts, MachineNumberQ], "Logarithmic times overflowed machine arithmetic."];
        samples = Length[ts]; steps = samples - 1,
        _, requireSweep[False, "timeGridType must be Linear or Log."]
    ];
    (* Establish a writable directory BEFORE any diagonalization or evolution. *)
    requireSweep[StringQ[baseOutputDirectory] && DirectoryQ[baseOutputDirectory],
        "baseOutputDirectory must be an existing directory (normally NotebookDirectory[])."];
    dir = makeRunDirectory[baseOutputDirectory, L, J, hz, N[hxExponents]];
    If[FailureQ[dir], Throw[dir, "sweep"]];
    bases = reflectionSectorBases[L]; bp = N[bases["Even"]]; bm = N[bases["Odd"]];
    ortho = Max[
        Norm[Transpose[bp] . bp - IdentityMatrix[Dimensions[bp][[2]]], "Frobenius"],
        Norm[Transpose[bm] . bm - IdentityMatrix[Dimensions[bm][[2]]], "Frobenius"],
        Norm[Transpose[bp] . bm, "Frobenius"]
    ];
    requireSweep[ortho <= reconstructionTolerance, "Reflection basis validation failed."];
    states = Table[BlockRandom[
        SeedRandom[seeds[[r]]];
        Developer`ToPackedArray[N[RandomChainProductState[L]]]
    ], {r, count}];
    requireSweep[AllTrue[states, VectorQ[#, NumericQ] && Length[#] == 2^L &&
        TrueQ[Abs[Norm[#]^2 - 1.] <= validationNormTolerance] &],
        "The QMB initial states have invalid dimensions or normalization."];
    values = Developer`ToPackedArray[N[10^N[hxExponents]]];
    cfg = <|
        "L" -> L, "J" -> N[J], "hz" -> N[hz],
        "hxExponents" -> N[hxExponents], "hxList" -> values,
        "dA" -> 2^(L/2), "dB" -> 2^(L/2), "BpN" -> bp, "BmN" -> bm,
        "calculationMode" -> mode, "stateCount" -> count, "stateSeeds" -> seeds,
        "initialStates" -> states,
        "initialStatesP" -> (Developer`ToPackedArray[Transpose[bp] . #] & /@ states),
        "initialStatesM" -> (Developer`ToPackedArray[Transpose[bm] . #] & /@ states),
        "timeGridType" -> timeGridType, "tmin" -> N[tmin], "tmax" -> N[tmax], "dt" -> N[dt],
        "nSteps" -> steps, "nSamples" -> samples, "timeGrid" -> ts,
        "logTMin" -> N[logTMin], "logTMax" -> N[logTMax],
        "pointsPerDecade" -> pointsPerDecade, "includeZeroInLogGrid" -> includeZeroInLogGrid,
        "batchSize" -> batchSize, "resetEvery" -> resetEvery,
        "requestedKernels" -> requestedKernels, "useParallel" -> useParallel,
        "outputDir" -> dir, "qmbInitPath" -> qmbInitPath,
        "hermiticityTolerance" -> hermiticityTolerance, "reflectionTolerance" -> reflectionTolerance,
        "reconstructionTolerance" -> reconstructionTolerance,
        "eigensystemResidualTolerance" -> eigensystemResidualTolerance,
        "validationEntropyTolerance" -> validationEntropyTolerance,
        "validationNormTolerance" -> validationNormTolerance,
        "basisOrthogonalityError" -> ortho, "WolframVersion" -> $Version
    |>;
    cfg
], "sweep"];

(* Reuse live kernels; all numerical inputs travel as a fresh literal snapshot. *)
initializeSweepWorkers[cfg_Association] := Module[{missing, status, visibility},
    If[!TrueQ[cfg["useParallel"]], Return[True]];
    missing = Max[0, cfg["requestedKernels"] - Length[Kernels[]]];
    If[missing > 0, LaunchKernels[missing]];
    If[Length[Kernels[]] == 0, Return[Failure["NoKernels", <||>]]];
    status = With[{path = cfg["qmbInitPath"]},
        ParallelEvaluate[Quiet[Check[Get[path]; True, False]]]
    ];
    If[!AllTrue[status, TrueQ], Return[Failure["WorkerPackageLoad", <|"status" -> status|>]]];
    DistributeDefinitions[
        runOneHx, checkedExport, validEntropyVectorQ, sortedEigensystem,
        phaseStepC, phaseBatchC, phaseAtTimesBatchC, reseedC, entropyFromSVC,
        entropyFast, entropyBatchFromSectors, linearEntropyEvolutionSector,
        arbitraryEntropyEvolutionSector, directEntropyAndNormSector, sampledEigensystemResidual
    ];
    (* Verify each local worker can write here before expensive work begins. *)
    visibility = With[{directory = cfg["outputDir"]}, ParallelEvaluate[
        Module[{p, result},
            p = FileNameJoin[{directory, "worker_" <> ToString[$KernelID] <> ".tmp"}];
            result = checkedExport[p, "write test", "Text"];
            If[!FailureQ[result], DeleteFile[p]];
            !FailureQ[result]
        ]
    ]];
    If[!AllTrue[visibility, TrueQ], Return[Failure["WorkerCannotWrite", <|"status" -> visibility|>]]];
    Print["Using ", Length[Kernels[]], " existing/available kernels; no kernels are closed."];
    True
];

ClearAll[writeSweepDescription];
writeSweepDescription[cfg_Association] := Module[{gridText, description, metadata},
    gridText = If[cfg["timeGridType"] === "Linear",
        "Linear grid: t(i) = tmin + i dt, i = 0,...,nSteps (both endpoints).\n" <>
        "tmin = " <> ToString[cfg["tmin"], InputForm] <> "; tmax = " <>
        ToString[cfg["tmax"], InputForm] <> "; dt = " <> ToString[cfg["dt"], InputForm],
        "Log grid: positive times = 10^Subdivide[logTMin,logTMax,nIntervals].\n" <>
        "logTMin = " <> ToString[cfg["logTMin"], InputForm] <> "; logTMax = " <>
        ToString[cfg["logTMax"], InputForm] <> "; pointsPerDecade = " <>
        ToString[cfg["pointsPerDecade"]] <> "; includeZero = " <>
        ToString[cfg["includeZeroInLogGrid"]] <>
        ".\nTime increments are nonuniform; the exact timeGrid is stored in every result WXF."
    ];
    description = StringRiffle[{
        "Prethermal entropy sweep", "Mode: " <> cfg["calculationMode"],
        "L = " <> ToString[cfg["L"]] <> "; J = " <> ToString[cfg["J"], InputForm] <>
            "; hz = " <> ToString[cfg["hz"], InputForm],
        "hxExponents = " <> ToString[cfg["hxExponents"], InputForm],
        "hx values = " <> ToString[cfg["hxList"], InputForm],
        "Initial states: QMB RandomChainProductState[L], one product of local spin states per seed.",
        "State count = " <> ToString[cfg["stateCount"]],
        "Seeds = " <> ToString[cfg["stateSeeds"], InputForm],
        "Each state is generated with BlockRandom[SeedRandom[stateSeed]; RandomChainProductState[L]].",
        "Exactly the same initial state list is used for every hx in this sweep.",
        "Initial vectors are retained in metadata.wxf for exact reuse; no individual entropy trajectories are saved.",
        gridText, "nSteps = " <> ToString[cfg["nSteps"]] <>
            "; nSamples = " <> ToString[cfg["nSamples"]],
        "Batch size = " <> ToString[cfg["batchSize"]] <>
            "; recurrence reseed interval = " <> ToString[cfg["resetEvery"]],
        "Observable: arithmetic mean of individual pure-state half-chain entanglement entropies, in nats.",
        "The ensemble average is taken before any plotting or moving average.",
        "The saved entropy is unsmoothed. No member trajectories or density-matrix averages are exported.",
        "Wolfram version: " <> $Version,
        "QMB init path: " <> cfg["qmbInitPath"]
    }, "\n\n"];
    checkedExport[FileNameJoin[{cfg["outputDir"], "description.txt"}], description, "Text"]
];

runEntropySweep[mode_String] := Module[{cfg, startup, summaries, seconds, failed, result, metadata},
    cfg = buildSweepConfig[mode];
    If[FailureQ[cfg], Print[cfg]; Return[cfg]];
    outputDir = cfg["outputDir"]; (* convenience for the existing plotting workflow *)
    Print["Saving to: ", outputDir];
    result = writeSweepDescription[cfg];
    If[FailureQ[result], Print[result]; Return[result]];
    metadata = KeyDrop[cfg, {"BpN", "BmN", "initialStatesP", "initialStatesM"}];
    result = checkedExport[FileNameJoin[{outputDir, "metadata.wxf"}], metadata, "WXF"];
    If[FailureQ[result], Print[result]; Return[result]];
    startup = initializeSweepWorkers[cfg];
    If[FailureQ[startup], Print[startup]; Return[startup]];
    seconds = First[AbsoluteTiming[
        summaries = If[TrueQ[cfg["useParallel"]],
            With[{snapshot = cfg}, ParallelMap[
                runOneHx[snapshot, #] &, Range[Length[snapshot["hxList"]]],
                Method -> "CoarsestGrained", DistributedContexts -> None
            ]],
            Map[runOneHx[cfg, #] &, Range[Length[cfg["hxList"]]]]
        ];
    ]];
    failed = If[ListQ[summaries],
        Select[summaries, !AssociationQ[#] || !TrueQ[Lookup[#, "success", False]] &],
        {summaries}];
    result = checkedExport[FileNameJoin[{outputDir, "run_summary.wxf"}], summaries, "WXF"];
    If[FailureQ[result], Print[result]; Return[result]];
    If[failed =!= {},
        checkedExport[FileNameJoin[{outputDir, "failed_runs.wxf"}], failed, "WXF"];
        Print["Some fields failed; their errors are saved. No incomplete ensemble is marked successful."];
        Return[Failure["SweepFailed", <|"outputDir" -> outputDir, "failedRuns" -> failed|>]]
    ];
    Print["Completed ", Length[summaries], " fields, ", cfg["stateCount"],
        " RPS per field; sweep wall time = ", seconds, " s."];
    <|"success" -> True, "outputDir" -> outputDir, "runSummary" -> summaries,
        "sweepSeconds" -> seconds, "stateCount" -> cfg["stateCount"]|>
];



(* ::Title::Closed:: *)
(* Replot saved results *)


ClearAll[startupMovingAverage, replotEntropySweep];
startupMovingAverage[values_List, window_Integer] := Module[{n = Length[values], prefix},
    If[n == 0 || window < 1, Return[{}]];
    If[n < window, Return[Accumulate[values]/Range[n]]];
    prefix = If[window == 1, {},
        Accumulate[Take[values, window - 1]]/Range[window - 1]];
    Developer`ToPackedArray[Join[prefix, MovingAverage[values, window]]]
];
replotEntropySweep[directory_String, window_Integer : 200] := Module[
    {files, runs, data, times, exps, colors, cf, lo, hi, legend, title,
     plotDirectory, paths = <||>, positive, averaged, averagedPositive,
     graphic, result, makePlot, samples, prefix, meta},
    If[window < 1, Return[Failure["InvalidAverageWindow", <||>]]];
    files = Sort[FileNames["run_" ~~ DigitCharacter .. ~~ ".wxf", directory]];
    If[files === {}, Return[Failure["NoSavedTrajectories", <|"directory" -> directory|>]]];
    runs = Import[#, "WXF"] & /@ files;
    If[!AllTrue[runs, AssociationQ[#] && KeyExistsQ[#, "entropy"] &],
        Return[Failure["InvalidResultFiles", <||>]]];
    runs = SortBy[runs, #["hxExponent"] &];
    meta = First[runs];
    data = Table[
        samples = Length[runs[[k]]["entropy"]];
        times = If[runs[[k]]["timeGridType"] === "Linear",
            N[runs[[k]]["tmin"] + Range[0, samples - 1]*runs[[k]]["dt"]],
            runs[[k]]["timeGrid"]];
        Transpose[{times, runs[[k]]["entropy"]}],
        {k, Length[runs]}
    ];
    exps = Lookup[runs, "hxExponent"];
    lo = Min[exps]; hi = Max[exps];
    cf = Blend[{Blue, Purple, Red}, #] &;
    colors = If[lo == hi, {cf[0.5]}, cf /@ Rescale[exps, {lo, hi}]];
    (* Continuous hx colorbar, with color position uniform in log10[hx]. *)
    legend = Module[{barRange, barTicks, barColors, tickExponents},
        barRange = If[lo == hi, {lo - 0.5, hi + 0.5}, {lo, hi}];
        barColors = If[lo == hi, Function[z, Blend[{Blue, Purple, Red}, 0.5]], cf];
        tickExponents = If[lo == hi, {lo}, Range[Ceiling[lo], Floor[hi]]];
        If[tickExponents === {}, tickExponents = {lo, hi}];
        barTicks = ({#, Style[Superscript["10", #], 24, Black]} & /@ tickExponents);
        BarLegend[
            {barColors, barRange},
            LegendLabel -> Style[Subscript["h", "x"], 30, Black],
            LabelStyle -> {Black, 20},
            ColorFunctionScaling -> True,
            LegendMarkerSize -> {12, 360},
            Ticks -> barTicks
        ]
    ];
    prefix = If[Lookup[meta, "stateCount", 1] == 1, "RPS", 
        "RPS ensemble mean (" <> ToString[meta["stateCount"]] <> " states)"];
    title = Row[{prefix, ", L = ", meta["L"], ", J = ", meta["J"], ", hz = ", meta["hz"]}];
    plotDirectory = FileNameJoin[{directory, "plots"}];
    If[!DirectoryQ[plotDirectory], Quiet[Check[CreateDirectory[plotDirectory], Null]]];
    If[!DirectoryQ[plotDirectory], Return[Failure["PlotDirectoryFailed", <||>]]];
    (* All requested samples are plotted; there is no automatic point-count skip. *)
    makePlot[coordinates_, logarithmic_, caption_] := Legended[
        If[TrueQ[logarithmic], ListLogLinearPlot, ListPlot][
            coordinates, Joined -> False,
            PlotStyle -> (Directive[#, Opacity[0.8], PointSize[0.002]] & /@ colors),
            PlotRange -> All, Frame -> True, Axes -> False, ImageSize -> 1000,
            FrameStyle -> Directive[Black, 22],
            PlotLabel -> Style[Column[{title, caption}, Alignment -> Center], 22, Black],
            FrameLabel -> {Style["t", 26, Black], Style["S(t)", 26, Black]},
            GridLines -> None
        ], Placed[legend, Right]
    ];
    graphic = makePlot[data, False, "Unsmooth entropy"];
    result = checkedExport[FileNameJoin[{plotDirectory, "entropy_hx_raw_linear.png"}], graphic, "PNG"];
    AssociateTo[paths, "rawLinear" -> result];
    positive = (Select[#, First[#] > 0. &] & /@ data);
    If[AnyTrue[positive, Length[#] > 0 &],
        graphic = makePlot[positive, True, "Unsmooth entropy"];
        result = checkedExport[FileNameJoin[{plotDirectory, "entropy_hx_raw_logtime.png"}], graphic, "PNG"];
        AssociateTo[paths, "rawLog" -> result]
    ];
    (* Average the SAVED ensemble mean only after the ensemble is complete. *)
    averaged = (Transpose[{#[[All, 1]], startupMovingAverage[#[[All, 2]], window]}] & /@ data);
    graphic = makePlot[averaged, False, "Moving average, window = " <> ToString[window]];
    result = checkedExport[FileNameJoin[{plotDirectory, "entropy_hx_moving_average_linear.png"}], graphic, "PNG"];
    AssociateTo[paths, "movingAverageLinear" -> result];
    averagedPositive = (Select[#, First[#] > 0. &] & /@ averaged);
    If[AnyTrue[averagedPositive, Length[#] > 0 &],
        graphic = makePlot[averagedPositive, True, "Moving average, window = " <> ToString[window]];
        result = checkedExport[FileNameJoin[{plotDirectory, "entropy_hx_moving_average_logtime.png"}], graphic, "PNG"];
        AssociateTo[paths, "movingAverageLog" -> result]
    ];
    If[AnyTrue[Values[paths], FailureQ], Print["One or more plot exports failed: ", paths],
        Print["Plots saved in: ", plotDirectory]];
    paths
];



(* ::Title::Closed:: *)
(* Parameters*)


L = 8;
J = 1.;
hz = 1.;
hxExponents = N[Range[-2, 0, 2/10]];

(* Linear or logarithmic grid; times are NEVER used in folder names. *)
timeGridType = "Linear";
tmin = 0.;
tmax = 10^5;
dt = 1.;
logTMin = -1.;
logTMax = 10.;
pointsPerDecade = 2000;
includeZeroInLogGrid = True;
resetEvery = 100000;
batchSize = Which[L <= 10, 64, L <= 12, 32, True, 16];

baseOutputDirectory = NotebookDirectory[];
useParallel = True;
requestedKernels = 8;
seed = 12345;

hermiticityTolerance = 10^-12;
reflectionTolerance = 10^-12;
reconstructionTolerance = 10^-10;
eigensystemResidualTolerance = 10^-10;
validationEntropyTolerance = 10^-8;
validationNormTolerance = 10^-10;

movingAverageWindow = 200;
exportPlots = True;


(* ::Title::Closed:: *)
(* Single RPS sweep*)


If[calculationMode === "Single",
    singleResult = runEntropySweep["Single"];
    If[AssociationQ[singleResult] && TrueQ[singleResult["success"]],
        singleOutputDir = singleResult["outputDir"];
        If[TrueQ[exportPlots], singlePlotFiles = replotEntropySweep[singleOutputDir, movingAverageWindow]]
    ];
];


singlePlotFiles = replotEntropySweep[singleOutputDir, movingAverageWindow];


(* ::Title::Closed:: *)
(*Ensemble of RPS*)


(* Choose which complete section runs when evaluating the entire file. *)
calculationMode = "Single"; (* change to "Ensemble" for the new RPS ensemble section *)
ensembleSize = 50;
ensembleSeed = 12345; (* member r uses ensembleSeed + r - 1 *)


ensembleResult = runEntropySweep["Ensemble"];

If[
    AssociationQ[ensembleResult] &&
    TrueQ[ensembleResult["success"]],

    ensembleOutputDir = ensembleResult["outputDir"];

    If[TrueQ[exportPlots],
        ensemblePlotFiles = replotEntropySweep[
            ensembleOutputDir,
            movingAverageWindow
        ];
    ];
];


L = 8;
J = 1.;
hz = 2.;
hxExponents = N[Range[-2, 0, 2/10]];

(* Linear or logarithmic grid; times are NEVER used in folder names. *)
timeGridType = "Linear";
tmin = 0.;
tmax = 10^10;
dt = 1.;
logTMin = -1.;
logTMax = 10.;
pointsPerDecade = 2000;
includeZeroInLogGrid = True;
resetEvery = 100000;
batchSize = Which[L <= 10, 64, L <= 12, 32, True, 16];

baseOutputDirectory = NotebookDirectory[];
useParallel = True;
requestedKernels = 8;
seed = 12345;

hermiticityTolerance = 10^-12;
reflectionTolerance = 10^-12;
reconstructionTolerance = 10^-10;
eigensystemResidualTolerance = 10^-10;
validationEntropyTolerance = 10^-8;
validationNormTolerance = 10^-10;

movingAverageWindow = 200;
exportPlots = True;


ensembleResult = runEntropySweep["Ensemble"];

If[
    AssociationQ[ensembleResult] &&
    TrueQ[ensembleResult["success"]],

    ensembleOutputDir = ensembleResult["outputDir"];

    If[TrueQ[exportPlots],
        ensemblePlotFiles = replotEntropySweep[
            ensembleOutputDir,
            movingAverageWindow
        ];
    ];
];


(* ::Title::Closed:: *)
(* Replot later from saved averages only *)


(* Set savedRunDirectory to a previous run folder, then evaluate:
   replotEntropySweep[savedRunDirectory, movingAverageWindow];
   You may also use ensembleOutputDir or singleOutputDir in the same session.
   No diagonalization, new random states, or time evolution is performed.
   Kernel shutdown is intentionally absent from this file. *)


(* ::Title::Closed:: *)
(*Only log sampling*)


(* Physical parameters *)
logOnlyL = 8;
logOnlyJ = 1.;
logOnlyHz = 2.;

(* hx = 10^alpha *)
logOnlyHxExponents = N[Range[-4, 0, 5/10]];
logOnlyHxList = Developer`ToPackedArray[
    N[10^logOnlyHxExponents]
];


(* ------------------------------------------------------------ *)
(* Logarithmic time grid                                        *)
(* ------------------------------------------------------------ *)

(* The grid is uniform in log10[t].

   With:
       logOnlyLogTMin = -1
       logOnlyLogTMax = 10
       logOnlyPointsPerDecade = 1000

   the code produces exactly

       (10 - (-1))*1000 = 11000

   positive samples, including both endpoints 10^-1 and 10^10.

   The parameter logOnlyPointsPerDecade is free to change.
*)

logOnlyLogTMin = -1.;
logOnlyLogTMax = 10.;
logOnlyPointsPerDecade = 2000;

logOnlyDecades = logOnlyLogTMax - logOnlyLogTMin;

logOnlyNSamples = Round[
    logOnlyDecades*logOnlyPointsPerDecade
];

If[
    logOnlyNSamples < 2,
    Print["ERROR: the logarithmic grid must contain at least two points."];
    Abort[];
];

logOnlyTimeGrid = Developer`ToPackedArray[
    N[
        10.^Subdivide[
            logOnlyLogTMin,
            logOnlyLogTMax,
            logOnlyNSamples - 1
        ]
    ]
];

If[
    Length[logOnlyTimeGrid] =!= logOnlyNSamples,
    Print["ERROR: unexpected logarithmic time-grid length."];
    Abort[];
];

Print[
    "Log-only grid: ",
    Length[logOnlyTimeGrid],
    " points from ",
    First[logOnlyTimeGrid],
    " to ",
    Last[logOnlyTimeGrid],
    "."
];


(* ------------------------------------------------------------ *)
(* Initial states                                                *)
(* ------------------------------------------------------------ *)

(* The first ensemble member is also the single RPS.
   Therefore the single-RPS curve and ensemble mean are obtained
   in the SAME sweep and from the SAME diagonalization for each hx. *)

logOnlyEnsembleSize = 50;
logOnlySeed = 12345;

logOnlyStateSeeds =
    logOnlySeed + Range[0, logOnlyEnsembleSize - 1];


(* ------------------------------------------------------------ *)
(* Numerical/performance parameters                              *)
(* ------------------------------------------------------------ *)

logOnlyBatchSize = Which[
    logOnlyL <= 10, 64,
    logOnlyL <= 12, 32,
    True, 16
];

logOnlyUseParallel = True;
logOnlyRequestedKernels = 8;

logOnlyHermiticityTolerance = 10^-12;
logOnlyReflectionTolerance = 10^-12;
logOnlyReconstructionTolerance = 10^-10;
logOnlyEigensystemResidualTolerance = 10^-10;
logOnlyValidationEntropyTolerance = 10^-8;
logOnlyValidationNormTolerance = 10^-10;

logOnlyMovingAverageWindow = 200;

logOnlyBaseOutputDirectory = NotebookDirectory[];


(* ::Title::Closed:: *)
(* Log-only setup *)


If[
    !IntegerQ[logOnlyL] ||
    !EvenQ[logOnlyL] ||
    logOnlyL <= 0,

    Print["ERROR: logOnlyL must be a positive even integer."];
    Abort[];
];

If[
    !IntegerQ[logOnlyEnsembleSize] ||
    logOnlyEnsembleSize < 1,

    Print["ERROR: logOnlyEnsembleSize must be a positive integer."];
    Abort[];
];

If[
    !IntegerQ[logOnlyRequestedKernels] ||
    logOnlyRequestedKernels < 1,

    Print["ERROR: logOnlyRequestedKernels must be a positive integer."];
    Abort[];
];

If[
    !IntegerQ[logOnlyBatchSize] ||
    logOnlyBatchSize < 1,

    Print["ERROR: logOnlyBatchSize must be a positive integer."];
    Abort[];
];


logOnlyDA = 2^(logOnlyL/2);
logOnlyDB = 2^(logOnlyL/2);


(* Reflection basis *)

logOnlyBases =
    reflectionSectorBases[logOnlyL];

logOnlyBpN =
    N[logOnlyBases["Even"]];

logOnlyBmN =
    N[logOnlyBases["Odd"]];


logOnlyBasisOrthogonalityError = Max[
    Norm[
        Transpose[logOnlyBpN] . logOnlyBpN -
        IdentityMatrix[
            Dimensions[logOnlyBpN][[2]]
        ],
        "Frobenius"
    ],

    Norm[
        Transpose[logOnlyBmN] . logOnlyBmN -
        IdentityMatrix[
            Dimensions[logOnlyBmN][[2]]
        ],
        "Frobenius"
    ],

    Norm[
        Transpose[logOnlyBpN] . logOnlyBmN,
        "Frobenius"
    ]
];


If[
    logOnlyBasisOrthogonalityError >
    logOnlyReconstructionTolerance,

    Print[
        "ERROR: reflection-basis validation failed: ",
        logOnlyBasisOrthogonalityError
    ];
    Abort[];
];


(* Generate the full ensemble once.
   State 1 is the single RPS used in the single-state plot. *)

logOnlyInitialStates = Table[

    BlockRandom[
        SeedRandom[
            logOnlyStateSeeds[[r]]
        ];

        Developer`ToPackedArray[
            N[
                RandomChainProductState[
                    logOnlyL
                ]
            ]
        ]
    ],

    {r, logOnlyEnsembleSize}
];


If[
    !AllTrue[
        logOnlyInitialStates,

        VectorQ[#, NumericQ] &&
        Length[#] == 2^logOnlyL &&
        TrueQ[
            Abs[
                Norm[#]^2 - 1.
            ] <= logOnlyValidationNormTolerance
        ] &
    ],

    Print[
        "ERROR: one or more initial RPS vectors have invalid ",
        "dimension or normalization."
    ];
    Abort[];
];


logOnlyInitialStatesP =
    (
        Developer`ToPackedArray[
            Transpose[logOnlyBpN] . #
        ] &
    ) /@ logOnlyInitialStates;


logOnlyInitialStatesM =
    (
        Developer`ToPackedArray[
            Transpose[logOnlyBmN] . #
        ] &
    ) /@ logOnlyInitialStates;



(* ::Title::Closed:: *)
(* Log-only output directory *)


ClearAll[logOnlyMakeRunDirectory];

logOnlyMakeRunDirectory[] := Module[
    {
        root,
        tag,
        configDir,
        runDir,
        probe,
        result,
        attempt = 0
    },

    root = FileNameJoin[
        {
            logOnlyBaseOutputDirectory,
            "entropy_log_only_sweeps"
        }
    ];

    If[
        !DirectoryQ[root],
        Quiet[
            Check[
                CreateDirectory[
                    root,
                    CreateIntermediateDirectories -> True
                ],
                Null
            ]
        ]
    ];

    If[
        !DirectoryQ[root],
        Return[
            Failure[
                "CreateDirectoryFailed",
                <|"path" -> root|>
            ]
        ]
    ];


    tag = StringRiffle[
        {
            "L" <> ToString[logOnlyL],
            "J" <> tagNumber[logOnlyJ],
            "hz" <> tagNumber[logOnlyHz],
            "hx" <>
                tagNumber[
                    Min[logOnlyHxExponents]
                ] <>
                "to" <>
                tagNumber[
                    Max[logOnlyHxExponents]
                ],
            "Nens" <> ToString[logOnlyEnsembleSize]
        },
        "_"
    ];


    configDir = FileNameJoin[
        {
            root,
            tag
        }
    ];

    If[
        !DirectoryQ[configDir],
        Quiet[
            Check[
                CreateDirectory[
                    configDir,
                    CreateIntermediateDirectories -> True
                ],
                Null
            ]
        ]
    ];

    If[
        !DirectoryQ[configDir],
        Return[
            Failure[
                "CreateDirectoryFailed",
                <|"path" -> configDir|>
            ]
        ]
    ];


    While[
        attempt < 10,

        attempt++;

        runDir = FileNameJoin[
            {
                configDir,
                "run_" <>
                StringTake[
                    CreateUUID[],
                    12
                ]
            }
        ];

        result = Quiet[
            Check[
                CreateDirectory[runDir],
                $Failed
            ]
        ];

        If[
            StringQ[result] &&
            DirectoryQ[runDir],
            Break[]
        ];
    ];


    If[
        !DirectoryQ[runDir],
        Return[
            Failure[
                "CreateRunFailed",
                <|"path" -> configDir|>
            ]
        ]
    ];


    probe = FileNameJoin[
        {
            runDir,
            "write_probe.tmp"
        }
    ];

    result = checkedExport[
        probe,
        "write test",
        "Text"
    ];

    If[
        FailureQ[result],
        Return[result]
    ];

    DeleteFile[probe];

    runDir
];


logOnlyOutputDir =
    logOnlyMakeRunDirectory[];


If[
    FailureQ[logOnlyOutputDir],
    Print[logOnlyOutputDir];
    Abort[];
];


Print[
    "Log-only output directory: ",
    logOnlyOutputDir
];



(* ::Title::Closed:: *)
(* Log-only configuration *)


logOnlyConfig = <|

    "L" -> logOnlyL,
    "J" -> N[logOnlyJ],
    "hz" -> N[logOnlyHz],

    "hxExponents" ->
        Developer`ToPackedArray[
            N[logOnlyHxExponents]
        ],

    "hxList" ->
        logOnlyHxList,

    "dA" -> logOnlyDA,
    "dB" -> logOnlyDB,

    "BpN" -> logOnlyBpN,
    "BmN" -> logOnlyBmN,

    "timeGrid" ->
        logOnlyTimeGrid,

    "nSamples" ->
        Length[logOnlyTimeGrid],

    "logTMin" ->
        N[logOnlyLogTMin],

    "logTMax" ->
        N[logOnlyLogTMax],

    "pointsPerDecade" ->
        logOnlyPointsPerDecade,

    "batchSize" ->
        logOnlyBatchSize,

    "stateCount" ->
        logOnlyEnsembleSize,

    "stateSeeds" ->
        logOnlyStateSeeds,

    "initialStates" ->
        logOnlyInitialStates,

    "initialStatesP" ->
        logOnlyInitialStatesP,

    "initialStatesM" ->
        logOnlyInitialStatesM,

    "outputDir" ->
        logOnlyOutputDir,

    "hermiticityTolerance" ->
        logOnlyHermiticityTolerance,

    "reflectionTolerance" ->
        logOnlyReflectionTolerance,

    "reconstructionTolerance" ->
        logOnlyReconstructionTolerance,

    "eigensystemResidualTolerance" ->
        logOnlyEigensystemResidualTolerance,

    "validationEntropyTolerance" ->
        logOnlyValidationEntropyTolerance,

    "validationNormTolerance" ->
        logOnlyValidationNormTolerance

|>;


logOnlyMetadata = <|

    "calculation" ->
        "Single RPS and RPS ensemble mean in one logarithmic-time sweep",

    "singleRPSIsEnsembleMember" ->
        1,

    "L" ->
        logOnlyL,

    "J" ->
        N[logOnlyJ],

    "hz" ->
        N[logOnlyHz],

    "hxExponents" ->
        logOnlyHxExponents,

    "hxList" ->
        logOnlyHxList,

    "logTMin" ->
        logOnlyLogTMin,

    "logTMax" ->
        logOnlyLogTMax,

    "pointsPerDecadeParameter" ->
        logOnlyPointsPerDecade,

    "numberOfTimeSamples" ->
        Length[logOnlyTimeGrid],

    "timeGrid" ->
        logOnlyTimeGrid,

    "ensembleSize" ->
        logOnlyEnsembleSize,

    "stateSeeds" ->
        logOnlyStateSeeds,

    "initialStates" ->
        logOnlyInitialStates,

    "batchSize" ->
        logOnlyBatchSize,

    "requestedKernels" ->
        logOnlyRequestedKernels,

    "basisOrthogonalityError" ->
        logOnlyBasisOrthogonalityError,

    "WolframVersion" ->
        $Version,

    "SystemID" ->
        $SystemID

|>;


logOnlyMetadataPath = FileNameJoin[
    {
        logOnlyOutputDir,
        "metadata.wxf"
    }
];


logOnlyMetadataExport = checkedExport[
    logOnlyMetadataPath,
    logOnlyMetadata,
    "WXF"
];


If[
    FailureQ[logOnlyMetadataExport],
    Print[logOnlyMetadataExport];
    Abort[];
];



(* ::Title::Closed:: *)
(* One hx: single RPS + ensemble mean *)


ClearAll[logOnlyRunOneHx];

logOnlyRunOneHx[
    cfg_Association,
    k_Integer
] := Catch[
    Module[
        {
            L = cfg["L"],
            J = cfg["J"],
            hz = cfg["hz"],
            hxExponents = cfg["hxExponents"],
            hxList = cfg["hxList"],
            dA = cfg["dA"],
            dB = cfg["dB"],
            BpN = cfg["BpN"],
            BmN = cfg["BmN"],
            times = cfg["timeGrid"],
            nSamples = cfg["nSamples"],
            batchSize = cfg["batchSize"],
            stateCount = cfg["stateCount"],
            initialStates = cfg["initialStates"],
            initialStatesP = cfg["initialStatesP"],
            initialStatesM = cfg["initialStatesM"],
            outputDir = cfg["outputDir"],

            hxExponent,
            hx,

            H,
            hScale,
            hermiticityError,
            crossBlock,
            reflectionBlockError,

            Hp,
            Hm,

            diagSeconds,
            eigensystems,
            evalsP,
            evecsP,
            evalsM,
            evecsM,
            eigResidualP,
            eigResidualM,

            WP,
            WM,
            maxImagV,

            member,
            psi0,
            psi0P,
            psi0M,
            coeffP,
            coeffM,
            reconstructedPsi0,
            reconstructionError,

            evolutionResult,
            entropyMember,
            evolutionNormError,

            validationIndices,
            validationTimes,
            validationResults,
            validationEntropyDifferences,
            validationNormErrors,

            maxValidationEntropyDifference = 0.,
            maxValidationNormError = 0.,
            maxEvolutionNormError = 0.,
            maxReconstructionError = 0.,

            singleEntropy,
            ensembleMean,

            runStart,
            diagWall,
            evolutionWall = 0.,
            computeSeconds,

            runFile,
            resultAssociation,
            exportSeconds,
            exportResult,
            totalSeconds
        },


        hxExponent =
            hxExponents[[k]];

        hx =
            hxList[[k]];

        runStart =
            AbsoluteTime[];


        (* ----------------------------------------------------- *)
        (* Hamiltonian and symmetry checks                       *)
        (* ----------------------------------------------------- *)

        H = N[
            IsingHamiltonian[
                hx,
                hz,
                J,
                L
            ]
        ];


        hScale = Max[
            1.,
            Norm[
                H,
                "Frobenius"
            ]
        ];


        hermiticityError =
            Norm[
                H -
                ConjugateTranspose[H],
                "Frobenius"
            ]/
            hScale;


        If[
            hermiticityError >
            cfg["hermiticityTolerance"],

            Throw[
                <|
                    "success" -> False,
                    "index" -> k,
                    "hx" -> hx,
                    "message" ->
                        "Hermiticity check failed.",
                    "hermiticityError" ->
                        hermiticityError
                |>,
                "logOnlyFailure"
            ]
        ];


        crossBlock =
            ConjugateTranspose[BpN] .
            H .
            BmN;


        reflectionBlockError =
            Norm[
                crossBlock,
                "Frobenius"
            ]/
            hScale;


        If[
            reflectionBlockError >
            cfg["reflectionTolerance"],

            Throw[
                <|
                    "success" -> False,
                    "index" -> k,
                    "hx" -> hx,
                    "message" ->
                        "Reflection block check failed.",
                    "reflectionBlockError" ->
                        reflectionBlockError
                |>,
                "logOnlyFailure"
            ]
        ];


        Hp =
            ConjugateTranspose[BpN] .
            H .
            BpN;


        Hm =
            ConjugateTranspose[BmN] .
            H .
            BmN;


        (* ----------------------------------------------------- *)
        (* Diagonalize once for this hx                          *)
        (* ----------------------------------------------------- *)

        diagWall = First[
            AbsoluteTiming[

                eigensystems = {
                    sortedEigensystem[Hp],
                    sortedEigensystem[Hm]
                };

            ]
        ];


        {
            {evalsP, evecsP},
            {evalsM, evecsM}
        } = eigensystems;


        eigResidualP =
            sampledEigensystemResidual[
                Hp,
                evalsP,
                evecsP
            ];


        eigResidualM =
            sampledEigensystemResidual[
                Hm,
                evalsM,
                evecsM
            ];


        If[
            Max[
                eigResidualP,
                eigResidualM
            ] >
            cfg[
                "eigensystemResidualTolerance"
            ],

            Throw[
                <|
                    "success" -> False,
                    "index" -> k,
                    "hx" -> hx,
                    "message" ->
                        "Sampled eigensystem residual check failed.",
                    "eigResidualP" ->
                        eigResidualP,
                    "eigResidualM" ->
                        eigResidualM
                |>,
                "logOnlyFailure"
            ]
        ];


        maxImagV = Max[
            Max[
                Abs[
                    Im[evecsP]
                ]
            ],
            Max[
                Abs[
                    Im[evecsM]
                ]
            ]
        ];


        If[
            maxImagV >
            10^-10,

            Throw[
                <|
                    "success" -> False,
                    "index" -> k,
                    "hx" -> hx,
                    "message" ->
                        "Sector eigenbasis is not numerically real.",
                    "maxImagV" ->
                        maxImagV
                |>,
                "logOnlyFailure"
            ]
        ];


        WP =
            Developer`ToPackedArray[
                Re[
                    Transpose[
                        evecsP
                    ]
                ]
            ];


        WM =
            Developer`ToPackedArray[
                Re[
                    Transpose[
                        evecsM
                    ]
                ]
            ];


        (* ----------------------------------------------------- *)
        (* Prepare output arrays                                 *)
        (* ----------------------------------------------------- *)

        ensembleMean =
            Developer`ToPackedArray[
                ConstantArray[
                    0.,
                    nSamples
                ]
            ];


        validationIndices =
            DeleteDuplicates[
                {
                    1,
                    1 + Floor[
                        (nSamples - 1)/3
                    ],
                    1 + Floor[
                        2 (nSamples - 1)/3
                    ],
                    nSamples
                }
            ];


        validationTimes =
            N[
                times[[
                    validationIndices
                ]]
            ];


        (* ----------------------------------------------------- *)
        (* Ensemble members                                      *)
        (* ----------------------------------------------------- *)

        Do[

            psi0 =
                initialStates[[member]];

            psi0P =
                initialStatesP[[member]];

            psi0M =
                initialStatesM[[member]];


            coeffP =
                Conjugate[evecsP] .
                psi0P;


            coeffM =
                Conjugate[evecsM] .
                psi0M;


            reconstructedPsi0 =
                BpN .
                (WP . coeffP) +
                BmN .
                (WM . coeffM);


            reconstructionError =
                Norm[
                    reconstructedPsi0 -
                    psi0
                ];


            If[
                reconstructionError >
                cfg[
                    "reconstructionTolerance"
                ],

                Throw[
                    <|
                        "success" -> False,
                        "index" -> k,
                        "hx" -> hx,
                        "member" -> member,
                        "message" ->
                            "Initial-state reconstruction check failed.",
                        "reconstructionError" ->
                            reconstructionError
                    |>,
                    "logOnlyFailure"
                ]
            ];


            maxReconstructionError =
                Max[
                    maxReconstructionError,
                    reconstructionError
                ];


            (* Every logarithmic sample is generated independently
               from the original spectral coefficients. *)

            evolutionWall += First[
                AbsoluteTiming[

                    evolutionResult =
                        arbitraryEntropyEvolutionSector[
                            WP,
                            WM,
                            BpN,
                            BmN,
                            evalsP,
                            evalsM,
                            coeffP,
                            coeffM,
                            times,
                            batchSize,
                            dA,
                            dB
                        ];

                ]
            ];


            If[
                !AssociationQ[
                    evolutionResult
                ],

                Throw[
                    <|
                        "success" -> False,
                        "index" -> k,
                        "hx" -> hx,
                        "member" -> member,
                        "message" ->
                            "Logarithmic evolution did not return an association."
                    |>,
                    "logOnlyFailure"
                ]
            ];


            entropyMember =
                Lookup[
                    evolutionResult,
                    "entropy",
                    Missing[
                        "NotAvailable"
                    ]
                ];


            If[
                !validEntropyVectorQ[
                    entropyMember,
                    nSamples
                ],

                Throw[
                    <|
                        "success" -> False,
                        "index" -> k,
                        "hx" -> hx,
                        "member" -> member,
                        "message" ->
                            "Invalid entropy trajectory."
                    |>,
                    "logOnlyFailure"
                ]
            ];


            evolutionNormError =
                Lookup[
                    evolutionResult,
                    "maxNormError",
                    Infinity
                ];


            If[
                evolutionNormError >
                cfg[
                    "validationNormTolerance"
                ],

                Throw[
                    <|
                        "success" -> False,
                        "index" -> k,
                        "hx" -> hx,
                        "member" -> member,
                        "message" ->
                            "Norm error exceeded tolerance.",
                        "maxEvolutionNormError" ->
                            evolutionNormError
                    |>,
                    "logOnlyFailure"
                ]
            ];


            maxEvolutionNormError =
                Max[
                    maxEvolutionNormError,
                    evolutionNormError
                ];


            (* Direct-phase validation at four representative times. *)

            validationResults =
                Table[

                    directEntropyAndNormSector[
                        WP,
                        WM,
                        BpN,
                        BmN,
                        evalsP,
                        evalsM,
                        coeffP,
                        coeffM,
                        N[
                            validationTimes[[j]]
                        ],
                        dA,
                        dB
                    ],

                    {
                        j,
                        Length[
                            validationIndices
                        ]
                    }
                ];


            validationEntropyDifferences =
                Abs[
                    entropyMember[[
                        validationIndices
                    ]] -
                    validationResults[[All, 1]]
                ];


            validationNormErrors =
                validationResults[[All, 2]];


            If[
                Max[
                    validationEntropyDifferences
                ] >
                    cfg[
                        "validationEntropyTolerance"
                    ] ||

                Max[
                    validationNormErrors
                ] >
                    cfg[
                        "validationNormTolerance"
                    ],

                Throw[
                    <|
                        "success" -> False,
                        "index" -> k,
                        "hx" -> hx,
                        "member" -> member,
                        "message" ->
                            "Direct-phase validation failed.",
                        "maxEntropyDifference" ->
                            Max[
                                validationEntropyDifferences
                            ],
                        "maxNormError" ->
                            Max[
                                validationNormErrors
                            ]
                    |>,
                    "logOnlyFailure"
                ]
            ];


            maxValidationEntropyDifference =
                Max[
                    maxValidationEntropyDifference,
                    Max[
                        validationEntropyDifferences
                    ]
                ];


            maxValidationNormError =
                Max[
                    maxValidationNormError,
                    Max[
                        validationNormErrors
                    ]
                ];


            (* Member 1 is the single RPS. *)

            If[
                member == 1,

                singleEntropy =
                    Developer`ToPackedArray[
                        entropyMember
                    ];
            ];


            (* Online arithmetic ensemble mean. *)

            ensembleMean +=
                (
                    entropyMember -
                    ensembleMean
                )/
                member;


            Clear[
                evolutionResult,
                entropyMember
            ];

            ,
            {
                member,
                stateCount
            }
        ];


        ensembleMean =
            Developer`ToPackedArray[
                ensembleMean
            ];


        computeSeconds =
            N[
                AbsoluteTime[] -
                runStart
            ];


        (* ----------------------------------------------------- *)
        (* Save both observables in one file                     *)
        (* ----------------------------------------------------- *)

        runFile =
            FileNameJoin[
                {
                    outputDir,
                    "log_run_" <>
                    IntegerString[
                        k,
                        10,
                        2
                    ] <>
                    ".wxf"
                }
            ];


        resultAssociation = <|

            "success" -> True,

            "index" -> k,

            "L" -> L,

            "J" -> N[J],

            "hz" -> N[hz],

            "hxExponent" ->
                N[hxExponent],

            "hx" ->
                N[hx],

            "timeGridType" ->
                "Log",

            "timeGrid" ->
                times,

            "logTMin" ->
                cfg["logTMin"],

            "logTMax" ->
                cfg["logTMax"],

            "pointsPerDecadeParameter" ->
                cfg[
                    "pointsPerDecade"
                ],

            "nSamples" ->
                nSamples,

            "singleStateSeed" ->
                cfg["stateSeeds"][[1]],

            "ensembleStateSeeds" ->
                cfg["stateSeeds"],

            "ensembleSize" ->
                stateCount,

            "singleEntropy" ->
                singleEntropy,

            "ensembleEntropy" ->
                ensembleMean,

            "evalsP" ->
                evalsP,

            "evalsM" ->
                evalsM,

            "hermiticityError" ->
                hermiticityError,

            "reflectionBlockError" ->
                reflectionBlockError,

            "eigResidualP" ->
                eigResidualP,

            "eigResidualM" ->
                eigResidualM,

            "maxImagV" ->
                maxImagV,

            "maxReconstructionError" ->
                maxReconstructionError,

            "maxValidationEntropyDifference" ->
                maxValidationEntropyDifference,

            "maxValidationNormError" ->
                maxValidationNormError,

            "maxEvolutionNormError" ->
                maxEvolutionNormError,

            "diagonalizationSeconds" ->
                diagWall,

            "ensembleEvolutionSeconds" ->
                evolutionWall,

            "computeSecondsBeforeExport" ->
                computeSeconds

        |>;


        exportSeconds = First[
            AbsoluteTiming[

                exportResult =
                    checkedExport[
                        runFile,
                        resultAssociation,
                        "WXF"
                    ];

            ]
        ];


        If[
            FailureQ[exportResult],

            Throw[
                <|
                    "success" -> False,
                    "index" -> k,
                    "hx" -> hx,
                    "message" ->
                        "Export failed.",
                    "runFile" ->
                        runFile
                |>,
                "logOnlyFailure"
            ]
        ];


        totalSeconds =
            computeSeconds +
            exportSeconds;


        <|
            "success" -> True,

            "index" -> k,

            "hxExponent" ->
                hxExponent,

            "hx" ->
                hx,

            "diagonalizationSeconds" ->
                diagWall,

            "ensembleEvolutionSeconds" ->
                evolutionWall,

            "exportSeconds" ->
                exportSeconds,

            "totalSeconds" ->
                totalSeconds,

            "maxReconstructionError" ->
                maxReconstructionError,

            "maxValidationEntropyDifference" ->
                maxValidationEntropyDifference,

            "maxValidationNormError" ->
                maxValidationNormError,

            "runFile" ->
                runFile
        |>

    ],
    "logOnlyFailure"
];



(* ::Title::Closed:: *)
(* Log-only worker initialization *)


ClearAll[logOnlyInitializeWorkers];

logOnlyInitializeWorkers[] := Module[
    {
        missing,
        packageStatus,
        writeStatus
    },


    If[
        !TrueQ[
            logOnlyUseParallel
        ],
        Return[True]
    ];


    missing =
        Max[
            0,
            logOnlyRequestedKernels -
            Length[
                Kernels[]
            ]
        ];


    If[
        missing > 0,
        LaunchKernels[
            missing
        ]
    ];


    If[
        Length[
            Kernels[]
        ] == 0,

        Return[
            Failure[
                "NoParallelKernels",
                <||>
            ]
        ]
    ];


    packageStatus =
        With[
            {
                path =
                    qmbInitPath
            },

            ParallelEvaluate[
                Quiet[
                    Check[
                        Get[path];
                        True,
                        False
                    ]
                ]
            ]
        ];


    If[
        !AllTrue[
            packageStatus,
            TrueQ
        ],

        Return[
            Failure[
                "WorkerPackageLoadFailed",
                <|
                    "status" ->
                        packageStatus
                |>
            ]
        ]
    ];


    DistributeDefinitions[

        logOnlyRunOneHx,

        checkedExport,

        validEntropyVectorQ,

        sortedEigensystem,

        phaseAtTimesBatchC,

        reseedC,

        entropyFromSVC,

        entropyFast,

        entropyBatchFromSectors,

        arbitraryEntropyEvolutionSector,

        directEntropyAndNormSector,

        sampledEigensystemResidual
    ];


    writeStatus =
        With[
            {
                directory =
                    logOnlyOutputDir
            },

            ParallelEvaluate[

                Module[
                    {
                        p,
                        result
                    },

                    p =
                        FileNameJoin[
                            {
                                directory,
                                "worker_" <>
                                ToString[
                                    $KernelID
                                ] <>
                                ".tmp"
                            }
                        ];


                    result =
                        checkedExport[
                            p,
                            "write test",
                            "Text"
                        ];


                    If[
                        !FailureQ[result],
                        DeleteFile[p]
                    ];


                    !FailureQ[result]
                ]
            ]
        ];


    If[
        !AllTrue[
            writeStatus,
            TrueQ
        ],

        Return[
            Failure[
                "WorkerCannotWrite",
                <|
                    "status" ->
                        writeStatus
                |>
            ]
        ]
    ];


    Print[
        "Log-only sweep using ",
        Length[Kernels[]],
        " parallel kernels."
    ];


    True
];



(* ::Title::Closed:: *)
(* Run log-only sweep *)


logOnlyStartup =
    logOnlyInitializeWorkers[];


If[
    FailureQ[logOnlyStartup],
    Print[logOnlyStartup];
    Abort[];
];


logOnlySweepSeconds = First[
    AbsoluteTiming[

        logOnlyRunSummary = If[

            TrueQ[
                logOnlyUseParallel
            ],

            With[
                {
                    snapshot =
                        logOnlyConfig
                },

                ParallelMap[
                    logOnlyRunOneHx[
                        snapshot,
                        #
                    ] &,

                    Range[
                        Length[
                            snapshot[
                                "hxList"
                            ]
                        ]
                    ],

                    Method ->
                        "CoarsestGrained",

                    DistributedContexts ->
                        None
                ]
            ],

            Map[
                logOnlyRunOneHx[
                    logOnlyConfig,
                    #
                ] &,

                Range[
                    Length[
                        logOnlyConfig[
                            "hxList"
                        ]
                    ]
                ]
            ]
        ];

    ]
];


logOnlyFailedRuns =
    If[
        ListQ[
            logOnlyRunSummary
        ],

        Select[
            logOnlyRunSummary,

            !AssociationQ[#] ||
            !TrueQ[
                Lookup[
                    #,
                    "success",
                    False
                ]
            ] &
        ],

        {
            logOnlyRunSummary
        }
    ];


checkedExport[
    FileNameJoin[
        {
            logOnlyOutputDir,
            "run_summary.wxf"
        }
    ],
    logOnlyRunSummary,
    "WXF"
];


If[
    logOnlyFailedRuns =!= {},

    checkedExport[
        FileNameJoin[
            {
                logOnlyOutputDir,
                "failed_runs.wxf"
            }
        ],
        logOnlyFailedRuns,
        "WXF"
    ];

    Print[
        "ERROR: one or more logarithmic hx runs failed."
    ];

    Print[
        logOnlyFailedRuns
    ];

    Abort[];
];


Print[
    "Log-only sweep complete. Wall time = ",
    logOnlySweepSeconds,
    " s."
];



(* ::Title::Closed:: *)
(* Log-only plotting helpers *)


ClearAll[logOnlyPowerLabel];

logOnlyPowerLabel[
    exponent_?NumericQ,
    size_: 22
] :=
    Style[
        Superscript[
            "10",
            If[
                Abs[
                    exponent -
                    Round[exponent]
                ] <
                10^-10,

                Round[exponent],

                N[exponent]
            ]
        ],
        size,
        Black
    ];


ClearAll[logOnlyPlotResults];

logOnlyPlotResults[
    directory_String,
    window_Integer : 200
] := Module[
    {
        files,
        runs,
        times,
        exponents,
        lo,
        hi,
        cf,
        colors,
        barRange,
        barTicks,
        legend,
        decadeTicks,
        plotDirectory,

        singleRawData,
        ensembleRawData,

        singleAverageData,
        ensembleAverageData,

        baseTitle,
        singleTitle,
        ensembleTitle,

        makeLogPlot,

        graphic,
        paths = <||>,
        result
    },


    If[
        window < 1,

        Return[
            Failure[
                "InvalidMovingAverageWindow",
                <||>
            ]
        ]
    ];


    files =
        Sort[
            FileNames[
                "log_run_" ~~
                DigitCharacter .. ~~
                ".wxf",
                directory
            ]
        ];


    If[
        files === {},

        Return[
            Failure[
                "NoSavedLogTrajectories",
                <|
                    "directory" ->
                        directory
                |>
            ]
        ]
    ];


    runs =
        Import[
            #,
            "WXF"
        ] &
        /@
        files;


    If[
        !AllTrue[
            runs,

            AssociationQ[#] &&
            KeyExistsQ[
                #,
                "singleEntropy"
            ] &&
            KeyExistsQ[
                #,
                "ensembleEntropy"
            ] &&
            KeyExistsQ[
                #,
                "timeGrid"
            ] &
        ],

        Return[
            Failure[
                "InvalidLogResultFiles",
                <||>
            ]
        ]
    ];


    runs =
        SortBy[
            runs,
            #[
                "hxExponent"
            ] &
        ];


    times =
        First[runs][
            "timeGrid"
        ];


    exponents =
        Lookup[
            runs,
            "hxExponent"
        ];


    lo =
        Min[
            exponents
        ];


    hi =
        Max[
            exponents
        ];


    cf =
        Blend[
            {
                Blue,
                Purple,
                Red
            },
            #
        ] &;


    colors =
        If[
            lo == hi,

            ConstantArray[
                cf[0.5],
                Length[runs]
            ],

            cf /@
            Rescale[
                exponents,
                {
                    lo,
                    hi
                }
            ]
        ];


    barRange =
        If[
            lo == hi,
            {
                lo - 0.5,
                hi + 0.5
            },
            {
                lo,
                hi
            }
        ];


    barTicks =
        (
            {
                #,
                logOnlyPowerLabel[
                    #,
                    24
                ]
            } &
        ) /@
        Range[
            Ceiling[lo],
            Floor[hi]
        ];


    legend =
        BarLegend[
            {
                cf,
                barRange
            },

            LegendLabel ->
                Style[
                    Subscript[
                        "h",
                        "x"
                    ],
                    30,
                    Black
                ],

            LabelStyle ->
                {
                    Black,
                    20
                },

            ColorFunctionScaling ->
                True,

            LegendMarkerSize ->
                {
                    12,
                    360
                },

            Ticks ->
                barTicks
        ];


    decadeTicks =
        Table[
            {
                10.^x,
                logOnlyPowerLabel[
                    x,
                    22
                ]
            },

            {
                x,
                Ceiling[
                    First[runs][
                        "logTMin"
                    ]
                ],
                Floor[
                    First[runs][
                        "logTMax"
                    ]
                ],
                1
            }
        ];


    singleRawData =
        Table[
            Transpose[
                {
                    times,
                    runs[[k]][
                        "singleEntropy"
                    ]
                }
            ],
            {
                k,
                Length[runs]
            }
        ];


    ensembleRawData =
        Table[
            Transpose[
                {
                    times,
                    runs[[k]][
                        "ensembleEntropy"
                    ]
                }
            ],
            {
                k,
                Length[runs]
            }
        ];


    (* Moving average starts at the first saved logarithmic sample.
       Time coordinates are never averaged. *)

    singleAverageData =
        Table[
            Transpose[
                {
                    times,

                    startupMovingAverage[
                        runs[[k]][
                            "singleEntropy"
                        ],
                        window
                    ]
                }
            ],

            {
                k,
                Length[runs]
            }
        ];


    ensembleAverageData =
        Table[
            Transpose[
                {
                    times,

                    startupMovingAverage[
                        runs[[k]][
                            "ensembleEntropy"
                        ],
                        window
                    ]
                }
            ],

            {
                k,
                Length[runs]
            }
        ];


    baseTitle =
        Row[
            {
                "L = ",
                N[
                    First[runs][
                        "L"
                    ]
                ],

                ", J = ",
                N[
                    First[runs][
                        "J"
                    ]
                ],

                ", ",
                Subscript[
                    "h",
                    "z"
                ],
                " = ",
                N[
                    First[runs][
                        "hz"
                    ]
                ]
            }
        ];


    singleTitle =
        Row[
            {
                "RPS, ",
                baseTitle
            }
        ];


    ensembleTitle =
        Row[
            {
                "RPS ensemble mean (",
                First[runs][
                    "ensembleSize"
                ],
                " states), ",
                baseTitle
            }
        ];


    plotDirectory =
        FileNameJoin[
            {
                directory,
                "plots_log_only"
            }
        ];


    If[
        !DirectoryQ[
            plotDirectory
        ],

        Quiet[
            Check[
                CreateDirectory[
                    plotDirectory
                ],
                Null
            ]
        ]
    ];


    If[
        !DirectoryQ[
            plotDirectory
        ],

        Return[
            Failure[
                "PlotDirectoryFailed",
                <|
                    "directory" ->
                        plotDirectory
                |>
            ]
        ]
    ];


    makeLogPlot[
        coordinates_,
        title_,
        caption_
    ] :=
        Legended[

            ListLogLinearPlot[
                coordinates,

                Joined ->
                    False,

                PlotStyle ->
                    (
                        Directive[
                            #,
                            Opacity[0.82],
                            PointSize[0.002]
                        ] &
                    ) /@
                    colors,

                PlotRange ->
                    {
                        All,
                        {
                            0.,
                            All
                        }
                    },

                Frame ->
                    True,

                Axes ->
                    False,

                ImageSize ->
                    1000,

                PlotTheme ->
                    "Detailed",

                FrameStyle ->
                    Directive[
                        Black,
                        22
                    ],

                PlotLabel ->
                    Style[
                        Column[
                            {
                                title,
                                caption
                            },
                            Alignment ->
                                Center
                        ],
                        22,
                        Black
                    ],

                FrameLabel ->
                    {
                        Style[
                            "t",
                            26,
                            Black
                        ],

                        Style[
                            "S(t)",
                            26,
                            Black
                        ]
                    },

                FrameTicks ->
                    {
                        {
                            Automatic,
                            None
                        },

                        {
                            decadeTicks,
                            None
                        }
                    },

                GridLines ->
                    None
            ],

            Placed[
                legend,
                Right
            ]
        ];


    (* --------------------------------------------------------- *)
    (* Single RPS: raw log plot                                  *)
    (* --------------------------------------------------------- *)

    graphic =
        makeLogPlot[
            singleRawData,
            singleTitle,
            "Unsmooth entropy"
        ];


    result =
        checkedExport[
            FileNameJoin[
                {
                    plotDirectory,
                    "single_RPS_raw_logtime.png"
                }
            ],
            graphic,
            "PNG"
        ];


    AssociateTo[
        paths,
        "singleRawLog" ->
            result
    ];


    (* --------------------------------------------------------- *)
    (* Single RPS: moving-average log plot                       *)
    (* --------------------------------------------------------- *)

    graphic =
        makeLogPlot[
            singleAverageData,
            singleTitle,
            "Moving average, window = " <>
            ToString[window]
        ];


    result =
        checkedExport[
            FileNameJoin[
                {
                    plotDirectory,
                    "single_RPS_moving_average_logtime.png"
                }
            ],
            graphic,
            "PNG"
        ];


    AssociateTo[
        paths,
        "singleMovingAverageLog" ->
            result
    ];


    (* --------------------------------------------------------- *)
    (* Ensemble: raw log plot                                    *)
    (* --------------------------------------------------------- *)

    graphic =
        makeLogPlot[
            ensembleRawData,
            ensembleTitle,
            "Unsmooth entropy"
        ];


    result =
        checkedExport[
            FileNameJoin[
                {
                    plotDirectory,
                    "ensemble_RPS_raw_logtime.png"
                }
            ],
            graphic,
            "PNG"
        ];


    AssociateTo[
        paths,
        "ensembleRawLog" ->
            result
    ];


    (* --------------------------------------------------------- *)
    (* Ensemble: moving-average log plot                         *)
    (* --------------------------------------------------------- *)

    graphic =
        makeLogPlot[
            ensembleAverageData,
            ensembleTitle,
            "Moving average, window = " <>
            ToString[window]
        ];


    result =
        checkedExport[
            FileNameJoin[
                {
                    plotDirectory,
                    "ensemble_RPS_moving_average_logtime.png"
                }
            ],
            graphic,
            "PNG"
        ];


    AssociateTo[
        paths,
        "ensembleMovingAverageLog" ->
            result
    ];


    Clear[
        graphic
    ];


    If[
        AnyTrue[
            Values[paths],
            FailureQ
        ],

        Print[
            "One or more log-only plot exports failed: ",
            paths
        ],

        Print[
            "Log-only plots saved in: ",
            plotDirectory
        ]
    ];


    paths
];



(* ::Title::Closed:: *)
(* Export log-only plots *)


logOnlyPlotFiles =
    logOnlyPlotResults[
        logOnlyOutputDir,
        logOnlyMovingAverageWindow
    ];


Print[
    "Log-only section complete."
];


(* ----------------------------------------------------------------- *)
(* No CloseKernels[] is used here.  Existing kernels remain available *)
(* for the rest of the original file.                                *)
(* ----------------------------------------------------------------- *)
