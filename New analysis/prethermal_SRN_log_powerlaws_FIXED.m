(* ::Package:: *)

(* ::Title:: *)
(* SRN prethermalization*)


(* ::Input:: *)
(**)


(* Based on prethermal_entropy_only_LOG_analysis(2).m.
   Definitions and plot functions precede the editable parameters and execution.
   h in the power plots means hx; g means hz. Entropy is in natural-log units.
   Both S*=0.6 and S*=1.4 are measured for the single RPS and ensemble mean.
   No individual ensemble trajectories are retained, except member 1, which is
   also the requested single-RPS result. Each Hamiltonian is diagonalized once.

   After the first evaluation, edit/evaluate the PARAMETERS section and call
       logOnlyRunAnalysis[]
   again. No kernel restart is required. The last line runs this automatically.
   To refit/replot existing results without any dynamics, evaluate definitions
   only and use the examples in the final commented section.

   Smoothing exactly preserves the supplied LOG file: at sample i average
   S[Max[1,i-window+1] ... i], keeping time t_i unchanged. Startup windows grow
   until full length. No time-zero point is invented on a positive log grid.
   Crossing = first smoothed sample with S >= S*. No interpolation is used.
   A crossing at the first saved sample is left-censored and is not fitted.
   No crossing before tmax is right-censored and is not fitted.
   The ensemble crossing is taken AFTER averaging entropies, not by averaging
   member crossing times or taking the entropy of the mixed ensemble state.

   Fit: Log[T_M] = logA + b Log[hx], so T_M = Exp[logA] hx^b.
   Default fits use all uncensored fields; edit the fit ranges to restrict SRN.
   A good R-squared alone does not establish an asymptotic power law. Check
   fitted-field range, smoothing window, grid density and ensemble convergence.
   Standard errors are regression diagnostics, not ensemble error bars.
   Changing points/decade at fixed window changes its span in log time.
   Existing machine-precision phase evolution and tolerances are preserved;
   at very long times the same-eigensystem checks do not bound spectral phase
   error t*deltaE. This addition does not claim arbitrary-precision accuracy.

   Workflow repair: one complete logOnlyRunAnalysis definition; no notebook
   cell/section markers occur inside a function. One worker initializer lives
   with the definitions. The final execution cell explicitly calls the run.
   Parallel availability is tracked separately from logOnlyUseParallel, so a
   serial fallback does not change your preference for the next run.

   Validation of this repair: entire-file and all cell-boundary syntax checks
   passed in Wolfram. Sequential evaluation installed exactly one run definition.
   The actual execution cell completed two small serial test sweeps after hz
   and time-grid changes, saved data, and returned distinct output directories.
   Tests used an explicit test Ising model/product-state generator and WVM
   compilation; native C, your QMB package, real parallel workers and PNG/PDF
   rendering were not tested here. Plot objects were constructed successfully.

   Wolfram reference:
   https://reference.wolfram.com/language/ref/LinearModelFit.html
   https://reference.wolfram.com/language/ref/MovingAverage.html
*)



(* ::Section::Closed:: *)
(* Dependency setup: edit this path once for your QMB installation *)


(* ::Input:: *)
(**)


qmbInitPath = "C:\\Users\\Miguel\\Github\\libs\\QMB\\Kernel\\init.m";
If[!FileExistsQ[qmbInitPath], Print["ERROR: QMB init.m not found: ", qmbInitPath]; Abort[]];
Get[qmbInitPath];

(* Resolve an absolute default without changing the working directory. *)
logOnlyDefaultDirectory = If[StringQ[$InputFileName] && StringLength[$InputFileName] > 0,
    DirectoryName[ExpandFileName[$InputFileName]],
    Quiet[Check[NotebookDirectory[], Directory[]]]];
If[!StringQ[logOnlyDefaultDirectory], logOnlyDefaultDirectory = Directory[]];



(* ::Section::Closed:: *)
(* General numerical definitions: unchanged logarithmic evolution *)


(* ::Input:: *)
(**)


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

(* Batched arbitrary/logarithmic grid: every sample is generated directly
   from the original spectral coefficients, so recurrence drift is absent. *)
ClearAll[arbitraryEntropyEvolutionSector];
arbitraryEntropyEvolutionSector[
    WP_, WM_, BpN_, BmN_, evalsP_, evalsM_, coeffP_, coeffM_,
    times_, batchSize_Integer, dA_Integer, dB_Integer
] := Catch[Module[
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
            Throw[Failure["InvalidEntropyBatch", <|"startIndex" -> pos|>], "srnBatchFailure"]
        ];
        maxNormError = Max[maxNormError, batchResult["maxNormError"]];
        entropy[[pos ;; pos + len - 1]] = values;
        pos += len;
    ];
    <|"entropy" -> entropy, "maxNormError" -> maxNormError|>
], "srnBatchFailure"];

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
ClearAll[tagNumber, checkedExport];
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



(* ::Section::Closed:: *)
(* Directory, field-sweep and worker functions *)


(* ::Input:: *)
(**)


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
            "entropy_SRN"
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

ClearAll[logOnlyInitializeWorkers];

logOnlyInitializeWorkers[] := Module[
    {missing, workers, status},

    logOnlyParallelActive = False;

    If[!TrueQ[logOnlyUseParallel],
        Return[True]
    ];

    missing = Max[
        0,
        logOnlyRequestedKernels - Length[Kernels[]]
    ];

    If[missing > 0,
        LaunchKernels[missing]
    ];

    workers = Kernels[];

    If[workers === {},
        Print[
            "No parallel kernels available. ",
            "Running this sweep serially."
        ];
        logOnlyParallelActive = False;
        Return[True]
    ];

    (* Initialize one worker at a time.
       Reuse QMB if its Hamiltonian is already defined. *)

    status = Table[
        With[
            {
                worker = workers[[k]],
                path = ExpandFileName[qmbInitPath],
                directory = logOnlyOutputDir
            },

            First[
                ParallelEvaluate[
                    Module[
                        {
                            ready,
                            loaded = True,
                            probe,
                            written
                        },

                        ready =
                            Length[DownValues[IsingHamiltonian]] > 0;

                        If[!ready,
                            If[!FileExistsQ[path],
                                Return[
                                    <|
                                        "kernel" -> $KernelID,
                                        "success" -> False,
                                        "reason" -> "QMB file not found",
                                        "path" -> path
                                    |>
                                ]
                            ];

                            (* Keep package messages visible. *)
                            loaded = CheckAbort[
                                Get[path];
                                True,
                                False
                            ];

                            ready =
                                Length[
                                    DownValues[IsingHamiltonian]
                                ] > 0;
                        ];

                        If[!TrueQ[loaded] || !TrueQ[ready],
                            Return[
                                <|
                                    "kernel" -> $KernelID,
                                    "success" -> False,
                                    "reason" ->
                                        "QMB loading aborted or Hamiltonian undefined"
                                |>
                            ]
                        ];

                        (* Confirm this worker can save results. *)

                        probe = FileNameJoin[{
                            directory,
                            "worker_" <>
                                ToString[$KernelID] <>
                                ".tmp"
                        }];

                        written = Quiet[
                            Check[
                                Export[probe, "write test", "Text"],
                                $Failed
                            ]
                        ];

                        If[
                            !StringQ[written] ||
                            !FileExistsQ[probe],

                            Return[
                                <|
                                    "kernel" -> $KernelID,
                                    "success" -> False,
                                    "reason" ->
                                        "Cannot write to output directory",
                                    "path" -> directory
                                |>
                            ]
                        ];

                        DeleteFile[probe];

                        <|
                            "kernel" -> $KernelID,
                            "success" -> True
                        |>
                    ],
                    {worker}
                ]
            ]
        ],
        {k, Length[workers]}
    ];

    logOnlyWorkerStatus = status;

    If[
        !AllTrue[
            status,
            AssociationQ[#] &&
                TrueQ[Lookup[#, "success", False]] &
        ],

        Print["Worker initialization results: ", status];

        Print[
            "Using serial execution for this sweep ",
            "because a worker failed."
        ];

        logOnlyParallelActive = False;
        Return[True]
    ];

    (* Redistribute the current numerical definitions on every run. *)

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

    logOnlyParallelActive = True;

    Print[
        "Using ",
        Length[workers],
        " initialized parallel kernels."
    ];

    True
];



(* ::Section::Closed:: *)
(* Entropy plot definitions, including the hx colorbar *)


(* ::Input:: *)
(**)


ClearAll[startupMovingAverage];
startupMovingAverage[values_List, window_Integer] := Module[{n = Length[values], prefix},
    If[n == 0 || window < 1, Return[{}]];
    If[n < window, Return[Accumulate[values]/Range[n]]];
    prefix = If[window == 1, {},
        Accumulate[Take[values, window - 1]]/Range[window - 1]];
    Developer`ToPackedArray[Join[prefix, MovingAverage[values, window]]]
];
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
                "plots"
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



(* ::Section::Closed:: *)
(* Moving-average crossing times and power-law fits *)


(* ::Input:: *)
(**)


ClearAll[srnFiniteRealQ, srnRangeQ, srnLoadRuns, srnCrossing,
    srnFitCrossings, srnPowerPlot, srnAnalyzeSaved, srnDescription];
srnFiniteRealQ[x_] := NumberQ[x] && TrueQ[Im[x] == 0] &&
    FreeQ[x, Indeterminate | _DirectedInfinity];
srnRangeQ[r_] := r === Automatic ||
    (MatchQ[r, {_?srnFiniteRealQ, _?srnFiniteRealQ}] && TrueQ[0 < r[[1]] < r[[2]]]);

srnLoadRuns[directory_String] := Module[{files, runs, first},
    files = Sort[FileNames["log_run_" ~~ DigitCharacter .. ~~ ".wxf", directory]];
    If[files === {}, Return[Failure["NoSavedLogTrajectories", <|"directory" -> directory|>]]];
    runs = Quiet[Check[Import[#, "WXF"], $Failed]] & /@ files;
    If[!AllTrue[runs, Function[r, AssociationQ[r] &&
        AllTrue[{"success", "hx", "hxExponent", "L", "J", "hz", "timeGrid",
            "singleEntropy", "ensembleEntropy", "ensembleSize"}, KeyExistsQ[r, #] &]]],
        Return[Failure["InvalidLogFiles", <||>]]];
    first = First[runs];
    If[!AllTrue[runs, Function[r,
        TrueQ[r["success"]] && srnFiniteRealQ[r["hx"]] && TrueQ[r["hx"] > 0] &&
        VectorQ[r["timeGrid"], srnFiniteRealQ] && Length[r["timeGrid"]] >= 2 &&
        TrueQ[Min[r["timeGrid"]] > 0] && TrueQ[Min[Differences[r["timeGrid"]]] > 0] &&
        r["timeGrid"] === first["timeGrid"] &&
        Lookup[r, {"L", "J", "hz", "ensembleSize"}] ===
            Lookup[first, {"L", "J", "hz", "ensembleSize"}] &&
        AllTrue[{"singleEntropy", "ensembleEntropy"}, Function[key,
            VectorQ[r[key], srnFiniteRealQ] && Length[r[key]] == Length[r["timeGrid"]]]]]],
        Return[Failure["InconsistentLogData", <||>]]];
    If[Length[DeleteDuplicates[Lookup[runs, "hx"]]] != Length[runs],
        Return[Failure["DuplicateFields", <||>]]];
    SortBy[runs, #["hx"] &]
];

(* Input is already smoothed. The crossing is never extracted from raw entropy.
   Bounds bracket the first observed upward crossing on the sampled curve;
   they are sampling bounds, not statistical confidence intervals. *)
srnCrossing[times_List, smooth_List, threshold_?srnFiniteRealQ] := Module[{pos, i},
    If[Length[times] != Length[smooth] || Length[times] < 2 ||
       !VectorQ[times, srnFiniteRealQ] || !VectorQ[smooth, srnFiniteRealQ] ||
       !TrueQ[Min[times] > 0] || !TrueQ[Min[Differences[times]] > 0],
        Return[Failure["InvalidCrossingData", <||>]]];
    pos = FirstPosition[smooth, _?(TrueQ[# >= threshold] &), Missing["NoCrossing"]];
    If[MissingQ[pos], Return[<|"status" -> "NotReached", "time" -> Missing["NotReached"],
        "lowerTime" -> Last[times], "upperTime" -> Missing["BeyondTMax"],
        "sampleIndex" -> Missing["NotReached"]|>]];
    i = First[pos];
    If[i == 1, Return[<|"status" -> "AtOrBeforeFirstSample",
        "time" -> Missing["LeftCensored"], "lowerTime" -> 0., "upperTime" -> First[times],
        "sampleIndex" -> 1|>]];
    <|"status" -> "Crossed", "time" -> times[[i]],
      "lowerTime" -> times[[i - 1]], "upperTime" -> times[[i]], "sampleIndex" -> i|>
];

srnFitCrossings[rows_List, range_: Automatic] := Module[
    {selected, pairs, lm, z, logA, b, residual, ss, sst, n, stderr, rsq, marked},
    If[!srnRangeQ[range], Return[Failure["InvalidFitRange", <|"range" -> range|>]]];
    marked = Map[Function[row, Join[row, <|"selectedForFit" ->
        (row["status"] === "Crossed" &&
         (range === Automatic || TrueQ[range[[1]] <= row["hx"] <= range[[2]]]))|>]], rows];
    selected = Select[marked, TrueQ[#["selectedForFit"]] &];
    pairs = Lookup[#, {"hx", "time"}] & /@ selected;
    n = Length[pairs];
    If[n < 3 || Length[DeleteDuplicates[pairs[[All, 1]]]] < 3,
        Return[<|"status" -> "InsufficientCrossings", "n" -> n,
            "range" -> range, "rows" -> marked, "pairs" -> pairs|>]];
    lm = Quiet[LinearModelFit[Log[pairs], {1, z}, z]];
    If[Head[lm] =!= FittedModel, Return[Failure["FitFailed", <|"pairs" -> pairs|>]]];
    logA = N[Normal[lm] /. z -> 0];
    b = N[Coefficient[Normal[lm], z]];
    residual = Log[pairs[[All, 2]]] - (logA + b Log[pairs[[All, 1]]]);
    ss = Total[residual^2];
    sst = Total[(Log[pairs[[All, 2]]] - Mean[Log[pairs[[All, 2]]]])^2];
    rsq = If[TrueQ[sst > 0], 1. - ss/sst, Missing["ConstantCrossingTimes"]];
    stderr = Sqrt[ss/(n - 2)/Total[(Log[pairs[[All, 1]]] - Mean[Log[pairs[[All, 1]]]])^2]];
    <|"status" -> "Fitted", "n" -> n, "range" -> range,
      "rows" -> marked, "pairs" -> pairs, "logA" -> logA, "A" -> Exp[logA],
      "b" -> b, "bStandardError" -> stderr, "RSquaredLogSpace" -> rsq,
      "logResiduals" -> residual,
      "localSlopes" -> Transpose[{Sqrt[Most[pairs[[All, 1]]] Rest[pairs[[All, 1]]]],
          Differences[Log[pairs[[All, 2]]]]/Differences[Log[pairs[[All, 1]]]]}]|>
];

(* Same visual convention as the supplied power plot: black circles, red dashed
   fit, full black frame, white background, serif h/T_M labels and inset power.
   Each threshold gets its own plot; no individual-h legend is added. *)
srnPowerPlot[fit_Association, meta_Association, dataset_String, threshold_] := Module[
    {pairs = fit["pairs"], curve, x, key, label, title, fitted},
    fitted = fit["status"] === "Fitted";
    title = Column[{
        Style[Row[{"SRN,   L = ", meta["L"], ",   J = ", meta["J"],
            ",   g = ", meta["hz"]}], 34, Black, FontFamily -> "Times"],
        Style[Row[{If[dataset === "single", "Single RPS", "RPS ensemble mean"],
            ",   ", Subscript["S", "*"], " = ", threshold}], 20, Black]
        }, Alignment -> Center];
    If[pairs === {}, Return[Graphics[Text[Style["No uncensored crossings in the fit range", 22]],
        PlotLabel -> title, ImageSize -> 1000, Background -> White]]];
    curve = If[fitted,
        x = Exp[Subdivide[Log[Min[pairs[[All, 1]]]], Log[Max[pairs[[All, 1]]]], 300]];
        Transpose[{x, Exp[fit["logA"] + fit["b"] Log[x]]}], {}];
    key = Graphics[{Red, Dashed, AbsoluteThickness[2.5], Line[{{0, .5}, {1, .5}}]},
        PlotRange -> {{0, 1}, {0, 1}}, ImagePadding -> 2, ImageSize -> 48, AspectRatio -> 1/5];
    label = Framed[If[fitted,
        Row[{key, Spacer[8], Style[Row[{Subscript["T", "M"], " \[Proportional] ",
            Superscript["h", NumberForm[fit["b"], {8, 3}]]}], 28, Black, FontFamily -> "Times"]}],
        Style["At least 3 crossing fields required", 18]],
        Background -> White, FrameStyle -> GrayLevel[.65], RoundingRadius -> 6, FrameMargins -> 8];
    ListLogLogPlot[If[fitted, {pairs, curve}, {pairs}],
        Joined -> If[fitted, {False, True}, {False}],
        PlotStyle -> {Directive[Black, AbsolutePointSize[9]],
            Directive[Red, Dashed, AbsoluteThickness[2.5]]},
        PlotMarkers -> If[fitted, {{"\[FilledCircle]", 14}, None}, {{"\[FilledCircle]", 14}}],
        Frame -> True, Axes -> False, PlotTheme -> "Detailed", Background -> White,
        FrameStyle -> Directive[Black, 26], BaseStyle -> {FontFamily -> "Times"},
        FrameLabel -> {Style["h", Italic, 34], Style[Subscript["T", "M"], Italic, 34]},
        PlotLabel -> title, GridLines -> None, PlotRange -> All, ImageSize -> 1000,
        Epilog -> Inset[label, Scaled[{.88, .90}], {Right, Center}]]
];

(* Read saved single/mean curves and create a short, unique analysis directory.
   Repeated fits with different windows/ranges never overwrite earlier fits.
   fitRanges is Automatic or a list aligned with thresholds, e.g.
       {Automatic, {10^-4, 10^-1}}
   The same threshold-specific range is applied to single and ensemble data. *)
srnAnalyzeSaved[directory_String, window_Integer: 200,
    thresholds_List: {0.6, 1.4}, fitRanges_: Automatic] := Catch[Module[
    {runs, ranges, out, smoothed, results = {}, rows, fit, fields, times, r,
     dataset, key, j, k, crossing, item, path, exports = {}, summaryKeys,
     crossingKeys, csv, checked, analysis, graphic, prefix},
    If[window < 1 || thresholds === {} || !VectorQ[thresholds, srnFiniteRealQ] ||
        !TrueQ[Min[thresholds] > 0] || Length[DeleteDuplicates[thresholds]] != Length[thresholds],
        Throw[Failure["InvalidPowerSettings", <||>], "srnAnalysisFailure"]];
    ranges = If[fitRanges === Automatic, ConstantArray[Automatic, Length[thresholds]], fitRanges];
    If[!ListQ[ranges] || Length[ranges] != Length[thresholds] || !AllTrue[ranges, srnRangeQ],
        Throw[Failure["InvalidFitRanges", <||>], "srnAnalysisFailure"]];
    runs = srnLoadRuns[directory];
    If[FailureQ[runs], Throw[runs, "srnAnalysisFailure"]];
    out = FileNameJoin[{directory, "power_" <> StringTake[StringDelete[CreateUUID[], "-"], 8]}];
    Quiet[Check[CreateDirectory[out], Null]];
    If[!DirectoryQ[out], Throw[Failure["CreatePowerDirectoryFailed", <|"path" -> out|>], "srnAnalysisFailure"]];
    checked[name_, data_, format_] := Module[{result},
        result = checkedExport[FileNameJoin[{out, name}], data, format];
        If[FailureQ[result], Throw[result, "srnAnalysisFailure"]];
        AppendTo[exports, result]; result];
    csv[data_, keys_] := Prepend[
        (Map[If[MissingQ[#], ToString[#, InputForm], #] &, Lookup[#, keys, ""]] & /@ data), keys];
    times = First[runs]["timeGrid"];
    smoothed = Map[Function[run, <|"hx" -> run["hx"], "timeGrid" -> times,
        "single" -> startupMovingAverage[run["singleEntropy"], window],
        "ensemble" -> startupMovingAverage[run["ensembleEntropy"], window]|>], runs];
    (* Save exact plotted/fitted moving averages for later use. *)
    checked["moving_average.wxf", <|"window" -> window, "curves" -> smoothed,
        "convention" -> "Trailing sample mean with growing startup; original times"|>, "WXF"];
    Do[
        Do[
            rows = Table[
                crossing = srnCrossing[times, smoothed[[k]][dataset], thresholds[[j]]];
                If[FailureQ[crossing], Throw[crossing, "srnAnalysisFailure"]];
                Join[<|"dataset" -> dataset, "threshold" -> thresholds[[j]],
                    "hx" -> runs[[k]]["hx"]|>, crossing], {k, Length[runs]}];
            fit = srnFitCrossings[rows, ranges[[j]]];
            If[FailureQ[fit], Throw[fit, "srnAnalysisFailure"]];
            item = Join[<|"dataset" -> dataset, "threshold" -> thresholds[[j]]|>, fit];
            AppendTo[results, item];
            Print[dataset, ", S* = ", thresholds[[j]], ": ", fit["status"],
                ", ", fit["n"], " crossing fields",
                If[fit["status"] === "Fitted", Row[{", b = ", fit["b"],
                    ", R^2(log) = ", fit["RSquaredLogSpace"]}], ""]];
        , {j, Length[thresholds]}];
    , {dataset, {"single", "ensemble"}}];
    analysis = <|"sourceDirectory" -> directory, "outputDirectory" -> out,
        "parameters" -> KeyTake[First[runs], {"L", "J", "hz", "ensembleSize"}],
        "window" -> window, "thresholds" -> thresholds, "fitRanges" -> ranges,
        "crossingRule" -> "First smoothed sample >= threshold; censored cases excluded",
        "fitRule" -> "Unweighted least squares: Log[T_M] = logA + b Log[hx]",
        "results" -> results|>;
    (* Data and diagnostics are written before rendering any graphics. *)
    checked["power_laws.wxf", analysis, "WXF"];
    summaryKeys = {"dataset", "threshold", "status", "n", "A", "b", "bStandardError", "RSquaredLogSpace"};
    crossingKeys = {"dataset", "threshold", "hx", "status", "time", "lowerTime", "upperTime", "sampleIndex", "selectedForFit"};
    checked["fits.csv", csv[results, summaryKeys], "CSV"];
    checked["crossings.csv", csv[Flatten[Lookup[results, "rows"], 1], crossingKeys], "CSV"];
    checked["analysis.txt", StringRiffle[{
        "SRN threshold power laws; h = hx, g = hz; entropy in nats.",
        "Parameters: " <> ToString[analysis["parameters"], InputForm],
        "Moving-average window: " <> ToString[window] <> " samples; growing startup, original times.",
        "Thresholds: " <> ToString[thresholds, InputForm],
        "hx fit ranges (aligned with thresholds): " <> ToString[ranges, InputForm],
        "First smoothed sample >= threshold, without interpolation. Censored fields excluded.",
        "Fit is Log[T] = logA + b Log[hx]; at least 3 distinct crossing fields required.",
        "Regression standard errors do not quantify ensemble, smoothing or sampling uncertainty.",
        "See crossings.csv for all fields, censoring, sampling brackets and fit selection.",
        "See power_laws.wxf for fitted pairs, local slopes and log residuals."
        }, "\n"], "Text"];
    Do[
        item = results[[r]];
        graphic = srnPowerPlot[item, analysis["parameters"], item["dataset"], item["threshold"]];
        prefix = item["dataset"] <> "_S" <> tagNumber[item["threshold"]];
        checked[prefix <> ".png", graphic, "PNG"];
        checked[prefix <> ".pdf", graphic, "PDF"];
    , {r, Length[results]}];
    Print["Power-law data and plots saved in: ", out];
    Join[analysis, <|"exports" -> exports|>]
], "srnAnalysisFailure"];

srnDescription[meta_Association, window_Integer, thresholds_List, ranges_] := StringRiffle[{
    "SRN logarithmic prethermalization analysis. h = hx; g = hz; entropy in nats.",
    "Hamiltonian: QMB IsingHamiltonian[hx,hz,J,L], unchanged from the supplied LOG workflow.",
    "Parameters: " <> ToString[KeyTake[meta, {"L", "J", "hz", "hxExponents", "hxList"}], InputForm],
    "RPS: QMB RandomChainProductState[L], generated once per seed and reused at every hx.",
    "Single RPS = ensemble member 1. Ensemble mean = arithmetic mean of member entropies.",
    "Ensemble size: " <> ToString[meta["ensembleSize"]] <> "; seeds: " <> ToString[meta["stateSeeds"], InputForm],
    "Initial state vectors and exact sampled time array are saved in metadata.wxf.",
    "Positive logarithmic grid including both endpoints, with no added t=0 point.",
    "log10(tmin), log10(tmax): " <> ToString[{meta["logTMin"], meta["logTMax"]}, InputForm],
    "pointsPerDecade parameter: " <> ToString[meta["pointsPerDecadeParameter"]],
    "N = Round[(logTMax-logTMin)*pointsPerDecade] = " <> ToString[meta["numberOfTimeSamples"]],
    "t_i = 10^(logTMin+(i-1)*(logTMax-logTMin)/(N-1)), i=1,...,N.",
    "Delta log10(t): " <> ToString[N[(meta["logTMax"]-meta["logTMin"])/(meta["numberOfTimeSamples"]-1)], InputForm],
    "Physical time step is variable: delta t_i = t_i*(10^deltaLog10t-1).",
    "Moving average: trailing " <> ToString[window] <> " samples, with growing startup windows; time coordinates unchanged.",
    "Crossing thresholds: " <> ToString[thresholds, InputForm] <> "; fit ranges: " <> ToString[ranges, InputForm],
    "Dynamics uses machine precision, direct absolute phases and the original entropy cutoff 1e-14.",
    "Changing sampling density at a fixed sample window changes the smoothing span in log time.",
    "Existing kernels are reused; each run receives a fresh immutable configuration and unique directory.",
    "Wolfram version: " <> meta["WolframVersion"]
}, "\n"];



(* ::Section::Closed:: *)
(* Complete execution function: ONE complete input cell *)


(* ::Input:: *)
(**)


ClearAll[logOnlyRunAnalysis];
logOnlyRunAnalysis[] := Module[{},
    Print["Starting SRN logarithmic analysis: L = ", logOnlyL,
        ", hz = ", logOnlyHz, ", ", Length[logOnlyHxExponents],
        " fields, ", logOnlyEnsembleSize, " RPS per field."];
    logOnlyParallelActive = False;
    Clear[logOnlyPlotFiles, logOnlyPowerResults, logOnlyRunSummary, logOnlyOutputDir];
    If[!IntegerQ[logOnlyL] || !EvenQ[logOnlyL] || logOnlyL < 2 ||
       !IntegerQ[logOnlyEnsembleSize] || logOnlyEnsembleSize < 1 ||
       !IntegerQ[logOnlySeed] ||
       !VectorQ[logOnlyHxExponents, srnFiniteRealQ] || logOnlyHxExponents === {} ||
       Length[DeleteDuplicates[logOnlyHxExponents]] != Length[logOnlyHxExponents] ||
       !AllTrue[{logOnlyJ, logOnlyHz, logOnlyLogTMin, logOnlyLogTMax,
           logOnlyPointsPerDecade}, srnFiniteRealQ] ||
       !TrueQ[logOnlyLogTMax > logOnlyLogTMin] || !TrueQ[logOnlyPointsPerDecade > 0] ||
       !IntegerQ[logOnlyMovingAverageWindow] || logOnlyMovingAverageWindow < 1 ||
       !VectorQ[logOnlyThresholds, srnFiniteRealQ] || logOnlyThresholds === {} ||
       !TrueQ[Min[logOnlyThresholds] > 0] ||
       Length[DeleteDuplicates[logOnlyThresholds]] != Length[logOnlyThresholds] ||
       !(logOnlyFitRanges === Automatic || (ListQ[logOnlyFitRanges] &&
           Length[logOnlyFitRanges] == Length[logOnlyThresholds] && AllTrue[logOnlyFitRanges, srnRangeQ])),
       Print["ERROR: invalid physical, time-grid, ensemble or fit parameters."]; Abort[]];
    If[!StringQ[logOnlyBaseOutputDirectory],
       Print["ERROR: set logOnlyBaseOutputDirectory to a directory string."]; Abort[]];
    logOnlyBaseOutputDirectory = ExpandFileName[logOnlyBaseOutputDirectory];
    logOnlyHxList = Developer`ToPackedArray[N[10^logOnlyHxExponents]];
    If[!VectorQ[logOnlyHxList, srnFiniteRealQ] || !TrueQ[Min[logOnlyHxList] > 0],
        Print["ERROR: hx values must be finite and positive."]; Abort[]];
    logOnlyStateSeeds = logOnlySeed + Range[0, logOnlyEnsembleSize - 1];
    logOnlyBatchSize = If[logOnlyBatchSizeSetting === Automatic,
        Which[logOnlyL <= 10, 64, logOnlyL <= 12, 32, True, 16], logOnlyBatchSizeSetting];
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

    If[!VectorQ[logOnlyTimeGrid, srnFiniteRealQ] || !TrueQ[Min[logOnlyTimeGrid] > 0] ||
       !TrueQ[Min[Differences[logOnlyTimeGrid]] > 0],
       Print["ERROR: logarithmic times overflow, underflow or are not strictly increasing."]; Abort[]];
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

(* Build the fresh configuration and metadata for this run. *)

    logOnlyConfig = <|
        "L" -> logOnlyL,
        "J" -> N[logOnlyJ],
        "hz" -> N[logOnlyHz],
        "hxExponents" ->
            Developer`ToPackedArray[N[logOnlyHxExponents]],
        "hxList" -> logOnlyHxList,
        "dA" -> logOnlyDA,
        "dB" -> logOnlyDB,
        "BpN" -> logOnlyBpN,
        "BmN" -> logOnlyBmN,
        "timeGrid" -> logOnlyTimeGrid,
        "nSamples" -> Length[logOnlyTimeGrid],
        "logTMin" -> N[logOnlyLogTMin],
        "logTMax" -> N[logOnlyLogTMax],
        "pointsPerDecade" -> logOnlyPointsPerDecade,
        "batchSize" -> logOnlyBatchSize,
        "stateCount" -> logOnlyEnsembleSize,
        "stateSeeds" -> logOnlyStateSeeds,
        "initialStates" -> logOnlyInitialStates,
        "initialStatesP" -> logOnlyInitialStatesP,
        "initialStatesM" -> logOnlyInitialStatesM,
        "outputDir" -> logOnlyOutputDir,
        "hermiticityTolerance" -> logOnlyHermiticityTolerance,
        "reflectionTolerance" -> logOnlyReflectionTolerance,
        "reconstructionTolerance" -> logOnlyReconstructionTolerance,
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
        "singleRPSIsEnsembleMember" -> 1,
        "L" -> logOnlyL,
        "J" -> N[logOnlyJ],
        "hz" -> N[logOnlyHz],
        "hxExponents" -> logOnlyHxExponents,
        "hxList" -> logOnlyHxList,
        "logTMin" -> logOnlyLogTMin,
        "logTMax" -> logOnlyLogTMax,
        "pointsPerDecadeParameter" -> logOnlyPointsPerDecade,
        "numberOfTimeSamples" -> Length[logOnlyTimeGrid],
        "timeGrid" -> logOnlyTimeGrid,
        "ensembleSize" -> logOnlyEnsembleSize,
        "stateSeeds" -> logOnlyStateSeeds,
        "initialStates" -> logOnlyInitialStates,
        "batchSize" -> logOnlyBatchSize,
        "requestedKernels" -> logOnlyRequestedKernels,
        "basisOrthogonalityError" ->
            logOnlyBasisOrthogonalityError,
        "WolframVersion" -> $Version,
        "SystemID" -> $SystemID,
        "movingAverageWindow" -> logOnlyMovingAverageWindow,
        "thresholds" -> logOnlyThresholds,
        "fitRanges" -> logOnlyFitRanges
    |>;

    (* Save metadata and the run description. *)

    logOnlyMetadataPath = FileNameJoin[{
        logOnlyOutputDir,
        "metadata.wxf"
    }];

    logOnlyMetadataExport = checkedExport[
        logOnlyMetadataPath,
        logOnlyMetadata,
        "WXF"
    ];

    If[FailureQ[logOnlyMetadataExport],
        Print[logOnlyMetadataExport];
        Abort[];
    ];

    logOnlyDescriptionExport = checkedExport[
        FileNameJoin[{logOnlyOutputDir, "description.txt"}],
        srnDescription[
            logOnlyMetadata,
            logOnlyMovingAverageWindow,
            logOnlyThresholds,
            logOnlyFitRanges
        ],
        "Text"
    ];

    If[FailureQ[logOnlyDescriptionExport],
        Print[logOnlyDescriptionExport];
        Abort[];
    ];

    (* Initialize or reuse the parallel kernels. *)

    logOnlyStartup = logOnlyInitializeWorkers[];

    If[FailureQ[logOnlyStartup],
        Print[logOnlyStartup];
        Abort[];
    ];

    (* Compute the single-RPS curve and ensemble mean at each hx. *)

    logOnlySweepSeconds = First[
        AbsoluteTiming[
            logOnlyRunSummary = With[
                {snapshot = logOnlyConfig},
                If[TrueQ[logOnlyParallelActive],
                    ParallelMap[
                        logOnlyRunOneHx[snapshot, #] &,
                        Range[Length[snapshot["hxList"]]],
                        Method -> "CoarsestGrained",
                        DistributedContexts -> None
                    ],
                    Map[
                        logOnlyRunOneHx[snapshot, #] &,
                        Range[Length[snapshot["hxList"]]]
                    ]
                ]
            ];
        ]
    ];

    (* Check and save the sweep status. *)

    logOnlyFailedRuns = If[
        ListQ[logOnlyRunSummary],
        Select[
            logOnlyRunSummary,
            Function[result,
                If[AssociationQ[result],
                    !TrueQ[Lookup[result, "success", False]],
                    True
                ]
            ]
        ],
        {logOnlyRunSummary}
    ];

    logOnlySummaryExport = checkedExport[
        FileNameJoin[{logOnlyOutputDir, "run_summary.wxf"}],
        logOnlyRunSummary,
        "WXF"
    ];

    If[FailureQ[logOnlySummaryExport],
        Print[logOnlySummaryExport];
        Abort[];
    ];

    If[logOnlyFailedRuns =!= {},
        logOnlyFailureExport = checkedExport[
            FileNameJoin[{logOnlyOutputDir, "failed_runs.wxf"}],
            logOnlyFailedRuns,
            "WXF"
        ];

        If[FailureQ[logOnlyFailureExport],
            Print[logOnlyFailureExport];
        ];

        Print["ERROR: one or more logarithmic hx runs failed."];
        Print[logOnlyFailedRuns];
        Abort[];
    ];

    Print[
        "Log-only sweep complete. Wall time = ",
        logOnlySweepSeconds,
        " s."
    ];

    (* Extract power laws from moving-average data and export.
       srnAnalyzeSaved saves its numerical results before rendering. *)

    logOnlyPowerResults = srnAnalyzeSaved[
        logOnlyOutputDir,
        logOnlyMovingAverageWindow,
        logOnlyThresholds,
        logOnlyFitRanges
    ];

    If[FailureQ[logOnlyPowerResults],
        Print[logOnlyPowerResults];
    ];

    (* Export entropy plots, retaining the hx colorbars. *)

    logOnlyPlotFiles = logOnlyPlotResults[
        logOnlyOutputDir,
        logOnlyMovingAverageWindow
    ];

    If[FailureQ[logOnlyPlotFiles],
        Print[logOnlyPlotFiles];
    ];

    Print["Run directory: ", logOnlyOutputDir];

    logOnlyOutputDir
];



(* ::Section:: *)
(* PARAMETERS: edit and reevaluate this section for the next run *)


(* ::Input:: *)
(**)


logOnlyL = 8;
logOnlyJ = 1.;
logOnlyHz = 2.;

(* h_x = 10^alpha. Preserve your attached file's chosen sweep.
   For finer SRN sampling, e.g. N[Range[-4, -1, 1/5]]. *)
logOnlyHxExponents = N[Range[-8, 0, 5/10]];

(* These are BASE-10 EXPONENTS: logOnlyLogTMax = 10 means tmax = 10^10.
   A log grid has a variable physical timestep. Increase points/decade to
   decrease the step in log10(t), preserving the original endpoint convention. *)
logOnlyLogTMin = -1.;
logOnlyLogTMax = 10.;
logOnlyPointsPerDecade = 2000;

logOnlyEnsembleSize = 40;  (* use 1 for one RPS; single and mean then coincide *)
logOnlySeed = 12345;
logOnlyBatchSizeSetting = Automatic;
logOnlyUseParallel = True;
logOnlyRequestedKernels = 8;
logOnlyHermiticityTolerance = 10^-12;
logOnlyReflectionTolerance = 10^-12;
logOnlyReconstructionTolerance = 10^-10;
logOnlyEigensystemResidualTolerance = 10^-10;
logOnlyValidationEntropyTolerance = 10^-8;
logOnlyValidationNormTolerance = 10^-10;

logOnlyMovingAverageWindow = 200;
logOnlyThresholds = {0.6, 1.4};
(* Automatic fits all uncensored fields for each threshold. No data-dependent
   trimming and no forced exponent. For an SRN-only range at both thresholds:
       logOnlyFitRanges = {{10^-4, 10^-1}, {10^-4, 10^-1}};
   Adjust ranges after inspecting crossing plots and local slopes. *)
logOnlyFitRanges = Automatic;

(* Absolute default is the script/notebook folder. A short custom root such as
   "C:\\Entropy" is also valid. Names contain size/fields and a unique run ID,
   never time arrays, log limits or timesteps. *)
logOnlyBaseOutputDirectory = logOnlyDefaultDirectory;



(* ::Section::Closed:: *)
(* EXECUTE: rerun this line after editing/evaluating the parameters above *)


(* ::Input:: *)
(* (*This cell calls the complete function defined above. *)*)


If[Length[DownValues[logOnlyRunAnalysis]] == 0,
    Print["ERROR: run definition is missing. Evaluate the definitions above first."];
    Abort[];
];

logOnlyLastRun = logOnlyRunAnalysis[];



(* ::Section::Closed:: *)
(* Optional replot/refit from saved files, without recomputing the dynamics *)


(* ::Input:: *)
(**)


(*
   logOnlyPlotResults[logOnlyOutputDir, 200];
   logOnlyPowerResults = srnAnalyzeSaved[logOnlyOutputDir, 200, {0.6, 1.4}, Automatic];

   Or, in a later session, after loading the definitions only:
   savedRun = "C:\\Entropy\\entropy_SRN\\L8_J1_hz2_hxm4to0_Nens50\\run_...";
   logOnlyPlotResults[savedRun, 200];
   srnAnalyzeSaved[savedRun, 200, {0.6, 1.4}, {{10^-4, 10^-1}, {10^-4, 10^-1}}];

   Files in each run:
     description.txt, metadata.wxf, run_summary.wxf, log_run_01.wxf, ...
     plots/ : four entropy plots, including their hx colorbars.
     power_<id>/ : moving_average.wxf, power_laws.wxf, fits.csv,
       crossings.csv, analysis.txt, single_S0p6.png/pdf, single_S1p4.png/pdf,
       ensemble_S0p6.png/pdf, ensemble_S1p4.png/pdf.
   A new short power_<id> folder preserves each reanalysis independently.
   Censored fields and fields outside the selected fit range remain in CSV/WXF.
   Data files are saved even when too few crossings exist to fit a power law.
   No CloseKernels[], Quit[] or ClearAll["Global`*"] is used.
*)


(* ::Section:: *)
(* Single RPS: full spectra and four-panel collage *)


(* ::Subsection::Closed:: *)
(*defs*)


(* Paste this ENTIRE section after the definitions in the repaired SRN file.
   It uses an existing successful run folder, not the current physical parameters.
   The first saved RPS (ensemble member 1) is reused exactly; no new state is drawn.
   Existing dependencies: IsingHamiltonian, reflectionSectorBases,
   sortedEigensystem, phaseAtTimesBatchC, sampledEigensystemResidual,
   startupMovingAverage, checkedExport.

   lambda_i = s_i^2 are the Schmidt probabilities (rho_A eigenvalues), sorted
   decreasingly at EACH time. xi_i = -Log[lambda_i] = -2 Log[s_i].
   All 2^(L/2) levels are retained, with no singular-value truncation.
   The saved s_i also allow plotting Schmidt coefficients instead of probabilities.
   Level labels denote instantaneous rank, not tracked Schmidt eigenvectors.
   P(t) = Abs[Conjugate[psi0].psi(t)]^2 is survival under the same Hamiltonian.

   Every observable uses the SAME trailing moving average with growing startup
   windows and unchanged log-time coordinates, exactly as in the main file.
   Energies are averaged AFTER taking -Log, and survival AFTER taking Abs[...]^2.
   Smoothed entropy is recomputed alongside the spectra and smoothed directly; it is not
   recomputed as the entropy of an averaged spectrum or an averaged state.
   Smoothing uses the full saved time grid BEFORE applying any display zoom.

   True zero Schmidt coefficients have infinite entanglement energy. They are
   retained as Infinity in WXF; a window containing Infinity stays infinite.
   Infinite values cannot be drawn on a finite plot and are explicitly counted.
   Tiny nonzero machine-precision levels may be dominated by numerical roundoff;
   keeping all levels does not increase their accuracy. There is no artificial
   positive floor, no renormalization, no interpolation and no sampling reduction.
   See https://reference.wolfram.com/language/ref/SingularValueList.html

   One selected hx gives one collage with this layout:
      top left: entropy; top right: full Schmidt spectrum;
      bottom left: full entanglement spectrum; bottom right: survival probability.
   The panels share logarithmic x scale, x range/ticks, linear y scale, point style,
   frame, padding, font, aspect ratio and moving-average window. Y ranges reflect
   the different quantities. The spectra share a rank colorbar, not an hx legend.
   Each call creates a short unique spectra_<id> folder and saves raw data first.
   Validation: whole-file and cell syntax passed in Wolfram. A four-spin test
   agreed with direct MatrixExp evolution to 1.7e-13 for entropy and 1.0e-13
   for survival. A product-state/zero-H test preserved {1,0}, {0,Infinity} and
   P=1 and verified saved-data replotting with a different window. Plot objects
   were constructed; local QMB, native C and PNG rendering were not run here.
   The original saved entropy is retained as savedEntropy (raw and smoothed).
   A mismatch is printed, recorded in WXF/summary and noted in the figure title;
   it no longer prevents saving and plotting the new, internally consistent data.
   It does not prove the origin of the mismatch or long-time numerical accuracy.
   No existing definitions, parameter settings, kernels or result files are changed.
*)

ClearAll[srnSpecFiniteQ, srnSpecRequire, srnSpecDirectory, srnSpecExport,
    srnSpecAverage, srnSpecCompute, srnSpecSmooth, srnSpecPlots,
    srnSpecReplot, srnSingleSpectra];

srnSpecFiniteQ[x_] := NumberQ[x] && TrueQ[Im[x] == 0] &&
    FreeQ[x, Indeterminate | _DirectedInfinity];
srnSpecRequire[test_, tag_String, details_: <||>] :=
    If[!TrueQ[test], Throw[Failure[tag, details], "srnSpectraFailure"]];
srnSpecExport[path_String, data_, format_String] := Module[{result},
    result = checkedExport[path, data, format];
    If[FailureQ[result], Throw[result, "srnSpectraFailure"]]; result
];
srnSpecDirectory[parent_String, prefix_String] := Module[{out, result},
    out = FileNameJoin[{parent, prefix <> "_" <>
        StringTake[StringDelete[CreateUUID[], "-"], 8]}];
    result = Quiet[Check[CreateDirectory[out], $Failed]];
    srnSpecRequire[StringQ[result] && DirectoryQ[out], "CreateSpectraDirectoryFailed", <|"path" -> out|>];
    out
];

(* The finite case calls the original smoother without modification.
   Infinity handling preserves the arithmetic mean on each exact same window. *)
srnSpecAverage[values_List, window_Integer] := Module[{n, bad, counts, finiteMean},
    If[AllTrue[values, srnSpecFiniteQ], Return[startupMovingAverage[values, window]]];
    srnSpecRequire[AllTrue[values, srnSpecFiniteQ[#] || # === Infinity &], "InvalidSpectrumValues"];
    n = Length[values];
    bad = Accumulate[(If[# === Infinity, 1, 0] &) /@ values];
    counts = bad - PadLeft[Take[bad, Max[0, n - window]], n];
    finiteMean = startupMovingAverage[Replace[values, Infinity -> 0., {1}], window];
    MapThread[If[#2 > 0, Infinity, #1] &, {finiteMean, counts}]
];

(* Recompute ONE state, with the same reflection sectors and batched absolute
   phases as the successful entropy sweep. Only one time batch of states is held. *)
srnSpecCompute[run_Association, psi0_List, bases_Association, batch_Integer] := Catch[Module[
    {L = run["L"], times = run["timeGrid"], dA, n, bp, bm, h, scale,
     hp, hm, ep, em, vp, vm, wp, wm, cp, cm, residual, recon,
     p, m, re, im, sv, coefficients, probabilities, energies, survival, entropyCheck,
     pos, len, j, row, tchunk, normError = 0., traceError = 0., difference,
     start = AbsoluteTime[], tolerance = 10^-8},
    dA = 2^(L/2); n = Length[times];
    bp = N[bases["Even"]]; bm = N[bases["Odd"]];
    srnSpecRequire[Length[psi0] == 2^L && Abs[Norm[psi0]^2 - 1.] <= 10^-10,
        "InvalidSavedInitialState"];
    h = N[IsingHamiltonian[run["hx"], run["hz"], run["J"], L]];
    srnSpecRequire[Dimensions[h] === {2^L, 2^L}, "InvalidHamiltonian"];
    scale = Max[1., Norm[h, "Frobenius"]];
    srnSpecRequire[Norm[h - ConjugateTranspose[h], "Frobenius"]/scale <= 10^-12,
        "NonHermitianHamiltonian"];
    srnSpecRequire[Norm[ConjugateTranspose[bp] . h . bm, "Frobenius"]/scale <= 10^-12,
        "ReflectionSymmetryFailed"];
    hp = ConjugateTranspose[bp] . h . bp; hm = ConjugateTranspose[bm] . h . bm;
    {ep, vp} = sortedEigensystem[hp]; {em, vm} = sortedEigensystem[hm];
    residual = Max[sampledEigensystemResidual[hp, ep, vp], sampledEigensystemResidual[hm, em, vm]];
    srnSpecRequire[residual <= 10^-10, "EigensystemResidualFailed", <|"error" -> residual|>];
    srnSpecRequire[Max[Abs[Im[vp]]] <= 10^-10 && Max[Abs[Im[vm]]] <= 10^-10,
        "ComplexSectorEigenbasis"];
    wp = Developer`ToPackedArray[Re[Transpose[vp]]];
    wm = Developer`ToPackedArray[Re[Transpose[vm]]];
    cp = Conjugate[vp] . (Transpose[bp] . psi0);
    cm = Conjugate[vm] . (Transpose[bm] . psi0);
    recon = Norm[bp . (wp . cp) + bm . (wm . cm) - psi0];
    srnSpecRequire[recon <= 10^-10, "InitialStateReconstructionFailed", <|"error" -> recon|>];
    coefficients = ConstantArray[0., {n, dA}];
    probabilities = ConstantArray[0., {n, dA}];
    energies = ConstantArray[0., {n, dA}];
    survival = ConstantArray[0., n]; entropyCheck = ConstantArray[0., n];
    pos = 1;
    While[pos <= n,
        len = Min[batch, n - pos + 1];
        tchunk = Developer`ToPackedArray[N[times[[pos ;; pos + len - 1]]]];
        p = phaseAtTimesBatchC[Developer`ToPackedArray[Re[cp]], Developer`ToPackedArray[Im[cp]], ep, tchunk];
        m = phaseAtTimesBatchC[Developer`ToPackedArray[Re[cm]], Developer`ToPackedArray[Im[cm]], em, tchunk];
        re = bp . (wp . p[[1]]) + bm . (wm . m[[1]]);
        im = bp . (wp . p[[2]]) + bm . (wm . m[[2]]);
        normError = Max[normError, Max[Abs[Total[re re + im im] - 1.]]];
        Do[
            row = pos + j - 1;
            sv = SingularValueList[ArrayReshape[re[[All, j]] + I im[[All, j]], {dA, dA}], Tolerance -> 0];
            srnSpecRequire[Length[sv] == dA && VectorQ[sv, srnSpecFiniteQ] && Min[sv] >= 0,
                "IncompleteSchmidtSpectrum", <|"timeIndex" -> row|>];
            coefficients[[row]] = sv;
            probabilities[[row]] = sv^2;
            energies[[row]] = (If[# > 0., -2. Log[#], Infinity] &) /@ sv;
            survival[[row]] = Abs[Conjugate[psi0] . (re[[All, j]] + I im[[All, j]])]^2;
            (* Same probability cutoff as entropyFromSVC in the supplied workflow. *)
            entropyCheck[[row]] = -Total[(If[# > 1.*^-14, # Log[#], 0.] &) /@ probabilities[[row]]];
            traceError = Max[traceError, Abs[Total[probabilities[[row]]] - 1.]];
        , {j, len}];
        pos += len;
    ];
    srnSpecRequire[normError <= 10^-10 && traceError <= 10^-10,
        "SpectrumNormalizationFailed", <|"normError" -> normError, "traceError" -> traceError|>];
    (* A comparison between two separately diagonalized evolutions is a
       diagnostic, not an internal consistency test for the new spectrum.
       At long times, eigenvalue/eigenvector roundoff can affect this comparison.
       Retain the discrepancy without asserting its cause or relaxing norm checks. *)
    difference = Max[Abs[entropyCheck - run["singleEntropy"]]];
    If[difference > tolerance,
        Print["Saved entropy comparison: max |new - saved| = ", difference,
            ", at t = ", times[[First[Ordering[Abs[entropyCheck - run["singleEntropy"]], -1]]]],
            ". Saving both curves. The collage uses entropy recomputed with the spectra.",
            " This comparison does not establish the accuracy of the smallest levels."]
    ];
    <|"L" -> L, "J" -> run["J"], "hz" -> run["hz"], "hx" -> run["hx"],
      "hxExponent" -> run["hxExponent"], "singleStateSeed" -> Lookup[run, "singleStateSeed", Missing["NotSaved"]],
      "timeGrid" -> times, "initialState" -> psi0,
      "entropy" -> Developer`ToPackedArray[entropyCheck],
      "savedEntropy" -> run["singleEntropy"],
      "entropySource" -> "Recomputed from the same state evolution as the spectra and survival",
      "savedEntropyAgreement" -> TrueQ[difference <= tolerance],
      "savedEntropyComparisonTolerance" -> tolerance,
      "savedEntropyDifferenceTime" -> times[[First[Ordering[Abs[entropyCheck - run["singleEntropy"]], -1]]]],
      "schmidtProbabilities" -> Developer`ToPackedArray[probabilities],
      "schmidtCoefficients" -> Developer`ToPackedArray[coefficients],
      "entanglementEnergies" -> energies, "survival" -> Developer`ToPackedArray[survival],
      "levelCount" -> dA, "rankConvention" -> "Descending probability at each time",
      "maxNormError" -> normError, "maxTraceError" -> traceError,
      "savedEntropyMaxDifference" -> difference, "eigensystemResidual" -> residual,
      "computeSeconds" -> N[AbsoluteTime[] - start]|>
], "srnSpectraFailure"];

srnSpecSmooth[data_Association, window_Integer] := Module[{},
    srnSpecRequire[window >= 1, "InvalidSpectrumWindow"];
    <|"window" -> window, "entropy" -> srnSpecAverage[data["entropy"], window],
      "savedEntropy" -> srnSpecAverage[Lookup[data, "savedEntropy", data["entropy"]], window],
      "schmidtProbabilities" -> Transpose[srnSpecAverage[#, window] & /@ Transpose[data["schmidtProbabilities"]]],
      "schmidtCoefficients" -> Transpose[srnSpecAverage[#, window] & /@ Transpose[data["schmidtCoefficients"]]],
      "entanglementEnergies" -> Transpose[srnSpecAverage[#, window] & /@ Transpose[data["entanglementEnergies"]]],
      "survival" -> srnSpecAverage[data["survival"], window]|>
];



(* ::Subsection::Closed:: *)
(* Common plotting function: identical time coordinates and display settings *)


(* ::Input:: *)
(**)


srnSpecPlots[data_Association, smooth_Association, timeWindow_: Automatic,
    schmidtQuantity_String: "Probabilities"] := Catch[Module[
    {times = data["timeGrid"], bounds, idx, t, n = data["levelCount"], cf, colors,
     ticks, tickStep, common, make, matrixData, entropyPlot, schmidtPlot,
     energyPlot, survivalPlot, bar, title, collage, schmidtKey, schmidtLabel, infiniteCount},
    bounds = If[timeWindow === Automatic, {First[times], Last[times]}, N[timeWindow]];
    srnSpecRequire[MatchQ[bounds, {_?srnSpecFiniteQ, _?srnSpecFiniteQ}] &&
        TrueQ[0 < bounds[[1]] < bounds[[2]]] && bounds[[1]] >= First[times] && bounds[[2]] <= Last[times],
        "InvalidSpectrumTimeWindow", <|"available" -> {First[times], Last[times]}|>];
    srnSpecRequire[MemberQ[{"Probabilities", "Coefficients"}, schmidtQuantity], "InvalidSchmidtQuantity"];
    idx = Flatten[Position[times, _?(TrueQ[bounds[[1]] <= # <= bounds[[2]]] &)]];
    srnSpecRequire[Length[idx] >= 2, "TooFewVisibleTimeSamples"];
    t = times[[idx]];
    cf = Blend[{Blue, Purple, Red}, #] &;
    colors = cf /@ Subdivide[0., 1., n - 1];
    tickStep = Max[1, Ceiling[(Log10[bounds[[2]]] - Log10[bounds[[1]]])/5]];
    ticks = Table[{10.^q, Style[Superscript[10, q], 16, Black]},
        {q, Ceiling[Log10[bounds[[1]]]], Floor[Log10[bounds[[2]]]], tickStep}];
    If[ticks === {}, ticks = ({#, ScientificForm[#, 2]} &) /@ bounds];
    common = {Joined -> False, Frame -> True, Axes -> False,
        PlotRange -> {bounds, {0., All}}, PlotRangePadding -> {None, Scaled[.04]},
        FrameTicks -> {{Automatic, None}, {ticks, None}},
        FrameStyle -> Directive[Black, 16], BaseStyle -> {FontFamily -> "Times", FontSize -> 16},
        ImageSize -> 550, AspectRatio -> 2/3, ImagePadding -> {{72, 20}, {55, 42}},
        PlotTheme -> "Detailed", GridLines -> None, Background -> White};
    make[coords_, ylabel_, caption_, styles_] := ListLogLinearPlot[coords,
        Evaluate[Sequence @@ common], PlotStyle -> styles,
        FrameLabel -> {Style["t", 20, Italic], Style[ylabel, 20]},
        PlotLabel -> Style[caption, 18, Black]];
    matrixData[mat_] := Table[
        Select[Transpose[{t, mat[[idx, rank]]}], srnSpecFiniteQ[Last[#]] &], {rank, n}];
    entropyPlot = make[{Transpose[{t, smooth["entropy"][[idx]]}]}, "S(t)",
        "(a) Entanglement entropy", {Directive[Black, Opacity[.82], PointSize[.002]]}];
    schmidtKey = If[schmidtQuantity === "Probabilities", "schmidtProbabilities", "schmidtCoefficients"];
    schmidtLabel = If[schmidtQuantity === "Probabilities", Subscript["\[Lambda]", "i"], Subscript["s", "i"]];
    schmidtPlot = make[matrixData[smooth[schmidtKey]], schmidtLabel,
        "(b) Full Schmidt spectrum", (Directive[#, Opacity[.82], PointSize[.002]] &) /@ colors];
    infiniteCount = Count[smooth["entanglementEnergies"][[idx]], Infinity, {2}];
    energyPlot = make[matrixData[smooth["entanglementEnergies"]], Subscript["\[Xi]", "i"],
        If[infiniteCount == 0, "(c) Full entanglement spectrum",
            "(c) Entanglement spectrum (infinite values omitted)"],
        (Directive[#, Opacity[.82], PointSize[.002]] &) /@ colors];
    survivalPlot = make[{Transpose[{t, smooth["survival"][[idx]]}]}, "P(t)",
        "(d) Survival probability", {Directive[Black, Opacity[.82], PointSize[.002]]}];
    bar = BarLegend[{cf, {1, n}}, ColorFunctionScaling -> True,
        LegendLabel -> Style["Spectrum rank i", 16, Black],
        Ticks -> ({#, #} &) /@ DeleteDuplicates[Round[Subdivide[1, n, Min[5, n - 1]]]],
        LegendMarkerSize -> {12, 240}, LabelStyle -> {Black, 14}];
    title = Column[{
        Style[Row[{"Single RPS,  L = ", data["L"], ",  J = ", data["J"],
            ",  ", Subscript["h", "z"], " = ", data["hz"], ",  ",
            Subscript["h", "x"], " = ", Superscript[10, data["hxExponent"]]}], 22, Black],
        Style[Row[{"Moving average: ", smooth["window"], " samples; ", n, " spectrum levels"}], 17, Black],
        If[TrueQ[Lookup[data, "savedEntropyAgreement", True]], Nothing,
            Style[Row[{"Entropy recomputed with spectra; max difference from saved curve = ",
                ScientificForm[data["savedEntropyMaxDifference"], 3]}], 14, Black]]
        }, Alignment -> Center];
    collage = Labeled[Legended[GraphicsGrid[{{entropyPlot, schmidtPlot}, {energyPlot, survivalPlot}},
        Spacings -> {.35, .5}, ImageSize -> 1140, Background -> White], Placed[bar, Right]], title, Top];
    <|"collage" -> collage,
      "schmidt" -> Labeled[Legended[schmidtPlot, Placed[bar, Right]], title, Top],
      "entanglement" -> Labeled[Legended[energyPlot, Placed[bar, Right]], title, Top],
      "panels" -> {entropyPlot, schmidtPlot, energyPlot, survivalPlot},
      "timeWindow" -> bounds, "infiniteEnergyValuesInWindow" -> infiniteCount|>
], "srnSpectraFailure"];

(* Replot a saved spec_XX.wxf; no QMB, eigensystem or evolution is used here. *)
srnSpecReplot[file_String, window_Integer: 200, timeWindow_: Automatic,
    schmidtQuantity_String: "Probabilities", display_: True] := Catch[Module[
    {data, smooth, plots, out, files},
    data = Quiet[Check[Import[file, "WXF"], $Failed]];
    srnSpecRequire[AssociationQ[data] && AllTrue[{"timeGrid", "entropy", "schmidtProbabilities",
        "schmidtCoefficients", "entanglementEnergies", "survival", "levelCount"}, KeyExistsQ[data, #] &],
        "InvalidSavedSpectra", <|"path" -> file|>];
    smooth = srnSpecSmooth[data, window];
    plots = srnSpecPlots[data, smooth, timeWindow, schmidtQuantity];
    If[FailureQ[plots], Throw[plots, "srnSpectraFailure"]];
    out = srnSpecDirectory[DirectoryName[file], FileBaseName[file] <> "_plots"];
    srnSpecExport[FileNameJoin[{out, "moving_average.wxf"}],
        <|"sourceFile" -> file, "timeGrid" -> data["timeGrid"], "smoothed" -> smooth,
          "timeWindow" -> plots["timeWindow"], "schmidtQuantity" -> schmidtQuantity|>, "WXF"];
    files = Table[srnSpecExport[FileNameJoin[{out, key <> ".png"}], plots[key], "PNG"],
        {key, {"schmidt", "entanglement", "collage"}}];
    If[TrueQ[display], Print[plots["collage"]]];
    Print["Spectra plots saved in: ", out];
    <|"sourceFile" -> file, "plotDirectory" -> out, "plotFiles" -> files,
      "infiniteEnergyValuesInWindow" -> plots["infiniteEnergyValuesInWindow"]|>
], "srnSpectraFailure"];



(* ::Subsection::Closed:: *)
(* Complete callable computation and export function *)


(* ::Input:: *)
(**)


srnSingleSpectra[directory_String, fieldIndices_: {1}, window_Integer: 200,
    timeWindow_: Automatic, schmidtQuantity_String: "Probabilities", display_: True] := Catch[Module[
    {meta, files, ordered, indices, psi0, bases, out, batch, run, result, path,
     summaries = {}, plots, notes, k},
    srnSpecRequire[DirectoryQ[directory], "RunDirectoryNotFound", <|"path" -> directory|>];
    srnSpecRequire[window >= 1, "InvalidSpectrumWindow"];
    meta = Quiet[Check[Import[FileNameJoin[{directory, "metadata.wxf"}], "WXF"], $Failed]];
    srnSpecRequire[AssociationQ[meta] && AllTrue[{"L", "J", "hz", "initialStates", "hxList", "timeGrid"},
        KeyExistsQ[meta, #] &], "MissingRunMetadata"];
    srnSpecRequire[EvenQ[meta["L"]] && meta["L"] >= 2 && ListQ[meta["initialStates"]] &&
        Length[meta["initialStates"]] >= 1, "InvalidRunMetadata"];
    files = FileNames["log_run_" ~~ DigitCharacter .. ~~ ".wxf", directory];
    ordered = SortBy[files, ToExpression[StringReplace[FileBaseName[#], "log_run_" -> ""]] &];
    srnSpecRequire[Length[ordered] > 0, "NoSavedLogRuns"];
    indices = If[fieldIndices === All, Range[Length[ordered]], fieldIndices];
    srnSpecRequire[ListQ[indices] && indices =!= {} && VectorQ[indices, IntegerQ] &&
        Min[indices] >= 1 && Max[indices] <= Length[ordered], "InvalidFieldIndices",
        <|"availableCount" -> Length[ordered]|>];
    indices = DeleteDuplicates[indices];
    psi0 = Developer`ToPackedArray[N[First[meta["initialStates"]]]];
    bases = reflectionSectorBases[meta["L"]];
    batch = Lookup[meta, "batchSize", 32];
    srnSpecRequire[IntegerQ[batch] && batch >= 1, "InvalidSpectrumBatchSize"];
    out = srnSpecDirectory[directory, "spectra"];
    notes = StringRiffle[{
        "Single RPS spectra and survival: saved ensemble member 1, no new random state.",
        "Physical parameters and time grid are read from the saved run, not live globals.",
        "Schmidt probabilities lambda=s^2; entanglement energies xi=-2 Log[s].",
        "All levels retained with Tolerance->0, sorted by descending lambda at each time.",
        "Rank is not an eigenvector label. Tiny machine-precision levels may be inaccurate.",
        "Survival P(t)=Abs[Conjugate[psi0].psi(t)]^2 under the same Hamiltonian.",
        "Same trailing sample moving average and growing startup as the main entropy plots.",
        "Moving-average window: " <> ToString[window],
        "Requested display time window: " <> ToString[timeWindow, InputForm],
        "Smoothing occurs before display zoom; time coordinates are never averaged.",
        "Energy smoothing is Mean[-Log[lambda]], not -Log[Mean[lambda]].",
        "Zero Schmidt coefficients imply Infinity; affected energy averages stay infinite.",
        "Raw WXF data are saved before plots; infinity is omitted only for drawing.",
        "Entropy is recomputed alongside spectra; the existing saved entropy is retained separately.",
        "The saved-entropy comparison is recorded and reported, not used to abort the new computation.",
        "Norm, trace, Hermiticity, reflection symmetry and reconstruction checks remain enforced.",
        "Use srnSpecReplot on any spec_XX.wxf to change window/zoom without dynamics."
    }, "\n"];
    srnSpecExport[FileNameJoin[{out, "description.txt"}], notes, "Text"];
    Do[
        Print["Computing single-RPS spectra: selected field ", k, " of ", Length[ordered], "."];
        run = Quiet[Check[Import[ordered[[k]], "WXF"], $Failed]];
        srnSpecRequire[AssociationQ[run] && TrueQ[Lookup[run, "success", False]] &&
            AllTrue[{"L", "J", "hz", "hx", "hxExponent", "timeGrid", "singleEntropy"}, KeyExistsQ[run, #] &],
            "InvalidSavedSingleRun", <|"path" -> ordered[[k]]|>];
        srnSpecRequire[Lookup[run, {"L", "J", "hz"}] === Lookup[meta, {"L", "J", "hz"}] &&
            run["timeGrid"] === meta["timeGrid"] &&
            VectorQ[run["timeGrid"], srnSpecFiniteQ] && Length[run["timeGrid"]] >= 2 &&
            Min[run["timeGrid"]] > 0 && Min[Differences[run["timeGrid"]]] > 0 &&
            VectorQ[run["singleEntropy"], srnSpecFiniteQ] &&
            Length[run["singleEntropy"]] == Length[run["timeGrid"]], "InconsistentSavedRun"];
        result = srnSpecCompute[run, psi0, bases, batch];
        If[FailureQ[result], Throw[result, "srnSpectraFailure"]];
        path = FileNameJoin[{out, "spec_" <> IntegerString[k, 10, 2] <> ".wxf"}];
        srnSpecExport[path, result, "WXF"];
        plots = srnSpecReplot[path, window, timeWindow, schmidtQuantity, display];
        If[FailureQ[plots], Throw[plots, "srnSpectraFailure"]];
        AppendTo[summaries, Join[<|"index" -> k, "hx" -> result["hx"],
            "computeSeconds" -> result["computeSeconds"],
            "savedEntropyMaxDifference" -> result["savedEntropyMaxDifference"],
            "savedEntropyAgreement" -> result["savedEntropyAgreement"],
            "savedEntropyDifferenceTime" -> result["savedEntropyDifferenceTime"]|>, plots]];
        srnSpecExport[FileNameJoin[{out, "summary.wxf"}], summaries, "WXF"];
    , {k, indices}];
    Print["Single-RPS spectra complete. Data folder: ", out];
    <|"outputDirectory" -> out, "fields" -> summaries|>
], "srnSpectraFailure"];



(* ::Subsection:: *)
(* Parameters for THIS section only *)


(* ::Input:: *)
(**)


srnSpectrumSourceDirectory = logOnlyOutputDir;  (* or an existing run folder string *)
srnSpectrumFieldIndices = Range[17];                 (* {1, 5, 10}, or All for every saved hx *)
srnSpectrumWindow = logOnlyMovingAverageWindow;
srnSpectrumTimeWindow = Automatic;             (* full saved interval; or {10., 10.^6} *)
srnSpectrumSchmidtQuantity = "Probabilities";   (* or "Coefficients" to show s_i *)
srnSpectrumDisplay = True;                     (* also print the collage in the notebook *)



(* ::Subsection:: *)
(* Execute this section: existing entropy sweep is NOT run again *)


(* ::Input:: *)
(**)


If[Length[DownValues[srnSingleSpectra]] == 0,
    Print["ERROR: evaluate all definitions in this section first."]; Abort[]];
If[!StringQ[srnSpectrumSourceDirectory],
    Print["ERROR: set srnSpectrumSourceDirectory to a completed run folder."]; Abort[]];
srnSpectrumResult = srnSingleSpectra[srnSpectrumSourceDirectory,
    srnSpectrumFieldIndices, srnSpectrumWindow, srnSpectrumTimeWindow,
    srnSpectrumSchmidtQuantity, srnSpectrumDisplay];
If[FailureQ[srnSpectrumResult], Print[srnSpectrumResult]];

(* Replot later without evolving any state:
   srnSpecReplot["C:\\...\\spectra_abcd1234\\spec_01.wxf", 200, {1.,10.^6}];
   To use the first result immediately:
   srnSpecReplot[srnSpectrumResult["fields"][[1]]["sourceFile"],
       srnSpectrumWindow, srnSpectrumTimeWindow, "Probabilities", True];
*)


(* ::Section:: *)
(* Single RPS: entropy, Loschmidt echo, H0-basis IPR and survival collage *)


(* ::Subsection::Closed:: *)
(* Definitions *)


(* Paste this ENTIRE section after the general definitions in the repaired SRN file.

   It uses an EXISTING successful run folder. It does NOT rerun the ensemble sweep
   and it does NOT draw a new random state. The first saved RPS (ensemble member 1)
   is reused exactly.

   The perturbed trajectory is
       psi_h(t) = Exp[-I H(hx) t] psi(0)

   and the reference/integrable trajectory is
       psi_0(t) = Exp[-I H0 t] psi(0),
       H0 = IsingHamiltonian[0, hz, J, L].

   Loschmidt echo / trajectory fidelity:
       L(t) = |<psi_0(t)|psi_h(t)>|^2.

   IPR is evaluated in the eigenbasis of H0:
       IPR_0(t) = Sum_n |<E_n^(0)|psi_h(t)>|^4.
   Thus the IPR is constant for an H0 evolution, while the perturbation can change it.
   It is an eigenvector-basis quantity: inside an exactly degenerate H0 eigenspace
   the numerical eigenbasis returned by Eigensystem fixes the convention.

   Survival probability:
       P(t) = |<psi(0)|psi_h(t)>|^2.

   Entanglement entropy is recomputed from the SAME perturbed state used for the
   other three diagnostics and compared with the entropy already saved in log_run_XX.wxf.

   ALL four observables use exactly the same trailing moving average already used
   by the main file:
       average at sample i = Mean[data[[Max[1,i-window+1] ;; i]]].
   The time coordinates are never averaged. Smoothing is done on the complete
   saved logarithmic grid BEFORE an optional display zoom.

   The four panels use the same plotting convention as the previous spectra collage:
      top left     : entanglement entropy
      top right    : Loschmidt echo
      bottom left  : H0-eigenbasis IPR
      bottom right : survival probability

   All panels have the same logarithmic time axis, identical visible time window,
   ticks, point style, frame, image padding, aspect ratio, font and moving-average
   window. Their y ranges are linear and allowed to adapt to the observable.

   Existing dependencies from the main file:
       IsingHamiltonian, reflectionSectorBases, sortedEigensystem,
       phaseAtTimesBatchC, sampledEigensystemResidual, entropyFast,
       startupMovingAverage, checkedExport.

   Raw numerical data are saved before plotting. Replotting with a different
   moving-average window or time window does not rerun any dynamics.
*)

ClearAll[
    srnGlobalFiniteQ,
    srnGlobalRequire,
    srnGlobalExport,
    srnGlobalDirectory,
    srnGlobalCompute,
    srnGlobalSmooth,
    srnGlobalPlots,
    srnGlobalReplot,
    srnSingleGlobalDiagnostics
];

srnGlobalFiniteQ[x_] :=
    NumberQ[x] &&
    TrueQ[Im[x] == 0] &&
    FreeQ[x, Indeterminate | _DirectedInfinity];

srnGlobalRequire[
    test_,
    tag_String,
    details_: <||>
] :=
    If[
        !TrueQ[test],
        Throw[
            Failure[tag, details],
            "srnGlobalFailure"
        ]
    ];

srnGlobalExport[
    path_String,
    data_,
    format_String
] := Module[{result},

    result =
        checkedExport[
            path,
            data,
            format
        ];

    If[
        FailureQ[result],
        Throw[
            result,
            "srnGlobalFailure"
        ]
    ];

    result
];

srnGlobalDirectory[
    parent_String,
    prefix_String
] := Module[
    {out, result},

    out =
        FileNameJoin[
            {
                parent,
                prefix <> "_" <>
                    StringTake[
                        StringDelete[
                            CreateUUID[],
                            "-"
                        ],
                        8
                    ]
            }
        ];

    result =
        Quiet[
            Check[
                CreateDirectory[out],
                $Failed
            ]
        ];

    srnGlobalRequire[
        StringQ[result] &&
        DirectoryQ[out],
        "CreateGlobalDirectoryFailed",
        <|"path" -> out|>
    ];

    out
];


(* ::Subsection::Closed:: *)
(* Compute the four raw observables for one saved field *)


srnGlobalCompute[
    run_Association,
    psiInitial_List,
    bases_Association,
    batch_Integer
] := Catch[
    Module[
        {
            L = run["L"],
            times = run["timeGrid"],
            dA,
            n,

            bp,
            bm,

            h,
            h0,
            scaleH,
            scale0,

            hp,
            hm,
            h0p,
            h0m,

            ep,
            em,
            vp,
            vm,

            e0p,
            e0m,
            v0p,
            v0m,

            wp,
            wm,
            w0p,
            w0m,

            cp,
            cm,
            c0p,
            c0m,

            reconH,
            recon0,
            residualH,
            residual0,

            p,
            m,
            p0,
            m0,

            re,
            im,
            re0,
            im0,

            psiH,
            psiRef,

            sectorH0P,
            sectorH0M,
            ampH0P,
            ampH0M,

            entropy,
            loschmidt,
            iprH0,
            survival,

            pos,
            len,
            j,
            row,
            tchunk,

            normErrorH = 0.,
            normError0 = 0.,
            h0BasisNormError = 0.,
            probabilityBoundError = 0.,

            savedDifference,
            maxDifferenceIndex,
            tolerance = 10^-8,

            start = AbsoluteTime[]
        },

        srnGlobalRequire[
            EvenQ[L] &&
            L >= 2,
            "InvalidGlobalSystemSize"
        ];

        dA = 2^(L/2);
        n = Length[times];

        srnGlobalRequire[
            IntegerQ[batch] &&
            batch >= 1,
            "InvalidGlobalBatchSize"
        ];

        srnGlobalRequire[
            VectorQ[times, srnGlobalFiniteQ] &&
            n >= 2 &&
            Min[times] > 0 &&
            Min[Differences[times]] > 0,
            "InvalidGlobalTimeGrid"
        ];

        srnGlobalRequire[
            Length[psiInitial] == 2^L &&
            Abs[Norm[psiInitial]^2 - 1.] <= 10^-10,
            "InvalidSavedInitialState"
        ];

        bp = N[bases["Even"]];
        bm = N[bases["Odd"]];

        (* Actual perturbed Hamiltonian for this saved hx. *)

        h =
            N[
                IsingHamiltonian[
                    run["hx"],
                    run["hz"],
                    run["J"],
                    L
                ]
            ];

        (* Reference/integrable Hamiltonian: same J,hz,L with hx exactly zero. *)

        h0 =
            N[
                IsingHamiltonian[
                    0.,
                    run["hz"],
                    run["J"],
                    L
                ]
            ];

        srnGlobalRequire[
            Dimensions[h] === {2^L, 2^L} &&
            Dimensions[h0] === {2^L, 2^L},
            "InvalidGlobalHamiltonianDimension"
        ];

        scaleH =
            Max[
                1.,
                Norm[h, "Frobenius"]
            ];

        scale0 =
            Max[
                1.,
                Norm[h0, "Frobenius"]
            ];

        srnGlobalRequire[
            Norm[
                h - ConjugateTranspose[h],
                "Frobenius"
            ]/scaleH <= 10^-12,
            "PerturbedHamiltonianNotHermitian"
        ];

        srnGlobalRequire[
            Norm[
                h0 - ConjugateTranspose[h0],
                "Frobenius"
            ]/scale0 <= 10^-12,
            "ReferenceHamiltonianNotHermitian"
        ];

        srnGlobalRequire[
            Norm[
                ConjugateTranspose[bp] . h . bm,
                "Frobenius"
            ]/scaleH <= 10^-12,
            "PerturbedReflectionSymmetryFailed"
        ];

        srnGlobalRequire[
            Norm[
                ConjugateTranspose[bp] . h0 . bm,
                "Frobenius"
            ]/scale0 <= 10^-12,
            "ReferenceReflectionSymmetryFailed"
        ];

        hp =
            ConjugateTranspose[bp] .
            h .
            bp;

        hm =
            ConjugateTranspose[bm] .
            h .
            bm;

        h0p =
            ConjugateTranspose[bp] .
            h0 .
            bp;

        h0m =
            ConjugateTranspose[bm] .
            h0 .
            bm;

        (* Diagonalize H(hx) and H0 once each. *)

        {ep, vp} =
            sortedEigensystem[hp];

        {em, vm} =
            sortedEigensystem[hm];

        {e0p, v0p} =
            sortedEigensystem[h0p];

        {e0m, v0m} =
            sortedEigensystem[h0m];

        residualH =
            Max[
                sampledEigensystemResidual[
                    hp,
                    ep,
                    vp
                ],
                sampledEigensystemResidual[
                    hm,
                    em,
                    vm
                ]
            ];

        residual0 =
            Max[
                sampledEigensystemResidual[
                    h0p,
                    e0p,
                    v0p
                ],
                sampledEigensystemResidual[
                    h0m,
                    e0m,
                    v0m
                ]
            ];

        srnGlobalRequire[
            residualH <= 10^-10,
            "PerturbedEigensystemResidualFailed",
            <|"error" -> residualH|>
        ];

        srnGlobalRequire[
            residual0 <= 10^-10,
            "ReferenceEigensystemResidualFailed",
            <|"error" -> residual0|>
        ];

        (* The supplied Ising matrices are real in this workflow. *)

        srnGlobalRequire[
            Max[Abs[Im[vp]]] <= 10^-10 &&
            Max[Abs[Im[vm]]] <= 10^-10 &&
            Max[Abs[Im[v0p]]] <= 10^-10 &&
            Max[Abs[Im[v0m]]] <= 10^-10,
            "ComplexGlobalSectorEigenbasis"
        ];

        wp =
            Developer`ToPackedArray[
                Re[
                    Transpose[vp]
                ]
            ];

        wm =
            Developer`ToPackedArray[
                Re[
                    Transpose[vm]
                ]
            ];

        w0p =
            Developer`ToPackedArray[
                Re[
                    Transpose[v0p]
                ]
            ];

        w0m =
            Developer`ToPackedArray[
                Re[
                    Transpose[v0m]
                ]
            ];

        (* Spectral coefficients of the SAME initial RPS. *)

        cp =
            Conjugate[vp] .
            (
                Transpose[bp] .
                psiInitial
            );

        cm =
            Conjugate[vm] .
            (
                Transpose[bm] .
                psiInitial
            );

        c0p =
            Conjugate[v0p] .
            (
                Transpose[bp] .
                psiInitial
            );

        c0m =
            Conjugate[v0m] .
            (
                Transpose[bm] .
                psiInitial
            );

        reconH =
            Norm[
                bp . (wp . cp) +
                bm . (wm . cm) -
                psiInitial
            ];

        recon0 =
            Norm[
                bp . (w0p . c0p) +
                bm . (w0m . c0m) -
                psiInitial
            ];

        srnGlobalRequire[
            reconH <= 10^-10,
            "PerturbedInitialStateReconstructionFailed",
            <|"error" -> reconH|>
        ];

        srnGlobalRequire[
            recon0 <= 10^-10,
            "ReferenceInitialStateReconstructionFailed",
            <|"error" -> recon0|>
        ];

        entropy =
            ConstantArray[
                0.,
                n
            ];

        loschmidt =
            ConstantArray[
                0.,
                n
            ];

        iprH0 =
            ConstantArray[
                0.,
                n
            ];

        survival =
            ConstantArray[
                0.,
                n
            ];

        pos = 1;

        While[
            pos <= n,

            len =
                Min[
                    batch,
                    n - pos + 1
                ];

            tchunk =
                Developer`ToPackedArray[
                    N[
                        times[[
                            pos ;;
                            pos + len - 1
                        ]]
                    ]
                ];

            (* Perturbed trajectory. *)

            p =
                phaseAtTimesBatchC[
                    Developer`ToPackedArray[Re[cp]],
                    Developer`ToPackedArray[Im[cp]],
                    ep,
                    tchunk
                ];

            m =
                phaseAtTimesBatchC[
                    Developer`ToPackedArray[Re[cm]],
                    Developer`ToPackedArray[Im[cm]],
                    em,
                    tchunk
                ];

            re =
                bp . (wp . p[[1]]) +
                bm . (wm . m[[1]]);

            im =
                bp . (wp . p[[2]]) +
                bm . (wm . m[[2]]);

            (* H0 reference trajectory at exactly the same times. *)

            p0 =
                phaseAtTimesBatchC[
                    Developer`ToPackedArray[Re[c0p]],
                    Developer`ToPackedArray[Im[c0p]],
                    e0p,
                    tchunk
                ];

            m0 =
                phaseAtTimesBatchC[
                    Developer`ToPackedArray[Re[c0m]],
                    Developer`ToPackedArray[Im[c0m]],
                    e0m,
                    tchunk
                ];

            re0 =
                bp . (w0p . p0[[1]]) +
                bm . (w0m . m0[[1]]);

            im0 =
                bp . (w0p . p0[[2]]) +
                bm . (w0m . m0[[2]]);

            normErrorH =
                Max[
                    normErrorH,
                    Max[
                        Abs[
                            Total[
                                re re +
                                im im
                            ] -
                            1.
                        ]
                    ]
                ];

            normError0 =
                Max[
                    normError0,
                    Max[
                        Abs[
                            Total[
                                re0 re0 +
                                im0 im0
                            ] -
                            1.
                        ]
                    ]
                ];

            Do[
                row =
                    pos +
                    j -
                    1;

                psiH =
                    re[[All, j]] +
                    I im[[All, j]];

                psiRef =
                    re0[[All, j]] +
                    I im0[[All, j]];

                (* (a) Same half-chain entropy as before. *)

                entropy[[row]] =
                    entropyFast[
                        re[[All, j]],
                        im[[All, j]],
                        dA,
                        dA
                    ];

                (* (b) Loschmidt echo / trajectory fidelity relative to H0. *)

                loschmidt[[row]] =
                    Abs[
                        Conjugate[psiRef] .
                        psiH
                    ]^2;

                (* (c) IPR of the perturbed state in the H0 eigenbasis.
                   Even and odd reflection sectors together form the full basis. *)

                sectorH0P =
                    Transpose[bp] .
                    psiH;

                sectorH0M =
                    Transpose[bm] .
                    psiH;

                ampH0P =
                    Conjugate[v0p] .
                    sectorH0P;

                ampH0M =
                    Conjugate[v0m] .
                    sectorH0M;

                iprH0[[row]] =
                    Total[
                        Abs[ampH0P]^4
                    ] +
                    Total[
                        Abs[ampH0M]^4
                    ];

                h0BasisNormError =
                    Max[
                        h0BasisNormError,
                        Abs[
                            Total[
                                Abs[ampH0P]^2
                            ] +
                            Total[
                                Abs[ampH0M]^2
                            ] -
                            1.
                        ]
                    ];

                (* (d) Survival probability under H(hx). *)

                survival[[row]] =
                    Abs[
                        Conjugate[psiInitial] .
                        psiH
                    ]^2;

                probabilityBoundError =
                    Max[
                        probabilityBoundError,
                        Max[
                            0.,
                            loschmidt[[row]] - 1.,
                            -loschmidt[[row]],
                            survival[[row]] - 1.,
                            -survival[[row]],
                            iprH0[[row]] - 1.,
                            -iprH0[[row]]
                        ]
                    ];

                ,
                {
                    j,
                    len
                }
            ];

            pos += len;
        ];

        srnGlobalRequire[
            normErrorH <= 10^-10,
            "PerturbedNormValidationFailed",
            <|"error" -> normErrorH|>
        ];

        srnGlobalRequire[
            normError0 <= 10^-10,
            "ReferenceNormValidationFailed",
            <|"error" -> normError0|>
        ];

        srnGlobalRequire[
            h0BasisNormError <= 10^-10,
            "H0BasisNormalizationFailed",
            <|"error" -> h0BasisNormError|>
        ];

        srnGlobalRequire[
            probabilityBoundError <= 10^-9,
            "GlobalProbabilityBoundsFailed",
            <|"error" -> probabilityBoundError|>
        ];

        srnGlobalRequire[
            VectorQ[entropy, srnGlobalFiniteQ] &&
            VectorQ[loschmidt, srnGlobalFiniteQ] &&
            VectorQ[iprH0, srnGlobalFiniteQ] &&
            VectorQ[survival, srnGlobalFiniteQ],
            "InvalidGlobalObservableVector"
        ];

        (* Compare against the entropy saved by the original sweep.
           As in the spectra section, record a mismatch rather than aborting. *)

        savedDifference =
            Max[
                Abs[
                    entropy -
                    run["singleEntropy"]
                ]
            ];

        maxDifferenceIndex =
            First[
                Ordering[
                    Abs[
                        entropy -
                        run["singleEntropy"]
                    ],
                    -1
                ]
            ];

        If[
            savedDifference > tolerance,

            Print[
                "Saved entropy comparison: max |new - saved| = ",
                savedDifference,
                ", at t = ",
                times[[maxDifferenceIndex]],
                ". The new collage uses entropy recomputed from the same ",
                "perturbed trajectory as echo, IPR and survival."
            ];
        ];

        <|
            "L" -> L,
            "J" -> run["J"],
            "hz" -> run["hz"],
            "hx" -> run["hx"],
            "hxExponent" -> run["hxExponent"],

            "referenceHx" -> 0.,
            "referenceHamiltonian" ->
                "IsingHamiltonian[0,hz,J,L]",

            "iprBasis" ->
                "Eigenbasis of H0 = IsingHamiltonian[0,hz,J,L]",

            "iprDefinition" ->
                "Sum_n Abs[<E_n^(0)|psi_h(t)>]^4",

            "loschmidtDefinition" ->
                "Abs[<psi_0(t)|psi_h(t)>]^2",

            "survivalDefinition" ->
                "Abs[<psi(0)|psi_h(t)>]^2",

            "singleStateSeed" ->
                Lookup[
                    run,
                    "singleStateSeed",
                    Missing["NotSaved"]
                ],

            "timeGrid" ->
                times,

            "initialState" ->
                psiInitial,

            "entropy" ->
                Developer`ToPackedArray[
                    entropy
                ],

            "savedEntropy" ->
                run["singleEntropy"],

            "loschmidtEcho" ->
                Developer`ToPackedArray[
                    loschmidt
                ],

            "iprH0" ->
                Developer`ToPackedArray[
                    iprH0
                ],

            "survival" ->
                Developer`ToPackedArray[
                    survival
                ],

            "savedEntropyAgreement" ->
                TrueQ[
                    savedDifference <= tolerance
                ],

            "savedEntropyComparisonTolerance" ->
                tolerance,

            "savedEntropyMaxDifference" ->
                savedDifference,

            "savedEntropyDifferenceTime" ->
                times[[maxDifferenceIndex]],

            "maxPerturbedNormError" ->
                normErrorH,

            "maxReferenceNormError" ->
                normError0,

            "maxH0BasisNormError" ->
                h0BasisNormError,

            "probabilityBoundError" ->
                probabilityBoundError,

            "perturbedEigensystemResidual" ->
                residualH,

            "referenceEigensystemResidual" ->
                residual0,

            "perturbedReconstructionError" ->
                reconH,

            "referenceReconstructionError" ->
                recon0,

            "computeSeconds" ->
                N[
                    AbsoluteTime[] -
                    start
                ]
        |>
    ],
    "srnGlobalFailure"
];


(* ::Subsection::Closed:: *)
(* Identical moving average for all four observables *)


srnGlobalSmooth[
    data_Association,
    window_Integer
] := Module[{},
    srnGlobalRequire[
        window >= 1,
        "InvalidGlobalMovingAverageWindow"
    ];

    <|
        "window" ->
            window,

        "entropy" ->
            startupMovingAverage[
                data["entropy"],
                window
            ],

        "savedEntropy" ->
            startupMovingAverage[
                Lookup[
                    data,
                    "savedEntropy",
                    data["entropy"]
                ],
                window
            ],

        "loschmidtEcho" ->
            startupMovingAverage[
                data["loschmidtEcho"],
                window
            ],

        "iprH0" ->
            startupMovingAverage[
                data["iprH0"],
                window
            ],

        "survival" ->
            startupMovingAverage[
                data["survival"],
                window
            ]
    |>
];


(* ::Subsection::Closed:: *)
(* Common plotting function: same window and same settings in all panels *)


srnGlobalPlots[
    data_Association,
    smooth_Association,
    timeWindow_: Automatic
] := Catch[
    Module[
        {
            times = data["timeGrid"],
            bounds,
            idx,
            t,
            ticks,
            tickStep,
            common,
            make,

            entropyPlot,
            loschmidtPlot,
            iprPlot,
            survivalPlot,

            title,
            collage
        },

        bounds =
            If[
                timeWindow === Automatic,
                {
                    First[times],
                    Last[times]
                },
                N[timeWindow]
            ];

        srnGlobalRequire[
            MatchQ[
                bounds,
                {
                    _?srnGlobalFiniteQ,
                    _?srnGlobalFiniteQ
                }
            ] &&
            TrueQ[
                0 <
                bounds[[1]] <
                bounds[[2]]
            ] &&
            bounds[[1]] >= First[times] &&
            bounds[[2]] <= Last[times],
            "InvalidGlobalTimeWindow",
            <|
                "available" ->
                    {
                        First[times],
                        Last[times]
                    }
            |>
        ];

        idx =
            Flatten[
                Position[
                    times,
                    _?(
                        TrueQ[
                            bounds[[1]] <=
                            # <=
                            bounds[[2]]
                        ] &
                    )
                ]
            ];

        srnGlobalRequire[
            Length[idx] >= 2,
            "TooFewVisibleGlobalTimeSamples"
        ];

        t =
            times[[idx]];

        tickStep =
            Max[
                1,
                Ceiling[
                    (
                        Log10[bounds[[2]]] -
                        Log10[bounds[[1]]]
                    )/5
                ]
            ];

        ticks =
            Table[
                {
                    10.^q,
                    Style[
                        Superscript[
                            10,
                            q
                        ],
                        16,
                        Black
                    ]
                },
                {
                    q,
                    Ceiling[
                        Log10[
                            bounds[[1]]
                        ]
                    ],
                    Floor[
                        Log10[
                            bounds[[2]]
                        ]
                    ],
                    tickStep
                }
            ];

        If[
            ticks === {},

            ticks =
                (
                    {
                        #,
                        ScientificForm[
                            #,
                            2
                        ]
                    } &
                ) /@
                bounds;
        ];

        (* EXACTLY the same panel geometry/style used by srnSpecPlots. *)

        common = {
            Joined ->
                False,

            Frame ->
                True,

            Axes ->
                False,

            PlotRange ->
                {
                    bounds,
                    {
                        0.,
                        All
                    }
                },

            PlotRangePadding ->
                {
                    None,
                    Scaled[.04]
                },

            FrameTicks ->
                {
                    {
                        Automatic,
                        None
                    },
                    {
                        ticks,
                        None
                    }
                },

            FrameStyle ->
                Directive[
                    Black,
                    16
                ],

            BaseStyle ->
                {
                    FontFamily ->
                        "Times",

                    FontSize ->
                        16
                },

            ImageSize ->
                550,

            AspectRatio ->
                2/3,

            ImagePadding ->
                {
                    {
                        72,
                        20
                    },
                    {
                        55,
                        42
                    }
                },

            PlotTheme ->
                "Detailed",

            GridLines ->
                None,

            Background ->
                White
        };

        make[
            coords_,
            ylabel_,
            caption_
        ] :=
            ListLogLinearPlot[
                coords,

                Evaluate[
                    Sequence @@
                    common
                ],

                PlotStyle ->
                    {
                        Directive[
                            Black,
                            Opacity[.82],
                            PointSize[.002]
                        ]
                    },

                FrameLabel ->
                    {
                        Style[
                            "t",
                            20,
                            Italic
                        ],

                        Style[
                            ylabel,
                            20
                        ]
                    },

                PlotLabel ->
                    Style[
                        caption,
                        18,
                        Black
                    ]
            ];

        entropyPlot =
            make[
                {
                    Transpose[
                        {
                            t,
                            smooth["entropy"][[idx]]
                        }
                    ]
                },
                "S(t)",
                "(a) Entanglement entropy"
            ];

        loschmidtPlot =
            make[
                {
                    Transpose[
                        {
                            t,
                            smooth["loschmidtEcho"][[idx]]
                        }
                    ]
                },
                "\[ScriptCapitalL](t)",
                "(b) Loschmidt echo"
            ];

        iprPlot =
            make[
                {
                    Transpose[
                        {
                            t,
                            smooth["iprH0"][[idx]]
                        }
                    ]
                },
                Row[
                    {
                        "IPR",
                        Subscript["H", "0"],
                        "(t)"
                    }
                ],
                "(c) IPR in the H0 eigenbasis"
            ];

        survivalPlot =
            make[
                {
                    Transpose[
                        {
                            t,
                            smooth["survival"][[idx]]
                        }
                    ]
                },
                "P(t)",
                "(d) Survival probability"
            ];

        title =
            Column[
                {
                    Style[
                        Row[
                            {
                                "Single RPS,  L = ",
                                data["L"],
                                ",  J = ",
                                data["J"],
                                ",  ",
                                Subscript[
                                    "h",
                                    "z"
                                ],
                                " = ",
                                data["hz"],
                                ",  ",
                                Subscript[
                                    "h",
                                    "x"
                                ],
                                " = ",
                                Superscript[
                                    10,
                                    data["hxExponent"]
                                ]
                            }
                        ],
                        22,
                        Black
                    ],

                    Style[
                        Row[
                            {
                                "Moving average: ",
                                smooth["window"],
                                " samples;  Loschmidt reference ",
                                Subscript[
                                    "h",
                                    "x"
                                ],
                                " = 0;  IPR basis = eigenbasis of ",
                                Subscript[
                                    "H",
                                    "0"
                                ]
                            }
                        ],
                        16,
                        Black
                    ],

                    If[
                        TrueQ[
                            Lookup[
                                data,
                                "savedEntropyAgreement",
                                True
                            ]
                        ],
                        Nothing,
                        Style[
                            Row[
                                {
                                    "Entropy recomputed with diagnostics; max difference from saved curve = ",
                                    ScientificForm[
                                        data[
                                            "savedEntropyMaxDifference"
                                        ],
                                        3
                                    ]
                                }
                            ],
                            14,
                            Black
                        ]
                    ]
                },
                Alignment ->
                    Center
            ];

        collage =
            Labeled[
                GraphicsGrid[
                    {
                        {
                            entropyPlot,
                            loschmidtPlot
                        },
                        {
                            iprPlot,
                            survivalPlot
                        }
                    },

                    Spacings ->
                        {
                            .35,
                            .5
                        },

                    ImageSize ->
                        1140,

                    Background ->
                        White
                ],

                title,
                Top
            ];

        <|
            "collage" ->
                collage,

            "entropy" ->
                Labeled[
                    entropyPlot,
                    title,
                    Top
                ],

            "loschmidt" ->
                Labeled[
                    loschmidtPlot,
                    title,
                    Top
                ],

            "ipr" ->
                Labeled[
                    iprPlot,
                    title,
                    Top
                ],

            "survival" ->
                Labeled[
                    survivalPlot,
                    title,
                    Top
                ],

            "panels" ->
                {
                    entropyPlot,
                    loschmidtPlot,
                    iprPlot,
                    survivalPlot
                },

            "timeWindow" ->
                bounds
        |>
    ],
    "srnGlobalFailure"
];


(* ::Subsection::Closed:: *)
(* Replot saved raw diagnostics without recomputing dynamics *)


srnGlobalReplot[
    file_String,
    window_Integer : 200,
    timeWindow_: Automatic,
    display_: True
] := Catch[
    Module[
        {
            data,
            smooth,
            plots,
            out,
            files
        },

        data =
            Quiet[
                Check[
                    Import[
                        file,
                        "WXF"
                    ],
                    $Failed
                ]
            ];

        srnGlobalRequire[
            AssociationQ[data] &&
            AllTrue[
                {
                    "timeGrid",
                    "entropy",
                    "loschmidtEcho",
                    "iprH0",
                    "survival"
                },
                KeyExistsQ[
                    data,
                    #
                ] &
            ],
            "InvalidSavedGlobalDiagnostics",
            <|"path" -> file|>
        ];

        smooth =
            srnGlobalSmooth[
                data,
                window
            ];

        plots =
            srnGlobalPlots[
                data,
                smooth,
                timeWindow
            ];

        If[
            FailureQ[plots],
            Throw[
                plots,
                "srnGlobalFailure"
            ]
        ];

        out =
            srnGlobalDirectory[
                DirectoryName[file],
                FileBaseName[file] <>
                    "_plots"
            ];

        srnGlobalExport[
            FileNameJoin[
                {
                    out,
                    "moving_average.wxf"
                }
            ],
            <|
                "sourceFile" ->
                    file,

                "timeGrid" ->
                    data["timeGrid"],

                "smoothed" ->
                    smooth,

                "timeWindow" ->
                    plots["timeWindow"],

                "movingAverageConvention" ->
                    "Trailing sample mean with growing startup; original times unchanged"
            |>,
            "WXF"
        ];

        files = <|

            "collage" ->
                srnGlobalExport[
                    FileNameJoin[
                        {
                            out,
                            "collage.png"
                        }
                    ],
                    plots["collage"],
                    "PNG"
                ],

            "entropy" ->
                srnGlobalExport[
                    FileNameJoin[
                        {
                            out,
                            "entropy.png"
                        }
                    ],
                    plots["entropy"],
                    "PNG"
                ],

            "loschmidt" ->
                srnGlobalExport[
                    FileNameJoin[
                        {
                            out,
                            "loschmidt_echo.png"
                        }
                    ],
                    plots["loschmidt"],
                    "PNG"
                ],

            "ipr" ->
                srnGlobalExport[
                    FileNameJoin[
                        {
                            out,
                            "ipr_H0.png"
                        }
                    ],
                    plots["ipr"],
                    "PNG"
                ],

            "survival" ->
                srnGlobalExport[
                    FileNameJoin[
                        {
                            out,
                            "survival_probability.png"
                        }
                    ],
                    plots["survival"],
                    "PNG"
                ]
        |>;

        If[
            TrueQ[display],
            Print[
                plots["collage"]
            ]
        ];

        Print[
            "Global-diagnostics plots saved in: ",
            out
        ];

        <|
            "sourceFile" ->
                file,

            "plotDirectory" ->
                out,

            "plotFiles" ->
                files,

            "timeWindow" ->
                plots["timeWindow"]
        |>
    ],
    "srnGlobalFailure"
];


(* ::Subsection::Closed:: *)
(* Complete callable computation and export function *)


srnSingleGlobalDiagnostics[
    directory_String,
    fieldIndices_: {1},
    window_Integer : 200,
    timeWindow_: Automatic,
    display_: True
] := Catch[
    Module[
        {
            meta,
            files,
            ordered,
            indices,
            psiInitial,
            bases,
            out,
            batch,
            run,
            result,
            path,
            summaries = {},
            plots,
            notes,
            k
        },

        srnGlobalRequire[
            DirectoryQ[directory],
            "RunDirectoryNotFound",
            <|"path" -> directory|>
        ];

        srnGlobalRequire[
            window >= 1,
            "InvalidGlobalMovingAverageWindow"
        ];

        meta =
            Quiet[
                Check[
                    Import[
                        FileNameJoin[
                            {
                                directory,
                                "metadata.wxf"
                            }
                        ],
                        "WXF"
                    ],
                    $Failed
                ]
            ];

        srnGlobalRequire[
            AssociationQ[meta] &&
            AllTrue[
                {
                    "L",
                    "J",
                    "hz",
                    "initialStates",
                    "hxList",
                    "timeGrid"
                },
                KeyExistsQ[
                    meta,
                    #
                ] &
            ],
            "MissingGlobalRunMetadata"
        ];

        srnGlobalRequire[
            EvenQ[meta["L"]] &&
            meta["L"] >= 2 &&
            ListQ[meta["initialStates"]] &&
            Length[meta["initialStates"]] >= 1,
            "InvalidGlobalRunMetadata"
        ];

        files =
            FileNames[
                "log_run_" ~~
                DigitCharacter .. ~~
                ".wxf",
                directory
            ];

        ordered =
            SortBy[
                files,
                ToExpression[
                    StringReplace[
                        FileBaseName[#],
                        "log_run_" ->
                            ""
                    ]
                ] &
            ];

        srnGlobalRequire[
            Length[ordered] > 0,
            "NoSavedGlobalLogRuns"
        ];

        indices =
            If[
                fieldIndices === All,
                Range[
                    Length[ordered]
                ],
                fieldIndices
            ];

        srnGlobalRequire[
            ListQ[indices] &&
            indices =!= {} &&
            VectorQ[
                indices,
                IntegerQ
            ] &&
            Min[indices] >= 1 &&
            Max[indices] <= Length[ordered],
            "InvalidGlobalFieldIndices",
            <|
                "availableCount" ->
                    Length[ordered]
            |>
        ];

        indices =
            DeleteDuplicates[
                indices
            ];

        psiInitial =
            Developer`ToPackedArray[
                N[
                    First[
                        meta[
                            "initialStates"
                        ]
                    ]
                ]
            ];

        bases =
            reflectionSectorBases[
                meta["L"]
            ];

        batch =
            Lookup[
                meta,
                "batchSize",
                32
            ];

        srnGlobalRequire[
            IntegerQ[batch] &&
            batch >= 1,
            "InvalidGlobalBatchSize"
        ];

        out =
            srnGlobalDirectory[
                directory,
                "global_diagnostics"
            ];

        notes =
            StringRiffle[
                {
                    "Single-RPS global diagnostics; saved ensemble member 1, no new random state.",
                    "Physical parameters and exact logarithmic time grid are read from the saved run.",
                    "Perturbed trajectory: psi_h(t)=Exp[-I H(hx)t] psi(0).",
                    "Reference trajectory: psi_0(t)=Exp[-I H0 t] psi(0), H0=IsingHamiltonian[0,hz,J,L].",
                    "Loschmidt echo: Abs[<psi_0(t)|psi_h(t)>]^2.",
                    "IPR: Sum_n Abs[<E_n^(0)|psi_h(t)>]^4 in the numerical eigenbasis of H0.",
                    "For exact H0 degeneracies, IPR depends on the chosen eigenvectors inside that degenerate subspace.",
                    "Survival probability: Abs[<psi(0)|psi_h(t)>]^2.",
                    "Entropy is recomputed from the same perturbed state used for all other diagnostics.",
                    "Same trailing sample moving average with growing startup as the main entropy plots.",
                    "Moving-average window: " <> ToString[window],
                    "Requested display time window: " <> ToString[timeWindow, InputForm],
                    "Smoothing occurs before display zoom; time coordinates are never averaged.",
                    "All four panels share the same visible logarithmic time interval and plot geometry.",
                    "Raw WXF data are saved before plotting.",
                    "Use srnGlobalReplot on any global_XX.wxf to change window/zoom without dynamics."
                },
                "\n"
            ];

        srnGlobalExport[
            FileNameJoin[
                {
                    out,
                    "description.txt"
                }
            ],
            notes,
            "Text"
        ];

        Do[
            Print[
                "Computing single-RPS global diagnostics: selected field ",
                k,
                " of ",
                Length[ordered],
                "."
            ];

            run =
                Quiet[
                    Check[
                        Import[
                            ordered[[k]],
                            "WXF"
                        ],
                        $Failed
                    ]
                ];

            srnGlobalRequire[
                AssociationQ[run] &&
                TrueQ[
                    Lookup[
                        run,
                        "success",
                        False
                    ]
                ] &&
                AllTrue[
                    {
                        "L",
                        "J",
                        "hz",
                        "hx",
                        "hxExponent",
                        "timeGrid",
                        "singleEntropy"
                    },
                    KeyExistsQ[
                        run,
                        #
                    ] &
                ],
                "InvalidSavedGlobalSingleRun",
                <|
                    "path" ->
                        ordered[[k]]
                |>
            ];

            srnGlobalRequire[
                Lookup[
                    run,
                    {
                        "L",
                        "J",
                        "hz"
                    }
                ] ===
                Lookup[
                    meta,
                    {
                        "L",
                        "J",
                        "hz"
                    }
                ] &&

                run["timeGrid"] ===
                meta["timeGrid"] &&

                VectorQ[
                    run["timeGrid"],
                    srnGlobalFiniteQ
                ] &&

                Length[
                    run["timeGrid"]
                ] >= 2 &&

                Min[
                    run["timeGrid"]
                ] > 0 &&

                Min[
                    Differences[
                        run["timeGrid"]
                    ]
                ] > 0 &&

                VectorQ[
                    run["singleEntropy"],
                    srnGlobalFiniteQ
                ] &&

                Length[
                    run["singleEntropy"]
                ] ==
                Length[
                    run["timeGrid"]
                ],
                "InconsistentSavedGlobalRun"
            ];

            result =
                srnGlobalCompute[
                    run,
                    psiInitial,
                    bases,
                    batch
                ];

            If[
                FailureQ[result],
                Throw[
                    result,
                    "srnGlobalFailure"
                ]
            ];

            path =
                FileNameJoin[
                    {
                        out,
                        "global_" <>
                            IntegerString[
                                k,
                                10,
                                2
                            ] <>
                            ".wxf"
                    }
                ];

            (* Save raw data BEFORE any smoothing or plot rendering. *)

            srnGlobalExport[
                path,
                result,
                "WXF"
            ];

            plots =
                srnGlobalReplot[
                    path,
                    window,
                    timeWindow,
                    display
                ];

            If[
                FailureQ[plots],
                Throw[
                    plots,
                    "srnGlobalFailure"
                ]
            ];

            AppendTo[
                summaries,
                Join[
                    <|
                        "index" ->
                            k,

                        "hx" ->
                            result["hx"],

                        "hxExponent" ->
                            result[
                                "hxExponent"
                            ],

                        "computeSeconds" ->
                            result[
                                "computeSeconds"
                            ],

                        "savedEntropyMaxDifference" ->
                            result[
                                "savedEntropyMaxDifference"
                            ],

                        "savedEntropyAgreement" ->
                            result[
                                "savedEntropyAgreement"
                            ],

                        "savedEntropyDifferenceTime" ->
                            result[
                                "savedEntropyDifferenceTime"
                            ]
                    |>,
                    plots
                ]
            ];

            srnGlobalExport[
                FileNameJoin[
                    {
                        out,
                        "summary.wxf"
                    }
                ],
                summaries,
                "WXF"
            ];

            ,
            {
                k,
                indices
            }
        ];

        Print[
            "Single-RPS global diagnostics complete. Data folder: ",
            out
        ];

        <|
            "outputDirectory" ->
                out,

            "fields" ->
                summaries
        |>
    ],
    "srnGlobalFailure"
];


(* ::Subsection:: *)
(* Parameters for THIS section only *)


srnGlobalSourceDirectory =
    logOnlyOutputDir;
    (* Or replace by an existing completed run-folder string. *)

srnGlobalFieldIndices =
    Range[17];
    (* For the supplied sweep -8,-7.5,...,0, index 10 is hx = 10^-3.5.
       Change to {1}, {1,5,10}, or All as needed. *)

srnGlobalWindow =
    logOnlyMovingAverageWindow;

srnGlobalTimeWindow =
    Automatic;
    (* Full saved interval; or e.g. {10., 10.^8}. *)

srnGlobalDisplay =
    True;
    (* Print the collage in the notebook as well as exporting it. *)


(* ::Subsection:: *)
(* Execute this section: the existing entropy sweep is NOT run again *)


If[
    Length[
        DownValues[
            srnSingleGlobalDiagnostics
        ]
    ] == 0,

    Print[
        "ERROR: evaluate all definitions in this section first."
    ];

    Abort[];
];

If[
    !StringQ[
        srnGlobalSourceDirectory
    ],

    Print[
        "ERROR: set srnGlobalSourceDirectory to a completed run folder."
    ];

    Abort[];
];

srnGlobalResult =
    srnSingleGlobalDiagnostics[
        srnGlobalSourceDirectory,
        srnGlobalFieldIndices,
        srnGlobalWindow,
        srnGlobalTimeWindow,
        srnGlobalDisplay
    ];

If[
    FailureQ[
        srnGlobalResult
    ],
    Print[
        srnGlobalResult
    ]
];


(* Replot later WITHOUT evolving any state:

   srnGlobalReplot[
       "C:\\...\\global_diagnostics_abcd1234\\global_10.wxf",
       200,
       {10.^-1, 10.^10},
       True
   ];

   Or use the first result immediately:

   srnGlobalReplot[
       srnGlobalResult["fields"][[1]]["sourceFile"],
       srnGlobalWindow,
       srnGlobalTimeWindow,
       True
   ];
*)
