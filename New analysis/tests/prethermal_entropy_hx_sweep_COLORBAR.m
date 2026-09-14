(* ::Package:: *)

(* ::Title::Closed:: *)
(*Setup*)


qmbInitPath = "C:\\Users\\Miguel\\Github\\libs\\QMB\\Kernel\\init.m";
Get[qmbInitPath];
SetDirectory[NotebookDirectory[]];


(* ::Title::Closed:: *)
(*Definitions and compiled kernels*)


(* ::Input:: *)
(**)


ClearAll[reflectionSectorBases];

reflectionSectorBases[L_Integer] := Module[
    {dim, perm, plusRules = {}, minusRules = {}, np = 0, nm = 0, j, k},

    dim = 2^L;

    perm = Table[
        1 + FromDigits[
            Reverse[IntegerDigits[j - 1, 2, L]],
            2
        ],
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


(* ::Input:: *)
(**)


ClearAll[sortedEigensystem];

sortedEigensystem[M_] := Module[
    {evals, evecs, ord},

    {evals, evecs} = Eigensystem[Normal[M]];
    ord = Ordering[evals];

    {
        evals[[ord]],
        evecs[[ord]]
    }
];


(* ::Input:: *)
(**)


ClearAll[phaseStepC];

phaseStepC = Compile[
    {
        {aR, _Real, 1},
        {aI, _Real, 1},
        {cosStep, _Real, 1},
        {sinStep, _Real, 1}
    },

    Module[
        {newR, newI},

        newR = aR*cosStep + aI*sinStep;
        newI = aI*cosStep - aR*sinStep;

        {newR, newI}
    ],

    CompilationTarget -> "C",
    RuntimeOptions -> "Speed"
];


(* ::Input:: *)
(**)


ClearAll[entropyFromSVC];

entropyFromSVC = Compile[
    {{sv, _Real, 1}},

    Module[
        {sum = 0., p = 0., i, n},

        n = Length[sv];

        For[
            i = 1,
            i <= n,
            i++,

            p = sv[[i]]*sv[[i]];

            If[
                p > 1.*^-14,
                sum -= p*Log[p]
            ];
        ];

        sum
    ],

    CompilationTarget -> "C",
    RuntimeOptions -> "Speed"
];


(* ::Input:: *)
(**)


ClearAll[reseedC];

reseedC = Compile[
    {
        {cR, _Real, 1},
        {cI, _Real, 1},
        {energies, _Real, 1},
        {t, _Real}
    },

    Module[
        {ct, st},

        ct = Cos[energies*t];
        st = Sin[energies*t];

        {
            cR*ct + cI*st,
            cI*ct - cR*st
        }
    ],

    CompilationTarget -> "C",
    RuntimeOptions -> "Speed"
];


(* ::Input:: *)
(**)


ClearAll[entropyFast];

entropyFast[
    psiR_,
    psiI_,
    dA_Integer,
    dB_Integer
] := Module[
    {mat, sv},

    mat = ArrayReshape[
        psiR + I psiI,
        {dA, dB}
    ];

    sv = SingularValueList[mat];

    entropyFromSVC[sv]
];


(* ::Title::Closed:: *)
(*Parameters*)


L = 8;
J = 1.;
hz = 1.;


hxExponents = N[Range[-2, 0, 1/10]];
hxList = Developer`ToPackedArray[N[10^hxExponents]];


tmin = 0.;
tmax = 10^5;
dt = 1.;

nSteps = Round[(tmax - tmin)/dt];

If[
    Abs[tmin + nSteps*dt - tmax] > 10^-10,
    Print["ERROR: tmin, tmax and dt do not define an integer number of steps."];
    Abort[];
];


(* ::Input:: *)
(* (*Number of recurrence steps between exact phase reseeds. *)*)


resetEvery = 100000;


(* ::Input:: *)
(* (*Fix the RPS across the entire hx sweep. *)*)


seed = 12345;


LA = L/2;
dA = 2^LA;
dB = 2^(L - LA);


outputDir = FileNameJoin[
    {
        NotebookDirectory[],
        "prethermal_entropy_hx_sweep"
    }
];

If[
    !DirectoryQ[outputDir],
    CreateDirectory[outputDir]
];


(* ::Title::Closed:: *)
(*Fixed reflection bases and fixed RPS*)


(* ::Input:: *)
(**)


bases = reflectionSectorBases[L];

Bp = bases["Even"];
Bm = bases["Odd"];


(* ::Input:: *)
(**)


SeedRandom[seed];

psi0 = Developer`ToPackedArray[
    N[RandomChainProductState[L]]
];


(* ::Input:: *)
(* (*These projections depend only on L, the reflection basis and psi0. *)*)


psi0P = Developer`ToPackedArray[
    N[ConjugateTranspose[Bp] . psi0]
];

psi0M = Developer`ToPackedArray[
    N[ConjugateTranspose[Bm] . psi0]
];


(* ::Input:: *)
(**)


Print[
    "Fixed RPS sector norm = ",
    Norm[psi0P]^2 + Norm[psi0M]^2
];


(* ::Input:: *)
(**)


manifest = <|
    "L" -> L,
    "J" -> J,
    "hz" -> hz,
    "seed" -> seed,
    "psi0" -> psi0,
    "hxExponents" -> hxExponents,
    "hxList" -> hxList,
    "tmin" -> tmin,
    "tmax" -> tmax,
    "dt" -> dt,
    "nSteps" -> nSteps,
    "resetEvery" -> resetEvery,
    "sectorDimensions" -> {
        Dimensions[Bp][[2]],
        Dimensions[Bm][[2]]
    }
|>;

Export[
    FileNameJoin[{outputDir, "manifest.wxf"}],
    manifest,
    "WXF"
];


(* ::Title:: *)
(*Automatic hx sweep*)


(* ::Input:: *)
(**)


runSummary = {};


(* ::Input:: *)
(**)


Do[

    hxExponent = hxExponents[[k]];
    hx = hxList[[k]];

    Print[""];
    Print["============================================================"];
    Print[
        "Run ", k, "/", Length[hxList],
        "   exponent = ", hxExponent,
        "   hx = ", hx
    ];
    Print["============================================================"];

    runStart = AbsoluteTime[];


    (* --------------------------------------------------------- *)
    (* Hamiltonian and reflection-sector projection              *)
    (* --------------------------------------------------------- *)

    H = IsingHamiltonian[hx, hz, J, L];

    Hp = Chop[
        ConjugateTranspose[Bp] . H . Bp
    ];

    Hm = Chop[
        ConjugateTranspose[Bm] . H . Bm
    ];


    (* --------------------------------------------------------- *)
    (* Diagonalize each parity sector separately                 *)
    (* --------------------------------------------------------- *)

    {diagSeconds, eigensystems} = AbsoluteTiming[
        {
            sortedEigensystem[Hp],
            sortedEigensystem[Hm]
        }
    ];

    {
        {evalsP, evecsP},
        {evalsM, evecsM}
    } = eigensystems;


    (* --------------------------------------------------------- *)
    (* Same fixed RPS, expressed in the new eigenbases           *)
    (* --------------------------------------------------------- *)

    coeffP = Conjugate[evecsP] . psi0P;
    coeffM = Conjugate[evecsM] . psi0M;


    (* --------------------------------------------------------- *)
    (* Reconstruct full eigenbasis once for this hx              *)
    (* --------------------------------------------------------- *)

    VP = Bp . Transpose[evecsP];
    VM = Bm . Transpose[evecsM];

    Vfull = Developer`ToPackedArray[
        N[Join[VP, VM, 2]]
    ];

    evalsFull = Developer`ToPackedArray[
        N[Join[evalsP, evalsM]]
    ];

    coeffFull = Developer`ToPackedArray[
        N[Join[coeffP, coeffM]]
    ];


    (* --------------------------------------------------------- *)
    (* Sanity checks                                             *)
    (* --------------------------------------------------------- *)

    reconstructionError = Norm[
        Vfull . coeffFull - psi0
    ];

    maxImagV = Max[
        Abs[Im[Vfull]]
    ];

    Print[
        "Diagonalization time [s] = ",
        diagSeconds
    ];

    Print[
        "Initial-state reconstruction error = ",
        reconstructionError
    ];

    Print[
        "Max imaginary part of Vfull = ",
        maxImagV
    ];

    If[
        maxImagV > 10^-10,
        Print[
            "ERROR: Vfull is not numerically real. ",
            "The two-real-Dot optimization is not safe."
        ];
        Abort[];
    ];


    (* --------------------------------------------------------- *)
    (* Packed arrays required by the long-time kernel            *)
    (* --------------------------------------------------------- *)

    V = Developer`ToPackedArray[
        Re[Vfull]
    ];

    cR = Developer`ToPackedArray[
        Re[coeffFull]
    ];

    cI = Developer`ToPackedArray[
        Im[coeffFull]
    ];

    cosStep = Developer`ToPackedArray[
        Cos[evalsFull*dt]
    ];

    sinStep = Developer`ToPackedArray[
        Sin[evalsFull*dt]
    ];


    (* --------------------------------------------------------- *)
    (* Initialize exactly at tmin                                *)
    (* --------------------------------------------------------- *)

    {aR, aI} = reseedC[
        cR,
        cI,
        evalsFull,
        N[tmin]
    ];


    (* --------------------------------------------------------- *)
    (* Long-time sequential entropy evolution                    *)
    (* --------------------------------------------------------- *)

    {evolutionSeconds, entropy} = AbsoluteTiming[

        Table[

            (* Full state in the computational basis *)
            psiR = V . aR;
            psiI = V . aI;

            (* Half-chain von Neumann entropy *)
            SvN = entropyFast[
                psiR,
                psiI,
                dA,
                dB
            ];

            (* Advance one uniform time step *)
            If[
                step < nSteps,

                If[
                    Mod[step + 1, resetEvery] == 0,

                    nextTime = tmin + (step + 1)*dt;

                    {aR, aI} = reseedC[
                        cR,
                        cI,
                        evalsFull,
                        N[nextTime]
                    ],

                    {aR, aI} = phaseStepC[
                        aR,
                        aI,
                        cosStep,
                        sinStep
                    ]
                ]
            ];

            SvN,

            {step, 0, nSteps}
        ]
    ];

    entropy = Developer`ToPackedArray[
        entropy
    ];


    (* --------------------------------------------------------- *)
    (* Save immediately: do not retain all long trajectories     *)
    (* --------------------------------------------------------- *)

    runFile = FileNameJoin[
        {
            outputDir,
            "run_" <> IntegerString[k, 10, 2] <> ".wxf"
        }
    ];

    totalSeconds = N[
        AbsoluteTime[] - runStart
    ];

    Export[
        runFile,

        <|
            "index" -> k,
            "L" -> L,
            "J" -> J,
            "hz" -> hz,
            "hxExponent" -> hxExponent,
            "hx" -> hx,
            "tmin" -> tmin,
            "tmax" -> tmax,
            "dt" -> dt,
            "nSteps" -> nSteps,
            "entropy" -> entropy,
            "evalsP" -> Developer`ToPackedArray[N[evalsP]],
            "evalsM" -> Developer`ToPackedArray[N[evalsM]],
            "reconstructionError" -> reconstructionError,
            "maxImagV" -> maxImagV,
            "diagonalizationSeconds" -> diagSeconds,
            "evolutionSeconds" -> evolutionSeconds,
            "totalSeconds" -> totalSeconds
        |>,

        "WXF"
    ];


    AppendTo[
        runSummary,
        {
            k,
            hxExponent,
            hx,
            diagSeconds,
            evolutionSeconds,
            totalSeconds,
            reconstructionError,
            runFile
        }
    ];

    Print[
        "Evolution time [s] = ",
        evolutionSeconds
    ];

    Print[
        "Total run time [s] = ",
        totalSeconds
    ];

    Print[
        "Saved: ",
        runFile
    ];


    (* --------------------------------------------------------- *)
    (* Release hx-dependent large arrays before the next run     *)
    (* --------------------------------------------------------- *)

    Clear[
        H, Hp, Hm,
        eigensystems,
        evalsP, evecsP,
        evalsM, evecsM,
        coeffP, coeffM,
        VP, VM,
        Vfull, evalsFull, coeffFull,
        V, cR, cI,
        cosStep, sinStep,
        aR, aI,
        psiR, psiI,
        entropy
    ];

    ,
    {k, Length[hxList]}
];


(* ::Title:: *)
(*Sweep summary*)


(* ::Input:: *)
(**)


summaryHeader = {
    "index",
    "hxExponent",
    "hx",
    "diagonalizationSeconds",
    "evolutionSeconds",
    "totalSeconds",
    "reconstructionError",
    "file"
};


(* ::Input:: *)
(**)


Export[
    FileNameJoin[{outputDir, "run_summary.csv"}],
    Prepend[runSummary, summaryHeader]
];


(* ::Input:: *)
(**)


Export[
    FileNameJoin[{outputDir, "run_summary.wxf"}],
    runSummary,
    "WXF"
];


(* ::Input:: *)
(**)


Print[""];
Print["============================================================"];
Print["Sweep complete."];
Print["Output directory: ", outputDir];
Print["Number of hx values: ", Length[hxList]];
Print["============================================================"];


(* ::Title:: *)
(*Plotting and export*)


(* ============================================================ *)
(* Load all completed hx runs                                   *)
(* ============================================================ *)

dataDir = FileNameJoin[
    {
        NotebookDirectory[],
        "prethermal_entropy_hx_sweep"
    }
];

files = Sort[
    FileNames[
        "run_" ~~ DigitCharacter ~~ DigitCharacter ~~ ".wxf",
        dataDir
    ]
];

runs = Import[#, "WXF"] & /@ files;


(* ============================================================ *)
(* Plot settings                                                *)
(* ============================================================ *)

plotDir = FileNameJoin[
    {
        dataDir,
        "plots"
    }
];

If[
    !DirectoryQ[plotDir],
    CreateDirectory[plotDir]
];


(* Continuous blue -> purple -> red convention in log10[hx].     *)
cf = Blend[
    {
        {0., Blue},
        {0.5, Purple},
        {1., Red}
    },
    #
] &;

hxExponentMin = Min[runs[[All, "hxExponent"]]];
hxExponentMax = Max[runs[[All, "hxExponent"]]];

colors = cf /@ Rescale[
    runs[[All, "hxExponent"]],
    {hxExponentMin, hxExponentMax}
];


(* Vertical color bar: hx = 10^x. *)
barTicks = Table[
    {
        x,
        Style[
            TraditionalForm[10^x],
            24,
            Black
        ]
    },
    {x, Ceiling[hxExponentMin], Floor[hxExponentMax], 1}
];

bar = BarLegend[
    {
        cf,
        {hxExponentMin, hxExponentMax}
    },
    LegendLabel -> Style[
        Subscript["h", "x"],
        30,
        Black
    ],
    LabelStyle -> {
        Black,
        20
    },
    ColorFunctionScaling -> True,
    LegendMarkerSize -> {12, 360},
    Ticks -> barTicks
];


(* ============================================================ *)
(* Raw data: no averaging, smoothing, interpolation, or joining  *)
(* ============================================================ *)

rawData = Table[
    Transpose[
        {
            Range[
                runs[[k]]["tmin"],
                runs[[k]]["tmax"],
                runs[[k]]["dt"]
            ],
            runs[[k]]["entropy"]
        }
    ],
    {k, Length[runs]}
];


(* ============================================================ *)
(* RAW: linear-time plot                                        *)
(* Bare points only                                             *)
(* ============================================================ *)

rawLinearPlot = Legended[
    ListPlot[
        rawData,
        Joined -> False,
        PlotStyle -> (
            Directive[
                #,
                Opacity[0.75],
                PointSize[0.002]
            ] & /@ colors
        ),
        PlotRange -> All,
        Frame -> True,
        Axes -> False,
        ImageSize -> 1000,
        PlotTheme -> "Detailed",
        FrameStyle -> Directive[
            Black,
            26
        ],
        PlotLabel -> Style[
            "RPS, L = " <> ToString[L] <>
            ", J = " <> ToString[J] <>
            ", hz = " <> ToString[hz],
            34,
            Black
        ],
        FrameLabel -> {
            Style["t", 30, Black],
            Style["S(t)", 30, Black]
        },
        GridLines -> None
    ],
    Placed[
        bar,
        Right
    ]
];

Export[
    FileNameJoin[
        {
            plotDir,
            "entropy_hx_raw_linear.png"
        }
    ],
    rawLinearPlot,
    ImageResolution -> 200
];


(* ============================================================ *)
(* RAW: logarithmic-time plot                                   *)
(* Bare points only                                             *)
(* ============================================================ *)

(* t = 0 cannot be represented on a logarithmic horizontal axis.
   The zero-time point is omitted only from this visualization;
   no entropy value is altered or processed. *)

rawDataPositiveTime = (
    Select[
        #,
        First[#] > 0. &
    ] &
) /@ rawData;

rawLogPlot = Legended[
    ListLogLinearPlot[
        rawDataPositiveTime,
        Joined -> False,
        PlotStyle -> (
            Directive[
                #,
                Opacity[0.75],
                PointSize[0.002]
            ] & /@ colors
        ),
        PlotRange -> All,
        Frame -> True,
        Axes -> False,
        ImageSize -> 1000,
        PlotTheme -> "Detailed",
        FrameStyle -> Directive[
            Black,
            26
        ],
        PlotLabel -> Style[
            "RPS, L = " <> ToString[L] <>
            ", J = " <> ToString[J] <>
            ", hz = " <> ToString[hz],
            34,
            Black
        ],
        FrameLabel -> {
            Style["t", 30, Black],
            Style["S(t)", 30, Black]
        },
        GridLines -> None
    ],
    Placed[
        bar,
        Right
    ]
];

Export[
    FileNameJoin[
        {
            plotDir,
            "entropy_hx_raw_logtime.png"
        }
    ],
    rawLogPlot,
    ImageResolution -> 200
];


(* ============================================================ *)
(* Moving-average visualization                                 *)
(* Logarithmic-time export only                                 *)
(* ============================================================ *)

movingAverageWindow = 200;

movingAverageData = (
    MovingAverage[
        #,
        movingAverageWindow
    ] &
) /@ rawData;

movingAverageLogPlot = Legended[
    ListLogLinearPlot[
        movingAverageData,
        Joined -> False,
        PlotStyle -> (
            Directive[
                #,
                Opacity[0.90],
                PointSize[0.0025]
            ] & /@ colors
        ),
        PlotRange -> All,
        Frame -> True,
        Axes -> False,
        ImageSize -> 1000,
        PlotTheme -> "Detailed",
        FrameStyle -> Directive[
            Black,
            26
        ],
        PlotLabel -> Style[
            "RPS, L = " <> ToString[L] <>
            ", J = " <> ToString[J] <>
            ", hz = " <> ToString[hz] <>
            ", moving average = " <> ToString[movingAverageWindow],
            34,
            Black
        ],
        FrameLabel -> {
            Style["t", 30, Black],
            Style["S(t)", 30, Black]
        },
        GridLines -> None
    ],
    Placed[
        bar,
        Right
    ]
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


(* ============================================================ *)
(* Do not display plots in the script                           *)
(* ============================================================ *)

Clear[
    rawLinearPlot,
    rawLogPlot,
    movingAverageLogPlot
];
