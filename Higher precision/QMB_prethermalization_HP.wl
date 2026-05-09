(* ::Package:: *)

BeginPackage["QMBPrecision`"];

(* ================================================================ *)
(*  QMBPrecision: Arbitrary-precision package for entanglement       *)
(*  entropy of the mixed-field Ising model (OBC).                    *)
(*                                                                    *)
(*  All numerical functions accept or use the global variable         *)
(*  $WorkingPrecision. Expansion matrices use exact Sqrt[2].          *)
(*  No Compile, no ToPackedArray — pure arbitrary-precision.          *)
(*  Parity sectors via M^dag.H.M (basis-consistent).                 *)
(* ================================================================ *)

$WorkingPrecision::usage = "Global variable controlling numerical precision. Set before calling any function.";

(* --- Hamiltonian --- *)
Pauli::usage = "Pauli[i] or Pauli[{i1,...}].";
IsingHamiltonianExact::usage = "IsingHamiltonianExact[hx,hz,J,L] returns the Ising Hamiltonian as an exact symbolic sparse matrix (no N[]).";
BoundaryConditions::usage = "Option for IsingHamiltonianExact.";

(* --- Parity sectors --- *)
revIndex::usage = "revIndex[i,L].";
buildMaps::usage = "buildMaps[L].";
makeExpansionMatrixExact::usage = "makeExpansionMatrixExact[map,L,type] returns the isometry with exact 1/Sqrt[2] entries.";
SectorHamiltoniansHP::usage = "SectorHamiltoniansHP[H,L,wp] returns {Heven,Hodd} at precision wp.";
DiagonalizeSectorHP::usage = "DiagonalizeSectorHP[Hsector,wp] returns {vals,vecs} at precision wp. Rows of vecs are eigenvectors.";
ProjectStateHP::usage = "ProjectStateHP[state,map,type,L] projects state into a sector using exact isometry.";
ExpandStateHP::usage = "ExpandStateHP[v,map,type,L] expands sector vector to full space.";

(* --- Evolution --- *)
EvolveSectorHP::usage = "EvolveSectorHP[eigvecs,eigvals,ck,t] evolves at time t using precomputed ck.";
ComputeOverlapsHP::usage = "ComputeOverlapsHP[eigvecs,sectorState] computes ck = Conj[V].psi once.";

(* --- Entropy --- *)
EntanglementEntropySVDHP::usage = "EntanglementEntropySVDHP[psi,dA,dB,wp] computes S_vN at precision wp.";

(* --- States --- *)
HaarStateHP::usage = "HaarStateHP[dim,wp] returns a Haar-random state at precision wp.";
RandomChainProductStateHP::usage = "RandomChainProductStateHP[L,wp] returns a random product state at precision wp.";

(* --- Reference --- *)
PageEntropy::usage = "PageEntropy[La,Lb].";

(* --- Precision tracking --- *)
ReportPrecision::usage = "ReportPrecision[label,expr] prints min precision of expr.";


Begin["`Private`"];

(* ---------------------------------------------------------------- *)
(* Precision tracking                                                *)
(* ---------------------------------------------------------------- *)

ReportPrecision[label_String, expr_] := Module[{p},
  p = Min[Precision[expr]];
  Print["  PREC | ", label, ": ", If[p === Infinity, "Exact",
    If[p === MachinePrecision, "Machine (~15.9)",
      ToString[NumberForm[p, {5, 1}]]]]];
  p
];


(* ---------------------------------------------------------------- *)
(* Pauli matrices (exact)                                            *)
(* ---------------------------------------------------------------- *)

Pauli[0] = Pauli[{0}] = SparseArray[{{1, 0}, {0, 1}}];
Pauli[1] = Pauli[{1}] = SparseArray[{{0, 1}, {1, 0}}];
Pauli[2] = Pauli[{2}] = SparseArray[{{0, -I}, {I, 0}}];
Pauli[3] = Pauli[{3}] = SparseArray[{{1, 0}, {0, -1}}];
Pauli[indices_List] := KroneckerProduct @@ (Pauli /@ indices)


(* ---------------------------------------------------------------- *)
(* Ising Hamiltonian — exact symbolic (no N[])                       *)
(* ---------------------------------------------------------------- *)

Options[IsingHamiltonianExact] = {BoundaryConditions -> "Open"};
IsingHamiltonianExact::badBC = "Invalid BoundaryConditions `1`.";

IsingHamiltonianExact[hx_, hz_, J_, L_, opts : OptionsPattern[]] :=
  Module[{NNIndices},
    NNIndices = Switch[OptionValue[BoundaryConditions],
      "Open",
        Normal[SparseArray[Thread[{#, # + 1} -> 3], {L}] & /@ Range[L - 1]],
      "Periodic",
        Normal[SparseArray[Thread[{#, Mod[# + 1, L, 1]} -> 3], {L}] & /@ Range[L]],
      _,
        Message[IsingHamiltonianExact::badBC, OptionValue[BoundaryConditions]];
        Return[$Failed]
    ];
    (* Return exact sparse matrix — caller applies SetPrecision *)
    Total[{hx Pauli[#] + hz Pauli[3 #] & /@ IdentityMatrix[L],
           J (Pauli /@ NNIndices)}, 2]
  ]


(* ---------------------------------------------------------------- *)
(* Parity sector maps and EXACT isometries                           *)
(* ---------------------------------------------------------------- *)

revIndex[i_, l_] := FromDigits[Reverse[IntegerDigits[i, 2, l]], 2]

buildMaps[l_] :=
  Module[{visited, evenMap = {}, oddMap = {}, j},
    visited = ConstantArray[False, 2^l];
    Do[
      If[!visited[[i + 1]],
        j = revIndex[i, l];
        visited[[i + 1]] = True;
        visited[[j + 1]] = True;
        If[i == j,
          AppendTo[evenMap, {i}],
          AppendTo[evenMap, {i, j}];
          AppendTo[oddMap, {i, j}]
        ]
      ],
      {i, 0, 2^l - 1}
    ];
    {evenMap, oddMap}
  ]

(* Exact entries: 1 and 1/Sqrt[2], no floating point *)
makeExpansionMatrixExact[map_, l_, type_] :=
  SparseArray[
    Flatten[
      Table[
        If[Length[pair] == 1,
          {{pair[[1]] + 1, k} -> 1},
          Which[
            type === "even",
              {{pair[[1]] + 1, k} -> 1/Sqrt[2],
               {pair[[2]] + 1, k} -> 1/Sqrt[2]},
            type === "odd",
              {{pair[[1]] + 1, k} -> 1/Sqrt[2],
               {pair[[2]] + 1, k} -> -1/Sqrt[2]}
          ]
        ],
        {k, Length[map]},
        {pair, {map[[k]]}}
      ], 2
    ],
    {2^l, Length[map]}
  ]


(* ---------------------------------------------------------------- *)
(* Sector Hamiltonians at controlled precision                       *)
(* ---------------------------------------------------------------- *)

SectorHamiltoniansHP[H_, L_, wp_] :=
  Module[{evenMap, oddMap, Me, Mo, Heven, Hodd, Hn},
    {evenMap, oddMap} = buildMaps[L];
    Me = makeExpansionMatrixExact[evenMap, L, "even"];
    Mo = makeExpansionMatrixExact[oddMap, L, "odd"];
    (* Convert H to numerical at working precision *)
    Hn = SetPrecision[Normal[H], wp];
    Heven = SetPrecision[Normal[ConjugateTranspose[Me] . Hn . Me], wp];
    Hodd  = SetPrecision[Normal[ConjugateTranspose[Mo] . Hn . Mo], wp];
    (* Force Hermiticity *)
    Heven = (Heven + ConjugateTranspose[Heven]) / 2;
    Hodd  = (Hodd + ConjugateTranspose[Hodd]) / 2;
    {Heven, Hodd}
  ]


(* ---------------------------------------------------------------- *)
(* Diagonalization at controlled precision                           *)
(* ---------------------------------------------------------------- *)

DiagonalizeSectorHP[Hsector_, wp_] :=
  Module[{vals, vecs, perm},
    {vals, vecs} = Eigensystem[SetPrecision[Hsector, wp]];
    perm = Ordering[Re[vals]];
    vals = vals[[perm]];
    vecs = vecs[[perm]];
    (* Rows of vecs are eigenvectors *)
    {vals, vecs}
  ]


(* ---------------------------------------------------------------- *)
(* State projection / expansion using exact isometries               *)
(* ---------------------------------------------------------------- *)

ProjectStateHP[state_, map_, type_, l_] :=
  ConjugateTranspose[makeExpansionMatrixExact[map, l, type]] . state

ExpandStateHP[v_, map_, type_, l_] :=
  makeExpansionMatrixExact[map, l, type] . v


(* ---------------------------------------------------------------- *)
(* Time evolution at arbitrary precision                              *)
(*                                                                   *)
(* Pre-compute ck once, then evolve for each t.                      *)
(* No Compile, no BLAS — pure Mathematica arbitrary precision.        *)
(* ---------------------------------------------------------------- *)

ComputeOverlapsHP[eigvecs_, sectorState_] :=
  Conjugate[eigvecs] . sectorState

EvolveSectorHP[eigvecs_, eigvals_, ck_, t_] :=
  Transpose[eigvecs] . (ck * Exp[-I eigvals t])


(* ---------------------------------------------------------------- *)
(* Entanglement entropy via SVD at controlled precision              *)
(* ---------------------------------------------------------------- *)

EntanglementEntropySVDHP[psi_, dA_, dB_, wp_] :=
  Module[{mat, svs, sv2, thresh},
    mat = ArrayReshape[psi, {dA, dB}];
    svs = SingularValueList[mat];
    sv2 = svs^2;
    (* Threshold scales with working precision *)
    thresh = 10^(-wp + 2);
    sv2 = Select[sv2, Re[#] > thresh &];
    sv2 = sv2 / Total[sv2];
    Re[-sv2 . Log[sv2]]
  ]


(* ---------------------------------------------------------------- *)
(* Random states at controlled precision                             *)
(* ---------------------------------------------------------------- *)

HaarStateHP[dim_, wp_] :=
  Module[{re, im, state},
    re = Table[SetPrecision[RandomVariate[NormalDistribution[]], wp], {dim}];
    im = Table[SetPrecision[RandomVariate[NormalDistribution[]], wp], {dim}];
    state = re + I im;
    state / Sqrt[Conjugate[state] . state]
  ]

RandomChainProductStateHP[L_, wp_] :=
  Module[{qubits},
    qubits = Table[
      Module[{x, y, z, th, ph},
        {x, y, z} = SetPrecision[RandomPoint[Sphere[]], wp];
        th = ArcCos[z];
        ph = Sign[y] ArcCos[x / Sqrt[x^2 + y^2]];
        {Cos[th/2], Exp[I ph] Sin[th/2]}
      ],
      {L}
    ];
    Flatten[KroneckerProduct @@ qubits]
  ]


(* ---------------------------------------------------------------- *)
(* Page entropy                                                      *)
(* ---------------------------------------------------------------- *)

PageEntropy[La_, Lb_] :=
  N[PolyGamma[0, 2^La * 2^Lb + 1] - PolyGamma[0, 2^Lb + 1] -
    (2^La - 1) / (2 * 2^Lb), 20]


End[];
EndPackage[];
