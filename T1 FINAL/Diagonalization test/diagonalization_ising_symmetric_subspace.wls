#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Diagonalization Ising*)


(* ::Chapter::Closed:: *)
(*Even parity subspace*)


(*SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files/wls"];
Import["../../Mathematica_packages/QMB.wl"];*)


(*PositiveParitySubspaceBasis[L_]:=DeleteDuplicatesBy[Tuples[{0,1},L],First[Sort[{#,Reverse[#]}]]&]*)


(*time1=AbsoluteTime[];*)


(*L = 16;

base = PositiveParitySubspaceBasis[L];
tagsbase = Tag /@ base;

(*Esto seguro se puede optimizar, pero queda para despu\[EAcute]s*)
Hx = Module[{m, dimSubspace = Length[tagsbase], q, l, tmp1},
m = SparseArray[Prepend[
	Table[
		(* Encontrar estados \sum_j \dyad{0}{1}_j \ket{k} *)
		tmp1 = Select[Join[ConstantArray[j,L]-IdentityMatrix[L],If[PalindromeQ[j],{},ConstantArray[Reverse[j],L]-IdentityMatrix[L]]], And@@NonNegative[#]&];
		q = ReplaceAll[PalindromeQ[j],{True->0,False->1}];
		(* Encontrar la posici\[OAcute]n en la base y asignar un valor *)
		Total[
			Table[
				l=ReplaceAll[PalindromeQ[k], {True->0, False->1}];
				SparseArray[Flatten[Position[tagsbase, If[MemberQ[tagsbase, Tag[k]], Tag[k], Tag[Reverse[k]]]]]-> 2.^(-(q+l)/2), dimSubspace]
		, {k, tmp1}]
		]
	,{j,base[[2;;]]}]
	,SparseArray[{},dimSubspace]
	]
	];
Transpose[m]+m
];*)


(*time2=AbsoluteTime[];*)


(*Print["Construcci\[OAcute]n de Hx tom\[OAcute] "<>ToString[NumberForm[time2-time1,{Infinity,2}]]<>" s"];*)


(*Hz=DiagonalMatrix[Total[(-1)^#]&/@base, TargetStructure->"Sparse"];
JJ=DiagonalMatrix[Total[(-1)^Differences[#]]&/@base, TargetStructure->"Sparse"];
Hamiltonian[hx_,hz_,J_]:=hx*Hx+hz*Hz-J*JJ;*)


(*points = Tuples[Join[Flatten[Table[Subdivide[#[[i]],#[[i+1]],3][[1;;3]],{i,9}]],{2.95}],2]&[Subdivide[0.05,2.95,9]];

time1=AbsoluteTime[];

LaunchKernels[35];
Print["Launched kernels : ", Length[ParallelKernels[]]];
eigenvalues = 
	ParallelTable[
		Module[{hz, J, H},
		{hz, J} = i;
		H = Normal[Hamiltonian[1., hz, J]];
		Flatten[{{hz, J}, Eigenvalues[H]}]
		]
	, {i, points},
	DistributedContexts->Full, Method->"CoarsestGrained"];
CloseKernels[];

time2=AbsoluteTime[];
Print["Tiempo diagonalizaciones: "<>ToString[NumberForm[time2-time1,{Infinity,2}]]<>" s"];

Export["data_eigenvalues_wisniacki_L_"<>ToString[L]<>"_hx_1.csv",Prepend[eigenvalues,{"hz,J,E_1,E_2,...,E_{2^L}"}],"CSV"];*)


(* ::Subsection::Closed:: *)
(*Nueva idea*)


(*RepeatedTiming[ConstantArray[0.,{2^#,2^#}];]&/@Range[3,10]*)


(*base2=DeleteDuplicates[If[PalindromeQ[#],{#},Sort[{#,Reverse[#]}]]&/@Tuples[{0,1},L]]*)


(*Map[Identity[3]-ConstantArray[#, 3]&,base2,{2}]*)


(*(*tags*)
L = 3;

base = PositiveParitySubspaceBasis[L];
tagsbase = Tag /@ base;

AbsoluteTiming[
(*Esto seguro se puede optimizar, pero queda para despu\[EAcute]s*)
Hx = Module[{m ,dimSubspace = Length[tagsbase], q, l, tmp1},
m = SparseArray[Prepend[
	Table[
		(* Encontrar estados \sum_j \dyad{0}{1}_j \ket{k} *)
		tmp1 = Select[Join[ConstantArray[j,L]-IdentityMatrix[L],If[PalindromeQ[j],{},ConstantArray[Reverse[j],L]-IdentityMatrix[L]]], And@@NonNegative[#]&];
		q = ReplaceAll[PalindromeQ[j],{True->0,False->1}];
		(* Encontrar la posici\[OAcute]n en la base y asignar un valor *)
		Total[
			Table[
				l=ReplaceAll[PalindromeQ[k], {True->0, False->1}];
				SparseArray[Flatten[Position[tagsbase, If[MemberQ[tagsbase, Tag[k]], Tag[k], Tag[Reverse[k]]]]]-> 2.^(-(q+l)/2), dimSubspace]
		, {k, tmp1}]
		]
	,{j,base[[2;;]]}]
	,SparseArray[{},dimSubspace]
	]
	];
Transpose[m]+m
]
]*)


(*Position[Map[Tag,base2,{2}], Tag@{0,0,1}][[1,1]]*)


(* ::Chapter:: *)
(*Odd parity subspace*)


(*SetDirectory["/home/jadeleon/Documents/chaos_meets_channels/ja_files/wls"];*)
SetDirectory["/home/miguel.gonzalez/Documents/chaos_meet_channels/ja_files/wls"];(*emmy*)
Import["../../Mathematica_packages/QMB.wl"];


L = 16; (*hx = 1.*)


(* ::Text:: *)
(*Funciones para este wls:*)


(*No se si se utiliza en alguna parte:*)
(*Parity[d_, list_]:=Fold[BitFlip[#1 - 1, #2]+1&, d, list]*)

HxTotal[vec_, L_] := BitFlip[vec, #]& /@ Range[0,L-1]


t1 = AbsoluteTime[];


TotalSigmaX = Total[Pauli /@ IdentityMatrix[L]];

(* Elemento representativo de los elementos de la base impar. Por ej, {0,1} para el estado 1/Sqrt[2](|01\[RightAngleBracket]+|10\[RightAngleBracket]) *)
repElementOddBasis = DeleteDuplicatesBy[Select[Tuples[{0, 1}, L], # != Reverse[#]&], First[Sort[{#,Reverse[#]}]]&];

(* Reglas de asignacion, por terminar de pulir *)
rulesCompToSymm = Join[Thread[Flatten[{FromDigits[#, 2],FromDigits[Reverse[#], 2]}&/@repElementOddBasis]-> Flatten[Transpose[ConstantArray[Range[Length[repElementOddBasis]], 2]]]],Thread[Complement[Range[0,2^L-1],Flatten[{FromDigits[#, 2],FromDigits[Reverse[#], 2]}&/@repElementOddBasis]]->Nothing]];
rulesSymmToComp = Thread[Range[Length[repElementOddBasis]]->({FromDigits[#,2],FromDigits[Reverse[#],2]}&/@repElementOddBasis)];
rules = Thread[Complement[Range[0,2^L-1],Flatten[{FromDigits[#, 2],FromDigits[Reverse[#], 2]}&/@repElementOddBasis]]->Nothing];
dimHodd = Length[repElementOddBasis];
rules2 = Join[Flatten[{FromDigits[#,2] -> {FromDigits[#,2],FromDigits[Reverse[#],2]}, FromDigits[Reverse[#],2] -> {FromDigits[#,2],FromDigits[Reverse[#],2]}}&/@repElementOddBasis,1],rules];


TotalSigmaXOdd = SparseArray[
Flatten[
Table[
	With[{vec = i /. rulesSymmToComp},
		{ReplaceAll[{vec[[1]], #[[1]]}, rulesCompToSymm] -> (1/2)*Total[{1., -1., -1., 1.}*(TotalSigmaX[[#1, #2]] &@@@ Tuples[{vec, #} + 1])]}& /@ ReplaceAll[Select[HxTotal[vec[[1]], L], ReplaceAll[#, rulesCompToSymm] > i&], rules2]
	]
, {i, dimHodd - 1}]
, 2]
, {Length[repElementOddBasis], Length[repElementOddBasis]}];
TotalSigmaXOdd = TotalSigmaXOdd + ConjugateTranspose[TotalSigmaXOdd];


t2 = AbsoluteTime[];
Print["Construcci\[OAcute]n de TotalSigmaXOdd tom\[OAcute] "<>ToString[NumberForm[t2 - t1, {Infinity, 2}]]<>" s"];


(* ::Text:: *)
(*Guardar TotalSigmaXOdd en un CSV:*)


Export["../data/eigenenergies/ising/TotalSigmaXOdd_L_"<>ToString[L]<>".mx", TotalSigmaXOdd];


TotalSigmaZOdd = DiagonalMatrix[Total[(-1)^#]&/@repElementOddBasis, TargetStructure->"Sparse"];
TotalSigmaZNNodd = DiagonalMatrix[Total[(-1)^Differences[#]]&/@repElementOddBasis, TargetStructure->"Sparse"];
Hamiltonian[hx_, hz_, J_] := hx*TotalSigmaXOdd + hz*TotalSigmaZOdd - J*TotalSigmaZNNodd;


(*points = Tuples[Join[Flatten[Table[Subdivide[#[[i]],#[[i+1]],3][[1;;3]],{i,9}]],{2.95}],2]&[Subdivide[0.05,2.95,9]];*)
points = Table[{i,1},{i,0,0.5,0.01}];
faltan = Length[points]; 

LaunchKernels[30];
Print["Launched kernels : ", Length[ParallelKernels[]]];
DistributeDefinitions[Hamiltonian];
SetSharedVariable[faltan];
ParallelDo[
	Module[{hz, J, H, eigenvals},
	{hz, J} = i;
	t1 = AbsoluteTime[];
	H = Normal[Hamiltonian[1., hz, J]];
	eigenvals = Chop[Eigenvalues[H]];
	Export["../data/eigenenergies/ising/L_"<>ToString[L]<>"/odd/hx_1_hz_"<>ToString[NumberForm[hz, {Infinity, 2}]]<>"_J_"<>ToString[NumberForm[J, {Infinity, 2}]]<>".csv", Sort[eigenvals], "CSV"];
	t2 = AbsoluteTime[];
	Print["hz = "<>ToString[NumberForm[hz, {Infinity,2}]]<>", J ="<>ToString[NumberForm[J, {Infinity,2}]]<>": "<>ToString[NumberForm[t2 - t1, {Infinity, 0}]]<>" s. Faltan "<>ToString[faltan-=1]];
	]
, {i, points}
, Method->"CoarsestGrained"];
CloseKernels[];


Print["Termin\[EAcute]"];


(*{#,784/#}&/@Range[56,64]*)


(* ::Title::Closed:: *)
(*Apendice*)


(* ::Chapter:: *)
(*Revisando los otros t\[EAcute]rminos*)


(*(Pauli[{3,0}]+Pauli[{0,3}]).{0,1,1,0}*)


(*L=4;base = PositiveParitySubspaceBasis[L];rep = DeleteDuplicatesBy[Select[Tuples[{0, 1}, L], # != Reverse[#]&], First[Sort[{#,Reverse[#]}]]&];
MatrixForm[DiagonalMatrix[Total[(-1)^#]&/@ rep, TargetStructure->"Sparse"]]*)


(*DiagonalMatrix[Total[(-1)^Differences[#]]&/@rep, TargetStructure->"Sparse"]//MatrixForm*)


(*rep*)


(* ::Section:: *)
(*Vieja nueva*)


(*L=5;Hx = Total[Pauli/@IdentityMatrix[L]];
rep = DeleteDuplicatesBy[Select[Tuples[{0, 1}, L], # != Reverse[#]&], First[Sort[{#,Reverse[#]}]]&];
oddbasis = Normalize[SparseArray[{FromDigits[#, 2], FromDigits[Reverse[#], 2]}+1 -> {1., -1.}, 2^L]]&/@ rep;
MatrixForm[HxOdd1=Table[i.Hx.j, {i, oddbasis}, {j, oddbasis}] // Chop]*)


(* ::Chapter::Closed:: *)
(*Nueva implementaci\[OAcute]n*)


(*Parity[d_, list_]:=Fold[BitFlip[#1 - 1, #2]+1&, d, list]*)


(*(*L = 3; *)
MatrixForm[Hx = Total[Pauli/@IdentityMatrix[L]]];
rep = DeleteDuplicatesBy[Select[Tuples[{0, 1}, L], # != Reverse[#]&], First[Sort[{#,Reverse[#]}]]&];
HxOdd = ConstantArray[0, {#, #}, SparseArray]&[Length[rep]];
rulesCompToSymm = Join[Thread[Flatten[{FromDigits[#, 2],FromDigits[Reverse[#], 2]}&/@rep]-> Flatten[Transpose[ConstantArray[Range[Length[rep]], 2]]]],Thread[Complement[Range[0,2^L-1],Flatten[{FromDigits[#, 2],FromDigits[Reverse[#], 2]}&/@rep]]->Nothing]];
rulesSymmToComp = Thread[Range[Length[rep]]->({FromDigits[#,2],FromDigits[Reverse[#],2]}&/@rep)];
rules = Thread[Complement[Range[0,2^L-1],Flatten[{FromDigits[#, 2],FromDigits[Reverse[#], 2]}&/@rep]]->Nothing];
dimHodd = Length[rep];
rules2 = Join[Flatten[{FromDigits[#,2] -> {FromDigits[#,2],FromDigits[Reverse[#],2]}, FromDigits[Reverse[#],2] -> {FromDigits[#,2],FromDigits[Reverse[#],2]}}&/@rep,1],rules];*)


(*ClearAll[HxTotal];
HxTotal[vec_, L_] := BitFlip[vec, #]& /@ Range[0,L-1]

AbsoluteTiming[
HxOdd = SparseArray[
Flatten[
Table[
	With[{vec = i /. rulesSymmToComp},
		{ReplaceAll[{vec[[1]], #[[1]]}, rulesCompToSymm] -> (1/2)*Total[{1., -1., -1., 1.}*(Hx[[#1, #2]] &@@@ Tuples[{vec, #} + 1])]}& /@ ReplaceAll[Select[HxTotal[vec[[1]], L], ReplaceAll[#, rulesCompToSymm] > i&], rules2]
	]
, {i, dimHodd - 1}]
, 2]
, {Length[rep], Length[rep]}];
MatrixForm[HxOdd = HxOdd + ConjugateTranspose[HxOdd]];
]*)


(*HxOdd//MatrixForm*)


(*HxOdd==HxOdd1*)


(* ::Chapter::Closed:: *)
(*Recycle bin*)


(*AbsoluteTiming[
MatrixForm[
HxOdd = SparseArray[
Flatten[
Table[
	With[{vec = i /. rulesSymmToComp},
		{ReplaceAll[{vec[[1]], #[[1]]}, rulesCompToSymm] -> (1/2)*Total[{1., -1., -1., 1.}*(Hx[[#1, #2]] &@@@ Tuples[{vec, #} + 1])]}&/@ReplaceAll[HxTotal[vec[[1]], L], rules2]
		(*{ReplaceAll[{vec[[1]], #[[1]]}, rulesCompToSymm] -> (1/2)*Total[{1., -1., -1., 1.}*(Hx[[#1, #2]] &@@@ Tuples[{vec, #} + 1])]}& /@ ReplaceAll[Select[HxTotal[vec[[1]], L], ReplaceAll[#, rulesCompToSymm] > i&], rules2]*)
	]
, {i, dimHodd - 1}]
, 2]
, {Length[rep], Length[rep]}]
];
]*)
