(* ::Package:: *)

BeginPackage["QMB`"];


(* ::Section::Closed:: *)
(*Notas*)


(* ::Text:: *)
(*Hay cosas en IAA_model.nb que 1) hay que migrar para ac\[AAcute] y 2) que hay que revisar si deber\[IAcute]a de poner ac\[AAcute]*)


(* ::Text:: *)
(*Hay cosas en los cuadernos del caometro donde hay rutinas para la secci\[OAcute]n de quantum chaos, como el unfolding etc*)


(* ::Text:: *)
(*Hay cosas de Heisenberg meets fuzzy que tambi\[EAcute]n tengo que pasar para ac\[AAcute]*)


(* ::Section::Closed:: *)
(*Usage definitions*)


(* ::Subsection::Closed:: *)
(*General quantum mechanics*)


DensityMatrix::usage = "DensityMatrix[\[Psi]] returns \!\(\*TemplateBox[{\"\[Psi]\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"\[Psi]\"},\n\"Bra\"]\).";


Pauli::usage= "Pauli[0-3] gives the Pauli matrices. 
Pauli[ { \!\(\*SubscriptBox[\(i\), \(1\)]\),\[Ellipsis],\!\(\*SubscriptBox[\(i\), \(n\)]\) } ] gives Pauli[\!\(\*SubscriptBox[\(i\), \(1\)]\)] \[CircleTimes] \[Ellipsis] \[CircleTimes] Pauli[\!\(\*SubscriptBox[\(i\), \(n\)]\)]."


MatrixPartialTrace::usage = "MatrixPartialTrace[mat, n, d] calculates the partial trace of mat over the nth subspace, where all subspaces have dimension d.
MatrixPartialTrace[mat, n, {\!\(\*SubscriptBox[\(d\), \(1\)]\),\!\(\*SubscriptBox[\(d\), \(1\)]\),\[Ellipsis]}] calculates the partial trace of matrix mat over the nth subspace, where mat is assumed to lie in a space constructed as a tensor product of subspaces with dimensions {d1,d2,\[Ellipsis]}.";


VectorFromKetInComputationalBasis::usage = "VectorFromKetInComputationalBasis[ket] returns the matrix representation of ket.";


KetInComputationalBasisFromVector::usage = "KetInComputationalBasisFromVector[vector] returns the ket representation in computational basis of vector.";


RandomQubitState::usage = "RandomQubitState[] returns a random qubit state.";


RandomChainProductState::usage = "RandomChainProductState[L] returns a random product state of L qubits.";


Dyad::usage = "Dyad[a] returns \!\(\*TemplateBox[{\"a\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"a\"},\n\"Bra\"]\).
Dyad[a,b] returns \!\(\*TemplateBox[{\"a\"},\n\"Ket\"]\)\!\(\*TemplateBox[{\"b\"},\n\"Bra\"]\).";


Commutator::usage="Commutator[A,B] returns A.B-B.A";
CommutationQ::usage="CommutationQ[A,B] yields True if A and B commute, and False otherwise.";
MutuallyCommutingSetQ::usage="MutuallyCommutingSetQ[ListOfMatrices] yields True if all matrices in the list mutually commute, and False otherwise.";


Braket::usage = "Braket[a,b] gives \!\(\*TemplateBox[{RowBox[{\"a\", \" \"}], RowBox[{\" \", \"b\"}]},\n\"BraKet\"]\).";


FixCkForStateEvoultion::usage = "FixCkForStateEvoultion[\!\(\*SubscriptBox[\(\[Psi]\), \(0\)]\), { \!\(\*TemplateBox[{SubscriptBox[\"E\", \"k\"]},\n\"Ket\"]\) }] fixes \!\(\*SubscriptBox[\(c\), \(k\)]\) = \!\(\*TemplateBox[{RowBox[{SubscriptBox[\"E\", \"k\"], \" \"}], RowBox[{\" \", SubscriptBox[\"\[Psi]\", \"0\"]}]},\n\"BraKet\"]\) for StateEvolution[]";


StateEvolution::usage = "StateEvolution[t, \!\(\*SubscriptBox[\(\[Psi]\), \(0\)]\), {E_i}, {\!\(\*TemplateBox[{\"E_i\"},\n\"Ket\"]\)} ] returns \!\(\*TemplateBox[{RowBox[{\"\[Psi]\", RowBox[{\"(\", \"t\", \")\"}]}]},\n\"Ket\"]\) = \!\(\*SubscriptBox[\(\[Sum]\), \(\(\\ \)\(i\)\)]\) \!\(\*SuperscriptBox[\(\[ExponentialE]\), \(\(-\[ImaginaryI]\)\\  \*SubscriptBox[\(E\), \(i\)]\\  t\)]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"], RowBox[{\" \", SubscriptBox[\"\[Psi]\", \"0\"]}]},\n\"BraKet\"]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"]},\n\"Ket\"]\).
StateEvolution[t, {\!\(\*SubscriptBox[\(E\), \(k\)]\)}] calculates \!\(\*TemplateBox[{RowBox[{\"\[Psi]\", RowBox[{\"(\", \"t\", \")\"}]}]},\"Ket\"]\) = \!\(\*SubscriptBox[\(\[Sum]\), \(\(\\\\\)\(i\)\)]\) \!\(\*SuperscriptBox[\(\[ExponentialE]\), \(\(-\[ImaginaryI]\)\\\\\*SubscriptBox[\(E\), \(i\)]\\\\t\)]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"], RowBox[{\" \", SubscriptBox[\"\[Psi]\", \"0\"]}]},\"BraKet\"]\)\!\(\*TemplateBox[{SubscriptBox[\"E\", \"i\"]},\"Ket\"]\) having fixed the \!\(\*SubscriptBox[\(c\), \(k\)]\)'s with FixCkForStateEvoultion[\!\(\*SubscriptBox[\(\[Psi]\), \(0\)]\), { \!\(\*TemplateBox[{SubscriptBox[\"E\", \"k\"]},\n\"Ket\"]\) }].";


BlochVector::usage = "BlochVector[\[Rho]] calculates \!\(\*SubscriptBox[\(r\), \(i\)]\) = Tr(\!\(\*SubscriptBox[\(\[Sigma]\), \(i\)]\)\[Rho]).";


KroneckerVectorProduct::usage = "KroneckerVectorProduct[a,b] calculates \!\(\*TemplateBox[{\"a\"},\n\"Ket\"]\)\[CircleTimes]\!\(\*TemplateBox[{\"b\"},\n\"Ket\"]\).";


Purity::usage = "Purity[\[Rho]] calculates the purity of \[Rho].";


qubit::usage = "Generates a state with the parametrization of the Bloch sphere (\[Theta],\[Phi])"


coherentstate::usage = "coherentstate[state,L] Generates a spin coherent state of L spins given a general single qubit state"


(* ::Subsection::Closed:: *)
(*Quantum chaos*)


(*buscar la rutina del unfolding para meterla aqu\[IAcute]. Quiz\[AAcute]s tambi\[EAcute]n las cosas de wigner dyson y poisson*)


MeanLevelSpacingRatio::usage = "MeanLevelSpacingRatio[\!\(\*
StyleBox[\"eigenvalues\",\nFontSlant->\"Italic\"]\)] gives \[LeftAngleBracket]\!\(\*SubscriptBox[\(r\), \(n\)]\)\[RightAngleBracket] of \!\(\*
StyleBox[\"eigenvalues\",\nFontSlant->\"Italic\"]\).";


(* ::Subsection::Closed:: *)
(*Quantum channels*)


Reshuffle::usage = "Reshuffle[m] applies the reshuffle transformation to the matrix m with dimension \!\(\*SuperscriptBox[\(d\), \(2\)]\)\[Times]\!\(\*SuperscriptBox[\(d\), \(2\)]\)."


(* ::Subsection::Closed:: *)
(*Spins chais*)


(* ::Subsubsection:: *)
(*Symmetries*)


SpinParityEigenvectors::usage = "SpinParityEigenvectors[L] gives a list of {even, odd} eigenvectors of the L-spin system parity operator P; P\!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"1\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"L\"]}]},\n\"Ket\"]\) = \!\(\*TemplateBox[{RowBox[{SubscriptBox[\"k\", \"L\"], \",\", \"\[Ellipsis]\", \",\", SubscriptBox[\"k\", \"1\"]}]},\n\"Ket\"]\), \!\(\*SubscriptBox[\(k\), \(i\)]\)=0,1.";


TranslationEigenvectorRepresentatives::usage = FormatUsage[
  "TranslationEigenvectorRepresentatives[L] returns a list of sublists, each containing:\
  the decimal representation of a bit-string representative eigenvector,\
  its pseudomomentum k, and the length of its translation orbit,\
  for a system of ```L``` qubits."
];


BlockDiagonalize::usage = FormatUsage[
	"BlockDiagonalize[matrix,opts] returns ```matrix``` in block-diagonal form. Only \
	option is '''Symmetry'''."
];


Transform::usage = FormatUsage[
	"BlockDiagonalize[matrix,opts] returns ```matrix``` in block-diagonal form. Only \
	option is '''Symmetry'''."
];


Symmetry::usage = FormatUsage[
	"Symmetry is an option for '''BlockDiagonalize''' to specify the symmetry according to which \
	```matrix``` in '''BlockDiagonalize'''[```matrix```] is block diagonalized. It takes the values: \
	\"Translation\". Default option is \"Translation\"."
];


(* ::Subsubsection::Closed:: *)
(*Spin chains*)


IsingHamiltonian::usage = FormatUsage[
	"IsingHamiltonian[h_x,h_z,J,L,opts] returns the Hamiltonian \ 
	H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - ```J``` \[Sum]_{*i=1*}^{*L-1*} \[Sigma]^z_i \[Sigma]^z_{*i+1*} \
	with boundary conditions specified by option BoundaryConditions (default is \"Open\")."
];


BoundaryConditions::usage = FormatUsage[
	"BoundaryConditions is an option for '''IsingHamiltonian''' to specify the boundary conditions. It \
	takes the values \"Open\" or \"Periodic\". Default option is \"Open\"."
];


IsingNNOpenHamiltonian::replaced = "Function `1` has been replaced by `2`.";


(*Quiet[
IsingNNOpenHamiltonian::usage = FormatUsage["IsingNNOpenHamiltonian[h_x,h_z,J,L] returns the Hamiltonian H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - ```J```\[Sum]_{*i=1*}^{*L-1*} \[Sigma]^z_i \[Sigma]^z_{*i+1*}.
IsingNNOpenHamiltonian[h_x,h_z,{J_1,...,J_L},L] returns the Hamiltonian H = \[Sum]_{*i=1*}^L (```h_x```\[Sigma]_i^x + ```h_z```\[Sigma]_i^z) - \[Sum]_{*i=1*}^{*L-1*} ```J_i``` \[Sigma]^z_i \[Sigma]^z_{*i+1*}."];
, {FrontEndObject::notavail, First::normal}];*)


IsingNNClosedHamiltonian::replaced = "Function `1` has been replaced by `2`.";


(*IsingNNClosedHamiltonian::usage = "IsingNNClosedHamiltonian[\!\(\*
StyleBox[SubscriptBox[\"h\", \"x\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[SubscriptBox[\"h\", \"z\"],\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"J\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)] returns the Hamiltonian \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(L\)]\)(\!\(\*SubscriptBox[\(h\), \(x\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(x\)]\) + \!\(\*SubscriptBox[\(h\), \(z\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\)) + \!\(\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), L]\) \!\(\*SubscriptBox[\(J\), \(i\)]\) \!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i\), \(z\)]\)\!\(\*SubsuperscriptBox[\(\[Sigma]\), \(i + 1\), \(z\)]\) with \!\(\*SubscriptBox[\(\[Sigma]\), \(L + 1\)]\) = \!\(\*SubscriptBox[\(\[Sigma]\), \(1\)]\).";*)


ClosedXXZHamiltonian::usage = "ClosedXXZHamiltonian[\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[CapitalDelta]\",\nFontSlant->\"Italic\"]\)] returns the closed XXZ 1/2-spin chain as in appendix A.1 of Quantum 8, 1510 (2024).";


OpenXXZHamiltonian::usage= "OpenXXZHamiltonian[\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[CapitalDelta]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"h1\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"h2\",\nFontSlant->\"Italic\"]\)] returns the open XXZ 1/2-spin chain as in appendix A.2 of Quantum 8, 1510 (2024).";


LeaSpinChainHamiltonian::usage = "LeaSpinChainHamiltonian[\!\(\*
StyleBox[\"Jxy\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"Jz\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Omega]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Epsilon]d\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"L\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\",\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"d\",\nFontSlant->\"Italic\"]\)] returns the spin-1/2 chain of eq. (1) in Am. J. Phys. 80, 246\[Dash]251 (2012).";


(* ::Section::Closed:: *)
(*Beginning of Package*)


Begin["`Private`"];


(* ::Section:: *)
(*Routine definitions*)


(*no poner los nombres de funciones p\[UAcute]blicas porque se joden la definici\[OAcute]n de uso*)
ClearAll[SigmaPlusSigmaMinus,SigmaMinusSigmaPlus,SigmaPlusSigmaMinus2,SigmaMinusSigmaPlus2]


(* ::Subsection::Closed:: *)
(*General quantum mechanics*)


DensityMatrix[\[Psi]_] := Outer[Times, \[Psi], Conjugate[\[Psi]]]


Pauli[0]=Pauli[{0}]=SparseArray[{{1,0}, {0,1}}]; 
Pauli[1]=Pauli[{1}]=SparseArray[{{0,1}, {1,0}}]; 
Pauli[2]=Pauli[{2}]=SparseArray[{{0,-I},{I,0}}]; 
Pauli[3]=Pauli[{3}]=SparseArray[{{1,0}, {0,-1}}];
Pauli[Indices_List] := KroneckerProduct @@ (Pauli /@ Indices)


VectorFromKetInComputationalBasis[ket_]:=Normal[SparseArray[FromDigits[ket,2]+1->1,Power[2,Length[ket]]]]


KetInComputationalBasisFromVector[vector_]:=IntegerDigits[Position[vector,1][[1,1]]-1,2,Log[2,Length[vector]]]


RandomQubitState[] := 
Module[{x,y,z,\[Theta],\[Phi]},
	{x,y,z} = RandomPoint[Sphere[]];
	{\[Theta],\[Phi]}={ArcCos[z],Sign[y]ArcCos[x/Sqrt[x^2+y^2]]};
	{Cos[\[Theta]/2],Exp[I \[Phi]]Sin[\[Theta]/2]}
]


RandomChainProductState[0] := {1}
RandomChainProductState[1] := RandomQubitState[]
RandomChainProductState[L_] := Flatten[KroneckerProduct@@Table[RandomQubitState[],L]]


Dyad[a_]:=Outer[Times,a,Conjugate[a]]
Dyad[a_,b_]:=Outer[Times,a,Conjugate[b]]


Commutator[A_,B_]:=A . B-B . A


ZeroMatrix[d_]:=ConstantArray[0,{d,d}]


CommutationQ[A_,B_]:=Commutator[A,B]==ZeroMatrix[Length[A]]


MutuallyCommutingSetQ[ListOfMatrices_]:=Module[{SetLength=Length[ListOfMatrices]},
AllTrue[Table[CommutationQ@@ListOfMatrices[[{i,j}]],{i,SetLength-1},{j,i+1,SetLength}],TrueQ,2]
]


Braket[a_,b_]:=Conjugate[a] . b


StateEvolution[t_,psi0_List,eigenvals_List,eigenvecs_List]:=
(*|\[Psi](t)\[RightAngleBracket] = Underscript[\[Sum], k] Subscript[c, k]\[ExponentialE]^(-Subscript[\[ImaginaryI]E, k]t)|Subscript[E, k]\[RightAngleBracket], Subscript[c, k]=\[LeftAngleBracket]Subscript[E, k]\[VerticalSeparator] Subscript[\[Psi], 0]\[RightAngleBracket]*)
	Module[{ck=Conjugate[eigenvecs] . psi0, newevals=Chop[eigenvals]},
		N[Total[ ck * Exp[-I*newevals*t] * eigenvecs]]
	]


FixCkForStateEvoultion[\[Psi]0_, eigenvecs_] :=
	Module[{},
		ck = N[ Chop[ Conjugate[eigenvecs] . \[Psi]0 ] ];
		Heigenvecs = eigenvecs;
	]


StateEvolution[t_,eigenvals_List]:=
(*|\[Psi](t)\[RightAngleBracket] = Underscript[\[Sum], k] Subscript[c, k]\[ExponentialE]^(-Subscript[\[ImaginaryI]E, k]t)|Subscript[E, k]\[RightAngleBracket], Subscript[c, k]=\[LeftAngleBracket]Subscript[E, k]\[VerticalSeparator] Subscript[\[Psi], 0]\[RightAngleBracket]*)
	N[Chop[Total[ ck * Exp[-I*eigenvals*N[t]] * Heigenvecs]]]


BlochVector[\[Rho]_]:=Chop[Tr[Pauli[#] . \[Rho]]&/@Range[3]]


KroneckerVectorProduct[a_,b_]:=Flatten[KroneckerProduct[a,b]]


Purity[\[Rho]_]:=Tr[\[Rho] . \[Rho]] 


qubit[\[Theta]_,\[Phi]_]:=FullSimplify[Normalize[{Cos[\[Theta]/2],Exp[I \[Phi]]Sin[\[Theta]/2]}]]


coherentstate[state_,L_]:=Flatten[KroneckerProduct@@Table[state,L]]


(* ::Subsection::Closed:: *)
(*Matrix partial trace*)


(* ::Input::Initialization:: *)
ClearAll[MatrixPartialTrace]
SyntaxInformation[MatrixPartialTrace] = {"ArgumentsPattern" -> {_, _, __}};
Options[MatrixPartialTrace]={Method->Automatic,"Verbose"->False};
(*Attributes[MatrixPartialTrace]={};*)
MatrixPartialTrace::usage="Compute a partial trace of a matrix.";


(* ::Input::Initialization:: *)
MatrixPartialTrace::notmat="The first argument '``' has to be a square matrix.";
MatrixPartialTrace::invmeth="Invalid value '``' of the Method option.";
MatrixPartialTrace::invdim="The third argument '``' must be a positive integer or a list of positive integers.";
MatrixPartialTrace::invdimint="The value of the third argument '``', representing the input dimension, is not consistent with the dimensions '``' of the matrix in the first argument.";
MatrixPartialTrace::incdim="The product '``' of subspace dimensions '``' in the third argument must be equal to the row/column dimension '``' of the matrix in the first argument.";
MatrixPartialTrace::invidx="The second argument '``' must be a nonzero integer or a list of nonzero integers or tokens such as All, Except or Span.";
MatrixPartialTrace::largeidx="The individual indices '``' in the second argument cannot exceed the total number '``' of subspace dimensions given in the third argument.";
MatrixPartialTrace::repidx="The indices '``' in the second argument cannot repeat.";
MatrixPartialTrace::invspan="Invalid Span specification '``' in the second argument.";
MatrixPartialTrace::invexc="Invalid Except specification '``' in the second argument.";


(* ::Input::Initialization:: *)
preprocessArguments[mat_,idx_,dims_]:=Module[{lmat=mat,lidx=idx,ldims=dims,matdims,subnum},

(* --- preprocess mat --- *)

(*if mat is not a matrix and/or is not a square matrix, throw an error*)
If[(!MatrixQ[lmat])||Unequal@@(matdims=Dimensions[lmat]),
ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::notmat,Shallow[lmat]];Return[$Failed]
];

(* --- preprocess dims --- *)

(*if dims is an (unevaluated) symbol, throw an error*)
If[Head[ldims]==Symbol,
ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::invdim,ldims];Return[$Failed]
];

(*if dims is a positive integer, try to create a list of subdimensions; it the resulting list is not valid, throw an error*)
If[IntegerQ[ldims],
If[Positive[ldims],
subnum=Log[ldims,First[matdims]];
If[!IntegerQ[subnum],
ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::invdimint,ldims,matdims];Return[$Failed];,
ldims=Table[ldims,subnum];
],
ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::invdim,ldims];Return[$Failed];
];
];

(*if dims is a list of positive integers (or the preprocessing above turned it into one), check whether the product of subdimensions gives the dimension of the input matrix and if not, throw an error; if dims is not a list of positive integers, throw an error*)
If[Head[ldims]==List&&TrueQ@AllTrue[ldims,IntegerQ[#]&&Positive[#]&],
If[Times@@ldims!=First[matdims],
ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::incdim,Times@@ldims,ldims,First@matdims];Return[$Failed];,
subnum=Length[ldims];
];
,ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::invdim,ldims];Return[$Failed];
];

(* --- preprocess idx --- *)

(*if idx is token All, expand it into an explicit list*)
If[lidx===All,lidx=Range[subnum]];

(*if idx is an integer, turn it into a list*)
If[IntegerQ[lidx],lidx={lidx}];

(*if idx is token Except, convert it into a list; if the resulting list is not valid, throw an error*)
If[Head[lidx]===Except,
lidx=List@@lidx;
If[MatchQ[lidx,{_Integer}|{{__Integer}}],
lidx=Select[Flatten[lidx],Abs[#]<=subnum&];
lidx=Complement[Range[subnum],Mod[lidx,subnum+1]];
,ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::invexc,lidx];Return[$Failed];
];
];

(*if idx is token Span, convert it into a list; if the resulting list is not valid, throw an error*)
If[Head[lidx]===Span,
lidx=List@@lidx;
lidx=Quiet@Check[Range@@lidx,$Failed];
If[lidx===$Failed,ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::invspan,lidx];Return[$Failed]];
];

(*if idx is not a list or the preprocessing above hasn't turned it into one, throw an error; otherwise do further preprocessing*)
If[Head[lidx]===List,

(*if idx is not empty and contains only nonzero integers, do further preprocessing; otherwise throw an error*)
If[lidx!={},
If[TrueQ@AllTrue[lidx,IntegerQ[#]&&(#!=0)&],

(*if idx contains integers that are out of range (both in positive and negative direction), throw an error; otherwise turn negative integers/indices into positive ones in range {1,...,number of subspaces} and sort them *)
If[!AllTrue[lidx,(Abs[#]<=subnum)&],
ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::largeidx,lidx,subnum];Return[$Failed]
];
lidx=Sort@Mod[lidx,subnum+1];

(*if some indices in idx repeat, throw an error*)
If[!DuplicateFreeQ[lidx],
ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::repidx,lidx];Return[$Failed]
];
,ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::invidx,lidx];Return[$Failed]
]];
,
ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::invidx,lidx];Return[$Failed]
];

(* --- return preprocessed arguments --- *)

{lmat,lidx,ldims,matdims,subnum}
]


(* ::Input::Initialization:: *)
ptraceSum[mat_,idx_,dims_]:=Module[{aux,lmat=mat,ldims=dims,pre,d,post},

(*partial trace over individual subsystems*)
Do[
(*for a given subsystem determine the dimension etc.*)
pre=Times@@ldims[[;;subidx-1]];
d=ldims[[subidx]];
post=Times@@ldims[[subidx+1;;]];

(*the resulting matrix is computed block-wise,
each block is a sum of appropriate blocks of the original matrix*)
aux=Table[
Sum[
lmat[[
post (d(is-1)+k-1)+1;;post(d(is-1)+k),
post (d(js-1)+k-1)+1;; post(d(js-1)+k)
]]
,{k,d}]
,{is,pre},{js,pre}
];

(*from block representation into 2D matrix representation*)
lmat=ArrayFlatten[aux];

(*drop the dimension over which we just traced over as a preparation for the next step*)
ldims=Drop[ldims,{subidx}];

,{subidx,Reverse[idx]}
];

lmat
]


(* ::Input::Initialization:: *)
ptraceTensor[mat_,idx_,dims_]:=Module[{lmat=mat,lidx=idx,resdim},

(*from 2D matrix representation into high.-dim. tensor representation*)
lmat=ArrayReshape[lmat,Flatten[{dims,dims}]];

(*the actual partial trace*)
lidx=Transpose[{idx,idx+Length[dims]}];
lmat=TensorContract[lmat,lidx];

(*in the special case, when partial trace is an actual trace that returns a number, turn the number into 2D array*)
If[Length[idx]==Length[dims],lmat=If[Head[mat]===SparseArray,SparseArray[{{lmat}}],{{lmat}}]];

(*from high.-dim. representation back into the 2D matrix representation*)
resdim=(Times@@dims)/(Times@@dims[[idx]]);
lmat=ArrayReshape[lmat,{resdim,resdim}];

(*return result*)
lmat
]


(* ::Input::Initialization:: *)
MatrixPartialTrace[mat_,idx_,dims_,OptionsPattern[]]:=Module[{method=OptionValue[Method],verbose=OptionValue["Verbose"],fun,lmat=mat,lidx=idx,ldims=dims,matdims,subnum,args},
(*
inputs:
 mat - input matrix; must be a square matrix
idx - an integer or a list of integers that index subspaces to be traced over; also tokens All, Except, and Span are allowed
dims - dimensions of individual subspaces entered either as a list or as a single integer; in the latter case all subspaces are assumed to have the same dimension and their number is deduced from the matrix mat dimensions;

outputs:
 a square matrix that is equal to the partial trace of the input matrix;

options:
 Method - possible choices: "TensorContract", "Sum", Automatic; "TensorContract" uses TensorContract internally and is fast for numerical matrices; "Sum" uses Sum internally and is fast for symbolic matrices, Automatic chooses "TensorContract" when the input matric is numerical and "Sum" when it is symbolic;
"Verbose" - possible choices: True, False; If True, then the summary of preprocessed input parameters is printed (together with the result);
*)

(* --- preprocess arguments --- *)

args=preprocessArguments[mat,idx,dims];
If[args===$Failed,Return[$Failed],{lmat,lidx,ldims,matdims,subnum}=args];

(* --- preprocess options --- *)

(*choose a particular method*)
method=If[method===Automatic,
If[MatrixQ[lmat,NumberQ],"TensorContract","Sum"],
method];

(*if "Verbose" is true, print the summary of input parameters*)
If[verbose,Print[Grid[
{{Style["Summary of input parameters",Bold],SpanFromLeft},
{"matrix dimensions",matdims},
{"numeric matrix",MatrixQ[lmat,NumberQ]},
{"method",Which[lidx=={},None,lidx==Range[subnum],"Tr",True,method]},
{"number of subspaces",subnum},
{"subspace dimensions",ldims},
{"subspace indices",lidx}},
Frame->All
]]];

(* --- choose the actual function and apply it --- *)

(*when no tracing is necessary, return the input matrix*)
If[lidx=={},Return[lmat]];

(*when tracing over all subspaces, skip successive partial tracing and apply instead directly Tr*)
If[lidx==Range[subnum],Return[If[Head[mat]===SparseArray,SparseArray,Identity]@{{Tr[lmat]}}]];

(*in all other cases, choose one of the two low-level functions and apply it*)
fun=Switch[method,
"TensorContract",ptraceTensor,
"Sum",ptraceSum,
_,ResourceFunction["ResourceFunctionMessage"][MatrixPartialTrace::invmeth,method];Return[$Failed]
];
fun[lmat,lidx,ldims]
]


(* ::Subsection::Closed:: *)
(*Quantum chaos*)


MeanLevelSpacingRatio[eigenvalues_]:=Mean[Min/@Transpose[{#,1/#}]&[Ratios[Differences[Sort[eigenvalues]]]]]


(* ::Subsection::Closed:: *)
(*Quantum channels*)


Reshuffle[m_] := ArrayFlatten[ArrayFlatten/@Partition[Partition[ArrayReshape[#,{Sqrt[Dimensions[m][[1]]],Sqrt[Dimensions[m][[1]]]}]&/@m,Sqrt[Dimensions[m][[1]]]],Sqrt[Dimensions[m][[1]]]],1];


(* ::Subsection:: *)
(*Spins*)


(* ::Subsubsection:: *)
(*Symmetries*)


SpinParityEigenvectors[L_]:=Module[{tuples,nonPalindromes,palindromes},
tuples=Tuples[{0,1},L];
nonPalindromes=Select[tuples,#!=Reverse[#]&];
palindromes=Complement[tuples,nonPalindromes];
nonPalindromes=DeleteDuplicatesBy[nonPalindromes,Sort[{#,Reverse[#]}]&];
Normal[
{
Join[SparseArray[FromDigits[#,2]+1->1.,2^L]&/@palindromes,Normalize[SparseArray[{FromDigits[#,2]+1->1.,FromDigits[Reverse[#],2]+1->1.},2^L]]&/@nonPalindromes],
Normalize[SparseArray[{FromDigits[#,2]+1->-1.,FromDigits[Reverse[#],2]+1->1.},2^L]]&/@nonPalindromes
}]
]


SpinParityEigenvectors[L_]:=Module[{tuples,nonPalindromes,palindromes},
tuples=Tuples[{0,1},L];
nonPalindromes=Select[tuples,#!=Reverse[#]&];
palindromes=Complement[tuples,nonPalindromes];
nonPalindromes=DeleteDuplicatesBy[nonPalindromes,Sort[{#,Reverse[#]}]&];
Normal[
{
Join[SparseArray[FromDigits[#,2]+1->1.,2^L]&/@palindromes,Normalize[SparseArray[{FromDigits[#,2]+1->1.,FromDigits[Reverse[#],2]+1->1.},2^L]]&/@nonPalindromes],
Normalize[SparseArray[{FromDigits[#,2]+1->-1.,FromDigits[Reverse[#],2]+1->1.},2^L]]&/@nonPalindromes
}]
]


Translation[state_, L_] := BitShiftRight[state] + BitAnd[state, 1]*2^(L-1)
BinaryNecklaces[L_Integer] := Module[{tuples=Tuples[{0,1},L]},
	Union[Table[First[Sort[NestList[RotateLeft, t, L-1]]], {t,tuples}]]
]
\[Omega][L_,k_] := Exp[2Pi*k*I/L]

TranslationEigenvectorRepresentatives[L_Integer] := Module[
	{
	necklaces = FromDigits[#,2]& /@ BinaryNecklaces[L],
	orbits
	},
	
	orbits = DeleteDuplicates[
	Sort[NestWhileList[Translation[#, L]&, #, UnsameQ[##]&, All][[;;-2]]]& /@ necklaces
	 ];
	Catenate[
		Outer[
			If[Mod[Length[#1]*#2, L]==0, {First[#1], #2, Length[#1]}, Nothing]&, 
			orbits, 
			Range[0, L-1], 
			1
		]
	]
]


(*RepresentativesOddBasis[basis_]:=DeleteDuplicatesBy[Discard[basis,PalindromeQ],Sort[{#,Reverse[#]}]&]*)
RepresentativesOddBasis[basis_]:=DeleteDuplicatesBy[Select[basis, Not@*PalindromeQ],Sort[{#,Reverse[#]}]&]
RepresentativesEvenBasis[basis_]:=DeleteDuplicatesBy[basis,Sort[{#,Reverse[#]}]&]


Options[BlockDiagonalize] = {
  Symmetry -> "Translation" (*coming soon: \"Parity\"*)
};

BlockDiagonalize[A_, opts:OptionsPattern[]]:= 
Switch[OptionValue[Symmetry],
	"Translation",
		Module[
			{
				L = Log[2, Length[A]],
				repsgathered,
				P
			},
				(*Gather translation eigenvectors by their pseudomomentum k*)
				repsgathered = GatherBy[TranslationEigenvectorRepresentatives[L], #[[2]]&];
				
				(*change of basis matrix*)
				P = SparseArray[
						Catenate[
							Apply[
								Normalize[SparseArray[
									Thread[
										FoldList[Translation[#, L]&, #1, Range[#3 - 1]] + 1 -> 
										Power[\[Omega][L, #2], Range[0., #3 - 1]]
									]
										, 2^L
								]]&,
								repsgathered,
								{2}
							]
						]
					];
					
				BlockDiagonalMatrix[Chop[Conjugate[P] . A . Transpose[P]]]
		],
	"Parity",
		Module[
		{
			L = Log[2, Length[A]],
			basis, reps, rules, Heven, Hodd
		},
			basis = Tuples[{0, 1}, L];
			rules=AssociationThread[basis->Range[Length[basis]]];
			(*reps=Comap[{RepresentativesEvenBasis, RepresentativesOddBasis},basis];*)
			reps=Through[{RepresentativesEvenBasis,RepresentativesOddBasis}[basis]];
			Heven=1/2 (A[[#1,#1]] + A[[#1,#2]] + A[[#2,#1]] + A[[#2,#2]] & @@ Map[rules, {#, Reverse/@#}&[reps[[1]]], {2}]);
			Hodd=1/2 (A[[#1,#1]] - A[[#1,#2]] - A[[#2,#1]] + A[[#2,#2]] & @@ Map[rules, {#, Reverse/@#}&[reps[[2]]], {2}]);
			(*Heven = # . Heven . #&[DiagonalMatrix[ReplacePart[ConstantArray[1.,Length[reps[[1]]]],Thread[Catenate[Position[reps[[1]],_?PalindromeQ,1]]->1/Sqrt[2.]]],TargetStructure->"Sparse"]];*)
			Heven = # . Heven . #&[SparseArray[DiagonalMatrix[ReplacePart[ConstantArray[1.,Length[reps[[1]]]],Thread[Catenate[Position[reps[[1]],_?PalindromeQ,1]]->1/Sqrt[2.]]]]]];
			Hodd=1/2 (A[[#1,#1]] - A[[#1,#2]] - A[[#2,#1]] + A[[#2,#2]] & @@ Map[rules, {#, Reverse/@#}&[reps[[2]]], {2}]);
			
			{Heven, Hodd}
		],
	_,
		Message[BlockDiagonalize::badSymmetry, OptionValue[Symmetry]];
		Return[$Failed];
]

(*Mensaje de error si la opci\[OAcute]n es inv\[AAcute]lida*)
BlockDiagonalize::badSymmetry = 
  "Option badSymmetry -> `1` is not valid. Valid options: \"Translation\".";


Options[Transform] = {
  Symmetry -> "Translation" (*coming soon: \"Parity\"*)
};

Transform[A_, opts:OptionsPattern[]]:= 
Switch[OptionValue[Symmetry],
	"Translation",
		Module[
			{
				L = Log[2, Length[A]],
				repsgathered,
				P
			},
				(*Gather translation eigenvectors by their pseudomomentum k*)
				repsgathered = GatherBy[TranslationEigenvectorRepresentatives[L], #[[2]]&];
				
				(*change of basis matrix*)
				P = SparseArray[
						Catenate[
							Apply[
								Normalize[SparseArray[
									Thread[
										FoldList[Translation[#, L]&, #1, Range[#3 - 1]] + 1 -> 
										Power[\[Omega][L, #2], Range[0., #3 - 1]]
									]
										, 2^L
								]]&,
								repsgathered,
								{2}
							]
						]
					];
					
				P
		],
	"Parity",
		Module[
		{
			L = Log[2, Length[A]],
			basis, reps, rules, Heven, Hodd
		},
			basis = Tuples[{0, 1}, L];
			rules=AssociationThread[basis->Range[Length[basis]]];
			(*reps=Comap[{RepresentativesEvenBasis, RepresentativesOddBasis},basis];*)
			reps=Through[{RepresentativesEvenBasis,RepresentativesOddBasis}[basis]];
			Heven=1/2 (A[[#1,#1]] + A[[#1,#2]] + A[[#2,#1]] + A[[#2,#2]] & @@ Map[rules, {#, Reverse/@#}&[reps[[1]]], {2}]);
			Hodd=1/2 (A[[#1,#1]] - A[[#1,#2]] - A[[#2,#1]] + A[[#2,#2]] & @@ Map[rules, {#, Reverse/@#}&[reps[[2]]], {2}]);
			(*Heven = # . Heven . #&[DiagonalMatrix[ReplacePart[ConstantArray[1.,Length[reps[[1]]]],Thread[Catenate[Position[reps[[1]],_?PalindromeQ,1]]->1/Sqrt[2.]]],TargetStructure->"Sparse"]];*)
			Heven = # . Heven . #&[SparseArray[DiagonalMatrix[ReplacePart[ConstantArray[1.,Length[reps[[1]]]],Thread[Catenate[Position[reps[[1]],_?PalindromeQ,1]]->1/Sqrt[2.]]]]]];
			Hodd=1/2 (A[[#1,#1]] - A[[#1,#2]] - A[[#2,#1]] + A[[#2,#2]] & @@ Map[rules, {#, Reverse/@#}&[reps[[2]]], {2}]);
			
			{Heven, Hodd}
		],
	_,
		Message[BlockDiagonalize::badSymmetry, OptionValue[Symmetry]];
		Return[$Failed];
]

(*Mensaje de error si la opci\[OAcute]n es inv\[AAcute]lida*)
BlockDiagonalize::badSymmetry = 
  "Option badSymmetry -> `1` is not valid. Valid options: \"Translation\".";


(* ::Subsubsection::Closed:: *)
(*Spin chains*)


Options[IsingHamiltonian] = {
  BoundaryConditions -> "Open"
};

IsingHamiltonian[hx_, hz_, J_, L_, opts:OptionsPattern[]] := Module[
	{NNIndices},
	NNIndices = Switch[OptionValue[BoundaryConditions],
		"Open",
			Normal[SparseArray[Thread[{#,#+1}->3],{L}]&/@Range[L-1]],
		"Periodic",
			Normal[SparseArray[Thread[{#,Mod[#+1,L,1]}->3],{L}]&/@Range[L]],
		_,
			Message[
                IsingHamiltonian::badBoundaryCondition, 
                OptionValue[BoundaryConditions]
            ];
            Return[$Failed];
	];
	
	Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]
]

(*Mensaje de error si la opci\[OAcute]n es inv\[AAcute]lida*)
IsingHamiltonian::badBoundaryCondition = 
  "Option BoundaryConditions -> `1` not valid. Valid options: \"Open\" o \"Periodic\".";


IsingNNOpenHamiltonian[args___] := Message[
	IsingNNOpenHamiltonian::replaced, "IsingNNOpenHamiltonian", "IsingHamiltonian"];
	
(*IsingNNOpenHamiltonian[hx_,hz_,J_,L_] := Module[{NNIndices},
	NNIndices=Normal[SparseArray[Thread[{#,#+1}->3],{L}]&/@Range[L-1]];
	N[Normal[Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]]]]*)


IsingNNClosedHamiltonian[args___] := Message[
	IsingNNClosedHamiltonian::replaced, "IsingNNClosedHamiltonian", "IsingHamiltonian"];

(*IsingNNClosedHamiltonian[hx_,hz_,J_,L_] := 
Module[{NNIndices},
	NNIndices=Normal[SparseArray[Thread[{#,Mod[#+1,L,1]}->3],{L}]&/@Range[L]];
	Total[{hx*Pauli[#]+hz*Pauli[3#]&/@IdentityMatrix[L],-J*(Pauli/@NNIndices)},2]
]*)


ClosedXXZHamiltonian[L_,\[CapitalDelta]_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L] ];
		N[Normal[-1/2*Total[Join[Pauli/@NNindices,Pauli/@(2NNindices),\[CapitalDelta] (Pauli[#]-IdentityMatrix[2^L])&/@(3NNindices)]]]]
	]


OpenXXZHamiltonian[L_,\[CapitalDelta]_,h1_,h2_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L-1] ];
		N[Normal[-1/2*Total[Join[Pauli/@NNindices,Pauli/@(2NNindices),\[CapitalDelta] (Pauli[#]-IdentityMatrix[2^L])&/@(3NNindices)]]  
		- 1/2*(h1 Pauli[Join[{1},ConstantArray[0,L-1]]] + h2*Pauli[Join[ConstantArray[0,L-1],{1}]])+ 1/2*(h1 + h2)IdentityMatrix[2^L]]]
	]


HamiltonianNN[Jxy_,Jz_,L_]:=
	Module[{NNindices},
		NNindices = Normal[ SparseArray[Thread[{#, Mod[# + 1, L, 1]}->1], {L}] &/@ Range[L-1] ];
		N[Normal[(1/4)*Total[Join[Jxy*(Pauli/@NNindices),Jxy*(Pauli/@(2NNindices)),Jz*(Pauli[#]&/@(3NNindices))]]]]
	]

HamiltonianZ[\[Omega]_,\[Epsilon]d_,L_,d_]:=N[(1/2)*(\[Omega]*Total[Pauli/@(3*IdentityMatrix[L])]+\[Epsilon]d*Pauli[Normal[SparseArray[d->3,L]]])]

LeaSpinChainHamiltonian[Jxy_,Jz_,\[Omega]_,\[Epsilon]d_,L_,d_]:=HamiltonianNN[Jxy,Jz,L]+HamiltonianZ[\[Omega],\[Epsilon]d,L,d]


(* ::Section::Closed:: *)
(*End of Package*)


End[];


EndPackage[];
