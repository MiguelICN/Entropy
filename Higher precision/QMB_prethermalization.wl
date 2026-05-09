(* ::Package:: *)

BeginPackage["QMBPrethermalization`"];


(* ================================================================ *)
(*  QMBPrethermalization: Minimal package for entanglement entropy   *)
(*  of the mixed-field Ising model with open boundary conditions.    *)
(*                                                                    *)
(*  Parity sectors are built exclusively via buildMaps/               *)
(*  makeExpansionMatrix, guaranteeing basis consistency.              *)
(*  BlockDiagonalize has been removed — it used an incompatible       *)
(*  basis ordering.                                                   *)
(* ================================================================ *)


(* ::Section:: *)
(* Public API *)

Pauli::usage = "Pauli[i] gives Pauli matrix i (0=Id,1=X,2=Y,3=Z). Pauli[{i1,...,iL}] gives the tensor product.";

IsingHamiltonian::usage = "IsingHamiltonian[hx,hz,J,L,BoundaryConditions->\"Open\"] returns H = Sum_i(hx*Sx_i + hz*Sz_i) + J*Sum_{<i,j>} Sz_i Sz_j.";
BoundaryConditions::usage = "Option for IsingHamiltonian. Values: \"Open\" (default) or \"Periodic\".";

revIndex::usage = "revIndex[i,L] returns the bit-reversal of integer i for L qubits.";
buildMaps::usage = "buildMaps[L] returns {evenMap, oddMap} for spatial-reflection parity sectors.";
makeExpansionMatrix::usage = "makeExpansionMatrix[map,L,type] returns the sparse isometry M mapping sector basis to full Hilbert space. type is \"even\" or \"odd\".";
ProjectBlockVector::usage = "ProjectBlockVector[state,map,type,L] projects a full-space state into a parity sector.";
ExpandBlockVector::usage = "ExpandBlockVector[v,map,type,L] expands a sector vector back to full Hilbert space.";

SectorHamiltonians::usage = "SectorHamiltonians[H,L] returns {Heven, Hodd} built via M^dag.H.M, guaranteeing basis consistency with buildMaps.";
DiagonalizeSector::usage = "DiagonalizeSector[Hsector] returns {eigenvalues, eigenvectorMatrix} sorted by eigenvalue. Rows of the matrix are eigenvectors.";

EvolveSector::usage = "EvolveSector[eigvecs,eigvals,sectorState,t] returns the time-evolved sector state at time t.";
EvolveSectorBatch::usage = "EvolveSectorBatch[eigvecs,eigvals,sectorState,tlist] returns a list of evolved sector states for each t in tlist. Uses batched BLAS multiply.";

EntanglementEntropySVD::usage = "EntanglementEntropySVD[psi,dA,dB] computes the von Neumann entropy of the reduced state of subsystem A via SVD.";

RandomQubitState::usage = "RandomQubitState[] returns a Haar-random single qubit state.";
RandomChainProductState::usage = "RandomChainProductState[L] returns a random product state of L qubits.";

PageEntropy::usage = "PageEntropy[La,Lb] returns the Page value for a random state in a Hilbert space of dimension 2^La * 2^Lb.";


(* ::Section:: *)
(* Implementation *)

Begin["`Private`"];


(* ---------------------------------------------------------------- *)
(* Pauli matrices and tensor products                                *)
(* ---------------------------------------------------------------- *)

Pauli[0] = Pauli[{0}] = SparseArray[{{1, 0}, {0, 1}}];
Pauli[1] = Pauli[{1}] = SparseArray[{{0, 1}, {1, 0}}];
Pauli[2] = Pauli[{2}] = SparseArray[{{0, -I}, {I, 0}}];
Pauli[3] = Pauli[{3}] = SparseArray[{{1, 0}, {0, -1}}];
Pauli[indices_List] := KroneckerProduct @@ (Pauli /@ indices)


(* ---------------------------------------------------------------- *)
(* Ising Hamiltonian                                                 *)
(* ---------------------------------------------------------------- *)

Options[IsingHamiltonian] = {BoundaryConditions -> "Open"};

IsingHamiltonian::badBC = "BoundaryConditions -> `1` is not valid. Use \"Open\" or \"Periodic\".";

IsingHamiltonian[hx_, hz_, J_, L_, opts : OptionsPattern[]] :=
  Module[{NNIndices},
    NNIndices = Switch[OptionValue[BoundaryConditions],
      "Open",
        Normal[SparseArray[Thread[{#, # + 1} -> 3], {L}] & /@ Range[L - 1]],
      "Periodic",
        Normal[SparseArray[Thread[{#, Mod[# + 1, L, 1]} -> 3], {L}] & /@ Range[L]],
      _,
        Message[IsingHamiltonian::badBC, OptionValue[BoundaryConditions]];
        Return[$Failed]
    ];
    Total[{hx Pauli[#] + hz Pauli[3 #] & /@ IdentityMatrix[L],
           J (Pauli /@ NNIndices)}, 2]
  ]


(* ---------------------------------------------------------------- *)
(* Spatial-reflection parity: sector maps and isometries             *)
(*                                                                   *)
(* The reflection operator P acts as:                                *)
(*   P |k_1, ..., k_L> = |k_L, ..., k_1>                           *)
(* which at the level of integers is bit-reversal.                   *)
(*                                                                   *)
(* Even sector: palindromes |i> and (|i>+|j>)/sqrt(2)               *)
(* Odd sector:  (|i>-|j>)/sqrt(2) for non-palindromic pairs         *)
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

makeExpansionMatrix[map_, l_, type_] :=
  SparseArray[
    Flatten[
      Table[
        If[Length[pair] == 1,
          (* Palindrome: single basis vector *)
          {{pair[[1]] + 1, k} -> 1.},
          (* Non-palindrome pair *)
          Which[
            type === "even",
              {{pair[[1]] + 1, k} -> 1./Sqrt[2.],
               {pair[[2]] + 1, k} -> 1./Sqrt[2.]},
            type === "odd",
              {{pair[[1]] + 1, k} -> 1./Sqrt[2.],
               {pair[[2]] + 1, k} -> -1./Sqrt[2.]}
          ]
        ],
        {k, Length[map]},
        {pair, {map[[k]]}}
      ], 2
    ],
    {2^l, Length[map]}
  ]

ExpandBlockVector[v_, map_, type_, l_] :=
  makeExpansionMatrix[map, l, type] . v

ProjectBlockVector[state_, map_, type_, l_] :=
  ConjugateTranspose[makeExpansionMatrix[map, l, type]] . state


(* ---------------------------------------------------------------- *)
(* Sector Hamiltonians via direct projection (replaces               *)
(* BlockDiagonalize to guarantee basis consistency)                  *)
(* ---------------------------------------------------------------- *)

SectorHamiltonians[H_, L_] :=
  Module[{evenMap, oddMap, Me, Mo, Heven, Hodd},
    {evenMap, oddMap} = buildMaps[L];
    Me = makeExpansionMatrix[evenMap, L, "even"];
    Mo = makeExpansionMatrix[oddMap, L, "odd"];
    Heven = ConjugateTranspose[Me] . H . Me;
    Hodd  = ConjugateTranspose[Mo] . H . Mo;
    (* Force Hermiticity to eliminate O(eps) asymmetry *)
    Heven = (Heven + ConjugateTranspose[Heven]) / 2;
    Hodd  = (Hodd + ConjugateTranspose[Hodd]) / 2;
    {Heven, Hodd}
  ]


(* ---------------------------------------------------------------- *)
(* Diagonalization with sorting and packing                          *)
(* ---------------------------------------------------------------- *)

DiagonalizeSector[Hsector_] :=
  Module[{vals, vecs, perm},
    {vals, vecs} = Eigensystem[N[Hsector]];
    perm = Ordering[vals];
    vals = Developer`ToPackedArray[N[vals[[perm]]]];
    vecs = Developer`ToPackedArray[N[vecs[[perm]]]];
    (* Rows of vecs are the eigenvectors *)
    {vals, vecs}
  ]


(* ---------------------------------------------------------------- *)
(* Time evolution in a parity sector                                 *)
(*                                                                   *)
(* |psi(t)> = Sum_k  c_k  e^{-i E_k t} |E_k>                       *)
(* where c_k = <E_k|psi_0> = Conjugate[vecs[[k]]] . psi_0           *)
(* and |psi(t)> in sector basis = Transpose[vecs] . (ck * phases)    *)
(* ---------------------------------------------------------------- *)

EvolveSector[eigvecs_, eigvals_, sectorState_, t_] :=
  Module[{ck},
    ck = Conjugate[eigvecs] . sectorState;
    Transpose[eigvecs] . (ck * Exp[-I eigvals t])
  ]

(* Batched version: evolve for a list of times using a single        *)
(* matrix-matrix multiply (DGEMM), much faster than looping.         *)
(* Returns a (Nt x dim_sector) matrix; row i is the state at t_i.   *)

EvolveSectorBatch[eigvecs_, eigvals_, sectorState_, tlist_] :=
  Module[{ck, phaseMat, evolvedSector},
    ck = Conjugate[eigvecs] . sectorState;
    (* phaseMat: dim_sector x Nt *)
    phaseMat = ck * Exp[-I Outer[Times, eigvals, tlist]];
    (* evolvedSector: dim_sector x Nt *)
    evolvedSector = Transpose[eigvecs] . phaseMat;
    (* Return as list of column vectors: Nt entries, each of length dim_sector *)
    Transpose[evolvedSector]
  ]


(* ---------------------------------------------------------------- *)
(* Entanglement entropy via SVD                                      *)
(* ---------------------------------------------------------------- *)

EntanglementEntropySVD[psi_, dA_, dB_] :=
  Module[{mat, sv2},
    mat = ArrayReshape[psi, {dA, dB}];
    sv2 = SingularValueList[mat]^2;
    (* Hard threshold at machine epsilon level *)
    sv2 = Select[sv2, # > 1.*^-14 &];
    (* Renormalize to absorb discarded numerical noise *)
    sv2 = sv2 / Total[sv2];
    -sv2 . Log[sv2]
  ]


(* ---------------------------------------------------------------- *)
(* Random product states                                             *)
(* ---------------------------------------------------------------- *)

RandomQubitState[] :=
  Module[{x, y, z, \[Theta], \[Phi]},
    {x, y, z} = RandomPoint[Sphere[]];
    {\[Theta], \[Phi]} = {ArcCos[z], Sign[y] ArcCos[x / Sqrt[x^2 + y^2]]};
    {Cos[\[Theta] / 2], Exp[I \[Phi]] Sin[\[Theta] / 2]}
  ]

RandomChainProductState[0] := {1}
RandomChainProductState[1] := RandomQubitState[]
RandomChainProductState[L_] := Flatten[KroneckerProduct @@ Table[RandomQubitState[], L]]


(* ---------------------------------------------------------------- *)
(* Page entropy                                                      *)
(* ---------------------------------------------------------------- *)

PageEntropy[La_, Lb_] :=
  N[PolyGamma[0, 2^La * 2^Lb + 1] - PolyGamma[0, 2^Lb + 1] -
    (2^La - 1) / (2 * 2^Lb)]


(* ---------------------------------------------------------------- *)
End[];
EndPackage[];
