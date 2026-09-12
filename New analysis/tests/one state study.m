(* ::Package:: *)

(* ::Title::Closed:: *)
(*Setup*)


(* ::Input:: *)
(*qmbInitPath = "C:\\Users\\Miguel\\Github\\libs\\QMB\\Kernel\\init.m";*)
(*Get[qmbInitPath];*)
(*SetDirectory[NotebookDirectory[]];*)


(* ::Input:: *)
(*LaunchKernels[12];*)


(* ::Input:: *)
(*Names["QMB`*"]*)


(* ::Title::Closed:: *)
(*Hamiltonian*)


(* ::Input:: *)
(*L=8;*)
(*J=1.;*)
(*hx=1.05;*)
(*hz=0.5;*)
(*H=IsingHamiltonian[hx,hz,J,L];*)
(*evals=Eigenvalues[Normal[H]];*)


(* ::Input:: *)
(*Histogram[evals,FrameLabel->{"E","DOS"}]*)


(* ::Input:: *)
(*newevals=Unfold[evals];*)
(*Histogram[newevals["UnfoldedLevels"],FrameLabel->{"E","DOS"}]*)


(* ::Title::Closed:: *)
(*State*)


(* ::Input:: *)
(*RandomChainProductState[3]*)


(* ::Title::Closed:: *)
(*Entanglement entropy*)


(* ::Input:: *)
(*Clear[pageEntropy];*)
(*pageEntropy[La_,Lb_]:=PolyGamma[0,(2^La)(2^Lb)+1]-PolyGamma[0,(2^Lb)+1]-((2^La)-1)/(2*(2^Lb))*)


(* ::Input:: *)
(*ClearAll[entanglement];*)
(*entanglement[psi_, LA_, L_] := Module[{mat, s, p},*)
(*    mat = ArrayReshape[psi, {2^LA, 2^(L - LA)}];*)
(*    s = SingularValueList[mat];*)
(*    p = s^2;*)
(*    p = Select[p, # > 10^-14 &];*)
(*    -Chop[Total[p Log[p]]]*)
(*];*)


(* ::Input:: *)
(*L=4;*)
(*ini=RandomChainProductState[L];*)


(* ::Input:: *)
(*entanglement[ini,1,L]*)


(* ::Title::Closed:: *)
(*Time evolution*)


(* ::Input:: *)
(*L=6;*)


(* ::Input:: *)
(*J=1.;*)
(*hx=1.05;*)
(*hz=0.5;*)


(* ::Input:: *)
(*H=IsingHamiltonian[hx,hz,J,L];*)


(* ::Input:: *)
(*{evals,evecs}=Transpose[Sort[Transpose[Eigensystem[Normal[H]]]]];*)


(* ::Input:: *)
(*psi0=RandomChainProductState[L];*)


(* ::Input:: *)
(*coeff=Conjugate[evecs] . psi0;*)


(* ::Input:: *)
(*Clear[psit];*)
(*psit[t_]:=Transpose[evecs] . (Exp[-I evals t] coeff);*)


(* ::Input:: *)
(*tlist=Table[i,{i,0,50,0.1}];*)


(* ::Input:: *)
(*data=ParallelTable[*)
(*With[{psi=Transpose[evecs] . (Exp[-I evals t] coeff)},*)
(*{t,entanglement[psi,L/2,L]}],*)
(*{t,tlist}];*)


(* ::Input:: *)
(*pagevalue=Plot[pageEntropy[L/2,L/2],{x,tlist[[1]],tlist[[-1]]},PlotStyle->Directive[Red,Dashed]];*)


(* ::Input:: *)
(*Show[ListPlot[data,PlotRange->All],pagevalue]*)


(* ::Title::Closed:: *)
(*Reflection operator*)


(* ::Input:: *)
(*ClearAll[reflectionOperator];*)
(*reflectionOperator[L_Integer]:=Module[{dim,perm},dim=2^L;*)
(*perm=Table[1+FromDigits[Reverse[IntegerDigits[j-1,2,L]],2],{j,dim}];*)
(*SparseArray[(#->1)&/@Transpose[{perm,Range[dim]}],{dim,dim}]];*)


(* ::Input:: *)
(*ClearAll[reflectionSectorBases];*)
(*reflectionSectorBases[L_Integer]:=Module[{dim,perm,plusRules={},minusRules={},np=0,nm=0,j,k},*)
(*dim=2^L;*)
(*perm=Table[1+FromDigits[Reverse[IntegerDigits[j-1,2,L]],2],{j,dim}];*)
(*Do[k=perm[[j]];*)
(*Which[j<k,np++;*)
(*AppendTo[plusRules,{j,np}->1/Sqrt[2]];*)
(*AppendTo[plusRules,{k,np}->1/Sqrt[2]];*)
(*nm++;*)
(*AppendTo[minusRules,{j,nm}->1/Sqrt[2]];*)
(*AppendTo[minusRules,{k,nm}->-1/Sqrt[2]],j==k,np++;*)
(*AppendTo[plusRules,{j,np}->1]],{j,dim}];*)
(*<|"Even"->SparseArray[plusRules,{dim,np}],"Odd"->SparseArray[minusRules,{dim,nm}]|>];*)


(* ::Input:: *)
(*J=1.;*)
(*hx=1.05;*)
(*hz=0.5;*)


(* ::Input:: *)
(*H=IsingHamiltonian[hx,hz,J,L];*)


(* ::Input:: *)
(*R=reflectionOperator[L];*)


(* ::Input:: *)
(*R//Dimensions*)


(* ::Input:: *)
(*Norm[H . R-R . H]*)


(* ::Input:: *)
(*bases=reflectionSectorBases[L];*)
(*Bp=bases["Even"];*)
(*Bm=bases["Odd"];*)


(* ::Input:: *)
(*Chop[ConjugateTranspose[Bp] . Bp]*)
(*Chop[ConjugateTranspose[Bm] . Bm]*)
(*Chop[ConjugateTranspose[Bp] . Bm]*)


(* ::Input:: *)
(*Hp=Chop[ConjugateTranspose[Bp] . H . Bp];*)
(*Hm=Chop[ConjugateTranspose[Bm] . H . Bm];*)


(* ::Input:: *)
(*Dimensions[H]*)
(*Dimensions[Hp]*)
(*Dimensions[Hm]*)


(* ::Input:: *)
(*ClearAll[sortedEigensystem];*)
(*sortedEigensystem[M_]:=Module[{evals,evecs,ord},*)
(*{evals,evecs}=Eigensystem[Normal[M]];*)
(*ord=Ordering[evals];*)
(*{evals[[ord]],evecs[[ord]]}];*)


(* ::Input:: *)
(*{evalsP,evecsP}=sortedEigensystem[Hp];*)
(*{evalsM,evecsM}=sortedEigensystem[Hm];*)


(* ::Input:: *)
(*(*Lift eigenvectors back to the full Hilbert space*)*)
(*fullEvecsP=(Bp . #)&/@evecsP;*)
(*fullEvecsM=(Bm . #)&/@evecsM;*)


(* ::Input:: *)
(*parityResidualP=Norm[R . #-#]&/@fullEvecsP;*)
(*parityResidualM=Norm[R . #+#]&/@fullEvecsM;*)


(* ::Input:: *)
(*Max[parityResidualP]*)
(*Max[parityResidualM]*)


(* ::Input:: *)
(*parityP=Chop[Conjugate[#] . (R . #)&/@fullEvecsP];*)
(*parityM=Chop[Conjugate[#] . (R . #)&/@fullEvecsM];*)


(* ::Input:: *)
(*Union[parityP]*)
(*Union[parityM]*)


(* ::Title::Closed:: *)
(*Full entropy evolution*)


(* ::Input:: *)
(*ClearAll[sortedEigensystem];*)
(*sortedEigensystem[M_]:=Module[{evals,evecs,ord},*)
(*{evals,evecs}=Eigensystem[Normal[M]];*)
(*ord=Ordering[evals];*)
(*{evals[[ord]],evecs[[ord]]}];*)


(* ::Input:: *)
(*L=8;*)
(*J=1.;*)
(*hz=1;*)


(* ::Input:: *)
(*hx=0.5;*)
(*H=IsingHamiltonian[hx,hz,J,L];*)


(* ::Input:: *)
(*Hp=Chop[ConjugateTranspose[Bp] . H . Bp];*)
(*Hm=Chop[ConjugateTranspose[Bm] . H . Bm];*)


(* ::Input:: *)
(*{evalsP,evecsP}=sortedEigensystem[Hp];*)
(*{evalsM,evecsM}=sortedEigensystem[Hm];*)


(* ::Input:: *)
(*psi0=RandomChainProductState[L];*)


(* ::Input:: *)
(*psi0P=ConjugateTranspose[Bp] . psi0;*)
(*psi0M=ConjugateTranspose[Bm] . psi0;*)


(* ::Input:: *)
(*Norm[psi0P]^2+Norm[psi0M]^2*)


(* ::Input:: *)
(*coeffP=Conjugate[evecsP] . psi0P;*)
(*coeffM=Conjugate[evecsM] . psi0M;*)


(* ::Input:: *)
(*ClearAll[psit];*)
(*psit[t_]:=Bp . (Transpose[evecsP] . (Exp[-I evalsP t] coeffP))+Bm . (Transpose[evecsM] . (Exp[-I evalsM t] coeffM));*)


(* ::Input:: *)
(*tlist=Table[i,{i,0,10000,1.}];*)
(*tlist//Length*)


(* ::Input:: *)
(*data=ParallelTable[With[{psi=psit[t]},{t,entanglement[psi,L/2,L]}],{t,tlist}];*)


(* ::Input:: *)
(*pagevalue=Plot[pageEntropy[L/2,L/2],{x,tlist[[1]],tlist[[-1]]},PlotStyle->Directive[Red,Dashed]];*)


(* ::Input:: *)
(*Show[ListPlot[data,PlotRange->All],pagevalue,PlotRange->All]*)


(* ::Input:: *)
(*Show[ListPlot[data,PlotRange->All],pagevalue,PlotRange->All]*)


(* ::Input:: *)
(*evals=Sort[Join[evalsP,evalsM]];*)
(*spectrum=SortBy[Join[({#,+1}&/@evalsP),({#,-1}&/@evalsM)],First];*)


(* ::Title::Closed:: *)
(*plots*)


(* ::Input:: *)
(*ListLogLinearPlot[{Transpose[{tlist,entropy[[q]]}],MovingAverage[Transpose[{tlist,entropy[[q]]}],200]},PlotStyle->{Directive[Gray,Opacity[0.10]],Directive[reds[[h]],Opacity[1]]},PlotRange->{{tlist[[1]],tlist[[-1]]},All},Joined->{True,False},ImageSize->800,PlotTheme->"Detailed",PlotTheme->"Detailed",FrameStyle->Directive[Black,20],PlotLabel->Style["PBC, L = "<>ToString[L]<>", J = "<>ToString[J]<>", h = 1",35,Black],FrameLabel->{Style["t",25,Black],Style["Svn(t)",25,Black]}]*)


(* ::Title:: *)
(*Full entropy evolution packed*)


(* ::Input:: *)
(*ClearAll[sortedEigensystem];*)
(*sortedEigensystem[M_]:=Module[{evals,evecs,ord},*)
(*{evals,evecs}=Eigensystem[Normal[M]];*)
(*ord=Ordering[evals];*)
(*{evals[[ord]],evecs[[ord]]}];*)


(* ::Input:: *)
(*L=8;*)
(*J=1.;*)
(*hz=1;*)


(* ::Input:: *)
(*hx=0.01;*)
(*H=IsingHamiltonian[hx,hz,J,L];*)


(* ::Input:: *)
(*Hp=Chop[ConjugateTranspose[Bp] . H . Bp];*)
(*Hm=Chop[ConjugateTranspose[Bm] . H . Bm];*)


(* ::Input:: *)
(*{evalsP,evecsP}=sortedEigensystem[Hp];*)
(*{evalsM,evecsM}=sortedEigensystem[Hm];*)


(* ::Input:: *)
(*psi0=RandomChainProductState[L];*)


(* ::Input:: *)
(*psi0P=ConjugateTranspose[Bp] . psi0;*)
(*psi0M=ConjugateTranspose[Bm] . psi0;*)


(* ::Input:: *)
(*coeffP=Conjugate[evecsP] . psi0P;*)
(*coeffM=Conjugate[evecsM] . psi0M;*)


(* ::Input:: *)
(*VP=Bp . Transpose[evecsP];*)
(*VM=Bm . Transpose[evecsM];*)


(* ::Input:: *)
(*Vfull=Join[VP,VM,2];*)
(*evalsFull=Join[evalsP,evalsM];*)
(*coeffFull=Join[coeffP,coeffM];*)


(* ::Input:: *)
(*Norm[Vfull . coeffFull-psi0]*)


(* ::Input:: *)
(*Vfull=Developer`ToPackedArray[N[Vfull]];*)
(*evalsFull=Developer`ToPackedArray[N[evalsFull]];*)
(*coeffFull=Developer`ToPackedArray[N[coeffFull]];*)


(* ::Input:: *)
(*Developer`PackedArrayQ/@{Vfull,evalsFull,coeffFull}*)


(* ::Input:: *)
(*V=Developer`ToPackedArray[Re[Vfull]];*)


(* ::Input:: *)
(*cR=Developer`ToPackedArray[Re[coeffFull]];*)
(*cI=Developer`ToPackedArray[Im[coeffFull]];*)
(**)
(*cosE=Developer`ToPackedArray[Cos[evalsFull]];*)
(*sinE=Developer`ToPackedArray[Sin[evalsFull]];*)


(* ::Input:: *)
(*ClearAll[phaseStepC];*)
(*phaseStepC=Compile[{{aR,_Real,1},{aI,_Real,1},{cosE,_Real,1},{sinE,_Real,1}},Module[{newR,newI},*)
(*newR=aR*cosE+aI*sinE;*)
(*newI=aI*cosE-aR*sinE;*)
(*{newR,newI}],CompilationTarget->"C",RuntimeOptions->"Speed"];*)


(* ::Input:: *)
(*ClearAll[entropyFromSVC];*)
(*entropyFromSVC=Compile[{{sv,_Real,1}},Module[{sum=0.,p=0.,i,n},*)
(*n=Length[sv];*)
(*For[i=1,i<=n,i++,p=sv[[i]]*sv[[i]];*)
(*If[p>1.*^-14,sum-=p*Log[p]];];*)
(*sum],CompilationTarget->"C",RuntimeOptions->"Speed"];*)


(* ::Input:: *)
(*L=8;*)


(* ::Input:: *)
(*(* ============================================================*)*)
(*(*Optimized single-state long-time entropy evolution*)*)
(*(* ============================================================*)*)
(*LA=L/2;*)
(*dA=2^LA;*)
(*dB=2^(L-LA);*)


(* ::Input:: *)
(*(*------------------------------------------------------------*)*)
(*(*Combine the already symmetry-resolved eigenbases*)*)
(*(*------------------------------------------------------------*)*)
(*VP=Bp . Transpose[evecsP];*)
(*VM=Bm . Transpose[evecsM];*)
(**)
(*Vfull=Developer`ToPackedArray[N@Join[VP,VM,2]];*)
(**)
(*evalsFull=Developer`ToPackedArray[N@Join[evalsP,evalsM]];*)
(**)
(*coeffFull=Developer`ToPackedArray[N@Join[coeffP,coeffM]];*)
(**)
(*(*Check reconstruction at t=0*)*)
(*Print["Initial-state reconstruction error = ",Norm[Vfull . coeffFull-psi0]];*)


(* ::Input:: *)
(*(*------------------------------------------------------------*)*)
(*(*This model should have a real eigenbasis*)*)
(*(*------------------------------------------------------------*)*)
(*Print["Max imaginary part of Vfull = ",Max[Abs[Im[Vfull]]]];*)
(*V=Developer`ToPackedArray[Re[Vfull]];*)


(* ::Input:: *)
(*(*------------------------------------------------------------*)*)
(*(*Initial coefficients and one-step phases*)*)
(*(*------------------------------------------------------------*)*)
(*cR=Developer`ToPackedArray[Re[coeffFull]];*)
(*cI=Developer`ToPackedArray[Im[coeffFull]];*)
(*cosE=Developer`ToPackedArray[Cos[evalsFull]];*)
(*sinE=Developer`ToPackedArray[Sin[evalsFull]];*)


(* ::Input:: *)
(*(*------------------------------------------------------------*)*)
(*(*C-compiled one-step eigenphase evolution*)*)
(*(*------------------------------------------------------------*)*)
(*ClearAll[phaseStepC];*)
(*phaseStepC=Compile[{{aR,_Real,1},{aI,_Real,1},{cosE,_Real,1},{sinE,_Real,1}},Module[{newR,newI},newR=aR*cosE+aI*sinE;*)
(*newI=aI*cosE-aR*sinE;*)
(*{newR,newI}],CompilationTarget->"C",RuntimeOptions->"Speed"];*)
(*(*------------------------------------------------------------*)*)
(*(*C-compiled entropy reduction*)*)
(*(*------------------------------------------------------------*)*)
(*ClearAll[entropyFromSVC];*)
(*entropyFromSVC=Compile[{{sv,_Real,1}},Module[{sum=0.,p=0.,i,n},n=Length[sv];*)
(*For[i=1,i<=n,i++,p=sv[[i]]*sv[[i]];*)
(*If[p>1.*^-14,sum-=p*Log[p]];];*)
(*sum],CompilationTarget->"C",RuntimeOptions->"Speed"];*)
(*(*------------------------------------------------------------*)*)
(*(*Periodic exact reseeding*)*)
(*(*------------------------------------------------------------*)*)
(*ClearAll[reseedC];*)
(*reseedC=Compile[{{cR,_Real,1},{cI,_Real,1},{energies,_Real,1},{t,_Real}},Module[{ct,st},ct=Cos[energies*t];*)
(*st=Sin[energies*t];*)
(*{cR*ct+cI*st,cI*ct-cR*st}],CompilationTarget->"C",RuntimeOptions->"Speed"];*)
(*(*------------------------------------------------------------*)*)
(*(*Entropy:SVD stays in optimized Wolfram numerical LA*)*)
(*(*------------------------------------------------------------*)*)
(*ClearAll[entropyFast];*)
(*entropyFast[psiR_,psiI_]:=Module[{mat,sv},mat=ArrayReshape[psiR+I psiI,{dA,dB}];*)
(*sv=SingularValueList[mat];*)
(*entropyFromSVC[sv]];*)


(* ::Input:: *)
(*(*------------------------------------------------------------*)*)
(*(*Long-time sequential evolution*)*)
(*(*------------------------------------------------------------*)*)
(*tmax=10^5;*)
(*resetEvery=100000;*)
(**)
(*aR=cR;*)
(*aI=cI;*)


(* ::Input:: *)
(*entropy=Table[*)
(*(*Computational-basis state*)*)
(*psiR=V . aR;*)
(*psiI=V . aI;*)
(*(*Half-chain von Neumann entropy*)*)
(*SvN=entropyFast[psiR,psiI];*)
(*(*Advance to t+1*)*)
(*If[t<tmax,If[Mod[t+1,resetEvery]==0,{aR,aI}=reseedC[cR,cI,evalsFull,N[t+1]],{aR,aI}=phaseStepC[aR,aI,cosE,sinE]]];*)
(*SvN,{t,0,tmax}];*)
(*entropy=Developer`ToPackedArray[entropy];*)


(* ::Input:: *)
(*ListPlot[Transpose[{Range[0,tmax],entropy}],PlotRange->All,PlotTheme->"Detailed",ImageSize->600]*)


(* ::Input:: *)
(*ListPlot[Transpose[{Range[0,tmax],entropy}],PlotRange->All,PlotTheme->"Detailed"]*)


(* ::Input:: *)
(*ListPlot[Transpose[{Range[0,tmax],entropy}],PlotRange->All,PlotTheme->"Detailed"]*)


(* ::Input:: *)
(*ListPlot[Transpose[{Range[0,tmax],entropy}],PlotRange->{{0,1000},{0,0.6}},PlotTheme->"Detailed",ImageSize->600]*)


(* ::Input:: *)
(*ListPlot[Transpose[{Range[0,tmax],entropy}],PlotRange->{{0,20},{0,0.6}},PlotTheme->"Detailed",ImageSize->600]*)
