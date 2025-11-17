(* ::Package:: *)

(* ::Section::Closed:: *)
(*Tabulation of thermal averages*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*InterThermalAverages[]:=*)
(*Module[{*)
(*dataC0k1=Import["files/k^1/dataC0,k^1.dat"],*)
(*dataC0k2=Import["files/k^2/dataC0,k^2.dat"],*)
(*dataC1k1=Import["files/k^1/dataC1,k^1.dat"],*)
(*dataC2k1=Import["files/k^1/dataC2,k^1.dat"],*)
(*dataC1k2=Import["files/k^2/dataC1,k^2.dat"],*)
(*dataC2k2=Import["files/k^2/dataC2,k^2.dat"],*)
(*dataC0k3=Import["files/k^3/dataC0,k^3.dat"],*)
(*dataC1k3=Import["files/k^3/dataC1,k^3.dat"],*)
(*dataC2k3=Import["files/k^3/dataC2,k^3.dat"],*)
(*m0k1,b0k1,m0k1NR,b0k1NR,m1k1,b1k1,m1k1NR,b1k1NR,m2k1,b2k1,m2k1NR,b2k1NR,*)
(*m0k2,b0k2,m1k2,b1k2,m2k2,b2k2,m0k2NR,b0k2NR,m1k2NR,b1k2NR,m2k2NR,b2k2NR,*)
(*m0k3,b0k3,m1k3,b1k3,m2k3,b2k3,m0k3NR,b0k3NR,m1k3NR,b1k3NR,m2k3NR,b2k3NR,*)
(*facC0 =  (9 Sqrt[3] )/(4096  (\[Pi]^2) ),facC1=(9 Sqrt[3] )/(8192  (\[Pi]^2) ),facC2=(9 Sqrt[3] )/(4096 (\[Pi]^2) )*)
(*},*)
(*Do[dataC0k2[[j]][[2]]=-Re[ToExpression[dataC0k2[[j]][[2]]]],{j,Length[dataC0k2]}];*)
(*Do[dataC1k2[[j]][[2]]=-Re[ToExpression[dataC1k2[[j]][[2]]]],{j,Length[dataC1k2]}];*)
(*Do[dataC2k2[[j]][[2]]=-Re[ToExpression[dataC2k2[[j]][[2]]]],{j,Length[dataC2k2]}];*)
(**)
(*Do[*)
(*dataC0k1[[j]][[2]]=10^(2Log10[dataC0k1[[j]][[1]]]-3Log10[BesselK[2,dataC0k1[[j]][[1]]]] +Log10[dataC0k1[[j]][[2]]]),*)
(*{j,Length[dataC0k1]}];*)
(*Do[*)
(*dataC0k2[[j]][[2]]=10^(2Log10[dataC0k2[[j]][[1]]]-3Log10[BesselK[2,dataC0k2[[j]][[1]]]] +Log10[dataC0k2[[j]][[2]]]),*)
(*{j,Length[dataC0k2]}];*)
(*Do[*)
(*dataC0k3[[j]][[2]]=10^(2Log10[dataC0k3[[j]][[1]]]-3Log10[BesselK[2,dataC0k3[[j]][[1]]]] +Log10[dataC0k3[[j]][[2]]]),*)
(*{j,Length[dataC0k3]}];*)
(**)
(*Do[*)
(*dataC1k1[[j]][[2]]=10^(4Log10[dataC1k1[[j]][[1]]]-3Log10[BesselK[2,dataC1k1[[j]][[1]]]] +Log10[dataC1k1[[j]][[2]]]),*)
(*{j,Length[dataC1k1]}];*)
(*Do[*)
(*dataC1k2[[j]][[2]]=10^(4Log10[dataC1k2[[j]][[1]]]-3Log10[BesselK[2,dataC1k2[[j]][[1]]]] +Log10[dataC1k2[[j]][[2]]]),*)
(*{j,Length[dataC1k2]}];*)
(*Do[*)
(*dataC1k3[[j]][[2]]=10^(4Log10[dataC1k3[[j]][[1]]]-3Log10[BesselK[2,dataC1k3[[j]][[1]]]] +Log10[dataC1k3[[j]][[2]]]),*)
(*{j,Length[dataC1k3]}];*)
(**)
(*Do[*)
(*dataC2k1[[j]][[2]]=10^(4Log10[dataC2k1[[j]][[1]]]-3Log10[BesselK[2,dataC2k1[[j]][[1]]]] +Log10[dataC2k1[[j]][[2]]]),*)
(*{j,Length[dataC2k1]}];*)
(*Do[*)
(*dataC2k2[[j]][[2]]=10^(4Log10[dataC2k2[[j]][[1]]]-3Log10[BesselK[2,dataC2k2[[j]][[1]]]] +Log10[dataC2k2[[j]][[2]]]),*)
(*{j,Length[dataC2k2]}];*)
(*Do[*)
(*dataC2k3[[j]][[2]]=10^(4Log10[dataC2k3[[j]][[1]]]-3Log10[BesselK[2,dataC2k3[[j]][[1]]]] +Log10[dataC2k3[[j]][[2]]]),*)
(*{j,Length[dataC2k3]}];*)
(**)
(*averagedC0k1=Interpolation[dataC0k1/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)
(*averagedC1k1=Interpolation[dataC1k1/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)
(*averagedC2k1=Interpolation[dataC2k1/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)
(*averagedC0k2=Interpolation[dataC0k2/.{x_,z_}->{x,Log10[Re[ToExpression[z]]]},Method->"Spline",InterpolationOrder->1];*)
(*averagedC1k2=Interpolation[dataC1k2/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)
(*averagedC2k2=Interpolation[dataC2k2/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)
(*averagedC0k3=Interpolation[dataC0k3/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)
(*averagedC1k3=Interpolation[dataC1k3/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)
(*averagedC2k3=Interpolation[dataC2k3/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)
(*m0k1=(averagedC0k1[dataC0k1[[1]][[1]]]-averagedC0k1[dataC0k1[[2]][[1]]])/(Log10[dataC0k1[[1]][[1]]]-Log10[dataC0k1[[2]][[1]]]);*)
(*b0k1=10^(averagedC0k1[dataC0k1[[1]][[1]]]-m0k1 Log10[dataC0k1[[1]][[1]]]);*)
(*m1k1=(averagedC1k1[dataC1k1[[1]][[1]]]-averagedC1k1[dataC1k1[[2]][[1]]])/(Log10[dataC1k1[[1]][[1]]]-Log10[dataC1k1[[2]][[1]]]);*)
(*b1k1=10^(averagedC1k1[dataC1k1[[1]][[1]]]-m1k1 Log10[dataC1k1[[1]][[1]]]);*)
(*m2k1=(averagedC2k1[dataC2k1[[1]][[1]]]-averagedC2k1[dataC2k1[[2]][[1]]])/(Log10[dataC2k1[[1]][[1]]]-Log10[dataC2k1[[2]][[1]]]);*)
(*b2k1=10^(averagedC2k1[dataC2k1[[1]][[1]]]-m2k1 Log10[dataC2k1[[1]][[1]]]);*)
(**)
(*(*m0k1NR=(averagedC0k1[dataC0k1[[-1]][[1]]]-averagedC0k1[dataC0k1[[-2]][[1]]])/(Log10[dataC0k1[[-1]][[1]]]-Log10[dataC0k1[[-2]][[1]]]);*)
(*b0k1NR=10^(averagedC0k1[dataC0k1[[-1]][[1]]]-m0k1NR Log10[dataC0k1[[-1]][[1]]]);*)
(*m1k1NR=(averagedC1k1[dataC1k1[[-1]][[1]]]-averagedC1k1[dataC1k1[[-2]][[1]]])/(Log10[dataC1k1[[-1]][[1]]]-Log10[dataC1k1[[-2]][[1]]]);*)
(*b1k1NR=10^(averagedC1k1[dataC1k1[[-1]][[1]]]-m1k1NR Log10[dataC1k1[[-1]][[1]]]);*)
(*m2k1NR=(averagedC2k1[dataC2k1[[-1]][[1]]]-averagedC2k1[dataC2k1[[-2]][[1]]])/(Log10[dataC2k1[[-1]][[1]]]-Log10[dataC2k1[[-2]][[1]]]);*)
(*b2k1NR=10^(averagedC2k1[dataC2k1[[-1]][[1]]]-m2k1NR Log10[dataC2k1[[-1]][[1]]]);*)*)
(*(*-----------*)*)
(*m0k2=(averagedC0k2[dataC0k2[[1]][[1]]]-averagedC0k2[dataC0k2[[2]][[1]]])/(Log10[dataC0k2[[1]][[1]]]-Log10[dataC0k2[[2]][[1]]]);*)
(*b0k2=10^(averagedC0k2[dataC0k2[[1]][[1]]]-m0k2 Log10[dataC0k2[[1]][[1]]]);*)
(*m1k2=(averagedC1k2[dataC1k2[[1]][[1]]]-averagedC1k2[dataC1k2[[2]][[1]]])/(Log10[dataC1k2[[1]][[1]]]-Log10[dataC1k2[[2]][[1]]]);*)
(*b1k2=10^(averagedC1k2[dataC1k2[[1]][[1]]]-m1k2 Log10[dataC1k2[[1]][[1]]]);*)
(*m2k2=(averagedC2k2[dataC2k2[[1]][[1]]]-averagedC2k2[dataC2k2[[2]][[1]]])/(Log10[dataC2k2[[1]][[1]]]-Log10[dataC2k2[[2]][[1]]]);*)
(*b2k2=10^(averagedC2k2[dataC2k2[[1]][[1]]]-m2k2 Log10[dataC2k2[[1]][[1]]]);*)
(**)
(*(*m0k2NR=(averagedC0k2[dataC0k2[[-1]][[1]]]-averagedC0k2[dataC0k2[[-2]][[1]]])/(Log10[dataC0k2[[-1]][[1]]]-Log10[dataC0k2[[-2]][[1]]]);*)
(*b0k2NR=10^(averagedC0k2[dataC0k2[[-1]][[1]]]-m0k2NR Log10[dataC0k2[[-1]][[1]]]);*)
(*m1k2NR=(averagedC1k2[dataC1k2[[-1]][[1]]]-averagedC1k2[dataC1k2[[-2]][[1]]])/(Log10[dataC1k2[[-1]][[1]]]-Log10[dataC1k2[[-2]][[1]]]);*)
(*b1k2NR=10^(averagedC1k2[dataC1k2[[-1]][[1]]]-m1k2NR Log10[dataC1k2[[-1]][[1]]]);*)
(*m2k2NR=(averagedC2k2[dataC2k2[[-1]][[1]]]-averagedC2k2[dataC2k2[[-2]][[1]]])/(Log10[dataC2k2[[-1]][[1]]]-Log10[dataC2k2[[-2]][[1]]]);*)
(*b2k2NR=10^(averagedC2k2[dataC2k2[[-1]][[1]]]-m2k2NR Log10[dataC2k2[[-1]][[1]]]);*)*)
(*(*-----------*)*)
(*m0k3=(averagedC0k3[dataC0k3[[1]][[1]]]-averagedC0k3[dataC0k3[[2]][[1]]])/(Log10[dataC0k3[[1]][[1]]]-Log10[dataC0k3[[2]][[1]]]);*)
(*b0k3=10^(averagedC0k3[dataC0k3[[1]][[1]]]-m0k3 Log10[dataC0k3[[1]][[1]]]);*)
(*m1k3=(averagedC1k3[dataC1k3[[1]][[1]]]-averagedC1k3[dataC1k3[[2]][[1]]])/(Log10[dataC1k3[[1]][[1]]]-Log10[dataC1k3[[2]][[1]]]);*)
(*b1k3=10^(averagedC1k3[dataC1k3[[1]][[1]]]-m1k3 Log10[dataC1k3[[1]][[1]]]);*)
(*m2k3=(averagedC2k3[dataC2k3[[1]][[1]]]-averagedC2k3[dataC2k3[[2]][[1]]])/(Log10[dataC2k3[[1]][[1]]]-Log10[dataC2k3[[2]][[1]]]);*)
(*b2k3=10^(averagedC2k3[dataC2k3[[1]][[1]]]-m2k3 Log10[dataC2k3[[1]][[1]]]);*)
(**)
(*(*m0k3NR=(averagedC0k3[dataC0k3[[-1]][[1]]]-averagedC0k3[dataC0k3[[-2]][[1]]])/(Log10[dataC0k3[[-1]][[1]]]-Log10[dataC0k3[[-2]][[1]]]);*)
(*b0k3NR=10^(averagedC0k3[dataC0k3[[-1]][[1]]]-m0k3NR Log10[dataC0k3[[-1]][[1]]]);*)
(*m1k3NR=(averagedC1k3[dataC1k3[[-1]][[1]]]-averagedC1k3[dataC1k3[[-2]][[1]]])/(Log10[dataC1k3[[-1]][[1]]]-Log10[dataC1k3[[-2]][[1]]]);*)
(*b1k3NR=10^(averagedC1k3[dataC1k3[[-1]][[1]]]-m1k3NR Log10[dataC1k3[[-1]][[1]]]);*)
(*m2k3NR=(averagedC2k3[dataC2k3[[-1]][[1]]]-averagedC2k3[dataC2k3[[-2]][[1]]])/(Log10[dataC2k3[[-1]][[1]]]-Log10[dataC2k3[[-2]][[1]]]);*)
(*b2k3NR=10^(averagedC2k3[dataC2k3[[-1]][[1]]]-m2k3NR Log10[dataC2k3[[-1]][[1]]]);*)*)
(**)
(*(*-------------*)*)
(*AveragedC0k1[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:= (facC0 k)/mphi^5 Piecewise[{{b0k1 xphi^m0k1,xphi<=dataC0k1[[1]][[1]]},{10^averagedC0k1[xphi],xphi>dataC0k1[[1]][[1]]}}];*)
(*AveragedC1k1[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:= (facC1 k)/mphi^5 (*xphi^4/BesselK[2,xphi]^3*)Piecewise[{{b1k1 xphi^m1k1,xphi<=dataC1k1[[1]][[1]]},{10^averagedC1k1[xphi],xphi>dataC1k1[[1]][[1]]}}];*)
(*AveragedC2k1[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:= (facC2 k)/mphi^5 (*xphi^4/BesselK[2,xphi]^3*)Piecewise[{{b2k1 xphi^m2k1,xphi<=dataC2k1[[1]][[1]]},{10^averagedC2k1[xphi],xphi>dataC2k1[[1]][[1]]}}];*)
(*(*---------------*)*)
(*AveragedC0k2[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:= (facC0 k^2)/mphi^5 (*xphi^2/BesselK[2,xphi]^3*)Piecewise[{{b0k2 xphi^m0k2,xphi<=dataC0k2[[1]][[1]]},{10^averagedC0k2[xphi],xphi>dataC0k2[[1]][[1]]}}];*)
(*AveragedC1k2[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:= (facC1 k^2)/mphi^5 (*xphi^4/BesselK[2,xphi]^3*)Piecewise[{{b1k2 xphi^m1k2,xphi<=dataC1k2[[1]][[1]]},{10^averagedC1k2[xphi],xphi>dataC1k2[[1]][[1]]}}];*)
(*AveragedC2k2[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:= (facC2 k^2)/mphi^5 (*xphi^4/BesselK[2,xphi]^3*)Piecewise[{{b2k2 xphi^m2k2,xphi<=dataC2k2[[1]][[1]]},{10^averagedC2k2[xphi],xphi>dataC2k2[[1]][[1]]}}];*)
(*(*---------------*)*)
(*AveragedC0k3[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:= (facC0 k^2)/mphi^5 (*xphi^2/BesselK[2,xphi]^3*)Piecewise[{{b0k3 xphi^m0k3,xphi<=dataC0k3[[1]][[1]]},{10^averagedC0k3[xphi],xphi>dataC0k3[[1]][[1]]}}];*)
(*AveragedC1k3[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:= (facC1 k^3)/mphi^5 (*xphi^4/BesselK[2,xphi]^3*)Piecewise[{{b1k3 xphi^m1k3,xphi<=dataC1k3[[1]][[1]]},{10^averagedC1k3[xphi],xphi>dataC1k3[[1]][[1]]}}];*)
(*AveragedC2k3[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:= (facC2 k^3)/mphi^5 (*xphi^4/BesselK[2,xphi]^3*)Piecewise[{{b2k3 xphi^m2k3,xphi<=dataC2k3[[1]][[1]]},{10^averagedC2k3[xphi],xphi>dataC2k3[[1]][[1]]}}];*)
(*averagedCif[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:=Re[AveragedC0k1[xphi,mphi,k]-AveragedC0k2[xphi,mphi,k]+AveragedC0k3[xphi,mphi,k]];*)
(*averagedC1if[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:=Re[AveragedC1k1[xphi,mphi,k]-AveragedC1k2[xphi,mphi,k]+AveragedC1k3[xphi,mphi,k]];*)
(*averagedC2if[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:=HeavisideTheta[xphi-0.5](Re[AveragedC2k1[xphi,mphi,k]-AveragedC2k2[xphi,mphi,k]+AveragedC2k3[xphi,mphi,k]]);*)
(*C2cannibal[xphi_?NumericQ,mphi_?NumericQ,k_?NumericQ]:=Re[2averagedC1if[xphi,mphi,k]-3averagedC2if[xphi,mphi,k]];*)
(*]*)


(* ::Input:: *)
(*InterThermalAverages[];//Quiet*)


(* ::Input:: *)
(**)


(* ::Section::Closed:: *)
(*Definition of H rate, entropy in equilibrium*)


(* ::Input:: *)
(*datag = Import["files/gEnergy.dat"];*)
(*datagentropy = Import["files/gEntropy.dat"];*)
(*For[i =0,i<=Length[datag]-1,i++;datag[[i]][[1]]=Log10[datag[[i]][[1]]]];*)
(*gg=Interpolation[datag,Method->"Spline"];*)
(*grel[T_]:=  gg[Log10[T]];*)
(*For[i =0,i<=Length[datagentropy]-1,i++;datagentropy[[i]][[1]]=Log10[datagentropy[[i]][[1]]]];*)
(*gentropy=Interpolation[datagentropy,InterpolationOrder->2,Method->"Spline"];*)
(*gent[xx_,mphi_]:=gentropy[Log10[mphi/xx]]*)
(*gentTemp[T_]:=gentropy[Log10[T]]*)
(*s[T_]:=(2*Pi^2*gentropy[Log10[T]]*T^3)/45;*)
(*ss[x_,mphi_]:=(2*Pi^2*gent[x,mphi])/45*mphi^3/x^3;*)
(*H[T_?NumericQ]:=(1.66*Sqrt[grel[T]](*10*))/(1.22089*10^19) T^2;*)
(*Hbar[T_?NumericQ]:=H[T]/(1+T/(3gentTemp[T])*gentTemp'[T]);*)
(*ssc[x_,mphi_]:=(2*Pi^2*106.69594335689291(mphi/x)^3)/45*)
(*averagedP4E3list=Import["files/intP4_over_E3.dat"];*)
(*averagedP4E3inter=Interpolation[averagedP4E3list/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)


(* ::Input:: *)
(*(*to calculate the relic density, remember that 0.12 = \[CapitalOmega]h^2= (mphi s Yobs)/\[Rho]c, Hence Yobs = (9.711808033846503`*^-48 gev^4)/(mphi(2*43*Pi^2/(11*45)((8.617333262145`*^-14)2.72548 gev)^3))*)*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*nphieqq[x_?NumericQ,mphi_?NumericQ]:=mphi^3 1/(x*2Pi^2)*BesselK[2,x];*)


(* ::Input:: *)
(*nphieq=Compile[{{T,_Real,0},{mphi,_Real,0}},*)
(*(mphi^2 T BesselK[2,mphi/T])/(2 \[Pi]^2),RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True];*)


(* ::Input:: *)
(*averagedP4E3ef[]:=Module[{m1,b1,m2,b2},*)
(*m1=(averagedP4E3inter[0.00015]-averagedP4E3inter[0.00013])/(Log10[0.00015]-Log10[0.00013]);*)
(*b1=averagedP4E3inter[0.00013]-m1 Log10[0.00013];*)
(*averagedP4E3e[xphi_?NumericQ(*,mphi_?NumericQ*)]:=*)
(*(*1/(2Pi^2*nphieq[xphi,mphi])*mphi^4*)Piecewise[*)
(*{*)
(*{10^averagedP4E3inter[xphi],0.00013<xphi},*)
(*{Abs[10^b1 xphi^m1],xphi<=0.00013}*)
(*}*)
(*];*)
(*]*)


(* ::Input:: *)
(*Y\[Gamma][T_]:=Zeta[3]/Pi^2*2*T^3/s[T]*)


(* ::Section::Closed:: *)
(*Definition of h->\[Phi]\[Phi] decay routines*)


(* ::Input:: *)
(*nphi=Compile[{{T,_Real,0},{m,_Real,0}},*)
(*(m^2 T BesselK[2,m/T])/(2 \[Pi]^2),RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True];*)


(* ::Input:: *)
(*\[CapitalGamma]htophiphii=Compile[{{T,_Real,0},{mphi,_Real,0},{Tc,_Real,0}},*)
(*Module[{mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},*)
(*(*Piecewise[{{*)mh/(16Pi^3) Sqrt[1-(4mphi^2)/mh^2]T BesselK[1,Re[mh/T]](*,T<=Tc},{0,T>Tc}}]*)*)
(*],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True*)
(*];*)
(*\[CapitalGamma]htophiphii2=Compile[{{T,_Real,0},{mphi,_Real,0},{Tc,_Real,0}},*)
(*Module[{mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},*)
(*Piecewise[{{Re[mh^2/(32Pi^3) Sqrt[1-(4mphi^2)/mh^2]  T BesselK[2,Re[mh/T]]],1<T<=Tc},{0,T>Tc||T<1}}]*)
(*]*)
(*,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True*)
(*];*)
(*\[CapitalGamma]htophiphii3=Compile[{{T,_Real,0},{mphi,_Real,0},{Tc,_Real,0}},*)
(*Module[{mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},*)
(*mh^2/(32Pi^3) Sqrt[1-(4mphi^2)/mh^2]T BesselK[2,mh/T]*)
(*]*)
(*,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True*)
(*];*)
(**)
(*\[CapitalGamma]htophiphi[T_?NumericQ,Tc_?NumericQ]:=\[CapitalGamma]htophiphii[T,0,Tc]*)
(*\[CapitalGamma]htophiphi2[T_?NumericQ,Tc_?NumericQ]:=\[CapitalGamma]htophiphii3[T,0,Tc]*)


(* ::Section::Closed:: *)
(*Definition of h->\[Phi]SS decay routines*)


(* ::Input:: *)
(*\[CapitalGamma]htophiSSC0[T_?NumericQ]:=Module[{mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},(mh^3 T BesselK[1,mh/T])/(1024 \[Pi]^5)]*)


(* ::Input:: *)
(*\[CapitalGamma]htophiSSC2[T_?NumericQ]:=Module[{mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},(mh^4 T BesselK[2,mh/T])/(3072 \[Pi]^5)]*)


(* ::Section::Closed:: *)
(*Definition of h->SS decay routines*)


(* ::Input:: *)
(*\[CapitalGamma]hToSS=Compile[{{T,_Real,0},{mS,_Real,0}},*)
(*Module[{mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},*)
(*mh/(16Pi^3) Sqrt[1-(4mS^2)/mh^2]T BesselK[1,mh/T]*)
(*],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True*)
(*];*)
(*\[CapitalGamma]hToSS2=Compile[{{T,_Real,0},{mS,_Real,0}},*)
(*Module[{mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},*)
(*mh^2/(32Pi^3) Sqrt[1-(4mS^2)/mh^2]  T BesselK[2,mh/T]*)
(*]*)
(*,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True*)
(*];*)
(**)


(* ::Section::Closed:: *)
(*Definition of hh->\[Phi]\[Phi] routine*)


(* ::Input:: *)
(*\[Sigma]annihilation:=Compile[{{s,_Real,0},{mphi,_Real,0},{lhphi,_Real,0},{v,_Real,0}},*)
(*Block[*)
(*{MH=125,*)
(*lH=0.13},*)
(*(Sqrt[-4 MH^2+s] (-MH^2+s+24 lH v^2)^2)/(64 \[Pi] (MH^2-s)^2 s Sqrt[-4 mphi^2+s]) lhphi^2]*)
(*,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True*)
(*];*)
(**)


(* ::Subsection:: *)
(*Definition annihilation averages*)


(* ::Input:: *)
(*C0aH[T_]:=Module[{(*MH=125*)MH=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},1/2 (MH^2 T^2 BesselK[1,Abs[MH]/T]^2)/(128 \[Pi]^5)];*)
(*C2aH[T_]:=Module[{MH=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},5 (E^(-((2 MH)/T)) MH^2 T^3)/(512 \[Pi]^4)];*)


(* ::Section::Closed:: *)
(*SS<->\[Phi]\[Phi]*)


(* ::Input:: *)
(*C0a[T_?NumericQ,m_?NumericQ]:=1/2! (m^2 T^2 BesselK[1,m/T]^2)/(128 \[Pi]^5);*)
(*C0av2[T_?NumericQ,m1_?NumericQ,m2_?NumericQ]:=1/2! (m1^2 T^2 BesselK[1,m1/T]^2)/(128 \[Pi]^5)-1/2 (E^(-((2 m1)/T)) m2^2 T^3)/(512 m1 \[Pi]^4);*)
(*C2a[xphi_?NumericQ,mphi_?NumericQ]:=mphi^5 c2ai[xphi];*)


(* ::Input:: *)
(*fAnn[T_?NumericQ,MH_?NumericQ,mphi_?NumericQ]:=1/2 T/(512 (\[Pi]^5) ) NIntegrate[(Sqrt[-4 MH^2+s] Sqrt[-4 mphi^2+s]  BesselK[1,Sqrt[s]/T])/Sqrt[s],{s,Max[4mphi^2,4MH^2],2Max[4mphi^2,4MH^2],3Max[4mphi^2,4MH^2],10Max[4mphi^2,4MH^2],Infinity}]*)
(*InterpolatingFann[mphi_?NumericQ,mS_?NumericQ]:=*)
(*Module[{C0ann,dat,dati,C0an},*)
(*dat=Table[{10^T//N,fAnn[10^T,mphi,mS]//N},{T,Log10[10^-5],Log10[10],(Log10[10]-Log10[10^-5])/800}];*)
(*dati=Interpolation[dat/.{x_,z_}->{x,Log10[z]},InterpolationOrder->6];*)
(*C0an[T_?NumericQ]:=10^dati[T];*)
(*C0ann[T_]:=Piecewise[{{C0an[T],T<10},{C0a[T,mphi],T>10}}];*)
(*Return[C0ann]*)
(*]*)


(* ::Section::Closed:: *)
(*Boson,Boson->Boson,\[Phi]*)


(* ::Input:: *)
(*\[CapitalGamma]boson0[T_?NumericQ]:=Module[{mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2],v=246,mZ = 91.2,mW=80.3},20 1/(Pi v^4) (mZ^2/12+mZ^2/36+mW^2/18+(mW^2 (8mW^2+mZ^2))/(18 mZ^2)+mW^2/36+2*(20mW^4-mW^2 mZ^2+mZ^4)/(36mZ^2))((E^(-((2 mh)/T)) T^6 (512 (mh/T)^(7/2)+(mh^(3/2) (992 mh+1225 T))/T^(5/2)))/(1024 \[Pi]^(7/2))-(mh^4 T^2 BesselK[2,(2 mh)/T])/\[Pi]^4)]*)


(* ::Input:: *)
(*\[CapitalGamma]boson2[T_?NumericQ]:=Module[{mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2],v=246,mZ = 91.2,mW=80.3},20 1/(Pi v^4) (mZ^2/12+mZ^2/36+mW^2/18+(mW^2 (8mW^2+mZ^2))/(18 mZ^2)+mW^2/36+2*(20mW^4-mW^2 mZ^2+mZ^4)/(36mZ^2))((E^(-((2 mh)/T)) (mh T)^(5/2) (512 mh^2+1632 mh T+2961 T^2))/(1024 \[Pi]^(7/2))-(mh^5 T^2 BesselK[3,(2 mh)/T])/\[Pi]^4)]*)


(* ::Section::Closed:: *)
(*Decay rates and lifetime*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)


(* ::Input:: *)
(*decaydata= Import["files/dec_rate.dat"];*)


(* ::Input:: *)
(*decaydataLinear = {};*)
(*Do[AppendTo[decaydataLinear, {10^ToExpression[StringDrop[j[[1]], -1]], 10^ToExpression[j[[2]]]}], {j, decaydata}];*)


(* ::Input:: *)
(*decaydataLinear = SortBy[decaydataLinear,First];*)
(*decaydataLinear =DeleteDuplicatesBy[decaydataLinear,First];*)


(* ::Input:: *)
(*decaydataLinearI = Interpolation[decaydataLinear,InterpolationOrder->1];*)


(* ::Input:: *)
(*\[CapitalGamma]\[Phi]toff[mphi_,\[Theta]_,mf_]:=Block[{v=246},*)
(*Piecewise[{{(mf^2 mphi)/(8Pi v^2) (1-((2mf)/mphi)^2)^(3/2) \[Theta]^2,mphi<decaydataLinear[[1]][[1]]},{\[Theta]^2 decaydataLinearI[mphi],mphi>decaydataLinear[[1]][[1]]}}]];*)


(* ::Input:: *)
(*ff[\[Beta]_]:=Piecewise[{{ArcSin[\[Beta]^(-1/2)]^2,\[Beta]>=1},{-1/4*(Log[(1+Sqrt[1+\[Beta]])/(1-Sqrt[1-\[Beta]])]-I*Pi)^2,\[Beta]<1}}]*)
(*FW[\[Beta]_]:=2+3\[Beta]+3\[Beta](2-\[Beta])*ff[\[Beta]];*)
(*Ff[\[Beta]_]:=-2\[Beta](1+(1-\[Beta])ff[\[Beta]]);*)
(*CC[mphi_]:=FW[(4*80.3^2)/mphi^2]+Ff[(4*0.10566^2)/mphi^2]+Ff[(4*1.7768^2)/mphi^2]*)
(*\[CapitalGamma]\[Phi]to\[Gamma]\[Gamma][mphi_,\[Theta]_]:=Block[{\[Alpha]=1/137,v=246},mphi/Pi*((mphi \[Alpha] \[Theta])/(16*Pi v))^2*Re[Conjugate[CC[mphi]]*CC[mphi]]];*)


(* ::Input:: *)
(*\[CapitalGamma]tot[mphi_,\[Theta]_]:=*)
(*Block[{mmu=0.10566,*)
(*me=1/2000},*)
(*Piecewise[{{\[CapitalGamma]\[Phi]toff[mphi,\[Theta],me]+\[CapitalGamma]\[Phi]toff[mphi,\[Theta],mmu]+\[CapitalGamma]\[Phi]to\[Gamma]\[Gamma][mphi,\[Theta]],2mmu<mphi},{\[CapitalGamma]\[Phi]toff[mphi,\[Theta],me]+\[CapitalGamma]\[Phi]to\[Gamma]\[Gamma][mphi,\[Theta]],2me<=mphi<2mmu},{\[CapitalGamma]\[Phi]to\[Gamma]\[Gamma][mphi,\[Theta]],mphi<2me}}]*)
(*]*)


(* ::Input:: *)
(*ltime[mphi_?NumericQ,\[Theta]_?NumericQ]:=(6.582*10^-25)/\[CapitalGamma]tot[mphi,\[Theta]](*1/dataElineari[mphi]*)(*1/31536000*1/(14*10^9)*);*)
(**)


(* ::Section::Closed:: *)
(*Definition of SS<->ff*)


(* ::Input:: *)
(*f[T_?NumericQ,mS_?NumericQ,mphi_?NumericQ,\[Theta]_?NumericQ,mf_]:=Block[{gf=2},gf T/(1024 \[Pi]^5) (mf/246)^2 NIntegrate[((-4 mf^2+s)^(3/2) Sqrt[-4 mS^2+s]  BesselK[1,Sqrt[s]/T])/( Sqrt[s] ((mphi^2-s)^2+mphi^2 \[CapitalGamma]tot[mphi,\[Theta]]^2)),{s,Max[4mS^2,4mf^2],1.2Max[4mS^2,4mf^2],2Max[4mS^2,4mf^2],10Max[4mS^2,4mf^2],20Max[4mS^2,4mf^2],100Max[4mS^2,4mf^2],Infinity}(*,WorkingPrecision->30,PrecisionGoal->25,AccuracyGoal->25*)]];*)


(* ::Input:: *)
(*TabulatingC0[mS_?NumericQ,mphi_?NumericQ,\[Theta]_?NumericQ]:=*)
(*Module[{AnnAverg,AnnAvergI,AnnAvergF},*)
(*AnnAverg=Table[{10^T//N,f[10^T,mS,mphi,\[Theta],1/2000]+f[10^T,mS,mphi,\[Theta],0.10566]+f[10^T,mS,mphi,\[Theta],1.77]},{T,Log10[mS/180],Log10[151],(Log10[151]-Log10[mS/180])/100}];*)
(*AnnAvergI=Interpolation[AnnAverg/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->2];*)
(*AnnAvergF[T_]:=10^AnnAvergI[T];*)
(*Return[AnnAvergF]*)
(*]*)


(* ::Input:: *)
(*f2[T_?NumericQ,mS_?NumericQ,mphi_?NumericQ,\[Theta]_?NumericQ,mf_]:=Block[{gf=2},gf  T/(2048 (\[Pi]^5) ) (mf/246)^2 NIntegrate[((-4 mf^2+s)^(3/2) Sqrt[-4 mS^2+s] BesselK[2,Sqrt[s]/T])/((mphi^2-s)^2+mphi^2 \[CapitalGamma]tot[mphi,\[Theta]]^2),{s,Max[4mS^2,4mf^2],1.2Max[4mS^2,4mf^2],2Max[4mS^2,4mf^2],10Max[4mS^2,4mf^2],20Max[4mS^2,4mf^2],100Max[4mS^2,4mf^2],Infinity}(*,WorkingPrecision->30,PrecisionGoal->25,AccuracyGoal->25*)]];*)


(* ::Input:: *)
(*f2NR[T_?NumericQ,mS_?NumericQ,mphi_?NumericQ,\[Theta]_?NumericQ,mf_]:=Module[{gf=2},gf T/(2048Pi^5) (mf/246)^2 1/( 2mS) NIntegrate[((-4 mf^2+s)^(3/2) Sqrt[-4 mS^2+s] (Sqrt[s] (-4 mS^2+s) BesselK[1,Sqrt[s]/T]+4 (-mS^2+s) T BesselK[2,Sqrt[s]/T]))/( s ((mphi^2-s)^2+mphi^2 \[CapitalGamma]tot[mphi,\[Theta]]^2)),{s,Max[4mS^2,4mf^2],1.2Max[4mS^2,4mf^2],2Max[4mS^2,4mf^2],10Max[4mS^2,4mf^2],20Max[4mS^2,4mf^2],100Max[4mS^2,4mf^2],Infinity}(*,WorkingPrecision->30,PrecisionGoal->25,AccuracyGoal->25*)]];*)


(* ::Input:: *)
(*TabulatingC2[mS_?NumericQ,mphi_?NumericQ,\[Theta]_?NumericQ]:=*)
(*Module[{AnnAverg2,AnnAvergI2,AnnAverg2NR,AnnAvergI2NR,AnnAvergF2f},*)
(*AnnAverg2=Table[{10^T//N,f2[10^T,mS,mphi,\[Theta],1/2000]+f2[10^T,mS,mphi,\[Theta],0.10566]+f2[10^T,mS,mphi,\[Theta],1.77]},{T,Log10[2mS],Log10[151],(Log10[151]-Log10[2mS])/80}];*)
(*AnnAvergI2=Interpolation[AnnAverg2/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)
(*AnnAverg2NR=Table[{10^T//N,f2NR[10^T,mS,mphi,\[Theta],1/2000]+f2NR[10^T,mS,mphi,\[Theta],0.10566]+f2NR[10^T,mS,mphi,\[Theta],1.77]},{T,Log10[mS/180],Log10[2mS],(Log10[2mS]-Log10[mS/180])/80}];*)
(*AnnAvergI2NR=Interpolation[AnnAverg2NR/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->1];*)
(*AnnAvergF2f[T_]:=Piecewise[{{10^AnnAvergI2[2mS/3]/10^AnnAvergI2NR[2mS/3] 10^AnnAvergI2NR[T],T<=2/3 mS},{10^AnnAvergI2[T],T>2/3 mS}}];*)
(*Return[AnnAvergF2f]*)
(*]*)


(* ::Section::Closed:: *)
(*Definition of annihilation and production functions for \[Phi]*)


(* ::Input:: *)
(*ProductionC0[xx_?NumericQ,mS_?NumericQ]:=\[CapitalGamma]hhToh\[Phi]C0[mS/xx]+C0aH[mS/xx]+\[CapitalGamma]htophiphi[mS/xx,150];*)


(* ::Input:: *)
(*AnnihilationC0[xphi_?NumericQ,mS_?NumericQ]:=(*Yi^2/Yieq^2*)(C0aH[mS/xphi]+\[CapitalGamma]htophiphi[mS/xphi,150]);*)


(* ::Input:: *)
(*ProductionC2[xx_?NumericQ,xphi_?NumericQ,mS_?NumericQ,mphi_?NumericQ]:=(*1/Yi*)((*-ProductionC0[xx,mS]*)+(xphi/(3mphi))(+C2aH[mS/xx]+\[CapitalGamma]htophiphi2[mS/xx,150]));*)


(* ::Input:: *)
(*AnnihilationC2[xx_?NumericQ,xphi_?NumericQ,mS_?NumericQ,mphi_?NumericQ]:=(*Yi/Yieq^2*)(- AnnihilationC0[xphi,mS]+xphi/(3mphi) (+C2aH[mS/xphi]+\[CapitalGamma]htophiphi2[mS/xphi,150]));*)


(* ::Section::Closed:: *)
(*Definition of annihilation and production functions for S*)


(* ::Input:: *)
(*ProductionC0S[xx_?NumericQ,mS_?NumericQ]:=\[CapitalGamma]hToSS[mS/xx];*)


(* ::Input:: *)
(*AnnihilationC0S[xphi_?NumericQ,mS_?NumericQ]:=(*Yi^2/Yieq^2*)(\[CapitalGamma]hToSS[mS/xphi]);*)


(* ::Input:: *)
(*ProductionC2S[xx_?NumericQ,xphi_?NumericQ,mS_?NumericQ,mphi_?NumericQ]:=(*1/Yi*)(-ProductionC0[xx,mS]+xphi/(3mphi) (\[CapitalGamma]hToSS2[mS/xx]));*)


(* ::Input:: *)
(*AnnihilationC2S[xx_?NumericQ,xphi_?NumericQ,mS_?NumericQ,mphi_?NumericQ]:=(*Yi/Yieq^2*)(- AnnihilationC0[xphi,mS]+xphi/(3mphi) (\[CapitalGamma]hToSS2[mS/xphi]));*)


(* ::Section::Closed:: *)
(*Scattering*)


(* ::Input:: *)
(*Cscatter[T1_?NumericQ,T2_?NumericQ,m1_?NumericQ,m2_?NumericQ]:=2 1/(128Pi^5) T1^2 T2^2 (T2-T1)E^(-(m1/T1)-m2/T2);*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Cscatterv2[T1_?NumericQ,T2_?NumericQ,m1_?NumericQ,m2_?NumericQ]:=2 1/(256Pi^5) m1 m2 T1 T2 (-m1 BesselK[1,m2/T2] BesselK[2,m1/T1]+m2 BesselK[1,m1/T1] BesselK[2,m2/T2])*)


(* ::Section::Closed:: *)
(*CMB SS->ee constraints*)


(* ::Input:: *)
(*\[Sigma]vNR[mS_?NumericQ,a_?NumericQ,mf_?NumericQ,AphiS_?NumericQ,\[Theta]_?NumericQ]:=*)
(*Module[{vh=246,GeVm2=(3.8240000000000003`*^-28)(*cm^2*),cons=29979245800 (*cm/s*)},(1-mf^2/mS^2)^(3/2)/(8 \[Pi] (16 (-1+a^2)^2 mS^4+4 a^2 mS^2 \[CapitalGamma]tot[2a mS,\[Theta]]^2))*AphiS^2*(mf/vh \[Theta])^2(*(AphiS^2 mf^2 (-mf^2+mS^2)^(3/2) \[Theta]^2)/(32 (-4+a^2)^2  mS^7 \[Pi] vh^2)*)GeVm2 (*cm^2*) cons(*cm/s*)]*)


(* ::Input:: *)
(*SetDirectory[NotebookDirectory[]];*)
(*dataCMBelectron= Import["files/electron_CMB_constraints.dat"];*)


(* ::Input:: *)
(*dataCMBelectronlinear = {};*)
(*Do[AppendTo[dataCMBelectronlinear, {10^ToExpression[StringDrop[j[[1]], -1]], 10^ToExpression[j[[2]]]}], {j,dataCMBelectron}];*)


(* ::Input:: *)
(*dataCMBelectronlinear = SortBy[dataCMBelectronlinear,First];*)


(* ::Input:: *)
(*dataCMBelectronlineari = Interpolation[dataCMBelectronlinear,InterpolationOrder->1];*)


(* ::Input:: *)
(*ratioE[mS_?NumericQ,a_?NumericQ,AphiS_?NumericQ,\[Theta]_?NumericQ]:=\[Sigma]vNR[mS,a,1/2000,AphiS,\[Theta]]/dataCMBelectronlineari[mS]*)


(* ::Input:: *)
(*\[Theta]maxff[mS_?NumericQ,a_?NumericQ,AphiS_?NumericQ]:=Block[{*)
(*min=-20,*)
(*max=4,*)
(*c,\[Epsilon]=0.0001,counter=0},*)
(*c=(min+max)/2;*)
(*While[Abs[ratioE[mS,a,AphiS,10^c]-1]>\[Epsilon],*)
(*If[(ratioE[mS,a,AphiS,10^c]-1)>0,*)
(*max=c,*)
(*min=c*)
(*];*)
(*c=(min+max)/2;*)
(*counter=counter+1;*)
(*If[counter>10000,*)
(*Break[]*)
(*]*)
(*];*)
(*Return[{10^c,ltime[2a mS,10^c]}//N]*)
(*];*)


(* ::Section::Closed:: *)
(*\[Sigma]^T/mS*)


(* ::Input:: *)
(*\[Sigma]TOverMs[mS_?NumericQ,a_?NumericQ,AphiS_?NumericQ,\[Theta]_?NumericQ,lS_?NumericQ,k_?NumericQ]:=Module[{v=10^-4},(-(1/(1024 mS^6 \[Pi] (1+v^2))) AphiS^4 (4/(v^2 (1-a^2+v^2))-(4 mS^2)/(4 mS^2 (1-a^2+v^2)^2+a^2 \[CapitalGamma]tot[2 a mS,\[Theta]]^2)+1/v^4 (Log[(a^2 (4 a^2 mS^2+\[CapitalGamma]tot[2 a mS,\[Theta]]^2))/(4 mS^2 (a^2+v^2)^2+a^2 \[CapitalGamma]tot[2 a mS,\[Theta]]^2)]+(4 a mS (ArcCot[(a \[CapitalGamma]tot[2 a mS,\[Theta]])/(2 a^2 mS+2 mS v^2)]-ArcTan[(2 a mS)/\[CapitalGamma]tot[2 a mS,\[Theta]]]))/\[CapitalGamma]tot[2 a mS,\[Theta]]+(a^2 Log[(a^2 (4 mS^2 (-1+a^2-v^2)-\[CapitalGamma]tot[2 a mS,\[Theta]]^2))/(4 mS^2 (-1+a^2-v^2) (a^2+v^2)-a^2 \[CapitalGamma]tot[2 a mS,\[Theta]]^2)] (4 mS^2 (1-a^2+v^2)+\[CapitalGamma]tot[2 a mS,\[Theta]]^2))/(mS^2 (1-a^2+v^2)^2)))+(lS^2 (4 v^4 (1+2 v^2) ((1+4 v^2) (3+4 v^2)^2+6 k (3+16 v^2+16 v^4)+9 k^2 (19+52 v^2+32 v^4))-3 k (3+16 v^2+16 v^4) (3+10 v^2+8 v^4+3 k (1+5 v^2+4 v^4)) Log[1/(1+4 v^2)]-3 k (3+16 v^2+16 v^4) (3+22 v^2+48 v^4+32 v^6+3 k (1+3 v^2+4 v^4)) Log[1+4 v^2]))/(512 mS^2 \[Pi] v^4 (1+v^2) (1+2 v^2) (1+4 v^2) (3+4 v^2)^2)) 0.00021839596186203026`(*1 GeV^-2=0.3894*10^-27 cm^2*) 1/mS(*11/GeV=1/(1.783*10^-24g)*)]*)


(* ::Section::Closed:: *)
(*Boltzmann equations*)


(* ::Input:: *)
(*Sol[mphi_?NumericQ]:=*)
(*Module[{},*)
(*sst[u_]:=ss[u,mphi];*)
(*averagedP4E3ef[];*)
(*];*)
(*FullsolND[mS_?NumericQ,mphi_?NumericQ,lhphi_?NumericQ,lphiS_?NumericQ,AphiS_?NumericQ,\[Theta]_?NumericQ,lphi_?NumericQ,k_?NumericQ,\[Xi]infDM_?NumericQ,\[Xi]infMediator_?NumericQ,xsmi_?NumericQ,xsmf_?NumericQ,Y0phi_?NumericQ,YS0_?NumericQ,f_,fc2_,wp_,ac_,pg_,ffi\[Theta]_,flagFermInt_,flag2_]:=*)
(*Module[{m1,m2,mh=125,v=246,lh=0.12909808976138543`},If[mphi>mS,m1=mphi;m2=mS,m1=mS;m2=mphi];*)
(*NDSolve[{Yphi'[xx]==1/(xx*Hbar[mS/xx]sst[xx]) (\[Theta]^2 \[CapitalGamma]boson0[mS/xx]+(-((mh^2 \[Theta]^2)/(2 v)) )^2 \[CapitalGamma]htophiphii[mS/xx,mphi,150]+(3/2 lh \[Theta]^2)^2 C0aH[mS/xx]+lh^2 \[Theta]^2 \[CapitalGamma]hhToh\[Phi]C0[mS/xx]+lphiS^2 \[Theta]^2 \[CapitalGamma]htophiSSC0[mS/xx]+f lphiS^2 (((YS[xx]sst[xx])/nphieq[mS/xS[xx],mS])^2 C0av2[mS/xS[xx],m1,m2]-((Yphi[xx]sst[xx])/nphieq[mphi/xphi[xx],mphi])^2 C0av2[mphi/xphi[xx],m1,m2])),xphi'[xx]==-xphi[xx]/(xx*Hbar[mS/xx]sst[xx]Yphi[xx]) ((xphi[xx]/(3mphi) ((3/2 lh \[Theta]^2)^2 C2aH[mS/xx]+(-((mh^2 \[Theta]^2)/(2 v)) \[Theta]^2)^2 \[CapitalGamma]htophiphi2[mS/xx,150]+lh^2 \[Theta]^2 \[CapitalGamma]hhToh\[Phi]C2[mS/xx]+lphiS^2 \[Theta]^2 \[CapitalGamma]htophiSSC2[mS/xx]+\[Theta]^2 \[CapitalGamma]boson2[mS/xx]))+((YS[xx]sst[xx])/nphieq[mS/xS[xx],mS])((Yphi[xx]sst[xx])/nphieq[mphi/xphi[xx],mphi])lphiS^2*xphi[xx]/(3mphi) (Cscatter[mphi/xphi[xx],mS/xS[xx],mphi,mS])+Yphi[xx]sst[xx] H[mS/xx]*xphi[xx]/(3*mphi) 1/(2Pi^2*nphieq[mphi/xphi[xx],mphi])*mphi^4 averagedP4E3e[xphi[xx]])-xphi[xx]((2sst'[xx])/(3sst[xx])-Yphi'[xx]/Yphi[xx]),*)
(*YS'[xx]==1/(xx*Hbar[mS/xx]sst[xx]) (ffi\[Theta] AphiS^2 \[Theta]^2 \[CapitalGamma]hToSS[mS/xx,mS]+lphiS^2 \[Theta]^2 \[CapitalGamma]htophiSSC0[mS/xx]+0.5lphi^3 averagedCif[xS[xx],mS,k]YS[xx]^2 sst[xx]^2 (-YS[xx]sst[xx]+nphieq[mS/xS[xx],mS])+f 0.5lphiS^2 (((Yphi[xx]sst[xx])/nphieq[mphi/xphi[xx],mphi])^2 C0av2[mphi/xphi[xx],m1,m2]-((YS[xx]sst[xx])/nphieq[mS/xS[xx],mS])^2 C0av2[mS/xS[xx],m1,m2])+flagFermInt AphiS^2 \[Theta]^2 (AnnAvergF[mS/xx]-((YS[xx]sst[xx])/nphieq[mS/xS[xx],mS])^2 AnnAvergF[mS/xS[xx]])+flag2 (Yphi[xx]sst[xx])/nphieq[mphi/xphi[xx],mphi] (lphiS \[Theta])^2 C0\[Phi]hToSS[mS/xx,mphi/xphi[xx],mphi]),*)
(*xS'[xx]==-xS[xx]/(xx*Hbar[mS/xx]sst[xx]YS[xx]) ((ffi\[Theta] xS[xx])/(3mS) AphiS^2 \[Theta]^2 \[CapitalGamma]hToSS2[mS/xx,mS]+xS[xx]/(3mS) lphiS^2 \[Theta]^2 \[CapitalGamma]htophiSSC2[mS/xx]+flag2 (lphiS \[Theta])^2 (Yphi[xx]sst[xx])/nphieq[mphi/xphi[xx],mphi] xS[xx]/(3mS) C2\[Phi]htoSSforS[mS/xx,mphi/xphi[xx],mphi]+((YS[xx]sst[xx])/nphieq[mS/xS[xx],mS])((Yphi[xx]sst[xx])/nphieq[mphi/xphi[xx],mphi])lphiS^2*xS[xx]/(3mS) (Cscatter[mS/xS[xx],mphi/xphi[xx],mS,mphi])+flagFermInt AphiS^2 \[Theta]^2 xS[xx]/(3mS) (AnnAvergF2ff[mS/xx]-((YS[xx]sst[xx])/nphieq[mS/xS[xx],mS])^2 AnnAvergF2ff[mS/xS[xx]])+sst[xx]^2*YS[xx]^2 (YS[xx]sst[xx]-nphieq[mS/xS[xx],mS])*0.5lphi^3 (C2cannibal[xS[xx],mS,k])+YS[xx]sst[xx] H[mS/xx]*xS[xx]/(3*mS)*1/(2Pi^2*nphieq[mS/xS[xx],mS])*mS^4 averagedP4E3e[xS[xx]])-xS[xx]((2sst'[xx])/(3sst[xx])-YS'[xx]/YS[xx])*)
(*,Yphi[xsmi]==Y0phi,xphi[xsmi]==(mphi xsmi)/(mS \[Xi]infMediator),YS[xsmi]==YS0,xS[xsmi]==xsmi/\[Xi]infDM},{xphi,Yphi,xS,YS},{xx,xsmi,xsmf},WorkingPrecision->wp,Method->{"DiscontinuityProcessing"->False},AccuracyGoal->ac,PrecisionGoal->pg*)
(*]];*)


(* ::Input:: *)
(*solND2[mS_?NumericQ,mphi_?NumericQ,lhphi_?NumericQ,lphiS_?NumericQ,AphiS_?NumericQ,\[Theta]_?NumericQ,lphi_?NumericQ,k_?NumericQ,\[Xi]infDM_?NumericQ,xsmi_?NumericQ,xsmf_?NumericQ,Y0phi_?NumericQ,YS0_?NumericQ,flag_?NumericQ,ffi\[Theta]_,flagFermInt_]:=*)
(*Module[{m1,m2,lh=0.12909808976138543`,mh=125,v=246},If[mphi>mS,m1=mphi;m2=mS,m1=mS;m2=mphi];*)
(*NDSolve[{Yphi'[xx]==1/(xx*Hbar[mS/xx]sst[xx]) (\[Theta]^2 \[CapitalGamma]boson0[mS/xx]+((*3lh \[Theta]^2*)-((mh^2 \[Theta]^2)/(2 v)))^2 \[CapitalGamma]htophiphii[mS/xx,mphi,150]+(3/2 lh \[Theta]^2)^2 C0aH[mS/xx]+lh^2 \[Theta]^2 \[CapitalGamma]hhToh\[Phi]C0[mS/xx]+lphiS^2 \[Theta]^2 \[CapitalGamma]htophiSSC0[mS/xx]+lphiS^2 (((YS[xx]sst[xx])/nphieq[mS/xS[xx],mS])^2 C0av2[mS/xS[xx],m1,m2]-flag ((Yphi[xx]sst[xx])/nphieq[mS/xS[xx],mphi])^2 C0av2[mS/xS[xx],m1,m2])),*)
(*YS'[xx]==1/(xx*Hbar[mS/xx]sst[xx]) (ffi\[Theta] AphiS^2 \[Theta]^2 \[CapitalGamma]hToSS[mS/xx,mS]+lphiS^2 \[Theta]^2 \[CapitalGamma]htophiSSC0[mS/xx]+0.5lphi^3 averagedCif[xS[xx],mS,k]YS[xx]^2 sst[xx]^2 (-YS[xx]sst[xx]+nphieq[mS/xS[xx],mS])+ lphiS^2 0.5(flag ((Yphi[xx]sst[xx])/nphieq[mS/xS[xx],mphi])^2 C0av2[mS/xS[xx],m1,m2]-((YS[xx]sst[xx])/nphieq[mS/xS[xx],mS])^2 C0av2[mS/xS[xx],m1,m2])+flagFermInt AphiS^2 \[Theta]^2 (AnnAvergF[mS/xx]-((YS[xx]sst[xx])/nphieq[mS/xS[xx],mS])^2 AnnAvergF[mS/xS[xx]])),*)
(*xS'[xx]==-xS[xx]/(xx*Hbar[mS/xx]sst[xx]YS[xx]) (ffi\[Theta] xS[xx]/(3mS) AphiS^2 \[Theta]^2 \[CapitalGamma]hToSS2[mS/xx,mS]+lphiS^2 \[Theta]^2 xS[xx]/(3mS) \[CapitalGamma]htophiSSC2[mS/xx]+sst[xx]^2*YS[xx]^2 (YS[xx]sst[xx]-nphieq[mS/xS[xx],mS])*0.5lphi^3 (C2cannibal[xS[xx],mS,k])+flagFermInt AphiS^2 \[Theta]^2 xS[xx]/(3mS) (AnnAvergF2ff[mS/xx]-((YS[xx]sst[xx])/nphieq[mS/xS[xx],mS])^2 AnnAvergF2ff[mS/xS[xx]])+YS[xx]sst[xx] H[mS/xx]*xS[xx]/(3*mS)*1/(2Pi^2*nphieq[mS/xS[xx],mS])*mS^4 averagedP4E3e[xS[xx]])-xS[xx]((2sst'[xx])/(3sst[xx])-YS'[xx]/YS[xx])*)
(*,Yphi[xsmi]==(*nphieq[mphi/((mphi xsmi)/(mS \[Xi]infDM)),mphi]/sst[xsmi]*)Y0phi,YS[xsmi]==(*nphieq[mS/(xsmi/\[Xi]infDM),mS]/sst[xsmi]*)YS0,xS[xsmi]==xsmi/\[Xi]infDM},{Yphi,xS,YS},{xx,xsmi,xsmf},WorkingPrecision->MachinePrecision(*,Method->{"DiscontinuityProcessing"->False}*)(*,AccuracyGoal->13,PrecisionGoal->13*)*)
(*]];*)


(* ::Input:: *)
(*FinalS1[mS_?NumericQ,mphi_?NumericQ,lhphi_?NumericQ,lphiS_?NumericQ,AphiS_?NumericQ,\[Theta]_?NumericQ,lS_?NumericQ,k_?NumericQ,\[Xi]infDM_?NumericQ,\[Xi]infMediator_?NumericQ,f_,fc2_,xf_,ffi\[Theta]_,flagFermInt_,flag2_]:=*)
(*Module[{x0,x02,RelDen,RelObs,ratRelics,mann,xp,xDM,Ymed,YDM,xff,psecond,pfirst,\[Xi]f,MedRelDen,x\[Epsilon],Tphip,Tp,\[Alpha]p,\[Xi]inf2,\[Xi]inf2med,Yphi02,YS02,Yphi0,YS0,FullsolNDD,solNDD},*)
(*Sol[mS];*)
(*(*If[flagFermInt==0,*)*)
(*FullsolNDD=FullsolND(*,*)
(*FullsolNDD=FullsolND*)
(*]*);*)
(*x0=mS/150;*)
(*Yphi0=(*10^-19*)nphieq[mphi/((mphi x0)/(mS \[Xi]infMediator)),mphi]/sst[x0];*)
(*YS0=(*10^-15*)nphieq[mS/(x0/\[Xi]infDM),mS]/sst[x0];*)
(*pfirst=FullsolNDD[mS,mphi,lhphi,lphiS,AphiS,\[Theta],lS,k,\[Xi]infDM,\[Xi]infMediator,x0,xf,Yphi0,YS0,f,fc2,MachinePrecision,Automatic,Automatic,ffi\[Theta],flagFermInt,flag2];*)
(*If[pfirst[[1]][[1]][[2]]["Domain"][[1]][[2]]>0.1,*)
(*x02=pfirst[[1]][[1]][[2]]["Domain"][[1]][[2]]-0.05,*)
(*x02=pfirst[[1]][[1]][[2]]["Domain"][[1]][[2]]*)
(*];*)
(*\[Xi]inf2=x02/(xS[x02]/.pfirst[[1]]);*)
(*\[Xi]inf2med=mphi/mS x02/(xphi[x02]/.pfirst[[1]]);*)
(*Yphi02=Yphi[x02]/.pfirst[[1]];*)
(*YS02=YS[x02]/.pfirst[[1]];*)
(*Tphip=mphi/xphi[x02]/.pfirst[[1]];*)
(*Tp = mS/x02;*)
(*\[Alpha]p=Tphip/Tp^2;*)
(*psecond=FullsolNDD[mS,mphi,0.,lphiS,AphiS,\[Theta],lS,k,\[Xi]inf2,\[Xi]inf2med,x02,180.,Yphi02,YS02,f,fc2,MachinePrecision,Automatic,Automatic,ffi\[Theta],flagFermInt,flag2];*)
(*xp[x_]:=Piecewise[{{xphi[x]/.pfirst[[1]],x<x02},{xphi[x]/.psecond[[1]],x>=x02}}];*)
(*Ymed[x_]:=Piecewise[{{Yphi[x]/.pfirst[[1]],x<x02},{Yphi[x]/.psecond[[1]],x>=x02}}];*)
(*xDM[x_]:=Piecewise[{{xS[x]/.pfirst[[1]],x<x02},{xS[x]/.psecond[[1]],x>=x02}}];*)
(*YDM[x_]:=Piecewise[{{YS[x]/.pfirst[[1]],x<x02},{YS[x]/.psecond[[1]],x>=x02}}];*)
(*xff=psecond[[1]][[1]][[2]]["Domain"][[1]][[2]];*)
(*x\[Epsilon]=10^-7;*)
(*\[Xi]f=mphi/mS (xff-x\[Epsilon])/xp[xff-x\[Epsilon]];*)
(*RelDen=2YDM[xff-x\[Epsilon]];*)
(*MedRelDen=Ymed[xff-x\[Epsilon]];*)
(*RelObs=4.371791351001918`*^-10/mS;*)
(*ratRelics=RelDen/RelObs;*)
(*If[pfirst===$Failed||psecond===$Failed,Return[{-1}],*)
(*Return[{xp,Ymed,xDM,YDM,xff-x\[Epsilon],\[Theta]//N,RelObs,RelDen,MedRelDen,\[Xi]f,ratRelics}]];*)
(*];*)


(* ::Input:: *)
(*FinalS2[mS_?NumericQ,mphi_?NumericQ,lhphi_?NumericQ,lphiS_?NumericQ,AphiS_?NumericQ,\[Theta]_?NumericQ,lphi_?NumericQ,k_?NumericQ,\[Xi]infDM_?NumericQ,\[Xi]infMediator_?NumericQ,f_,fc2_,xf_,ffi\[Theta]_,flagFermInt_,flag2_]:=*)
(*Module[{x0,x02,RelDen,RelObs,ratRelics,mann,xp,xDM,Ymed,YDM,psecond,pfirst,\[Xi]f,MedRelDen,x\[Epsilon],Tphip,Tp,\[Alpha]p,\[Xi]inf2,\[Xi]inf2med,Yphi02,YS02,Yphi0,YS0,xff,FullsolNDD,solNDD},*)
(*Sol[mS];*)
(*(*If[flagFermInt==0,*)*)
(*FullsolNDD=FullsolND;*)
(*solNDD=solND2(*,*)
(*FullsolNDD=FullsolND;*)
(*solNDD=solND*)
(*]*);*)
(*x0=mS/150;*)
(*Yphi0=(*10^-15*)nphieq[mphi/((mphi x0)/(mS \[Xi]infMediator)),mphi]/sst[x0];*)
(*YS0=(*10^-15*)nphieq[mS/(x0/\[Xi]infDM),mS]/sst[x0];*)
(*pfirst=FullsolNDD[mS,mphi,lhphi,lphiS,AphiS,\[Theta],lphi,k,\[Xi]infDM,\[Xi]infMediator,x0,xf,Yphi0,YS0,f,fc2,MachinePrecision,Automatic,Automatic,ffi\[Theta],flagFermInt,flag2];*)
(*If[pfirst[[1]][[1]][[2]]["Domain"][[1]][[2]]>0.1,*)
(*x02=pfirst[[1]][[1]][[2]]["Domain"][[1]][[2]]-0.05,*)
(*x02=pfirst[[1]][[1]][[2]]["Domain"][[1]][[2]]*)
(*];*)
(*\[Xi]inf2=x02/(xS[x02]/.pfirst[[1]]);*)
(*\[Xi]inf2med=mphi/mS x02/(xphi[x02]/.pfirst[[1]]);*)
(*Yphi02=Yphi[x02]/.pfirst[[1]];*)
(*YS02=YS[x02]/.pfirst[[1]];*)
(*Tphip=mphi/xphi[x02]/.pfirst[[1]];*)
(*Tp = mS/x02;*)
(*\[Alpha]p=Tphip/Tp^2;*)
(*psecond=solNDD[mS,mphi,0.,lphiS,AphiS,\[Theta],lphi,k,\[Xi]inf2,x02,180.,Yphi02,YS02,0.,ffi\[Theta],flagFermInt];*)
(*xp[x_]:=Piecewise[{{xphi[x]/.pfirst[[1]],x<x02},{mphi/(\[Alpha]p (mS/x)^2),x>=x02}}];*)
(*Ymed[x_]:=Piecewise[{{Yphi[x]/.pfirst[[1]],x<x02},{Yphi[x]/.psecond[[1]],x>=x02}}];*)
(*xDM[x_]:=Piecewise[{{xS[x]/.pfirst[[1]],x<x02},{xS[x]/.psecond[[1]],x>=x02}}];*)
(*YDM[x_]:=Piecewise[{{YS[x]/.pfirst[[1]],x<x02},{YS[x]/.psecond[[1]],x>=x02}}];*)
(*xff=psecond[[1]][[1]][[2]]["Domain"][[1]][[2]];*)
(*\[Xi]f=mphi/mS x02/xp[x02];*)
(*MedRelDen=Ymed[x02];*)
(*x\[Epsilon]=10^-7;*)
(*RelDen=2YDM[xff-x\[Epsilon]];*)
(*RelObs=4.371791351001918`*^-10/mS;*)
(*ratRelics=RelDen/RelObs;*)
(*If[pfirst===$Failed||psecond===$Failed,Return[{-1}],*)
(*Return[{xp,Ymed,xDM,YDM,xff-x\[Epsilon],\[Theta]//N,RelObs,RelDen,MedRelDen,\[Xi]f,ratRelics}]];*)
(*];*)


(* ::Input:: *)
(*FinalS[mS_?NumericQ,a_?NumericQ,lphiS_?NumericQ,b_?NumericQ,\[Theta]_?NumericQ,lS_?NumericQ,k_?NumericQ,\[Xi]infDM_?NumericQ]:=Module[{},*)
(*If[mS>2a mS||lphiS<=10^-5,*)
(*Return[FinalS2[mS,2a mS,0.,lphiS,2a mS b,\[Theta],lS,k,\[Xi]infDM,\[Xi]infDM,1.,1.,100.,1.,0.,0.]],*)
(*Return[FinalS1[mS,2a mS,0.,lphiS,2a mS b,\[Theta],lS,k,\[Xi]infDM,\[Xi]infDM,1.,1.,100.,1.,0.,0.]]*)
(*]*)
(*]*)


(* ::Section::Closed:: *)
(*Scan routine*)


(* ::Input:: *)
(*Fit\[CapitalOmega]obs[mS_?NumericQ,a_?NumericQ,lphiS_?NumericQ,b_?NumericQ,lS_?NumericQ,k_?NumericQ,\[Xi]infDM_?NumericQ,min_,max_,flag_]:=*)
(*Module[{\[Epsilon]=0.01,counter,r,minn,maxx,th*)
(*},*)
(*counter=1;*)
(*minn=min;*)
(*maxx=max;*)
(*th=(minn+maxx)/2;*)
(*r=FinalS[mS,a,lphiS,b,10^th,lS,k,\[Xi]infDM];*)
(*(*Print["first log10(\[Theta]) is: "<>ToString[th//N]<>" with \[CapitalOmega] ratio: "<>ToString[r[[-1]]//N]];*)*)
(*While[Abs[r[[-1]]-1]>\[Epsilon],*)
(*If[r[[-1]]-1<0,*)
(*minn=th,*)
(*maxx=th*)
(*];*)
(*th=(minn+maxx)/2;*)
(*counter=counter+1;*)
(*If[counter>15,*)
(*(*Print["max iterations achieved"];*)*)
(*Return[{-10}]*)
(*];*)
(*r=FinalS[mS,a,lphiS,b,10^th,lS,k,\[Xi]infDM];*)
(*(*Print["updated log10(\[Theta]) is: "<>ToString[th//N]<>" with \[CapitalOmega] ratio: "<>ToString[r[[-1]]//N]];*)*)
(*];*)
(*Return[{10^th//N,ltime[2a mS,10^th],r[[-2]](*final ratio of temperatures*),r[[-3]] (*final mediator abundance*),r[[5]](*final xf*)}];*)
(*Clear["`*"]*)
(*];*)


(* ::Input:: *)
(*fit[mS_?NumericQ,a_?NumericQ,lphiS_?NumericQ,b_?NumericQ,lS_?NumericQ,k_?NumericQ,\[Xi]infDM_?NumericQ,min_,max_]:=*)
(*Module[{temp,err=1.},*)
(*temp=Fit\[CapitalOmega]obs[mS,a,lphiS,b,lS,k,\[Xi]infDM,min,max,0.];*)
(*If[Length[temp]==5,*)
(*(*AnnAvergF=TabulatingC0[mS,2 a mS,temp[[1]]];*)
(*AnnAvergF2ff=TabulatingC2[mS,2a mS,temp[[1]]];*)
(*Return[temp(*(*Fit\[CapitalOmega]obs[mS,a,lphiS,b,lS,k,\[Xi]infDM,Log10[temp[[1]]]-err,Log10[temp[[1]]]+err,0.]*)*)];*)*)
(*Return[temp],*)
(*Return[{-10}]*)
(*]*)
(*]*)


(* ::Section::Closed:: *)
(*Scan*)


(* ::Input:: *)
(*Off[NIntegrate::izero];*)
(*Off[General::munfl];*)
(*Off[NIntegrate::ncvb];*)
(*Off[CompiledFunction::cfsa];*)
(*Off[NDSolve::nlnum];*)
(*Scan[\[Xi]inf_,mSmin_,mSmax_,bmin_,bmax_,lphismin_,lphismax_,k_,it_]:=*)
(*Module[{fitt,mS,a,b,lphis,\[Theta]temp,fitt2,lslist,dat,\[Theta]min=-13,\[Theta]max=-8.},*)
(*dat=Import["completescan3(mS,a,b,th,lphis,ls,ltime,rat,MedAbun,xf)/k="<>ToString[k//N]<>",v4.dat"];*)
(*lslist={-4,-3.5,-3,-2.5,-2,-1.5,-1};*)
(*Do[*)
(*mS=10^RandomReal[{Log10[mSmin],Log10[mSmax]}];*)
(*a=RandomReal[{0.,1}];*)
(*b=10^RandomReal[{Log10[bmin],Log10[bmax]}];*)
(*lphis=10^RandomReal[{Log10[lphismin],Log10[lphismax]}];*)
(*If[\[Sigma]TOverMs[mS,a,2 a b mS,10^-5,0.,k]>1,Continue[]];*)
(*Print["mS = "<>ToString[mS//N]<>", a = "<>ToString[a//N]<>", log10(b) = "<>ToString[Log10[b]//N]<>", log10(lphiS)="<>ToString[Log10[lphis]]<>", fitting \[Theta] to the correct relic..."];fitt=TimeConstrained[fit[mS,a,lphis,b,0.,k,\[Xi]inf,\[Theta]min,\[Theta]max],60*12,{-10}]//Quiet;*)
(*If[Length[fitt]==5,PrintTemporary["\[Theta] = "<>ToString[fitt[[1]]//N]<>"with ltime = "<>ToString[fitt[[2]]]<>" s, saving the results..."];*)
(*\[Theta]temp=fitt[[1]];*)
(*AppendTo[dat,{mS 1000,a,b,fitt[[1]] (*=\[Theta]*),lphis,0.,fitt[[2]](*=lifetime*),fitt[[3]](*ratio of med temp at dec.*),fitt[[4]](*final med abundance*),fitt[[5]] (*final xf*)}];*)
(*Export["completescan3(mS,a,b,th,lphis,ls,ltime,rat,MedAbun,xf)/k="<>ToString[k//N]<>",v4.dat",dat];*)
(*Do[Print["Including self interactions..."];*)
(*If[\[Sigma]TOverMs[mS,a,2 a b mS,\[Theta]temp,10^lS,k]>1,Continue[]];*)
(*fitt2=TimeConstrained[fit[mS,a,lphis,b,10^lS//N,k,\[Xi]inf,Log10[\[Theta]temp]-2,Log10[\[Theta]temp]+2],60*12,{-10}]//Quiet;*)
(*If[Length[fitt2]==1,Continue[]];*)
(*\[Theta]temp=fitt2[[1]];*)
(*AppendTo[dat,{mS 1000,a,b,fitt2[[1]] (*=\[Theta]*),lphis,10^lS//N,fitt2[[2]](*=lifetime*),fitt2[[3]](*ratio of med temp at dec.*),fitt2[[4]](*final med abundance*),fitt2[[5]] (*final xf*)}];*)
(*Export["completescan3(mS,a,b,th,lphis,ls,ltime,rat,MedAbun,xf)/k="<>ToString[k//N]<>",v4.dat",dat]*)
(**)
(*,{lS,lslist}]*)
(*,*)
(*PrintTemporary["max iterations achieved, carrying on wih the next point..."];*)
(*Continue[]];,{it}];*)
(*Return[dat]];*)


(* ::Input:: *)
(*\[Xi]inf=0.001;*)
(*mSmin=0.0005;*)
(*mSmax=2.;*)
(*bmin=10^-5;*)
(*bmax=1.;*)
(*lphiSmin=10^-8;*)
(*lphiSmax=10^-1;*)
(*it=10000;*)
(*k=0.5;*)


(* ::Input:: *)
(*Scan[\[Xi]inf,mSmin,mSmax,bmin,bmax,lphiSmin,lphiSmax,k,it]*)
