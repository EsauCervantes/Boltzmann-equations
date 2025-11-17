(* ::Package:: *)

(* ::Text:: *)
(*Tabulation of thermal averages*)


SetDirectory[NotebookDirectory[]];
InterThermalAverages[]:=
Module[{dataC0=Import["files/dataC0.dat"],
dataC1=Import["files/dataC1V2.dat"],
dataC2=Import["files/dataC2V2.dat"],
m0,b0,m1,b1,m2,b2,facC0=1/(2(2Pi)^4),fac1=1/(8(2Pi)^4) 2,factor=(81 Sqrt[3])/(16384 \[Pi]^8),
dataC00={},
dataC01={},
dataC02={}
},
Do[AppendTo[dataC00,{El[[1]],El[[2]](1/(2!*3!))facC0(*(El[[1]]*2Pi^2)^3/(El[[1]](BesselK[2,El[[1]]])^3)*)(*mphi(9mphi^2)/mphi^9(0.1/mphi)^23mphi^3*)}],{El,dataC0}];
Do[AppendTo[dataC01,{El[[1]],El[[2]](El[[1]]/3 *1/3!)fac1(*(El[[1]]*2Pi^2)^3/(BesselK[2,El[[1]]])^3*)(*1/mphi(9mphi^2)/mphi^9(0.1/mphi)^2(3mphi^2)(3mphi^3)*)}],{El,dataC1}];
Do[AppendTo[dataC02,{El[[1]],El[[2]](El[[1]]/3*1/2! 1/2!)factor(*(El[[1]]*2Pi^2)^3/( BesselK[2,El[[1]]])^3*)}],{El,dataC2}];
averagedCi=Interpolation[dataC00/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->2];
averagedC1i=Interpolation[dataC01/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->2];(*thermal average with p^2/E outside the cross section*)
averagedC2i=Interpolation[dataC02/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->2];
m0=(averagedCi[dataC0[[1]][[1]]]-averagedCi[dataC0[[2]][[1]]])/(Log10[dataC0[[1]][[1]]]-Log10[dataC0[[2]][[1]]]);
b0=10^(averagedCi[dataC0[[1]][[1]]]-m0 Log10[dataC0[[1]][[1]]]);
m2=(averagedC2i[dataC2[[1]][[1]]]-averagedC2i[dataC2[[2]][[1]]])/(Log10[dataC2[[1]][[1]]]-Log10[dataC2[[2]][[1]]]);
b2=10^(averagedC2i[dataC2[[1]][[1]]]-m2 Log10[dataC2[[1]][[1]]]);
m1=(averagedC1i[dataC1[[1]][[1]]]-averagedC1i[dataC1[[2]][[1]]])/(Log10[dataC1[[1]][[1]]]-Log10[dataC1[[2]][[1]]]);
b1=10^(averagedC1i[dataC1[[1]][[1]]]-m1 Log10[dataC1[[1]][[1]]]);
averagedCif[xphi_?NumericQ,mphi_?NumericQ]:=0.27/mphi^5 Piecewise[{{b0 xphi^(m0-1) (xphi*2Pi^2)^3/( BesselK[2,xphi])^3,xphi<10^-6},{(xphi*2Pi^2)^3/( BesselK[2,xphi])^3 10^averagedCi[xphi]/xphi ,xphi>=10^-6}}];
averagedC1if[xphi_?NumericQ,mphi_?NumericQ]:=0.81/mphi^5 Piecewise[{{(xphi*2Pi^2)^3/( BesselK[2,xphi])^3 b1 xphi^m1,xphi<10^-2},{(xphi*2Pi^2)^3/( BesselK[2,xphi])^3 10^averagedC1i[xphi] ,xphi>=10^-2}}];
averagedC2if[xphi_?NumericQ,mphi_?NumericQ]:=1/mphi^5 Piecewise[{{(xphi*2Pi^2)^3/( BesselK[2,xphi])^3 b2 xphi^m2 ,xphi<10^-2},{(xphi*2Pi^2)^3/( BesselK[2,xphi])^3 10^averagedC2i[xphi],10^-2<=xphi}}];
C2cannibal[xphi_,mphi_]:=HeavisideTheta[xphi-0.5](averagedC1if[xphi,mphi]-averagedC2if[xphi,mphi])];



(* ::Text:: *)
(*Definition of H rate, entropy in equilibrium*)


datag = Import["files/gEnergy.dat"];
datagentropy = Import["files/gEntropy.dat"];
For[i =0,i<=Length[datag]-1,i++;datag[[i]][[1]]=Log10[datag[[i]][[1]]]];
gg=Interpolation[datag,Method->"Spline"];
grel[T_]:=  gg[Log10[T]];
(*grelx[x_,m_]:=grel[m/x];*)
(*For[i =0,i<=Length[datagentropy]-1,i++;datagentropy[[i]][[1]]=Log10[datagentropy[[i]][[1]]]];*)
gentropy=Interpolation[datagentropy/.{x_,z_}->{Log10[x],z},Method->"Spline"];
gent[xx_,mphi_]:=gentropy[Log10[mphi/xx]];
gentTemp[T_]:=gentropy[Log10[T]];
Clear[s];
s[T_]:=Piecewise[{{(2*Pi^2 gentropy[Log10[T]]*T^3)/45,T<1000},{(2*Pi^2 gentropy[Log10[1000]]*T^3)/45,T>=1000}}];
Clear[ss];
ss[x_,mphi_]:=s[mphi/x];(*Piecewise[{{(2*(Pi^2)*gent[x,mphi](*106.69594335689291*))/45*mphi^3/x^3,x>=mphi/1000},{(2*Pi^2*gent[0.1/1000,mphi](mphi/x)^3)/45,x<mphi/1000}}];*)
\[Rho]eqq[T_?NumericQ]:=grel[T] Pi^2/30 T^4;
H[T_?NumericQ]:=Module[{Mpl=2.4*10^18},Sqrt[\[Rho]eqq[T]/(3Mpl^2)]];(*(1.66*Sqrt[grel[T]](*10*))/(1.22089*10^19)T^2;*)
Hbar[T_?NumericQ]:=H[T]/(1+T/(3gentTemp[T])*gentTemp'[T]);
ssc[x_,mphi_]:=(*Piecewise[{{*)(2*Pi^2*(*gent[x,mphi]*)106.69594335689291(mphi/x)^3)/45;(*,x>=mphi/1000},{(2*Pi^2*gent[0.1/1000,mphi](mphi/x)^3)/45,x<mphi/1000}}];*)
nphieq=Compile[{{x,_Real,0},{m,_Real,0}},
1/(x*2Pi^2)*m^3*BesselK[2,x],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True];
(*averagedP4E3f[xx_?NumericQ]:=(*1/(2Pi^2*nphieq[xx,mphi])*mphi^4*)NIntegrate[(1/Et^2-1)^(5/2)Exp[-xx/Et],{Et,0,10^-9,10^-8,10^-7,10^-6,10^-5,10^-4,10^-3,10^-2,10^-1,0.5,0.8,0.9,1},AccuracyGoal->35,PrecisionGoal->35,WorkingPrecision->MachinePrecision]
temp=Table[{10^xdm//N,averagedP4E3f[10^xdm]},{xdm,Log10[10^-4],Log10[1000],(Log10[1000]-Log10[10^-4])/1400}];
Export["files/intP4_over_E3.dat",temp];*)
averagedP4E3list=Import["files/intP4_over_E3.dat"];
averagedP4E3inter=Interpolation[averagedP4E3list,Method->"Spline"];
\[Rho]eq[T_,m_]:=1/(2Pi^2) m^3 T(BesselK[1,m/T]+3 T/m BesselK[2,m/T]);


D\[Rho]eq[T_,m_,mprime_]:=(m ((BesselK[0,m/T] m (12 T^2+m^2-T m mprime))/T+BesselK[1,m/T] (24 T^2+5 m^2-T m mprime)))/(2 \[Pi]^2);


(*to calculate the relic density, remember that 0.12 = \[CapitalOmega]h^2= (mphi s Yobs)/\[Rho]c, Hence Yobs = (9.711808033846503`*^-48 gev^4)/(mphi(2*43*Pi^2/(11*45)((8.617333262145`*^-14)2.72548 gev)^3))*)


(* ::Text:: *)
(*Implementing background eqs + cBE*)


Neq[T_?NumericQ,m_?NumericQ,a_?NumericQ]:=a^3 (m^2 T BesselK[2,m/T])/(2 \[Pi]^2);(*number in equilibrium.*)


Nobs[ms_]:=2 1.0153045443906683`*^25/ms;


EnergyDensities[]:=Module[{Mpl = 2.4*10^18(*GeV*),ai = 1,
af=10^17,arh,Trh=150,Hsol,TsoldatStandard,TStandardI,afinal,Tnonst,domain},
TsoldatStandard = Table[{ai (s[Trh]/s[10^T])^(1/3)(*arh 5/4 Sqrt[H[TNonStandardI[arh] ]/H[10^T]]*),10^T},{T,Log10[0.0000001],Log10[1.Trh],(Log10[1.Trh]-Log10[0.0000001])/100000}];
(*ainit=FindRoot[TsoldatStandard];*)
TStandardI=Interpolation[TsoldatStandard,InterpolationOrder->2,Method->"Spline"];
(*aEWPTStandard=FindRoot[TStandardI[a]==150,{a,10^7},AccuracyGoal->20,PrecisionGoal->20,WorkingPrecision->30][[1]][[2]];*)
domain= TStandardI(*[[1]][[1]][[2]]*)["Domain"](*[[1]][[2]]*);
afinal = af;
Hsol[a_]:=Sqrt[\[Rho]eqq[TStandardI[a]]/(3Mpl^2)](*Total Hubble rate*);
Return[{Hsol,TStandardI,ai,afinal,domain(*,TsoldatStandard*)}]
];


en=EnergyDensities[];


Tstandard=en[[2]];


af=en[[-1]][[-1]][[-1]]


Definition of h -> \[Phi]\[Phi] decay routines


Clear[\[CapitalGamma]htophiphii]
\[CapitalGamma]htophiphii=Compile[{{T,_Real,0},{mphi,_Real,0}},
Module[{v=246,mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2],Tc=150},
(*Piecewise[{{*)(v^2 mh)/(16Pi^3) Sqrt[1-(4mphi^2)/mh^2]T BesselK[1,mh/T](*,T<=Tc},{0,T>Tc}}]*)
](*,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True*)
];
Clear[\[CapitalGamma]htophiphii2]
\[CapitalGamma]htophiphii2=Compile[{{T,_Real,0},{mphi,_Real,0}},
Module[{v=246,mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2],Tc=150},
(*Piecewise[{{*)HeavisideTheta[1-(4mphi^2)/mh^2] (v^2 mh^2)/(32Pi^3) Sqrt[1-(4mphi^2)/mh^2]  (T BesselK[2,mh/T])(*,T<=Tc},{0,T>Tc}}]*)
](*,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True*)
];
\[CapitalGamma]htophiphi[T_?NumericQ,mphi_?NumericQ]:=\[CapitalGamma]htophiphii[T,mphi];
\[CapitalGamma]htophiphi2[T_?NumericQ,mphi_?NumericQ]:=\[CapitalGamma]htophiphii2[T,mphi];


\[CapitalGamma]htophiphi2Energydens=Compile[{{T,_Real,0},{mphi,_Real,0}},
Module[{v=246,(*mh=125*)mh=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2],Tc=150},
HeavisideTheta[1-(4mphi^2)/mh^2] (v^2 mh^2)/(32Pi^3) Sqrt[1-(4mphi^2)/mh^2]  (T BesselK[2,mh/T])
]
,RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


\[CapitalGamma]htophiphii2Energydens[T_?NumericQ,mphi_?NumericQ]:=\[CapitalGamma]htophiphi2Energydens[T,mphi];


(* ::Text:: *)
(*Definition of hh -> \[Phi]\[Phi] routine*)


Tabulation of annihilation averages


C0a[T_]:=Module[{(*MH=125*)MH=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},(MH^2 T^2 BesselK[1,MH/T]^2)/(256 \[Pi]^5)];


C2a[T_]:= Module[{(*MH=125*)MH=Sqrt[125^2+(1/2 0.14+1/16 0.65^2+3/16 0.35^2+1/4 (176/246)^2) T^2]},5 (E^(-((2 MH)/T)) MH^2 T^3)/(512 \[Pi]^4)];


Definition of \[Phi]->\[Chi]\[Chi]


Clear[part2]
part2[]:=Block[{int,dat,part22i},
(*int[x_?NumericQ]:=NIntegrate[Sqrt[1-1/Ent^2]Exp[-Ent x],{Ent,1,1.00000001,1.0001,1.001,1.01,1.03,1.05,1.1,1.5,1.8,2,5,10,100,1000,10^4,10^5,10^6,10^7,10^8,10^9,Infinity},WorkingPrecision->50,PrecisionGoal->30,AccuracyGoal->30];*)
dat=Import["files/phi_to_2chi.dat"];(*Table[{10^x//N,int[10^x]},{x,-7,7,14/4000}];*)
(*Export["files/phi_to_2chi.dat",dat];*)
part22i=Interpolation[dat/.{x_,z_}->{x,Log10[z]},Method->"Spline",InterpolationOrder->2];
Return[{dat,part22i}]
];


p2=part2[];


Clear[\[CapitalGamma]\[Phi]to2chi]
\[CapitalGamma]\[Phi]to2chi=Compile[{{T,_Real,0},{mphi,_Real,0},{mchi,_Real,0}},
Module[{Tc=150,coupl=2  (mphi^2-4mchi^2)},
(*Piecewise[{{*)coupl mphi/(16Pi^3) Sqrt[1-(4mchi^2)/mphi^2]HeavisideTheta[1-(4mchi^2)/mphi^2]T BesselK[1,mphi/T](*,T<=Tc},{0,T>Tc}}]*)
],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


Clear[\[CapitalGamma]\[Phi]to2chi2];
\[CapitalGamma]\[Phi]to2chi2=Compile[{{T,_Real,0},{mphi,_Real,0},{mchi,_Real,0}},
Module[{Tc=150,coupl=2  (mphi^2-4mchi^2)},
(*Piecewise[{{*)coupl 1/(32Pi^3) Sqrt[1-(4mchi^2)/mphi^2]HeavisideTheta[1-(4mchi^2)/mphi^2]mphi^2 (T BesselK[2,mphi/T]-E^(-(mphi/T)) Sqrt[1/mphi] Sqrt[\[Pi]/2] T^(3/2)(*-mphi 10^p2[[2]][mphi/T]*))(*,T<=Tc},{0,T>Tc}}]*)
],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


Clear[\[CapitalGamma]\[Phi]to2chiEnergy];
\[CapitalGamma]\[Phi]to2chiEnergy=Compile[{{T,_Real,0},{mphi,_Real,0},{mchi,_Real,0}},
Module[{Tc=150,coupl=2  (mphi^2-4mchi^2)},
coupl 1/(32Pi^3) Sqrt[1-(4mchi^2)/mphi^2]HeavisideTheta[1-(4mchi^2)/mphi^2]mphi^2 (T BesselK[2,mphi/T])
],RuntimeOptions->"EvaluateSymbolically"->False,RuntimeOptions->"Speed",CompilationTarget->"C",Parallelization->True
];


\[CapitalGamma]\[Phi]to2\[Chi][T_?NumericQ,mphi_?NumericQ,mchi_?NumericQ]:=\[CapitalGamma]\[Phi]to2chi[T,mphi,mchi];


\[CapitalGamma]\[Phi]to2\[Chi]C2[T_?NumericQ,mphi_?NumericQ,mchi_?NumericQ]:=\[CapitalGamma]\[Phi]to2chi2[T,mphi,mchi];


\[CapitalGamma]\[Phi]to2chiEner[T_?NumericQ,mphi_?NumericQ,mchi_?NumericQ]:=\[CapitalGamma]\[Phi]to2chiEnergy[T,mphi,mchi];


\[CapitalGamma]\[Phi]to2\[Chi][10,mphi,mchi];


(* ::Text:: *)
(*Solver Routines*)


InterThermalAverages[];
nphieqq[x_?NumericQ,mphi_?NumericQ]:=mphi^3 1/(x*2Pi^2)*BesselK[2,x];
averagedP4E3ef[mphi_?NumericQ]:=Module[{m1,b1,m2,b2},
m1=(Log10[averagedP4E3inter[0.00015]]-Log10[averagedP4E3inter[0.00013]])/(Log10[0.00015]-Log10[0.00013]);
b1=Log10[averagedP4E3inter[0.00013]]-m1 Log10[0.00013];
m2=(Log10[averagedP4E3inter[99.9]]-Log10[averagedP4E3inter[100]])/(Log10[99.9]-Log10[100]);
b2=Log10[averagedP4E3inter[100]]-m2 Log10[100];
averagedP4E3e[xphi_?NumericQ]:=
Piecewise[
{
{1/(2Pi^2*nphieqq[xphi,mphi])*mphi^4 averagedP4E3inter[xphi],0.00013<xphi<=200},
{1/(2Pi^2*nphieqq[xphi,mphi])*mphi^4 Abs[10^b1 xphi^m1],xphi<=0.00013}
}
]
];
Sol[mphi_?NumericQ]:=
Module[{},
sst[u_]:=ss[u,mphi];
averagedP4E3ef[mphi]
]
Clear[solND]
solND[mphi_?NumericQ,\[Lambda]hphi_?NumericQ,lphi_?NumericQ,\[Xi]inf_?NumericQ,mchi_?NumericQ,\[Kappa]_?NumericQ,xsmi_?NumericQ,xsmf_?NumericQ,pg_,ac_,MSF_]:=
Module[{xf,sol},
Sol[mphi];
sol=NDSolve[{Yphi'[xx]==1/(xx*Hbar[mphi/xx]sst[xx]) (\[Lambda]hphi^2 (\[CapitalGamma]htophiphi[mphi/xx,mphi](*-((Yphi[xx]sst[xx])/nphieqq[xphi[xx],mphi])^2\[CapitalGamma]htophiphi[mphi/xphi[xx],mphi]*))-\[Kappa]^2 (Yphi[xx]sst[xx])/nphieqq[xphi[xx],mphi] \[CapitalGamma]\[Phi]to2\[Chi][mphi/xphi[xx],mphi,mchi]+lphi^3 averagedCif[xphi[xx],mphi]Yphi[xx]^2 ss[xx,mphi]^2 (-Yphi[xx]ss[xx,mphi]+nphieqq[xphi[xx],mphi])+\[Lambda]hphi^2 (C0a[mphi/xx](*-((Yphi[xx]sst[xx])/nphieqq[xphi[xx],mphi])^2C0a[mphi/xphi[xx]]*))),xphi'[xx]==-xphi[xx]/(xx*Hbar[mphi/xx]sst[xx]) (\[Lambda]hphi^2/Yphi[xx] (xphi[xx]/(3mphi) (\[CapitalGamma]htophiphi2[mphi/xx,mphi](*-((Yphi[xx]sst[xx])/nphieqq[xphi[xx],mphi])^2\[CapitalGamma]htophiphi2[mphi/xphi[xx],mphi]*)))-\[Kappa]^2/1(*Yphi[xx]*) (*Yphi[xx]*)sst[xx]/nphieq[xphi[xx],mphi] (xphi[xx]/(3mphi)*\[CapitalGamma]\[Phi]to2\[Chi]C2[mphi/xphi[xx],mphi,mchi])+\[Lambda]hphi^2/Yphi[xx] (xphi[xx]/(3mphi) (C2a[mphi/xx](*-((Yphi[xx]sst[xx])/nphieqq[xphi[xx],mphi])^2C2a[mphi/xphi[xx]]*)))+ss[xx,mphi]^2*Yphi[xx](Yphi[xx]ss[xx,mphi]-nphieqq[xphi[xx],mphi])*lphi^3 (C2cannibal[xphi[xx],mphi])+ss[xx,mphi] H[mphi/xx]*xphi[xx]/(3*mphi)*averagedP4E3e[xphi[xx]])+xphi[xx] Yphi'[xx]/Yphi[xx]-xphi[xx] (2sst'[xx])/(3ss[xx,mphi]),Yphi[xsmi]==nphieqq[xsmi/\[Xi]inf,mphi]/sst[xsmi](*Y0*),xphi[xsmi]==xsmi/\[Xi]inf},{xphi,Yphi},{xx,xsmi,xsmf},WorkingPrecision->MachinePrecision,(*Method->{"DiscontinuityProcessing"->False},*)PrecisionGoal->pg,AccuracyGoal->ac(*,MaxStepFraction\[Rule]MSF*)];
xf= sol[[1]][[1]][[2]]["Domain"][[1]][[2]];
Return[{sol[[1]][[1]][[2]](*=xphi*),sol[[1]][[2]][[2]](*=Yphi*),xf}]
];


Clear[solNDv2]
solNDv2[mphi_?NumericQ,\[Lambda]hphi_?NumericQ,lphi_?NumericQ,\[Xi]inf_?NumericQ,mchi_?NumericQ,\[Kappa]_?NumericQ,xsmi_?NumericQ,xint_?NumericQ,xsmf_?NumericQ,pg_,ac_,MSF_]:=
Module[{xf,sol1,sol2,\[Xi]2,xphi,Yphi,xfinal},
sol1=solND[mphi,\[Lambda]hphi,lphi,\[Xi]inf,mchi,\[Kappa],xsmi,xint(*xsmf*),7,15,MSF];
xf=xint (*sol1[[-1]](1-1/100)*);
\[Xi]2=xf/sol1[[1]][xf];
sol2=solND[mphi,\[Lambda]hphi,lphi,\[Xi]2,mchi,\[Kappa],xf,xsmf,pg,ac,MSF];
xfinal=sol2[[-1]];
xphi[x_]:=Piecewise[{{sol1[[1]][x],x<xf},{sol2[[1]][x],x>=xf}}];
Yphi[x_]:=Piecewise[{{sol1[[2]][x],x<xf},{sol2[[2]][x],x>=xf}}];
Return[{xphi,Yphi,xfinal}]
];


Clear[cBEwithA];
cBEwithA[ms_,ls_,lhs_,\[Xi]inf_,mchi_?NumericQ,\[Kappa]_?NumericQ,ai_,af_(*,Ni_*)]:=
Module[{solStandard,aFstandard,TsStand,NsS,Yf,datTemps,datNes,TsstandardAsFuncOfTds,TsstandardAsFuncOfNds,xff},
(*Note the factor 1/(1+\[Alpha]) in the temperature evolution. One has to include this term if one wants to account for thermal corrections-> <(n/(2E))>((dm^2)/dt)*)
solStandard=NDSolve[{Nn'[a]==(ls^3 averagedCif[ms/Ts[a],ms]/(a^7 H[Tstandard[a]]) Nn[a]^2 (Neq[Ts[a],ms,a]-Nn[a])+(lhs^2 a^2)/H[Tstandard[a]] (\[CapitalGamma]htophiphi[Tstandard[a],ms])-(\[Kappa]^2 a^2)/H[Tstandard[a]] Nn[a]/Neq[Ts[a],ms,a] \[CapitalGamma]\[Phi]to2\[Chi][Ts[a],ms,mchi]),
Ts'[a]==(*1/(1+\[Alpha])*)(lhs^2 /(a H[Tstandard[a]]) \[CapitalGamma]htophiphii2Energydens[Tstandard[a],ms]-\[Kappa]^2/(a H[Tstandard[a]]) Nn[a]/Neq[Ts[a],ms,a] \[CapitalGamma]\[Phi]to2chiEner[Ts[a],ms,mchi]-3/a ( Nn[a]Ts[a]/a^3)-Nn'[a]/Neq[Ts[a],ms,a] \[Rho]eq[Ts[a],ms])/(Nn[a]/Neq[Ts[a],ms,a] D\[Rho]eq[Ts[a],ms,0]-a^3 (\[Rho]eq[Ts[a],ms]/(Ts[a]Neq[Ts[a],ms,a]))^2 Nn[a]),Nn[ai]==Neq[\[Xi]inf Tstandard[ai],ms,ai](*Ni*),Ts[ai]==\[Xi]inf Tstandard[ai]
},{Nn,Ts},{a,ai,af//N},PrecisionGoal->9,AccuracyGoal->14,MaxSteps->10^7(*,Method->{"DiscontinuityProcessing"->False}*),MaxSteps->10^7];
aFstandard = solStandard[[1]][[1]][[2]]["Domain"][[1]][[2]];
TsStand[a_]:=Ts[a]/.solStandard[[1]];
NsS[a_]:=Nn[a]/.solStandard[[1]];
Yf=NsS[aFstandard]/(aFstandard^3 s[Tstandard[aFstandard]]);
datTemps=Table[{ms/Tstandard[10^a],ms/TsStand[10^a]},{a,Log10[1],Log10[aFstandard],(Log10[aFstandard]-Log10[1])/1000}];
datNes=Table[{ms/Tstandard[10^a],NsS[10^a]/(10^(3a) s[Tstandard[10^a]])},{a,Log10[1],Log10[aFstandard],(Log10[aFstandard]-Log10[1])/1000}];
xff=datTemps[[-1]][[1]];
TsstandardAsFuncOfTds=Interpolation[datTemps,InterpolationOrder->1,Method->"Spline"];
TsstandardAsFuncOfNds=Interpolation[datNes,InterpolationOrder->1,Method->"Spline"];
Return[{TsstandardAsFuncOfTds(*=(x,xphi)*),TsstandardAsFuncOfNds(*(x,Yphi)*),TsStand(*=(a,Tphi)*),NsS(*=(a,Nphi)*),aFstandard,xff}]
];


(* ::Text:: *)
(*Test*)


(*mphi = 10^-1(*10^-4*) (*gev*);
lhphi = 10^-8; (*portal*);
lphi =  10^-2(*strength of self/cannibal reactions*);
mchi = 5 10.^-6;(*gev*)
\[Kappa]=  10^-10;(*trilinear coupling of \[Phi]\[Chi]\[Chi]*)
\[Xi]inf=0.001(*initial ratio of Tdm/Tsm*);
Ti = 150.(*Gev, initial temperature*);
xsmi = mphi/Ti;
xf=7;
pg =(*Automatic*) 7;
ac=(*Automatic*)15;
MSF=Automatic;*)


(*tryEn=cBEwithA[mphi,lphi,lhphi,\[Xi]inf,mchi,0\[Kappa],(*ai=*)1.0,af(*=(10^9), In case of failure, make it smaller. *)]*)


(*tryEn0=cBEwithA[mphi,0. lphi,lhphi,\[Xi]inf,mchi,0\[Kappa],(*ai=*)1.0,(*5*10^4*)49621.01656475034`]*)


(*try=solND[mphi,lhphi,lphi,\[Xi]inf,mchi,0\[Kappa],xsmi,100,6,13,Automatic]*)


(*aff=Solve[Tstandard[a]==ms/4,a][[1]][[1]][[2]]*)


(*try0=solND[mphi,lhphi,0lphi,\[Xi]inf,mchi,0\[Kappa],xsmi,100,6,13,Automatic]*)


(*LogLogPlot[{tryEn[[3]][a],tryEn0[[3]][a]},{a,1,10^5},Frame->True,FrameLabel->{"a/\!\(\*SubscriptBox[\(a\), \(ewpt\)]\)","Tphi [GeV]"}]*)


(*LogLogPlot[{tryEn[[4]][a],tryEn0[[4]][a]},{a,1,10^5},Frame->True,FrameLabel->{"a/\!\(\*SubscriptBox[\(a\), \(ewpt\)]\)","Nphi [\!\(\*SuperscriptBox[\(GeV\), \(3\)]\)]"}]*)


(*LogLogPlot[{mphi/tryEn[[1]][x],mphi/tryEn0[[1]][x]},{x,xsmi,tryEn[[-1]]},Frame->True,FrameLabel->{"x=mphi/T","Tphi [GeV]"}]*)


(*LogLogPlot[{tryEn[[2]][x],tryEn0[[2]][x]},{x,xsmi,tryEn[[-1]]},Frame->True,FrameLabel->{"x=mphi/T","Yphi"}]*)


(* ::Text:: *)
(*Rates*)


(*\[Phi]->2\[Chi]*)


(*Clear[RatedecayC0]
RatedecayC0[sol_,k_,mphi_,mchi_]:=Module[{r,ratec0,Y=sol[[2]],Ts},
Sol[mphi];
Ts[x_]:=mphi/sol[[1]][x];
ratec0[x_]:=k^2/(sst[x]Y[x]) \[CapitalGamma]\[Phi]to2\[Chi][Ts[x],mphi,mchi];
Return[ratec0]
];*)


(*Clear[RatedecayC2]
RatedecayC2[sol_,k_,mphi_,mchi_]:=Module[{ratec2,Y=sol[[2]],Ts},
Sol[mphi];
Ts[x_]:=mphi/sol[[1]][x];
ratec2[x_]:=k^2/(sst[x]Y[x]) Ts[x]/(3mphi) \[CapitalGamma]\[Phi]to2\[Chi]C2[Ts[x],mphi,mchi];
Return[ratec2]
]*)


(* ::Text:: *)
(*Cannibal*)


(*Clear[RateCannibal3to2]
RateCannibal3to2[sol_,ls_,mphi_,mchi_]:=Module[{r,rate3to2,Y=sol[[2]],Ts},
Sol[mphi];
Ts[x_]:=mphi/sol[[1]][x];
rate3to2[x_]:=ls^3 averagedCif[sol[[1]][x],mphi]sst[x]^2 Y[x]^2;
Return[rate3to2]
];*)


(*Clear[RateCannibal2to3]
RateCannibal2to3[sol_,ls_,mphi_,mchi_]:=Module[{rate2to3,Y=sol[[2]],Ts},
Sol[mphi];
Ts[x_]:=mphi/sol[[1]][x];
rate2to3[x_]:=ls^3 averagedCif[sol[[1]][x],mphi]sst[x]Y[x]nphieqq[sol[[1]][x],mphi] ;
Return[rate2to3]
]*)


(*r3to2=RateCannibal3to2[try1,lphi,mphi,mchi];*)


(*r2to3=RateCannibal2to3[try1,lphi,mphi,mchi];*)


(*rdecay=RatedecayC0[try1,\[Kappa],mphi,mchi];*)


(*LogLogPlot[{r3to2[x],r2to3[x],(*rdecay[x],*)H[mphi/x]},{x,xsmi,0.5},Frame->True,FrameLabel->{"x=ms/T","Rates [GeV]"},PlotLegends->{"3->2","2->3","H"}]*)


(*LogLogPlot[{try0[[2]][x],try1[[2]][x]},{x,xsmi,3try0[[3]]},(*PlotRange->All,*) Frame->True, FrameLabel->{"x=mphi/T","Y"}]*)


(*LogLogPlot[{mphi/try0[[1]][x],mphi/try1[[1]][x]},{x,xsmi,try0[[3]]}, Frame->True, FrameLabel->{"x=mphi/T","Tphi [GeV]"}]*)


(* ::Text:: *)
(*Energy density*)


(*\[Rho][xs_,ms_,z_]:=z (ms^4  (xs BesselK[1,xs]+3 BesselK[2,xs]))/(2 \[Pi]^2 xs^2)*)


(*LogLogPlot[{\[Rho][tryEn[[1]][x], mphi,tryEn[[2]][x]/nphieqq[tryEn[[1]][x],mphi] sst[x] ]},{x,xsmi,10},Frame->True,FrameLabel->{"x=mphi/T","\!\(\*SubscriptBox[\(\[Rho]\), \(\[Phi]\)]\) [\!\(\*SuperscriptBox[\(GeV\), \(4\)]\)]"}]*)


(*LogLogPlot[{\[Rho][tryEn[[1]][x], mphi,tryEn[[2]][x]/nphieqq[tryEn[[1]][x],mphi] sst[x] ]/\[Rho][tryEn0[[1]][x], mphi,tryEn0[[2]][x]/nphieqq[tryEn0[[1]][x],mphi] sst[x] ],\[Rho][try[[1]][x], mphi,try[[2]][x]/nphieqq[try[[1]][x],mphi] sst[x] ]/\[Rho][try0[[1]][x], mphi,try0[[2]][x]/nphieqq[try0[[1]][x],mphi] sst[x] ]},{x,0.001,10},Frame->True,FrameLabel->{"x=mphi/T","\!\(\*SubscriptBox[\(\[Rho]\), \(\[Phi]\)]\)/\!\(\*SubscriptBox[\(\[Rho]\), \(\[Phi]\)]\)(l\[Phi]=0)"},PlotRange->All,PlotLegends->{"new solver","old solver"}]*)
