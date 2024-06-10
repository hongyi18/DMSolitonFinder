(* ::Package:: *)

(* Mathematica Package *)
(* Title: DMSolitonFinder *)
(* Author: Hong-Yi Zhang *)
(* Website: https://hongyi18.github.io *)
(* Package Version: 1.0 *)
(* Mathematica Version: 13.3.1 *)
(* Date: Jun 10, 2024 *) 

BeginPackage["DMSolitonFinder`"];
Print["DMSolitonFinder 1.0 has been loaded."];
Print["Please cite ", Hyperlink["arXiv:2406.05031", "https://arxiv.org/abs/2406.05031"], " if you find this package helpful for your research."]


MassDensity::usage = "
MassDensity[f, r] gives the mass density \[Rho].
MassDensity[f, r, NGI\[Rule]\[Xi]] gives the mass density \!\(\*SubscriptBox[\(\[Rho]\), \(\[Xi]\)]\) with a nonminimal gravitational interaction \[Xi].";

Eqf::usage = "
Eqf[f, \[CapitalPsi], r] gives the Schroedinger equation for soliton profiles.
Eqf[f, \[CapitalPsi], r, SI\[Rule]\[Lambda], NGI\[Rule]\[Xi], Spin\[Rule]\[Sigma]] gives the Schroedinger equation with a self-interaction \[Lambda], nonminimal gravitational interaction \[Xi], and particle spin \[Sigma]";

EqPsi::usage = "
EqPsi[f, \[CapitalPsi], r] gives Poisson's equation for the gravitational potential.
EqPsi[f, \[CapitalPsi], r, NGI\[Rule]\[Xi]] gives Poisson's equation modified by a nonminimal gravitational interaction \[Xi].";

SolveFields::usage = "
SolveFields[fi, \[CapitalPsi]i, rf] returns a list of numerical solutions {f[r], \[CapitalPsi][r]} up to the radius rf, given the boundary conditions f[MinRadius]=fi and \[CapitalPsi][MinRadius]=\[CapitalPsi]i.
Possible options include {MinRadius\[Rule]10^-7} and those for Eqf, EqPsi, and NDSolveValue. MinRadius specifies the inner boundary and the default value is 10^-7.";

ShootFields::usage = "
ShootFields[fi] returns a list of numerical soliton solutions {{f[r], \[CapitalPsi][r]}, rf}, given the boundary condition f[MinRadius]=fi. 
The solutions are valid up to the radius rf, and MinRadius specifies the inner boundary with the default value given by 10^-7. 
Possible options include those for SolveFields and {AmpTolerance\[Rule]10^-4, AmpCheckRangeRatio\[Rule]90/100, FirstMinimumLocationRatio\[Rule]80/100, IfPrintProgress\[Rule]0, InitialValuePsi\[Rule]-1, InitialStepPsi\[Rule]1/10, InitialMaxRadius\[Rule]5, InitialStepRadius\[Rule]1, MaxIteration\[Rule]10^3}. 

If one aims for soliton profiles with better convergence at large radius, decreasing AmpTolerance might help.
If the shooting method fails to find soliton solutions, slightly changing AmpCheckRangeRatio and FirstMinimumLocationRatio might help.
If one would like to see intermediate values of {\[CapitalPsi]i, d\[CapitalPsi], rf, drf} during the shooting loops, e.g., for debugging, set IfPrintProgress=1.
If the shooting method is stopped at very few loops or even at the beginning, changing InitialValuePsi and InitialStepPsi might help.
If soliton profiles do not behave well at large radius, changing InitialMaxRadius and InitialStepRadius might help.
If soliton solutions are not found after the maximum loops of calculations, try increasing MaxIteration.
For high-precision calculations, try increasing WorkingPrecision, AccuracyGoal and PrecisionGoal.";

CalMass::usage = "
CalMass[f, rf] calculates the mass enclosed by the radius rf for a soliton solution f[r].
Possible options include those for MassDensity and NIntegrate and {MinRadius\[Rule]10^-7}. MinRadius specifies the radius where the integration begins.";

CalRadius::usage = "
CalRadius[f, rf] calculates the radius of a soliton solution f[r] that encloses 95 percent of its total mass, which is evaluated by integrating the mass density up to the radius rf.
Possible options include those for CalMass and {EnclosedMass\[Rule]95/100, InitialGuess\[Rule]25/100}, where EnclosedMass specifies the percentage of mass enclosed by the radius and InitialGuess is a guess for the location of the soliton radius in terms of rf.";

CalParticleNumber::usage = "
CalParticleNumber[f, rf] calculates the total particle number enclosed by the radius rf for a soliton solution f[r].
Possible options include those for NIntegrate and {MinRadius\[Rule]10^-7}. MinRadius specifies the radius where the integration begins.";

CalFrequency::usage = "
CalFrequency[\[CapitalPsi], rf] calculates the oscillating frequency (chemical potential) for a soliton, whose modified gravitational potential is \[CapitalPsi][r]. To calculate the frequency accurately, the soliton profile f[r] must be significantly suppressed at the radius rf.
Possible options include those for FindRoot.";

CalEnergy::usage = "
CalEnergy[f, \[CapitalPsi], \[Mu], rf] calculates the energy enclosed by the radius rf for a soliton, whose profiles are f[r] and \[CapitalPsi][r] and oscillating frequency is \[Mu].
Possible options include those for MassDensity and NIntegrate and {SI\[Rule]0, Spin\[Rule]0, MinRadius\[Rule]10^-7}. SI and Spin specify the self-interaction \[Lambda] and particle spin \[Sigma]. MinRadius specifies the radius where the integration begins.";


Begin["`Private`"];


(* Set up field equations. *)
Unprotect[MassDensity];
Options[MassDensity] = {NGI->0};
SyntaxInformation[MassDensity] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
MassDensity[f_, r_, OptionsPattern[]] := f[r]^2+2 OptionValue[NGI] (Derivative[1][f][r]^2+f[r] ((2 Derivative[1][f][r])/r+ Derivative[2][f][r]));
Protect[MassDensity];

Unprotect[Eqf];
Options[Eqf] = Join[{SI->0, Spin->0}, Options[MassDensity]];
SyntaxInformation[Eqf] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};
Eqf[f_, \[CapitalPsi]_, r_, opts:OptionsPattern[]] := -(f''[r]+2/r f'[r])/2 + \[CapitalPsi][r] f[r] + (1-Abs[OptionValue[Spin]]/3)/8 OptionValue[SI] f[r]^3 + OptionValue[NGI] MassDensity[f, r, Evaluate@FilterRules[{opts, Options[Eqf]}, Options[MassDensity]]] f[r];
Protect[Eqf];

Unprotect[EqPsi];
Options[EqPsi] = Options[MassDensity];
SyntaxInformation[EqPsi] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};
EqPsi[f_, \[CapitalPsi]_, r_, opts:OptionsPattern[]] := -\[CapitalPsi]''[r]-2/r \[CapitalPsi]'[r] + MassDensity[f, r, Evaluate@FilterRules[Join[{opts}, Options[EqPsi]], Options[MassDensity]]];
Protect[EqPsi];


(* Solve field equations. *)
Unprotect[SolveFields];
Options[SolveFields] = Join[{MinRadius->10^-7}, Options[Eqf], Options[EqPsi], Options[NDSolveValue]];
SetOptions[SolveFields, AccuracyGoal->8/10 MachinePrecision, PrecisionGoal->8/10 MachinePrecision, MaxSteps->10^7];
SyntaxInformation[SolveFields] = {"ArgumentsPattern" -> {_, _, _, OptionsPattern[]}};
SolveFields[fi_, \[CapitalPsi]i_, rf_, opts:OptionsPattern[]] := NDSolveValue[
	{Eqf[f, \[CapitalPsi], r, Evaluate@FilterRules[Join[{opts}, Options[SolveFields]], Options[Eqf]]]==0
	, EqPsi[f, \[CapitalPsi], r, Evaluate@FilterRules[Join[{opts}, Options[SolveFields]], Options[EqPsi]]]==0
	, f[OptionValue[MinRadius]]==fi, f'[OptionValue[MinRadius]]==0
	, \[CapitalPsi][OptionValue[MinRadius]]==\[CapitalPsi]i, \[CapitalPsi]'[OptionValue[MinRadius]]==0}
	, {f, \[CapitalPsi]}, {r, OptionValue[MinRadius], rf}
	, Evaluate@FilterRules[Join[{opts}, Options[SolveFields]], Options[NDSolveValue]]
	];
Protect[SolveFields];


(* Calculate the first local minimum of f[r]. *)
(* This function will be used to determine if \[CapitalPsi]i needs to be increased in a slower rate and if rf needs to be increased. *)
calFirstMinimum[f_, ri_, rf_, workingPrecision_] := FindMinimum[f[r], {r, rf/100, ri, rf}, AccuracyGoal->Infinity, WorkingPrecision->workingPrecision]//Quiet;

(* Calculate the maximum amplitude of f[r] at nPointLargeRadius points within AmpCheckRangeRatio*rf and rf. *)
(* This function will be used to verify if the solution satisfies the amplitude tolerance at large radius . *)
calAmpLargeRadius[f_, rf_, AmpCheckRangeRatio_, nPointLargeRadius_] := Max[ Abs@Array[f, nPointLargeRadius, {AmpCheckRangeRatio rf, rf}] ];

(* Solve soliton profiles by a dynamical numerical shooting method. *)
Unprotect[ShootFields];
Options[ShootFields] = Join[{AmpTolerance->10^-4, AmpCheckRangeRatio->90/100, FirstMinimumLocationRatio->80/100, IfPrintProgress->0
	, InitialValuePsi->-1, InitialStepPsi->1/10, InitialMaxRadius->5, InitialStepRadius->1, MaxIteration->10^3}, Options[SolveFields]];
SyntaxInformation[ShootFields] = {"ArgumentsPattern" -> {_, OptionsPattern[]}};
ShootFields[fi_, opts:OptionsPattern[]] := Block[
	(* Local parameters for internal calculations. Do not modify them. *)
	{\[CapitalPsi]i=OptionValue[InitialValuePsi], d\[CapitalPsi]=OptionValue[InitialStepPsi], rf=OptionValue[InitialMaxRadius], drf=OptionValue[InitialStepRadius]
	, ifSolsUpdate\[CapitalPsi]iError=0, ifSolsUpdaterfError=0, nPointLargeRadius=5, precision=10^(-OptionValue[WorkingPrecision]+1)
	, iffiReachCriticalValue=0, ifInitialError=0, ifInitial\[CapitalPsi]iTooBig=0, ifd\[CapitalPsi]TooSmall=0, ifdrfTooSmall=0, ifMaxIterationReached=0
	, sols, firstMinimum, ampLargeRadius, \[CapitalPsi]iTest, rfTest, nStop}
	
	(* Solve field equations. If error occurs in solving the equations, return the trivial solutions of the system. *)
	(* For the shooting algorithm to work, the initial \[CapitalPsi]i should be such small that the first local minimum of f[r] is negative. *)
	, If[OptionValue[NGI]!=0, If[fi>=1/2/Abs[OptionValue[NGI]], iffiReachCriticalValue=1]]
	; sols = Quiet@Check[SolveFields[fi, \[CapitalPsi]i, rf, Evaluate@FilterRules[Join[{opts}, Options[ShootFields]], Options[SolveFields]]]
		, ifInitialError=1; SolveFields[0, 0, rf, Evaluate@FilterRules[Join[{opts}, Options[ShootFields]], Options[SolveFields]]]]
	; If[iffiReachCriticalValue==0 && ifInitialError==0
		, firstMinimum = calFirstMinimum[sols[[1]], OptionValue[MinRadius], rf, OptionValue[WorkingPrecision]]
		; If[firstMinimum[[1]]>0 && N[firstMinimum[[2,1,2]],3]<N[rf,3], ifInitial\[CapitalPsi]iTooBig=1]]
	
	; If[iffiReachCriticalValue==0 && ifInitialError==0 && ifInitial\[CapitalPsi]iTooBig==0
		, If[OptionValue[IfPrintProgress]==1, Print["Values of {\[CapitalPsi]i, d\[CapitalPsi], rf, drf} are:"] ] 
		; Do[
			(* Try increasing \[CapitalPsi]i. *)
			\[CapitalPsi]iTest = \[CapitalPsi]i+d\[CapitalPsi]
			; sols = Quiet@Check[SolveFields[fi, \[CapitalPsi]iTest, rf, Evaluate@FilterRules[Join[{opts}, Options[ShootFields]], Options[SolveFields]]]
				, ifSolsUpdate\[CapitalPsi]iError=1; SolveFields[fi, \[CapitalPsi]i, rf, Evaluate@FilterRules[Join[{opts}, Options[ShootFields]], Options[SolveFields]]]]
			; If[ifSolsUpdate\[CapitalPsi]iError==0, \[CapitalPsi]i = \[CapitalPsi]iTest; firstMinimum = calFirstMinimum[sols[[1]], OptionValue[MinRadius], rf, OptionValue[WorkingPrecision]]]
			
			(* If the first local minimum of soliton profile is positive and near the boundary rf, then rf is too small. *)
			(* "Near the boundary" is determined by whether the first minimum is located within OptionValue[FirstMinimumLocationRatio]*rf. *)
			; If[firstMinimum[[1]]>0 && firstMinimum[[2,1,2]]>OptionValue[FirstMinimumLocationRatio] rf
				(* Try increasing rf. Increase rf with bigger steps if the first minimum is located at rf. *)
				, If[N[firstMinimum[[2,1,2]],3]<N[rf,3], rfTest = rf+drf, rfTest = rf+10 drf]
				; rfTest = If[drf==OptionValue[InitialStepRadius]
					, rfTest
					, Min[rfTest, Max[rf+drf, 2/(OptionValue[AmpCheckRangeRatio]+1) firstMinimum[[2,1,2]]]]] (* Avoid increasing rf too much. *)
				; sols = Quiet@Check[
					SolveFields[fi, \[CapitalPsi]i-d\[CapitalPsi], rfTest, Evaluate@FilterRules[Join[{opts}, Options[ShootFields]], Options[SolveFields]]]
					; SolveFields[fi, \[CapitalPsi]i, rfTest, Evaluate@FilterRules[Join[{opts}, Options[ShootFields]], Options[SolveFields]]]
					, ifSolsUpdaterfError=1; SolveFields[fi, \[CapitalPsi]i, rf, Evaluate@FilterRules[Join[{opts}, Options[ShootFields]], Options[SolveFields]]]]
				; If[ifSolsUpdaterfError==0, rf=rfTest]]
			
			; If[OptionValue[IfPrintProgress]==1, Print["For n=", n, ": ", {\[CapitalPsi]i, d\[CapitalPsi], rf, drf}] ] (* Usually used for debugging purposes. *)
			(* End the calculation if the desired amplitude is achieved. *)
			; If[firstMinimum[[1]]>0
				, ampLargeRadius=calAmpLargeRadius[sols[[1]], rf, OptionValue[AmpCheckRangeRatio], nPointLargeRadius]
				; If[ampLargeRadius<=OptionValue[AmpTolerance] fi, Break[]]
				; \[CapitalPsi]i=\[CapitalPsi]i-d\[CapitalPsi]]
			; If[(firstMinimum[[1]]>0 && firstMinimum[[2,1,2]]<=OptionValue[FirstMinimumLocationRatio] rf)||(ifSolsUpdate\[CapitalPsi]iError==1), d\[CapitalPsi]=d\[CapitalPsi]/2; ifSolsUpdate\[CapitalPsi]iError=0]
			; If[ifSolsUpdaterfError==1, drf=drf/2; ifSolsUpdaterfError=0]
			
			(* Break the loop if d\[CapitalPsi] is too small, drf is too small, or the maximum number of loops is reached. *)
			; If[Abs@d\[CapitalPsi]<precision, ifd\[CapitalPsi]TooSmall=1; nStop=n; Break[]]
			; If[Abs@drf<rf/10/OptionValue[MaxIteration], ifdrfTooSmall=1; nStop=n; Break[]]
			; If[n==OptionValue[MaxIteration], ifMaxIterationReached=1; nStop=n; Break[]]	
			, {n, OptionValue[MaxIteration]}]]
	; If[iffiReachCriticalValue==1, Print["For fi=", fi, ", the field amplitude has reached the critical value due to nonminiaml gravitational interactions and no soliton solutions are expected to exist."]]
	; If[ifInitialError==1
		, Print["For fi=", fi, ", initial guesses are bad or some other error messages are generated."]
		; SolveFields[fi, \[CapitalPsi]i, rf, Evaluate@FilterRules[Join[{opts}, Options[ShootFields]], Options[SolveFields]]] ]
	; If[ifInitial\[CapitalPsi]iTooBig==1, Print["For fi=", fi, ", the initial value \[CapitalPsi]i is too big. Try decreasing InitialValuePsi."]]
	; If[ampLargeRadius > OptionValue[AmpTolerance] fi
		, If[ifd\[CapitalPsi]TooSmall==1, Print["For fi=", fi, ", requested precision is reached after ", nStop, " calculations, relative amplitude at large radius is ", ampLargeRadius/fi, "."]]
		; If[ifdrfTooSmall==1, Print["For fi=", fi, ", increasing rf doesn't improve results after ", nStop, " calculations, relative amplitude at large radius is ", ampLargeRadius/fi, "."]]
		; If[ifMaxIterationReached==1, Print["For fi=", fi, ", desired solution is not found after ", nStop, " calculations, relative amplitude at large radius is ", ampLargeRadius/fi, ". Try increasing MaxIteration."]]]
	; {sols, rf}];
Protect[ShootFields];


(* Calculate physical properties of solitons. *)
Unprotect[CalMass];
Options[CalMass] = Join[{MinRadius->10^-7}, Options[MassDensity], Options[NIntegrate]];
SetOptions[CalMass, AccuracyGoal->8, PrecisionGoal->8, Method->{Automatic,"SymbolicProcessing"->0}];
SyntaxInformation[CalMass] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
CalMass[f_, rf_, opts:OptionsPattern[]] := 4 \[Pi] NIntegrate[r^2 MassDensity[f, r, Evaluate@FilterRules[Join[{opts}, Options[CalMass]], Options[MassDensity]]]
	, {r, OptionValue[MinRadius], rf}
	, Evaluate@FilterRules[Join[{opts}, Options[CalMass]], Options[NIntegrate]]];
Protect[CalMass];

Unprotect[CalRadius];
Options[CalRadius] = Join[{EnclosedMass->95/100, InitialGuess->25/100}, Options[CalMass]];
SyntaxInformation[CalRadius] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
CalRadius[f_, rf_, opts:OptionsPattern[]] := FindRoot[
	CalMass[f, R, Evaluate@FilterRules[Join[{opts}, Options[CalRadius]], Options[CalMass]]]==OptionValue[EnclosedMass] CalMass[f, rf, Evaluate@FilterRules[Join[{opts}, Options[CalRadius]], Options[CalMass]]]
	, {R, OptionValue[InitialGuess] rf, OptionValue[MinRadius], rf}
	, AccuracyGoal->OptionValue[AccuracyGoal], PrecisionGoal->OptionValue[PrecisionGoal]][[1,2]]//Quiet;
Protect[CalRadius];

Unprotect[CalParticleNumber];
Options[CalParticleNumber] = Join[{MinRadius->10^-7}, Options[NIntegrate]];
SetOptions[CalParticleNumber, AccuracyGoal->8, PrecisionGoal->8, Method->{Automatic,"SymbolicProcessing"->0}];
SyntaxInformation[CalParticleNumber] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
CalParticleNumber[f_, rf_, opts:OptionsPattern[]] := 4 \[Pi] NIntegrate[r^2 f[r]^2
	, {r, OptionValue[MinRadius], rf}
	, Evaluate@FilterRules[Join[{opts}, Options[CalParticleNumber]], Options[NIntegrate]]];
Protect[CalParticleNumber];

Unprotect[CalFrequency];
Options[CalFrequency] = Options[FindRoot];
SetOptions[CalFrequency, AccuracyGoal->8, PrecisionGoal->8];
SyntaxInformation[CalFrequency] = {"ArgumentsPattern" -> {_, _, OptionsPattern[]}};
CalFrequency[\[CapitalPsi]_, rf_, opts:OptionsPattern[]] := Block[{\[Mu]1, \[Mu]2}
	, \[Mu]1 = FindRoot[rf(\[CapitalPsi][rf]-\[Mu]) == 97/100 rf(\[CapitalPsi][97/100 rf]-\[Mu]), {\[Mu], 1/10}, Evaluate@FilterRules[Join[{opts}, Options[CalFrequency]], Options[FindRoot]]][[1,2]] //Quiet
	; \[Mu]2 = FindRoot[97/100 rf(\[CapitalPsi][97/100 rf]-\[Mu]) == 95/100 rf(\[CapitalPsi][95/100 rf]-\[Mu]), {\[Mu], 1/10}, Evaluate@FilterRules[Join[{opts}, Options[CalFrequency]], Options[FindRoot]]][[1,2]] //Quiet
	; If[\[Mu]1<0, Print["Warning, the frequency is negative."]]
	; If[Abs[(\[Mu]1-\[Mu]2)/\[Mu]1]>10^-3, Print["Frequency is inaccurate because the input soliton profile is not small enough at large radius."]]
	; (\[Mu]1+\[Mu]2)/2];
Protect[CalFrequency];

Unprotect[CalEnergy];
Options[CalEnergy] = Join[{SI->0, Spin->0, MinRadius->10^-7}, Options[MassDensity], Options[NIntegrate]];
SetOptions[CalEnergy, AccuracyGoal->8, PrecisionGoal->8, Method->{Automatic,"SymbolicProcessing"->0}];
CalEnergy[f_, \[CapitalPsi]_, \[Mu]_, rf_, opts:OptionsPattern[]] := 4 \[Pi] NIntegrate[r^2 (1/2 Derivative[1][f][r]^2 + 1/2 (\[CapitalPsi][r]-\[Mu]) MassDensity[f, r, Evaluate@FilterRules[Join[{opts}, Options[CalEnergy]], Options[MassDensity]]] + (1-Abs[OptionValue[Spin]]/3)/16 OptionValue[SI] f[r]^4)
	, {r, OptionValue[MinRadius], rf}
	, Evaluate@FilterRules[Join[{opts}, Options[CalEnergy]], Options[NIntegrate]]];
Protect[CalEnergy];


End[];
EndPackage[];
