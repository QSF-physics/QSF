BeginPackage["QSF`DataAnalysis`", {"cmdline`log`","cmdline`opt`"}];
OpenBin;
IntegerChop;
FourierAt;
GaussianFilterOpt;
Begin["`Private`"];
OpenBin:=CHA[OpenRead[#,BinaryFormat->True], "File "<>#<>" cannot be found or cannot be read as binary."] &;
IntegerChop = With[{r = Round[#]}, r + Chop[# - r]] &;

Options[GaussianFilterOpt] = {"GaussianBlurMult"->1,"GaussianBlurRadius" -> None};
GaussianFilterOpt[data_List,opt:OptionsPattern[]]:=
With[{rad=COptionValue[{opt,GaussianFilterOpt},"GaussianBlurRadius"]},
  If[rad===None, data, LOG["Applying GaussianBlur with radius ",rad, " mult: ",COptionValue[{opt,GaussianFilterOpt},"GaussianBlurMult"]];
  GaussianFilter[data,{rad COptionValue[{opt,GaussianFilterOpt},"GaussianBlurMult"]}]]
]; 

(* GaussianFilterOpt[ass_Association,mult_,opt:OptionsPattern[]]:=
Map[GaussianFilterOpt[#,mult,opt]&,ass]; *)
(* Makes a fourier transform at dimension d *)
FourierAt[tens_List, d_Integer, inv_ : False] := Module[
    {dim = Length[Dimensions[tens]]},

    Transpose[
        Map[If[inv, InverseFourier[#], Fourier[#] ] &, 
            Transpose[tens, d <-> dim], 
            {dim - 1}], 
    d <-> dim]
];

decorator[LOGF]@
FourierAt[tens_List, dims_List, inv_ : False] := Module[
    {pos = Flatten@Position[dims, 1], dim = Length[Dimensions[tens]], ret = tens},
    If[Length[dims] != dim, LOGE[dims," ",dim, " ",Dimensions[tens]]; Abort[]];
   
    If[AllTrue[dims, 1 == # &], 
        If[inv, InverseFourier[ret], Fourier[ret] ],
        Fold[FourierAt[#1, #2, inv] &, ret, pos] 
    ]
];

End[];
EndPackage[];