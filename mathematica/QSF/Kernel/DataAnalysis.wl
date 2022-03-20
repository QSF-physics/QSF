BeginPackage["QSF`DataAnalysis`", {"cmdline`log`"}];

FourierAt;

Begin["`Private`"];

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