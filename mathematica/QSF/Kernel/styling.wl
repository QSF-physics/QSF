#!/usr/bin/env wolframscript
BeginPackage["QSF`styling`"];

PrettyPlots;
FontFix;

Begin["`Private`"];
(* For making any plot pretty! *)
PrettyPlots[]:=
    Map[SetOptions[Print["Changing style of ", #]; #, BaseStyle -> 
        {FontFamily -> "Latin Modern Math", FontSize -> 12,FontColor->Black}] &,
        {Plot,ListLinePlot, ArrayPlot, "Graphics"}
    ];

SetAttributes[FontFix, Listable];
FontFix[labs_String] := Style[labs, FontFamily -> "Latin Modern Math", FontColor->Black];

End[];
EndPackage[];