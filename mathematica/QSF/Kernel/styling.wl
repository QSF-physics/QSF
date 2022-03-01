#!/usr/bin/env wolframscript
BeginPackage["QSF`styling`"];

(* For making any plot pretty! *)
PrettyPlots[]:=Block[{style={FontFamily -> "Latin Modern Math", FontSize -> 12}},
    Map[SetOptions[#, BaseStyle -> style]&,
    {Plot,ListLinePlot}]
];

SetAttributes[FontFix, Listable];
FontFix[str_String] := Style[str, FontFamily -> "Latin Modern Math"];


EndPackage[];