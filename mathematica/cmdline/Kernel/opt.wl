#!/usr/bin/env wolframscript
BeginPackage["cmdline`opt`",{"cmdline`log`"}]

CurrentColorList;
GetColor;
UpdateOpts;
AddOpts;
BackupOpts;
RestoreOpts;
OptString;
GetOpts;
options;

Begin["`Private`"];
options2=<||>;

(* options=<|"ColorList"->ColorData[97, "ColorList"],"LineList"->{Dashed, DotDashed, Dotted},"colorIndex" -> <||>|>; *)
options=<|"colorIndex" -> <||>|>;
CurrentColorList[] := "DefaultPlotStyle" /. (Method /. Charting`ResolvePlotTheme[$PlotTheme, ListPlot]);
(* CurrentColorList[] := ColorData[97, "ColorList"]; *)
GetColor[key_] := (If[MissingQ[options["colorIndex"][key]], AssociateTo[options["colorIndex"], key -> (Length[options["colorIndex"]] + 1)]]; Part[CurrentColorList[], options["colorIndex"][key]]);

SetAttributes[UpdateOpts, Listable];
UpdateOpts[rule_] := AssociateTo[options, rule];
SetAttributes[AddOpts, Listable];
AddOpts[rule_] := If[MissingQ[options[First[rule] ] ], AssociateTo[options, rule] ];
Options[QSFcmdline]:=Normal@options;

BackupOpts[]:=(options2=options)&;
RestoreOpts[]:=(options=options2)&;


OptString[args_:If[$Notebooks,"",Quiet@Check[$ScriptCommandLine[[3;;]],""]  ] ]:=
    StringRiffle[
        Select[args,
        (!(StringContainsQ[#,"<|"|"|>"]||StringContainsQ[#,"OptionGroup"]))& 
        ],
    "_"]<>"/";

GetOpts[args_:If[$Notebooks,"",Quiet@Check[$ScriptCommandLine[[3;;]],""] ]]:=Module[
    {keys=args[[1;; ;;2]],vals=args[[2;; ;;2]],opts},
    CHA[
    opts=AssociationThread[ToString/@keys,ToExpression/@vals];
    LOGS["Using options: " <> ToString[opts] ];
    AssociateTo[options,opts],
    AssociateTo[options,mainOpts->OptString[args]],
    "Invalid number of options. Options should come in pairs: key value"]
];



End[];
EndPackage[];