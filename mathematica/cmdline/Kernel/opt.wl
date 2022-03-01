#!/usr/bin/env wolframscript
Print["opt.wls"];

BeginPackage["cmdline`opt`",{"cmdline`log`"}]

Begin["`Private`"];
options2=<||>;
End[];
options=<|"defaultColors"->ColorData[97, "ColorList"]|>;

SetAttributes[UpdateOpts, Listable];
UpdateOpts[rule_] := AssociateTo[options, rule];
SetAttributes[AddOpts, Listable];
AddOpts[rule_] := If[MissingQ[options[First[rule] ] ], AssociateTo[options, rule] ];
Options[QSFcmdline]:=Normal@options;

BackupOpts[]:=(`Private`options2=options)&;
RestoreOpts[]:=(options=`Private`options2)&;


OptString[args_:If[$Notebooks,"",$ScriptCommandLine[[3;;]] ] ]:=
    StringRiffle[
        Select[args,
        (!(StringContainsQ[#,"<|"|"|>"]||StringContainsQ[#,"OptionGroup"]))& 
        ],
    "_"]<>"/";

GetOpts[args_:If[$Notebooks,"",$ScriptCommandLine[[3;;]] ]]:=Block[
    {keys=args[[1;; ;;2]],vals=args[[2;; ;;2]],opts},
    CHA[
    opts=AssociationThread[ToString/@keys,ToExpression/@vals];
    LOGS["Using options: " <> ToString[opts] ];
    AssociateTo[options,opts],
    AssociateTo[options,mainOpts->OptString[args]],
    "Invalid number of options. Options should come in pairs: key value"]
];




EndPackage[];