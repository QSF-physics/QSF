#!/usr/bin/env wolframscript
BeginPackage["cmdline`opt`",{"cmdline`log`"}]

QSFcmdline;
UpdateOpts;
AddOpts;
OptString;
GetOpts;
AppendToKey;
PopKey;

Begin["`Private`"];

options=<||>;
Options[QSFcmdline]:=Normal@options;

UpdateOpts[rule_] := AssociateTo[options, rule];
AddOpts[rule_] := If[MissingQ[options[First[rule]]], AssociateTo[options, rule] ];
AppendToKey[rule_] := If[MissingQ[options[First[rule]]], AssociateTo[options, rule/.{Key[x_]->x}];,

options=Merge[{options, <|rule/.{Key[x_]->x}|>},Flatten]; ];
PopKey[key_]:=(options= Delete[options,{key,-1}]);

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