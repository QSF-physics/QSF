#!/usr/bin/env wolframscript
BeginPackage["cmdline`opt`",{"cmdline`log`"}]

QSFcmdline;
UpdateOpts;
AddOpts;
OptString;
GetOpts;
AppendToKey;
PopKey;
COptionValue;
options;

Begin["`Private`"];

options=<||>;
Options[QSFcmdline]:=Normal@options;
COptionValue[a_,b_]:=OptionValue[{QSFcmdline,a},b];
COptionValue[a,b__]:=OptionValue[{QSFcmdline,a},{b}];
COptionValue[b_]:=OptionValue[QSFcmdline,b];

(*ListPlotDefault->Explicit Pass->CMDLine->Final*)
(* Options[TK7] = {"sad" -> 3};
TK7[k : OptionsPattern[]] := {COptionValue[{k, TK7, ListPlot}, "sad"],
COptionValue[{k, ListPlot}, PlotRange]}; *)
ClearOpts[]:=(options=<||>);
UpdateOpts[rule_]:=AssociateTo[options,rule];
AddOpts[rule_]:=AssociateTo[options,rule];
AppendToKey[key_String][val_]:=If[KeyExistsQ[options,key],
options=MapAt[Append[#,val]&,options,key],AddOpts[key->{val}]];
PopKey[key_]:=(options=Delete[options,{key,-1}]);

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