#!/usr/bin/env wolframscript
BeginPackage["cmdline`log`"];
Begin["`Private`"];

COLORS = <|
   "green" -> "\033[1;32m", 
   "yellow" -> "\033[1;33m",
   "blue" -> "\033[1;34m", 
   "magenta" -> "\033[1;35m",
   "red" -> "\033[1;31m",
   "black" -> "\033[1;30m",
   "reset" -> "\033[1;0m"
|>;

template = StringTemplate["``[``````=``````]"];

End[];
SetAttributes[VARVALUE, HoldAll];
VARVALUE[msg_] := `Private`template[
   `Private`COLORS["reset"],
   `Private`COLORS["blue"],
   ToString@Unevaluated@msg,
   `Private`COLORS["reset"],
   `Private`COLORS["green"],
   ToString@msg, 
   `Private`COLORS["reset"]]/; ValueQ[msg];
(* VARVALUE[msg_] := Sequence["V:", msg]; *)
VARVALUE[msg_] := Sequence[`Private`COLORS[[
      If[
         OptionValue[{QSFcmdline, LOG},"LC"]===None,
         Mod[OptionValue[{QSFcmdline, LOG},"LL"],Length[`Private`COLORS]-1,1],
         Key[OptionValue[{QSFcmdline, LOG},"LC"] ] 
      ]
   ]],msg];

LL[]:=1;
Options[LOG]={"LL":>LL[], "LC"->None};
LOG[msgs___,op:OptionsPattern[{QSFcmdline, LOG}] ]:=Print[
   If[OptionValue[{QSFcmdline, LOG},"LL"]==0,"",
   If[OptionValue[{QSFcmdline, LOG},"LL"]==1,"+","|"] ], 
   StringRepeat["-", Mod[OptionValue[{QSFcmdline, LOG},"LL"],50] ],
   `Private`COLORS[[
      If[
         OptionValue[{QSFcmdline, LOG},"LC"]===None,
         Mod[OptionValue[{QSFcmdline, LOG},"LL"],Length[`Private`COLORS]-1,1],
         Key[OptionValue[{QSFcmdline, LOG},"LC"] ] 
      ]
   ]],
   msgs, 
   `Private`COLORS["reset"] 
];

LOGF[f_, args_, rhs_] := Module[{res},
   LOG["// ", ToString[f]];   
   LL[]++;
   res = rhs;
   LL[]--;
   res
];
ClearAll[decorator];
decorator /: SetDelayed[decorator[dec_][f_[args___]], rhs_]:=
  f[a:PatternSequence[args]]:=
     dec[Unevaluated @ f, Unevaluated @ {a}, Unevaluated @ rhs];

LOGI[msgs___]:=LOG[msgs,"LC"->"blue","LL"->0];
LOGII[msgs__]:=LOG[msgs,"LC"->"magenta","LL"->0];
LOGW[msgs__]:=LOG[msgs,"LC"->"yellow","LL"->0];
LOGS[msgs__]:=LOG[msgs,"LC"->"green","LL"->0];
LOGE[msgs__]:=(LOG[msgs,"LC"->"red","LL"->0];Pause[0.2];Abort[];);
LOGQ[msgs__]:=LOG[msgs,"LC"->"black","LL"->0];
SetAttributes[LOGL,{Listable}];
LOGL[msgs__,op:OptionsPattern[{QSFcmdline, LOGL}] ]:=LOG[msgs,op];
SetAttributes[CHA, HoldFirst];
CHA[exp_, msgs__]:=Quiet@Check[exp, LOGE[msgs]];


Attributes[LOGV] = {HoldAll, SequenceHold,Listable};
LOGV[msgs__] := LOG[VARVALUE /@ Unevaluated[msgs]];
LOGV[msgs_] := LOG[VARVALUE @ Unevaluated[msgs]];

(* Debug Print *)
DP[a_List]:=(LOG["List ->"]; LL[]++; 
If[Length[a]>10, LOG["Large list, First: "]; DP[First[a]], Map[DP,a] ]; 
LL[]--;);
DP[a_Association]:=(LOG["Assoc ->"]; LL[]++; 
KeyValueMap[(LOG["key: ",#1]; DP[#2];)&,a]; LL[]--;);
DP[a_WF]:=(LOG["WF ->"]; LL[]++; DP[First@a]; DP[Last@a]; LL[]--;);
DP[a_Data]:=(LOG["Data ->"]; LL[]++; DP[First@a]; LL[]--;);
DP[a_Number]:=a/;!NumberQ[a];
DP[a]:=LOG[a];

EndPackage[];