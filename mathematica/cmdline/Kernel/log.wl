#!/usr/bin/env wolframscript
BeginPackage["cmdline`log`"];
LOG;
LL;
LOGF;
LOGI;
LOGJ;
LOGW;
LOGS;
LOGE;
LOGQ;
LOGL;
LOGV;
CHA;
decorator;


Begin["`Private`"];

CSting[x_]=If[$Notebooks, "",x];
CPrint=If[$Notebooks,PrintTemporary,Print];
COLORS = Map[CSting,<|
   "green"->"\033[1;32m", 
   "yellow"->"\033[1;33m",
   "blue"->"\033[1;34m", 
   "magenta"->"\033[1;35m",
   "red"->"\033[1;31m",
   "black"->"\033[1;30m",
   "reset"->"\033[1;0m"
|>];

template = StringTemplate["``[``````=``````]"];

SetAttributes[VARVALUE, HoldAll];
VARVALUE[msg_] := template[
   COLORS["reset"],
   COLORS["blue"],
   StringReplace[ToString[Unevaluated@msg,InputForm],"$"~~DigitCharacter..->""],
   COLORS["reset"],
   COLORS["green"],
   msg, 
   COLORS["reset"]]/; ValueQ[msg];
(* VARVALUE[msg_] := Sequence["V:", msg]; *)
VARVALUE[msg_] := Sequence[COLORS[[
      If[
         OptionValue[{QSFcmdline, LOG},"LC"]===None,
         Mod[OptionValue[{QSFcmdline, LOG},"LL"],Length[COLORS]-1,1],
         Key[OptionValue[{QSFcmdline, LOG},"LC"] ] 
      ]
   ]],msg];

LL[]:=1;
SetAttributes[ToStringList,{Listable}];
ToStringList[msg_]:=If[StringQ[msg],msg,ToString[msg,InputForm]];


LogColor[op:OptionsPattern[{QSFcmdline, LOG}]]:=With[{
   l=OptionValue[{QSFcmdline, LOG},"LL"],
   c=OptionValue[{QSFcmdline, LOG},"LC"]},
   COLORS[[If[c===None,Mod[l,Length[COLORS]-1,1],Key[c]]]]
];

LogIndent[op:OptionsPattern[{QSFcmdline, LOG}]]:=
   With[{l=OptionValue[{QSFcmdline, LOG},"LL"]},
   Switch[l,0,"",1,"+",_,StringTake["|-------------------",Mod[l,20]]]];
   
Options[LOG]={"LL":>LL[], "LC"->None};
LOG[msgs__,op:OptionsPattern[{QSFcmdline, LOG}] ]:=
CPrint[
   LogIndent[op],
   LogColor[op],

   StringJoin[ToStringList[If[ListQ[#],ToString[#,InputForm],Identity[#]]&/@{msgs}]], 
   COLORS["reset"] 
];

LOGF[f_, args_, rhs_] := Module[{res},   
   LOG["[", f,"]"];
   LL[]++;
   res = rhs;
   LL[]--;
   res
];

(* ClearAll[decorator]; *)
decorator /: SetDelayed[decorator[dec_][f_[args___]], rhs_]:=
  f[a:PatternSequence[args]]:=
     dec[Unevaluated @ f, Unevaluated @ {a}, Unevaluated @ rhs];



LOGI:=LOG[##,"LC"->"blue","LL"->0] &;
LOGJ:=LOG[##,"LC"->"magenta","LL"->0] &;
LOGW:=LOG[##,"LC"->"yellow","LL"->0] &;
LOGS:=LOG[##,"LC"->"green","LL"->0] &;
LOGE:=(Print[##] &; Pause[0.2] &;Abort[] &;);
LOGQ:=LOG[##,"LC"->"black","LL"->0] &;
SetAttributes[LOGL,{Listable}];
LOGL[msgs__,op:OptionsPattern[{QSFcmdline, LOGL}]]:=LOG[msgs,op];
SetAttributes[CHA, HoldFirst];
CHA[exp_, msgs__]:=Quiet@Check[exp, LOGE[msgs]];


Attributes[LOGV] = {HoldAll, SequenceHold,Listable};
LOGV[msgs__] := LOG[VARVALUE /@ Unevaluated[msgs]];
LOGV[msgs_] := LOG[VARVALUE @ Unevaluated[msgs]];

End[];
EndPackage[];