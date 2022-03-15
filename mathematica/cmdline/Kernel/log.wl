#!/usr/bin/env wolframscript
BeginPackage["cmdline`log`"];
LOG;
LL;
LOGF;
LOGI;
LOGII;
LOGI;
LOGII;
LOGW;
LOGS;
LOGE;
LOGQ;
LOGL;
LOGV;
CHA;
decorator;


Begin["`Private`"];
CmdLineStr[x_]=If[$Notebooks, "",x];
CmdLinePrint=If[$Notebooks,PrintTemporary,Print];
COLORS = <|
   "green" -> CmdLineStr@"\033[1;32m", 
   "yellow" -> CmdLineStr@"\033[1;33m",
   "blue" -> CmdLineStr@"\033[1;34m", 
   "magenta" -> CmdLineStr@"\033[1;35m",
   "red" -> CmdLineStr@"\033[1;31m",
   "black" -> CmdLineStr@"\033[1;30m",
   "reset" -> CmdLineStr@"\033[1;0m"
|>;

template = StringTemplate["``[``````=``````]"];

SetAttributes[VARVALUE, HoldAll];
VARVALUE[msg_] := template[
   COLORS["reset"],
   COLORS["blue"],
   ToString@Unevaluated@msg,
   COLORS["reset"],
   COLORS["green"],
   ToString@msg, 
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
Options[LOG]={"LL":>LL[], "LC"->None};
LOG[msgs___,op:OptionsPattern[{QSFcmdline, LOG}] ]:=
CmdLinePrint[
   If[OptionValue[{QSFcmdline, LOG},"LL"]==0,"",
   If[OptionValue[{QSFcmdline, LOG},"LL"]==1,"+","|"] ], 
   StringRepeat["-", Mod[OptionValue[{QSFcmdline, LOG},"LL"],50] ],
   COLORS[[
      If[
         OptionValue[{QSFcmdline, LOG},"LC"]===None,
         Mod[OptionValue[{QSFcmdline, LOG},"LL"],Length[COLORS]-1,1],
         Key[OptionValue[{QSFcmdline, LOG},"LC"] ] 
      ]
   ]],
   msgs, 
   COLORS["reset"] 
];

LOGF[f_, args_, rhs_] := Module[{res},   
   LOG["[", ToString[f],"]"];
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
LOGE[msgs__]:=(Print[msgs]; Pause[0.2];Abort[];);
LOGQ[msgs__]:=LOG[msgs,"LC"->"black","LL"->0];
SetAttributes[LOGL,{Listable}];
LOGL[msgs__,op:OptionsPattern[{QSFcmdline, LOGL}] ]:=LOG[msgs,op];
SetAttributes[CHA, HoldFirst];
CHA[exp_, msgs__]:=Quiet@Check[exp, LOGE[msgs]];


Attributes[LOGV] = {HoldAll, SequenceHold,Listable};
LOGV[msgs__] := LOG[VARVALUE /@ Unevaluated[msgs]];
LOGV[msgs_] := LOG[VARVALUE @ Unevaluated[msgs]];

End[];
EndPackage[];