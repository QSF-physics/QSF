#!/usr/bin/env wolframscript
<<"cmdline`opt`";
<<"filestruct`";
Needs["QSF`"];
Needs["QSF`wf`"];
Needs["ForScience`"];
GetOpts[];
fileStruct=ParsePattern[$ScriptCommandLine[[2]] ];

Level 

(* filterFun:=Select[#, StringEndsQ[#, "1"|"2"|"3"] &] &; *)
processed=StructProcess[fileStruct, "Operations"->
{2 ->Average@*RemoveBoundedPart@*WFLoad,
    {0->WFPlot@*MergeOrthants,StructPeak, WFExport}
}];

(* {2 ->((Print[Head[#]])&/@WFLoad[#]&),1->PrintHead}]; *)
(* WFMultiExport[processed]; *)