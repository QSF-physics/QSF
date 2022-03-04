#!/usr/bin/env wolframscript
<<"cmdline`opt`";
<<"filestruct`";
Needs["QSF`"];
Needs["QSF`wf`"];
Needs["ForScience`"];
GetOpts[];
fileStruct=ParsePattern[$ScriptCommandLine[[2]] ];


(* filterFun:=Select[#, StringEndsQ[#, "1"|"2"|"3"] &] &; *)
processed=StructProcess[fileStruct, "Operations"->
{
    2 ->Average@*RemoveBoundedPart@*WFLoad,
        {GridLines->Full, Axes->False, PlotRange->Full, 0->Combine@*TransverseDiagSum@*MergeOrthants,StructPeak, WFExport}
}];

(* {2 ->((Print[Head[#]])&/@WFLoad[#]&),1->PrintHead}]; *)
(* WFMultiExport[processed]; *)