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
    TrimMarginsPercent->0.2,
    0 ->TrimMargins@*Average@*RemoveBoundedPart@*WFLoad,StructPeak,
    {
        FrameLabel->{"X","Y"}, GridLines->Full, 
        Axes->False, PlotRange->{{0,2},Full}, PrettyPlots, 
        2->Combine@*GaussianBlur@*TransverseDiagSum@*MergeOrthants,
        StructPeak,
        3->Multicolumn,
        StructPeak, WFExport
    }
        (* {GridLines->Full, Axes->False, PlotRange->Full, 0->WFPlot,StructPeak, WFExport} *)
}];

(* {2 ->((Print[Head[#]])&/@WFLoad[#]&),1->PrintHead}]; *)
(* WFMultiExport[processed]; *)