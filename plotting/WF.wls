#!/usr/bin/env wolframscript
(* Plot single, merged probability distribution from many input files *)

(* You can use this script in two ways:
1) Run this script from a parent directory (relative to input files) 
in which you'd like your files to be located and pass string "path/to/XXXX.psi" as first arg. 
2) Run from a different directory tree - this will locate the output files next to the input files.

Both options support "*" wildcards for finding multiple files at once in the input pathnames.
Important!: Place your input paths in doublequotes "" *)

SetDirectory[DirectoryName[$InputFileName]]
<< "common/parsing.wls";
<< "common/wf.wls";
ResetDirectory[];
Needs["ForScience`"]
Print["\nGroupPlot"]
filterFun:=Select[#, StringEndsQ[#, "1"|"2"|"3"] &] &;
ParseInput["*.psi*", filterFun];

MergeAndBlur[paths_]:=Block[{all,hd,data},
    all=(Echo[#, "Processing: "]; WFLoad[#,Probability->False]) & /@ paths;
    hd=First[First[all]];
    data=Total[all[[All,2]] ];
    Print["Making merged plot..."];
    WFPlot[hd, Abs[data]^2, BlurRadius->blur, QuarterMerge->All, PlotRange->{{0,4},{0,4}, Full}]
];
outputLocation="blur_radius_"<>ToString[blur]<>"/all_" <> groupOutput;
processed=Map[MergeAndBlur, groupInput];
processed=KeySortBy[ToExpression[Last[#]]&][processed];
(* groupBins=DeleteDuplicates/@Transpose[Most/@Keys[processed] ]; *)
Echo[groupCount, "Exporing the final image..., groupCount"];
Export[
    outputLocation, 
    Multicolumn[KeyValueMap[Labeled[#2,StringRiffle[Most[#1],"/"] ]&, processed ],groupCount] 
];
Echo[outputLocation, "Output location:"];