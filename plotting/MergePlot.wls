#!/usr/bin/env wolframscript
(* Plot single, merged probability distribution from many input files *)

(* You can use this script in two ways:
1) Run this script from a parent directory (relative to input files) 
in which you'd like your files to be located and pass string "path/to/XXXX.psi" as first arg. 
2) Run from a different directory tree - this will locate the output files next to the input files.

Both options support "*" wildcards for finding multiple files at once in the input pathnames.
Important!: Place your input paths in doublequotes "" *)

SetDirectory[DirectoryName[$InputFileName]]
<< "common/wf.wls"
ResetDirectory[];
Needs["ForScience`"]
Print["MERGE PLOT"]
filterFun:=Select[#, StringEndsQ[#, "1"|"2"|"3"] &] &;
ParseInput[filterFun];
all=WFLoad /@ inPaths;
hd=First[Last[all]];
data=Total[all[[All,2]] ]/ Length[all];


outPath=Last[outPathsMerge];
Echo[outPath, "Output location: "];
Export[PrependName[outPath, "blur_radius_0/qm_none_"], WFPlot[hd, data, BlurRadius->0.0, QuarterMerge->None, PlotRange->{{0,4},{0,4}, Full}] ];
Export[PrependName[outPath, "blur_radius_0/qm_corr_"], WFPlot[hd, data, BlurRadius->0.0, QuarterMerge->Correlated, PlotRange->{{0,4},{0,4}, Full}] ];
Export[PrependName[outPath, "blur_radius_0/qm_all_"], WFPlot[hd, data, BlurRadius->0.0, QuarterMerge->All, PlotRange->{{0,4},{0,4}, Full}] ];
Export[PrependName[outPath, "blur_radius_0.1/qm_none_"], WFPlot[hd, data, BlurRadius->0.1, QuarterMerge->None, PlotRange->{{0,4},{0,4}, Full}] ];
Export[PrependName[outPath, "blur_radius_0.1/qm_corr_"], WFPlot[hd, data, BlurRadius->0.1, QuarterMerge->Correlated, PlotRange->{{0,4},{0,4}, Full}] ];
Export[PrependName[outPath, "blur_radius_0.1/qm_all_"], WFPlot[hd, data, BlurRadius->0.1, QuarterMerge->All, PlotRange->{{0,4},{0,4}, Full}] ];
