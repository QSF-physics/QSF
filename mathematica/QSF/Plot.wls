#!/usr/bin/env wolframscript
(* Plot probability distributions from input files *)

(* You can use this script in two ways:
1) Run this script from a parent directory (relative to input files) 
in which you'd like your files to be located and pass string "path/to/XXXX.psi" as first arg. 
2) Run from a different directory tree - this will locate the output files next to the input files.

Both options support "*" wildcards for finding multiple files at once in the input pathnames.
Important!: Place your input paths in doublequotes "" *)

SetDirectory[DirectoryName[$InputFileName]];
<< "common/parsing.wl"
<< "common/wf.wl"
ResetDirectory[];
Needs["ForScience`"];
Print["\nPlot"];
ParsePattern[];

ExportMultiple[path_,hd_,data_]:=Block[{tmpName},
    tmpName=PrependName[path, "blur_radius_"<>ToString[blur]<>"/qm_none_"];
    Export[tmpName, 
    WFPlot[hd, data, BlurRadius->0.0, "QuarterMerge"->None, PlotRange->{{0,4},{0,4}, Full}] ];
];

MapThread[(
    {hd,data}=WFLoad[#1, Probability->True]; 
    outF=PrependName[#2, "blur_radius_0/qm_"<>ToString[]<>"_"];
    Print[#1 <> " --> " #2];
    
    Echo[WFNorm[hd,data], "Norm: "];
    Export[PrependName[#2, "blur_radius_0/qm_none_"], WFPlot[hd, data, BlurRadius->0.0, "QuarterMerge"->None, PlotRange->{{0,4},{0,4}, Full}] ];
    Export[PrependName[#2, "blur_radius_0/qm_corr_"], WFPlot[hd, data, BlurRadius->0.0, "QuarterMerge"->"Correlated", PlotRange->{{0,4},{0,4}, Full}] ];
    Export[PrependName[#2, "blur_radius_0/qm_all_"], WFPlot[hd, data, BlurRadius->0.0, "QuarterMerge"->All, PlotRange->{{0,4},{0,4}, Full}] ];
    Export[PrependName[#2, "blur_radius_0.1/qm_none_"], WFPlot[hd, data, BlurRadius->0.1, "QuarterMerge"->None, PlotRange->{{0,4},{0,4}, Full}] ];
    Export[PrependName[#2, "blur_radius_0.1/qm_corr_"], WFPlot[hd, data, BlurRadius->0.1, "QuarterMerge"->"Correlated", PlotRange->{{0,4},{0,4}, Full}] ];
    Export[PrependName[#2, "blur_radius_0.1/qm_all_"], WFPlot[hd, data, BlurRadius->0.1, "QuarterMerge"->All, PlotRange->{{0,4},{0,4}, Full}] ];
    ) &, {inPaths,outPaths}];