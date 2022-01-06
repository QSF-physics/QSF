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
<< "common/flx.wls";
ResetDirectory[];
Needs["ForScience`"]
Print["\nGroupFLX"]
(* filterFun:=Select[#, StringEndsQ[#, "1"|"2"|"3"] &] &; *)
ParseInput["RE.dat"];


MergeAndBlur[paths_]:=Block[{all,labels,data},
    (* pk=GroupBy[paths,StringExtract[#, "phase_"->2,"/"->1]&]; *)

    all=(Echo[#, "Processing: "]; FLXLoad[#]) & /@ paths;
    wavelength=ToExpression[StringExtract[First[paths],"nm_"->2,"/"->1] ];
    model=StringExtract[First[paths],"/"->1];
    Print[model];
    labels=First[First[all]];
    data=Total[all[[All,2]] ]/ Length[all[[All,2]]];
    Print["Making merged plot..."];
    FLXPlot[data, labels, wavelength, blur, If[model==="es", PlotStyle->Red,PlotStyle->Blue]]
];
outputLocation="blur_radius_"<>ToString[blur]<>"/color_" <> groupOutput;
processed=Map[MergeAndBlur, groupInput];
(* Print[PlotRange /. AbsoluteOptions[First[Values[processed] ]]]; *)
(* Export[outputLocation, First[Values[processed] ] ]; *)
(* processed=KeySortBy[ToExpression[Last[#]]&][processed]; *)
processed=KeySortBy[First[#]&][processed];

Echo[groupCount, "Exporing the final image..., groupCount"];
Export[
    outputLocation, 
    Multicolumn[KeyValueMap[Labeled[#2,Style[#,FontFamily -> "Latin Modern Math"]&/@{"time [cycles]",Style[Rotate[StringRiffle[Most[#1],"/"], 90 Degree ],FontSize -> 20]},{Bottom,Left} ]&, processed ]//Transpose ,3, Appearance -> "Horizontal"] 
];
Echo[outputLocation, "Output location:"];