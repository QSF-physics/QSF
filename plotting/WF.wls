#!/usr/bin/env wolframscript
(* Plot merged probability distribution from many input files *)

(*Run this script from a parent directory (relative to input files) 
in which you'd like your files to be located and pass string "path/to/XXXX.psi" as first arg. 
Important!: Place your input paths in doublequotes "" *)

SetDirectory[DirectoryName[$InputFileName]]
<< "common/parsing.wls";
<< "common/wf.wls";
ResetDirectory[];
Needs["ForScience`"]
LOG[yellow, "\nGroupWF"];

filterFun:=Select[#, StringEndsQ[#, "1"|"2"|"3"] &] &;
ParseInput["*.psi*", filterFun];
subdir="wf"; 

outputLocation="blur_radius_"<>ToString[blur]<>"/"<>subdir<>"/" <> groupOutput;


ProcessAvgs[avgs_]:=Block[{all,data,hd},
    (* Loading *)
    all=Map[WFLoad[#,Probability->False]&, avgs];
    (* hd=First[First[all]]; *)
    (* data=Total[all[[All,2]] ]; *)
    data=Total@Map[Last, all];
    (* data=Merge[all[[All,2]], Total]; *)
    hd=First[First[all]];
    (* {hd, data} *)

    WFPlot[hd, Abs[data]^2, BlurRadius->blur, QuarterMerge->All, PlotRange->{{0,4},{0,4}, Full}]
]

ProcessJoins[joins_]:=(
    LOG[blue,"Processing JOINs..."]; 
    First[Map[ProcessAvgs, joins]]
);

processed=Map[ProcessJoins, groupInput];

LOG[yellow, "Exporting the final image..."];
Export[
    outputLocation, 
    Multicolumn[
        KeyValueMap[
            Labeled[#2,FontFix/@{#1},{Top} ]&,
        processed ]//Transpose ,
    3, Appearance -> "Horizontal"] 
];
LOG[green, "Output location: ", outputLocation];

(* 
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
Echo[outputLocation, "Output location:"]; *)