#!/usr/bin/env wolframscript
(* Plot FLUX information*)

(* 1) Run this script from a parent directory (relative to input files) 
in which you'd like your files to be located and pass string "path/to/XXXX.psi" as first arg.  *)

(* Important!: Place your input paths in doublequotes "" *) 

SetDirectory[DirectoryName[$InputFileName] ];
<< "common/parsing.wl";
<< "common/flx.wl";
ResetDirectory[];

BeginPackage["QSF`"];
If[scriptMode,ParsePattern[] ];
LOGW["\nFLX"];

(* filterFun:=Select[#, StringEndsQ[#, "1"|"2"|"3"] &] &; *)
outputLocation="flx/_"<>options["groupOutput"];
(* For filtering unnecessary data *)
labelFilt:=StringContainsQ[#,(StartOfString~~"A"~~EndOfString)|(StartOfString~~"F"~~EndOfString)|"N2S"|"N2D"|"S2D"]&;

ProcessJoins[joins_]:=(
    LOGI["Processing JOINs..."]; 
    FLXJoin[Map[ProcessAvgs, joins], ImageSize->500, LabelTransform->KeyCleanup]
);

processed=Map[ProcessJoins, groupInput];

LOGW["Exporting the final image..."];
Export[
    outputLocation, 
    Multicolumn[
        KeyValueMap[
            Labeled[#2,
            FontFix@{
                "time [cycles]",#1
                (* Style[Rotate[#1, 90 Degree ],FontSize -> 20] *)
            },
            {Bottom,Top} ]&,
        processed ]//Transpose ,
    3, Appearance -> "Horizontal"] 
];
LOGS["Output location: ", outputLocation];

EndPackage[];