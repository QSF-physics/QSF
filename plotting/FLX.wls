#!/usr/bin/env wolframscript -noprompt 
(* Plot FLUX information*)

(* 1) Run this script from a parent directory (relative to input files) 
in which you'd like your files to be located and pass string "path/to/XXXX.psi" as first arg.  *)

(* Important!: Place your input paths in doublequotes "" *) 
SetDirectory[DirectoryName[$InputFileName]];
<< "common/parsing.wls";
<< "common/flx.wls";
ResetDirectory[];
LOG[yellow, "\nGroupFLX"];

(* filterFun:=Select[#, StringEndsQ[#, "1"|"2"|"3"] &] &; *)
ParseInput["RE.dat"];
(* Whether to integrate the fluxes *)
integrateFlux=Quiet[ToExpression@Check[$ScriptCommandLine[[5]],False] ];
subdir="flx"; 
subdir=If[integrateFlux,"i",""]<>subdir;

outputLocation="blur_radius_"<>ToString[blur]<>"/"<>subdir<>"/" <> groupOutput;
(* For filtering unnecessary data *)
labelFilt:=StringContainsQ[#,"A"|"F"|"N2S"|"N2D"|"S2D"]&;
(* For extracting wavelength data *)
lambdaExtract:=ToExpression[StringExtract[#,"nm_"->2,"/"->1] ]&;


ProcessAvgs[avgs_]:=Block[{lambda,all,data,meta},
    (* Loading *)
    all=Map[(
        lambda=lambdaExtract[#]; 
        
        FLXLoad[#, labelFilt, lambda]
        ) &, avgs];
    (* TODO: check that the labels are identical *)
    (* TODO: check that the metadata are identical or at least similar*)
    data=Merge[all[[All,1]], Mean];
    meta=Last[Last[all]];
    (* Smooth the data before std dev gets calculated *)
    Map[
        (data[#]=
        (If[integrateFlux, meta["dt"] Accumulate[#], Identity[#] ] &)@
        GaussianFilter[data[#],blur meta["T"]/meta["dt"] ];) &, 
        Select[Keys[data],FLXFilter]
    ];
    {data, meta}
]

ProcessJoins[joins_]:=(LOG[blue,"Processing JOINs..."]; 



FLXJoin[Map[ProcessAvgs, joins], ImageSize->500, LabelTransform->KeyCleanup]);

processed=Map[ProcessJoins, groupInput];

LOG[yellow, "Exporting the final image..."];
Export[
    outputLocation, 
    Multicolumn[
        KeyValueMap[
            Labeled[#2,
            FontFix/@{
                "time [cycles]",#1
                (* Style[Rotate[#1, 90 Degree ],FontSize -> 20] *)
            },
            {Bottom,Top} ]&,
        processed ]//Transpose ,
    3, Appearance -> "Horizontal"] 
];
LOG[green, "Output location: ", outputLocation];