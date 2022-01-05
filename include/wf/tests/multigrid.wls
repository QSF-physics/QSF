#!/usr/bin/env wolframscript 
(* ::Package:: *)
(* -format PNG *)

(* Usage: pass path to the directory containing region *.csv files *)
colors={White,White,Yellow,Orange, Red};
w2={-0.1,0.0,1/6,1/2,1};
colorfn3=Evaluate[Blend[Transpose[{w2,colors}],#]]&;
Pl3D:=ListDensityPlot3D[#,PlotLegends->Automatic,ImageSize->Medium,ColorFunction->colorfn3,ColorFunctionScaling->False,OpacityFunction->None,AxesLabel->{"X","Y","Z"},InterpolationOrder->1]&;
Pl2D:=ListDensityPlot[#,PlotLegends->Automatic,ImageSize->Medium,ColorFunction->colorfn3,ColorFunctionScaling->False,FrameLabel->{"X","Y"},InterpolationOrder->0]&;
Pl:=If[#1==3,Pl3D[#2],Pl2D[#2]]&;
Echo[$ScriptCommandLine[[2]], "args"];

DebugPlot[basename_]:=Block[{fs, fsd, width, output},
    fs=SortBy[FileNames[FileNameJoin[{$ScriptCommandLine[[2]], basename<>"*.csv"}]], ToExpression[StringExtract[#,"rID"->2,"_"->1]]&];
    Echo[FileBaseName[#],"Found matching mask files:"]&/@fs;
    fs=GatherBy[fs,ToExpression@StringExtract[#,basename->2,"_"->1]&];
    Echo[Map[Length,fs], "Number of files in each group"];

    fsd=Map[Import, fs, {2}];
    fsd=Map[Flatten[#/.{{}}->Nothing,1] &, fsd];
    Echo[Dimensions/@fsd,"Dims by group:"];
    width=Last[Dimensions[fsd] ];
    output=FileNameJoin[{$ScriptCommandLine[[2]], basename<>".png"}];
    Echo[output,"Processing... output will be saved in:"];
    {output, width, fsd}
];

{output, width, fsd} = DebugPlot["multigrid_mask_region"];
Export[output, 
    Multicolumn[{Labeled[Multicolumn[Pl[width,#[[All,;;width-1]] ]& /@fsd,
    3,Appearance->"Horizontal"],"mask"]},1,Appearance->"Vertical"] ];

{output, width, fsd} = DebugPlot["multigrid_slice_mask_region"];
Export[output, 
    Multicolumn[{
        Labeled[Multicolumn[Pl[width-2,#[[All,;;width-1]] ]&/@fsd,3,Appearance->"Horizontal"],"raw_mask"],
        Labeled[Multicolumn[Pl[width-2,#[[All,Append[Range[width-2],width]] ] ]&/@fsd,3,Appearance->"Horizontal"],"correction"],
        Labeled[Multicolumn[Pl[width-2,Transpose[Append[Transpose[#[[All,;;width-2]] ] ,#[[All,width-1]]*#[[All,width]] ] ] ]&/@fsd,3,Appearance->"Horizontal"],"final_mask"]
        },1,Appearance->"Vertical"]
    ];
