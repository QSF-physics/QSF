#!/usr/bin/env wolframscript 
(* ::Package:: *)
(* -format PNG *)
colors={White,White,Yellow,Orange, Red};
w2={-0.1,0.0,1/6,1/2,1};
colorfn3=Evaluate[Blend[Transpose[{w2,colors}],#]]&;
Pl3D:=ListDensityPlot3D[#,PlotLegends->Automatic,ImageSize->Medium,ColorFunction->colorfn3,ColorFunctionScaling->False,OpacityFunction->None,AxesLabel->{"X","Y","Z"},InterpolationOrder->1]&;
Pl2D:=ListDensityPlot[#,PlotLegends->Automatic,ImageSize->Medium,ColorFunction->colorfn3,ColorFunctionScaling->False,FrameLabel->{"X","Y"},InterpolationOrder->0]&;
Pl:=If[#1==5,Pl3D[#2],Pl2D[#2]]&;
Echo[$ScriptCommandLine[[2]], "args"];
fs=SortBy[FileNames[FileNameJoin[{$ScriptCommandLine[[2]], "/multigrid_mask_region*_rID*.csv"}]], ToExpression[StringExtract[#,"rID"->2,"_"->1]]&];
Echo[FileBaseName[#],"Found matching file:"]&/@fs;
fs=GatherBy[fs,ToExpression@StringExtract[#,"multigrid_mask_region"->2,"_"->1]&];
Echo[Map[Length,fs], "Number of files in each group"];
fsd=Map[Import, fs, {2}];
fsd=Map[Flatten[#/.{{}}->Nothing,1] &, fsd];
Echo[Dimensions/@fsd,"Dims by group:"];
dim=Last[Dimensions[fsd]];
output=FileNameJoin[{$ScriptCommandLine[[2]], "multigrid_mask.pdf"}];
Echo[output,"Processing... output will be saved in:"]
Export[output, 
Multicolumn[{
    Labeled[Multicolumn[Pl[dim,#[[All,;;dim-1]]]&/@fsd,3,Appearance->"Horizontal"],"raw_mask"],
    Labeled[Multicolumn[Pl[dim,#[[All,Append[Range[dim-2],dim]]]]&/@fsd,3,Appearance->"Horizontal"],"correction"],
    Labeled[Multicolumn[Pl[dim,Transpose[Append[Transpose[#[[All,;;dim-2]]] ,#[[All,dim-1]]*#[[All,dim]] ] ] ]&/@fsd,3,Appearance->"Horizontal"],"final_mask"]
    },1,Appearance->"Vertical"]
];

        (* Multicolumn[Pl2D[#[[All,;;3]]]&/@fsd,2,Appearance->"Horizontal"],
        Multicolumn[Pl2D[#[[All,{1,2,4}]]]&/@fsd,2,Appearance->"Horizontal"],
        Multicolumn[Pl2D[Append[(#[[All,;;2]])\[Transpose] ,#[[All,3]]*#[[All,4]]]\[Transpose]]&/@fsd,2,Appearance->"Horizontal"] *)