#!/usr/bin/env wolframscript
(* Plot single, merged probability distribution from many input files *)

(* You can use this script in two ways:
1) Run this script from a parent directory (relative to input files) 
in which you'd like your files to be located and pass string "path/to/XXXX.psi" as first arg. 
2) Run from a different directory tree - this will locate the output files next to the input files.

Both options support "*" wildcards for finding multiple files at once in the input pathnames.
Important!: Place your input paths in doublequotes "" *)

SetDirectory[DirectoryName[$InputFileName]]
<< "common/wf.wl"
ResetDirectory[];
Needs["ForScience`"];

Print["\nMomenta Ionized"];
filterFun[paths_]:=Block[{hd, res, bt}, 
   res=Select[paths, StringContainsQ[#, "repP.psi"|"repX.psi"] &];
   bt=Total[hd["bounded"]];
   Select[res, (hd=WFHeader[#]; (hd["rep"]==2 && bt==0) || (hd["rep"] ==1 && bt>0)) &]
];
ParsePattern[filterFun];
hd0=None;

ratio=Quiet[ToExpression@Check[$ScriptCommandLine[[5]],0.5] ];
Echo[blur, "blurRadius:"];
Echo[ratio, "ratio:"];

all=MapThread[
   (
      {hd, data} = WFLoad[#1]; 
      pscale = (Times @@ hd["dxs"]) Sqrt[(Times @@ Dimensions[data])/Power[(2 \[Pi]), hd["dim"] ] ];
   
      data=If[hd["rep"]==2, data, 
      hd["rep"]=2; Fourier[Transpose[data] WFShiftTrick[hd] WFInvMask[hd, ratio] ] pscale];
      data=Abs[data]^2;
      hd0=hd;
      tmpname=PrependName[#2, "blur_radius_"<>ToString[blur]<>"/ratio_"<>ToString[ratio]<>"_qm_none_"];
      Echo[tmpname, "Processing:"];
      Export[
         tmpname,
         WFPlot[hd, data, BlurRadius->blur, QuarterMerge->None, PlotRange->Full]
      ];
      data
   ) &, 
   {inPaths,outPaths}]; 

all=Total[all];
Echo[hd0,"Selected header:"];
outPath=PrependName[Last[outPathsMerge], "blur_radius_"<>ToString[blur]<>"/ratio_"<>ToString[ratio]<>"_qm_none_"];
Echo[outPath, "Output location: "];
Export[outPath, WFPlot[hd0, all, BlurRadius->blur, QuarterMerge->None, PlotRange->Full] ];

outPath=PrependName[Last[outPathsMerge], "blur_radius_"<>ToString[blur]<>"/ratio_"<>ToString[ratio]<>"_qm_all_"];
Echo[outPath, "Output location: "];
Export[outPath, WFPlot[hd0, all, BlurRadius->blur, QuarterMerge->All, PlotRange->{{0,4},{0,4}, Full}] ];
