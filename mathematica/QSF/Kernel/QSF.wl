
BeginPackage["QSF`",{"cmdline`opt`", "QSF`DataAnalysis`", "QSF`wf`", "QSF`flx`"}];
GaussianBlur;
Average;
(* GaussianBlur2; *)

Begin["`Private`"];

(* Average[x_][ass_Association]:=Mean@Cases[ass,_x,\[Infinity]]; *)
(* FLX /:Downsample[FLX[hd_Association,data_Association],x_:2]:=
Map[Downsample *)
(* Options[GaussianBlur2]=Join[Options[GaussianFilterOpt]];
GaussianBlur2[data_List,opt:OptionsPattern[]]:=GaussianFilterOpt[data, opt]; *)



(* GaussianBlur[props_List][FLX[hd_Association, data_Association],opt:OptionsPattern[]]:=
FLX[hd,MapAt[GaussianFilterOpt[#, hd["T"]/hd["dt"],opt]&, data, props]]; *)

(* GaussianBlur[matchF_][FLX[hd_Association, data_Association],opt:OptionsPattern[]]:=
FLX[hd,MapIndexed[If[matchF[#2],GaussianFilterOpt[#1, hd["T"]/hd["dt"],opt],data]&, data]]; *)

SetAttributes[GaussianBlur,{Listable}];
Options[GaussianBlur]=Options[GaussianFilterOpt];
GaussianBlur[WF[hd_Association, data_List],opt:OptionsPattern[]]:=
WF[hd,GaussianFilterOpt[data,"DataStep"->WFDataStep[hd],opt]];

End[];

EndPackage[];
