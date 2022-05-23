
BeginPackage["QSF`",{"cmdline`log`","cmdline`opt`", "QSF`DataAnalysis`", "QSF`wf`", "QSF`flx`"}];
GaussianBlur;
(* Average; *)
Expo;
ExportableQ;
(* GaussianBlur2; *)

Begin["`Private`"];

ExportableQ[_Legended|_Graphics|_Grid]:=True;
ExportableQ[_]:=False;

Options[Expo]={"ExportPath"->"plots/", "TreePath"->{}, "LeafPath"->{},"FileFormat"->".png"};
decorator[LOGF]@
Expo[gr_?ExportableQ,opt:OptionsPattern[]]:=
(LOG@Export[
	StringJoin[{
		FileNameJoin[Flatten[{
			StringJoin[Flatten[
				{
					LOG["TreePath: ",COptionValue[Expo,"TreePath"]];
					LOG["LeafPath: ",COptionValue[{opt,Expo},"LeafPath"]];
					COptionValue[{opt,Expo},"ExportPath"] 
					,COptionValue[{opt,Expo},"TreePath"]
				}]] 
				,COptionValue[{opt,Expo},"LeafPath"]
		}]],
		StringPadLeft[COptionValue[Expo,"FileFormat"],4, "."]
	}]
,gr];gr)

Expo[ass_Association,opt:OptionsPattern[]]:=MapIndexed[Expo[#1,"LeafPath"->Flatten@Join[COptionValue[{opt,Expo},"LeafPath"],#2/.{Key[x_]->x}]]&, ass];

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
