
BeginPackage["QSF`flx`", {"cmdline`log`","cmdline`opt`", "QSF`StyleUtils`","QSF`DataAnalysis`","PlotGrid`", "QSF`Physics`"}];
LaserDataQ
FluxDataQ;
LatexDataLabels;
EnrichMetadata;

FLX;
FLXLoad;
(* Average; *)
DataPlots;
DataMerge;
DataPlotGrid;
DataFileNameQ;
DataMap;
GlobalOpt;

Begin["`Private`"];

SetAttributes[LookInto,Listable];
LookInto[FLX[hd_Association, data_Association],x_Symbol|x_Function ,opt:OptionsPattern[]]:=With[{act=COptionValue[{opt},"ActOn"], rl=Sequence@@FilterRules[{opt,"DataStep"->QuantityMagnitude@hd["dt"]}, Options@x]},
FLX[hd,
  If[act===System`All,
    Map[x[#,rl]&,data],
    MapAt[x[#,rl]&,data,
      Switch[act,
        _List,Intersection[act, Keys@data], 
        _String,Intersection[{act}, Keys@data],
        _,Map[List,Select[Keys@data,act]] (* function *)
        (* Select[Keys@data,COptionValue[{opt,DataMap},"ActOn"] ]) *)
      ]
    ]
  ]]
];

Options[DataMap]={"ActOn"->System`All};
DataMap[x_Symbol|x_Function,Optional[l_String|l_List|l_Symbol|l_Function,COptionValue[{DataMap},"ActOn"]]][ass_Association]:=LookInto[ass,x,"ActOn"->l];

Options[GlobalOpt]={"ActOn"->System`All};
GlobalOpt[op_String, x_Symbol|x_Function,Optional[l_String|l_List|l_Symbol|l_Function,COptionValue[{GlobalOpt},"ActOn"]]][ass_Association]:=With[
  {opv=x[Flatten[Cases[ass,FLX[_,data_]:>Values[KeySelect[data,l]],\[Infinity]]]]},
  UpdateOpts[op->opv];
  LOGV[op, opv];
  ass
];


Options[DataThead]={"ThreadOn"->System`All};
DataMerge[x_Symbol|x_Function,Optional[l_String|l_List|l_Symbol|l_Function,COptionValue[{DataThead},"ThreadOn"]]][ass_Association]:=
FLX[
  (* Print[Merge[Cases[ass,FLX[hd_,data_]:>KeySelect[data,l],\[Infinity]],First]]; *)
  Merge[Cases[ass,FLX[hd_,data_]:>hd,\[Infinity]],First],
  Join[
    Merge[Cases[ass,FLX[_,data_]:>KeySelect[data,Not@*l],\[Infinity]],First],
    Merge[Cases[ass,FLX[_,data_]:>KeySelect[data,l],\[Infinity]],x]
  ]
];

(* Average /: SetDelayed[Average[FLX][ass_Association], rhs_]:=
  f[a:PatternSequence[args]]:=
     dec[Unevaluated @ f, Unevaluated @ {a}, Unevaluated @ rhs]; *)
(* For extracting wavelength data from path *)


DataFileNameQ:=StringEndsQ[#,".dat"]&;

fluxPattern=(StartOfString~~ (LetterCharacter|PunctuationCharacter)~~"2"~~ (LetterCharacter|PunctuationCharacter|"_")..~~EndOfString);
laserPattern=((StartOfString~~"A")|(StartOfString~~"F"))~~(("_"|DigitCharacter)...)~~EndOfString;

FluxDataQ:=StringContainsQ[#,fluxPattern]&;
LaserDataQ:=StringContainsQ[#,laserPattern]&;

LatexDataLabels:=Which[
  FluxDataQ[#], 
  StringReplace[#,{"*"->"","2"->" \\rightarrow ","ASYM"->"{-}","SYM"->"{+}"}]<>"\,[au]",
  LaserDataQ[#],
  StringJoin[#,"(t)\,[au]"]
]&;

(* TODO: placeholder - should be loaded from binary header *)
FLXHeader[st_,opt:OptionsPattern[]]:=<|"labels"-><|
  "step"->1,"time"->2,"A"->3,"norm"->4,"eta"->5,
  "N2S"->6,"N2D_SYM"->7,"N2D_ASYM"->8,"S2D_SYM"->9,
  "S2D_ASYM"->10,"S2CAP"->11,"D2CAP"->12|>,
  "TimeUnit"->"[au]","DataUnit"->""
|>;
FLXHeader2[st_,opt:OptionsPattern[]]:=<|"labels"-><|
  "step"->1,"time"->2,"A"->3,"eta"->4,"E_kin"->5,"E_pot"->6,
  "N2S"->7,"N2D_SYM"->8,"N2D_ASYM"->9,"S2D_SYM"->10,
  "S2D_ASYM"->11,"S2CAP"->12,"D2CAP"->13|>,
  "TimeUnit"->"[au]","DataUnit"->""
|>;
(* proper for flux_quiver_ratio_1.0000 *)

(* ArtificialLabelQ:=StringContainsQ[#,"*"]&;
NonArtificialLabels[hd_]:=KeySelect[hd["labels"],Not@*ArtificialLabelQ]; *)
FLXDerivedData[data_]:={
  "N2*"->(data["N2S"]+data["N2D_SYM"]+data["N2D_ASYM"]),
  "*2S"->(data["N2S"]-data["S2D_SYM"]-data["S2D_ASYM"]),
  "*2D"->(data["N2D_SYM"]+data["N2D_ASYM"]+data["S2D_SYM"]+data["S2D_ASYM"]),
  "*2D_SYM"->(data["N2D_SYM"]+data["S2D_SYM"]),
  "*2D_ASYM"->(data["N2D_ASYM"]+data["S2D_ASYM"]),
  "*2CAP"->(data["S2CAP"]+data["D2CAP"])
  };

FLXData[st_,hd_,opt:OptionsPattern[]]:=Module[{data},
  data=Partition[BinaryReadList[st,"Real64"],Length[hd["labels"]]];
  data=Transpose[data];
  data=Map[Part[data,#]&,hd["labels"]];
  CHA[AssociateTo[data,COptionValue[{opt,FLXLoad},"DerivedData"][data]],"Cannot derived data"];
  (* labeled=KeySelect[labeled,!StringContainsQ[#,"CAP"]&]; *)
  LOG["Imported ",Length@Keys@data ," keys: ", Keys[data]];
  data
];

Options[FLXLoad]={EnrichMetadata->{},"DerivedData"->FLXDerivedData};

SetAttributes[FLXLoad, Listable];
decorator[LOGF]@
FLXLoad[path_,opt:OptionsPattern[]]:=Block[{st,hd,data},
  st=OpenBin[path];
  hd=If[StringContainsQ[path,"flux_quiver_ratio_2.000000"],FLXHeader2[st],FLXHeader[st]];
  hd["path"]=path;  
  data=FLXData[st,hd,opt];
  Close[st];
  AssociateTo[hd,
  {
    "tmax"->Quantity[data["time"][[-1]],"AtomicUnitOfTime"],
    "dt"->Quantity[data["time"][[2]],"AtomicUnitOfTime"]}];
  LOG["dt:", hd["dt"], " tmax:", hd["tmax"]];
  AssociateTo[hd,COptionValue[{opt,FLXLoad},EnrichMetadata][hd, data]];
  FLX[hd, data]
];


(* decorator[LOGF]@ *)
(* Average[matchF_][inp_, OptionsPattern[]]:=FLX[
    Merge[Cases[inp,FLX[hd_,data_]:>hd,\[Infinity]],First],
		Join[
      Merge[Cases[inp,FLX[_,data_]:>KeySelect[data,Not@*matchF],\[Infinity]],First],
      Merge[Cases[inp,FLX[_,data_]:>KeySelect[data,matchF],\[Infinity]],Mean]
    ]
]; *)

fg[min_,max_]:=DeleteDuplicates[Flatten[Table[
  Table[{i,Directive[Thickness[0.00025 lz],Dashed,GrayLevel[(5.0-lz)/5.0]]},
  {i,Ceiling[min,10^lz],Floor[max,10^lz],10^lz}],{lz,5,1,-1}],1]
  ,First[#1]==First[#2]&];

Options[DataPlots]={
  System`DataRange->{0,"tmax"}, 
  System`PlotRange->{System`Automatic,System`Full}, 
  FrameLabel->{
  {HoldForm["DataLabel""DataUnit"],HoldForm["DataLabel""DataUnit"]},{HoldForm["time [au]"],None}},
  "LegendPlacement"->Bottom,
  "LeafPath"->{}
};
DataPlots[props_List][FLX[hd_Association,data_Association],opt:OptionsPattern[]]:=
  KeyValueMap[
    ListLinePlot[
        (* Print[COptionValue[{Options[QSFcmdline],opt,DataPlots},GridLines]]; *)
      Legended[Placed[#2,COptionValue[{opt,DataPlots},"LegendPlacement"]],COptionValue[{opt,DataPlots},"LeafPath"]]
      ,GridLines->{Automatic, fg[-10000,10000]}
      ,Frame->True
      ,PlotStyle->GetColor[COptionValue[{opt,DataPlots},"LeafPath"]]
      (* ,GridLines->{Automatic,fg} *)
      ,ReleaseHold@FilterRules[
        {Options[QSFcmdline],opt,Options[DataPlot]}/.Rule["DataLabel",#1]/.hd, 
        Options[ListLinePlot]]
    ]&
    ,KeyTake[data,props]
  ];

(* decorator[LOGF]@ *)
DataPlots[props_String][ass_Association,opt:OptionsPattern[]]:=
DataPlots[{props}][ass,opt];

DataPlots[props_List][ass_Association,opt:OptionsPattern[]]:=
MapThread[
  (* Show[##,AbsoluteOptions[#1,PlotRange]]/.FixNestedLegends &, *)
  Show[##
    ,ReplacePart[First@AbsoluteOptions[#1,PlotRange],{-1,-1}->Full]
  ]/.FixNestedLegends &,
  Cases[
		MapIndexed[
			If[MatchQ[#1,_FLX],DataPlots[props][#1,"LeafPath" -> #2],#1]&
			,ass
			,ArrayDepth[ass, AllowedHeads -> Association]
		]
    , _List, ArrayDepth[ass, AllowedHeads -> Association]
	]
];

Options[DataPlotGrid]={"GridLabels"->{},"GridTranspose"->False};
DataPlotGrid[ass_Association,opt:OptionsPattern[]]:=DataPlotGrid[Transpose@GriddedLeaves[ass],"GridLabels"->GridKeys[ass]];
DataPlotGrid[x_List,opt:OptionsPattern[]]:=If[MatrixQ[x]
  ,Legended[PlotGrid1[RemoveLegends[x],opt],UnifyLegends[x]]
  (* ,Grid[{x},Spacings -> Scaled[-0.04] *)
  ,Legended[PlotGrid1[RemoveLegends[Transpose[{x}]],opt],UnifyLegends[Transpose[{x}]]]
];

(* Too general :-() *)
(* FLX /:x_[FLX[hd_Association,data_Association],c___,opt:OptionsPattern[]]:=
FLX[hd,
  If[COptionValue[{opt,x},"ActOn"]===All,
    Map[x[#,c]&,data], 
    MapAt[x[#,c]&,data,
      If[ListQ[COptionValue[{opt,x},"ActOn"]],
        Intersection[COptionValue[{opt, x}, "ActOn"], Keys@data], 
        KeySelect[data, COptionValue[{opt, x}, "ActOn"]]
      ]
    ]
  ]
]; *)
End[];
EndPackage[];
