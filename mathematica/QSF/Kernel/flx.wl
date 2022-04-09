
BeginPackage["QSF`flx`", {"cmdline`log`","cmdline`opt`", "QSF`StyleUtils`","QSF`DataAnalysis`","PlotGrid`", "QSF`Physics`"}];
LaserDataQ
FluxDataQ;
LatexDataLabels;
EnrichMetadata;

FLX;
FLXLoad;
Average;
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

fluxPattern=(StartOfString~~ (LetterCharacter|PunctuationCharacter)~~"2"~~ (LetterCharacter|"_")..~~EndOfString);
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
(* proper for flux_quiver_ratio_1.0000 *)

(* ArtificialLabelQ:=StringContainsQ[#,"*"]&;
NonArtificialLabels[hd_]:=KeySelect[hd["labels"],Not@*ArtificialLabelQ]; *)
FLXDerivedData[data_]:={"*2S"->(data["N2S"]-data["S2D_SYM"]-data["S2D_ASYM"]),
"*2D"->(data["N2D_SYM"]+data["N2D_ASYM"]+data["S2D_SYM"]+data["S2D_ASYM"])};

FLXData[st_,hd_,opt:OptionsPattern[]]:=Module[{data},
  data=Partition[BinaryReadList[st,"Real64"],Length[hd["labels"]]];
  data=Transpose[data];
  data=Map[Part[data,#]&,hd["labels"]];
  CHA[AssociateTo[data,COptionValue[{opt,FLXLoad},"DerivedData"][data]],"Cannot derived data"];
  (* labeled=KeySelect[labeled,!StringContainsQ[#,"CAP"]&]; *)
  LOG["Imported ",Length@Keys@data ," keys: ", Keys[data]];
  data
];

Options[FLXLoad]={EnrichMetadata->{}, "DerivedData"->FLXDerivedData};

SetAttributes[FLXLoad, Listable];
decorator[LOGF]@
FLXLoad[path_,opt:OptionsPattern[]]:=Block[{st,hd, data},
  st=OpenBin[path];
  hd=FLXHeader[st];
  hd["path"]=path;  
  data=FLXData[st,hd,opt];
  Close[st];
  AssociateTo[hd,
  {"tmax"->Quantity[data["time"][[-1]],"AtomicUnitOfTime"],
  "dt"->Quantity[data["time"][[2]],"AtomicUnitOfTime"]}];
  AssociateTo[hd,COptionValue[{opt,FLXLoad},EnrichMetadata][hd, data]];
  FLX[hd, data]
];


(* decorator[LOGF]@ *)
Average[matchF_][inp_, OptionsPattern[]]:=FLX[
    Merge[Cases[inp,FLX[hd_,data_]:>hd,\[Infinity]],First],
		Join[
      Merge[Cases[inp,FLX[_,data_]:>KeySelect[data,Not@*matchF],\[Infinity]],First],
      Merge[Cases[inp,FLX[_,data_]:>KeySelect[data,matchF],\[Infinity]],Mean]
    ]
];

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
      Legended[Placed[#2,COptionValue[{opt,DataPlots},"LegendPlacement"]],COptionValue[{opt,DataPlots},"LeafPath"]]
      ,Frame->True
      ,PlotStyle->GetColor[COptionValue[{opt,DataPlots},"LeafPath"]]
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
  Show[##,ReplacePart[First@AbsoluteOptions[#1,PlotRange],{-1,-1}->Full]]/.FixNestedLegends &,
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
  ,Grid[{x},Spacings -> Scaled[-0.04]]
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

(* FLXJoin[all_, opt :OptionsPattern[] ]:=Block[{legend,size,colors},
    
    (* legend for data on the same plot and colors*)
    legend=Keys[all];
    size=Length[legend];
    colors=AssociationThread[legend -> ColorData[97, "ColorList"][[;;size]] ];
    
    (* Normalize the FLX data *)
    AddOpts["Scale"-> 10/StandardDeviation[Flatten[Values[KeySelect[#,FluxDataQ] ]& /@ Values[all[[All, 1]] ] ] ] ];
    (* val = If[val===None,, val]; *)
    Labeled[
        FLXColumnPlot[all,
            FilterRules[{Options[QSFcmdline],opt}, Options[FLXColumnPlot] ],
            ColorAssoc->colors, Scale->OptionValue[{QSFcmdline, FLXJoin}, "Scale"]
        ],
        If[size>1,LineLegend[
            Values[colors], 
            FontFix@Keys[colors], 
            LegendLabel -> "",
            LegendLayout -> "Row"],""], 
        Top]
]; *)
(* Options[FLXColumnPlot] = Join[{ColorAssoc->None, PlotRangeFunction->FLXDefRangeFun, Scale->1, LabelTransform->Identity}, Options[ListLinePlot]];
FLXColumnPlot[all_, opt :OptionsPattern[]]:=Block[{dat,tmax,T, labels},
    (* find all keys *)
    
    labels = DeleteDuplicates[Flatten[Keys /@ Values[all[[All, 1]] ] ]];
    GraphicsColumn[
    Table[
        LOG[magenta,  "Making joined plot for...", label];
        Show[KeyValueMap[(
            Check[dat = First[#2][label];
                meta=Last[#2];
                tmax=meta["tmax"];
                T=meta["T"];
                (* Debug to check that field values are right *)
                (* If[label==="A",Print[2 \[Pi] Max[dat]/T] ]; *)
                ListLinePlot[
                    If[FluxDataQ@label,OptionValue[{QSFcmdline, FLXColumnPlot},Scale],1] dat,
                    FilterRules[{Options[QSFcmdline], opt}, Options[ListLinePlot] ],
                    PlotStyle -> OptionValue[{QSFcmdline, FLXColumnPlot},ColorAssoc][#1],
                    AspectRatio -> 1/2.5, 
                    GridLines -> {
                                    (* Join[ *)
                                        (* {#, Gray}&/@Range[IntegerPart[tmax/T] ], *)
                                        tmp=OptionValue[{QSFcmdline, FLXColumnPlot},PlotRangeFunction][meta];
                                        Table[{First[tmp]+i*(Last[tmp]-First[tmp])/8.0, Red},{i,1,8}]   
                                    (* ] *)
                                    ,Automatic
                                },
                    PlotLegends->None,
                    Frame->True,
                    PlotRange -> {OptionValue[{QSFcmdline, FLXColumnPlot},PlotRangeFunction][meta],Full},
                    FrameTicks->{Automatic,{If[Last@labels===label,Automatic,None],None}},
                    DataRange -> {0, tmax/T}
                ], 
                Nothing]) &, all], 
        (* PlotRange -> Full,  *)
        PlotRange -> {OptionValue[{QSFcmdline, FLXColumnPlot},PlotRangeFunction][meta],Full},
        FrameLabel->{{None,OptionValue[{QSFcmdline, FLXColumnPlot}, LabelTransform][label]},{None,None}},
        (* ImagePadding -> {{30, 20}, {30, 10}}, *)
        ImagePadding -> {{30, 20}, {30 , 10 }} OptionValue[{QSFcmdline, FLXColumnPlot},ImageSize]/400,
        (* ImageSize -> 400, *)
        Frame -> True], 
    {label, labels}],
    (* Alignment -> Right, *)
    Spacings->{-50.0,-45.5 }OptionValue[{QSFcmdline, FLXColumnPlot},ImageSize]/400,
    ImagePadding->{{0,0},{0,15}}
    (* ImagePadding->20 OptionValue[ImageSize]/400 *)
    ]
    
]; *)