
BeginPackage["QSF`flx`", {"cmdline`log`","cmdline`opt`", "QSF`StyleUtils`","QSF`DataAnalysis`","PlotGrid`"}];

FLX;
FLXLoad;
FluxDataQ;
Average;
DataPlots;
DataPlotGrid;
DataFileNameQ;
DataMap;
Begin["`Private`"];

SetAttributes[LookInto,Listable];
LookInto[FLX[hd_Association, data_Association],x_Symbol|x_Function ,opt:OptionsPattern[]]:=With[{act=COptionValue[{opt},"ActOn"], rl=Sequence@@FilterRules[{opt,"DataStep"->hd["dt"]}, Options@x]},
FLX[hd,
  If[act===System`All,Print["map"];
    Map[x[#,rl]&,data],
    MapAt[x[#,rl]&,data,
      Switch[act,
        _List,Print["L"];Intersection[act, Keys@data], 
        _String,Print["S"]; Intersection[{act}, Keys@data],
        _, Print["Any"]; Map[List,Select[Keys@data,act]] (* function *)
        (* Select[Keys@data,COptionValue[{opt,DataMap},"ActOn"] ]) *)
      ]
    ]
  ]]
];

Options[DataMap]={"ActOn"->System`All};
DataMap[x_Symbol|x_Function,Optional[l_String|l_List|l_Symbol|l_Function,COptionValue[{DataMap},"ActOn"]]][ass_Association]:=LookInto[ass,x,"ActOn"->l];


(* DataMap[x_Symbol,l_List][FLX[hd_Association, data_Association]]:= *)



(* Average /: SetDelayed[Average[FLX][ass_Association], rhs_]:=
  f[a:PatternSequence[args]]:=
     dec[Unevaluated @ f, Unevaluated @ {a}, Unevaluated @ rhs]; *)
(* For extracting wavelength data from path *)
LambdaFromPath:=ToExpression[StringExtract[#,"nm_"->2,"/"->1]]&;
FWHMFromPath:=ToExpression[StringExtract[#,"fwhm_cycles_"->2,"/"->1]]&;
DataFileNameQ:=StringEndsQ[#,".dat"]&;

FLXEnrichMetadata[hd_Association,data_Association]:={
  "tmax"->data["time"][[-1]],
  "dt"->data["time"][[2]],
  "T"->2*\[Pi]/(45.563352529/LambdaFromPath[hd["path"]]),
  "fwhm"->FWHMFromPath[hd["path"]]
};

FLXTimeStep[hd_]:=hd["dt"];

FluxDataQ:=StringContainsQ[#,(StartOfString~~ (LetterCharacter|PunctuationCharacter)~~"2"~~ (LetterCharacter|"_")..~~EndOfString)]&;

FLXFilterData:=StringContainsQ[#,(StartOfString~~"A"~~EndOfString)|(StartOfString~~"F"~~EndOfString)|(StartOfString~~ (LetterCharacter|PunctuationCharacter)~~"2"~~ (LetterCharacter|"_")..~~EndOfString)]&;

(* TODO: placeholder - should be loaded from binary header *)
FLXHeader[st_]:=
<|"labels"-><|"step"->1,"time"->2,"A"->3,"norm"->4,"eta"->5,"N2S"->6,"N2D_SYM"->7,
"N2D_ASYM"->8,"S2D_SYM"->9,"S2D_ASYM"->10,"S2CAP"->11,"D2CAP"->12|>|>;
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

Options[FLXLoad]={"EnrichMetadata"->FLXEnrichMetadata, 
                  "DerivedData"->FLXDerivedData};

SetAttributes[FLXLoad, Listable];
decorator[LOGF]@
FLXLoad[path_,opt:OptionsPattern[]]:=Block[{st,hd, data},
  st=OpenBin[path];
  hd=FLXHeader[st];
  hd["path"]=path;  
  data=FLXData[st,hd,opt];
  Close[st];
  AssociateTo[hd,COptionValue[{opt,FLXLoad},"EnrichMetadata"][hd, data]];
  FLX[hd, data]
];


(* FLX /: a_FLX[x__] + b_FLX[x__] := FLX[Extract[a, {x}] + Extract[b, {x}]]; *)
Options[Average] = {ShowPopulations->False};
(* decorator[LOGF]@ *)
Average[matchF_][inp_, OptionsPattern[]]:=FLX[
    Merge[Cases[inp,FLX[hd_,data_]:>hd,\[Infinity]],First],
		Join[
      Merge[Cases[inp,FLX[_,data_]:>KeySelect[data,Not@*matchF],\[Infinity]],First],
      Merge[Cases[inp,FLX[_,data_]:>KeySelect[data,matchF],\[Infinity]],Mean]
    ]
];


Options[FLXColumnPlot] = Join[{ColorAssoc->None, PlotRangeFunction->FLXDefRangeFun, Scale->1, LabelTransform->Identity}, Options[ListLinePlot]];
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
    
];


Options[DataPlots]={"PlotRange"->Full,"LegendLabels"->Identity
  ,"LegendPlacement"->Bottom,"LeafPath"->{}};


FLXPlot[FLX[hd,data], opt:OptionsPattern[]]:=ListLinePlot[data];

DataPlots[props_List][FLX[hd_Association,data_Association],opt:OptionsPattern[]]:=
With[{cycles=hd["tmax"]/hd["T"],fwhm=hd["fwhm"]},
  KeyValueMap[
    ListLinePlot[
      Print[{(cycles-4)/2-(fwhm/2), (cycles-4)/2+fwhm}];
      Legended[#2,COptionValue[{opt,DataPlots},"LeafPath"]]
      (* #2 *)
      ,DataRange -> {0, cycles}
      ,PlotRange->{{(cycles-4)/2-(fwhm/2), (cycles-4)/2+fwhm},Full}
      ,PlotStyle->GetColor[COptionValue[{opt,DataPlots},"LeafPath"]]
      (* ,PlotRange->COptionValue[{opt,DataPlots},"PlotRange"] *)
      ,Frame->True
      ,FrameLabel->{"time",#1}
    ]&
    ,KeySelect[data,MemberQ[props,#]&]
  ]
];

(* decorator[LOGF]@ *)
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

(* //.FixNestedLegends; *)

FLXJoin[all_, opt :OptionsPattern[] ]:=Block[{legend,size,colors},
    
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