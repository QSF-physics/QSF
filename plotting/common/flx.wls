
OpenBin:=Check[OpenRead[#,BinaryFormat->True],Abort[]] &;

FLXLoad[path_, labelFilter_:Identity, lambda_:1, every_:1]:=Block[{st,meta, labels, data},
    st=OpenBin[path];
    (* TODO: placeholder - should be loaded from binary header *)
    labels = <|"step"->1,"time"->2,"A"->3,"norm"->4,"eta"->5,"N2S"->6,"N2D_SYM"->7,"N2D_ASYM"->8,"S2D_SYM"->9,"S2D_ASYM"->10,"S2CAP" -> 11, "D2CAP" -> 12|>;
    data = BinaryReadList[st, "Real64"];
    Close[st];
    data=Partition[data, 12];
    (* optional sparsification *)
    data=If[every>1, data[[;; ;;every]], data];
    (* meta= {tmax, dt, T}*)
    meta =<|"tmax"->data[[-1, 2]],"dt"->data[[2, 2]], "T"->2*\[Pi]/(45.563352529/lambda)|>;
    (* data time series is accessible by label index *)
    labels=KeySelect[labels, If[labelFilter===Identity, True &,labelFilter] ];
    data=Transpose[data];
    
    {Map[data[[#]]&, labels], meta}
];
FLXFilter:=StringContainsQ[#,"N2S"|"N2D"|"S2D"]&;

(* DefRangeFun[tmax_,T_]:={Round[0.4*tmax/T,2],Round[tmax/T]-Round[0.2*tmax/T,2]}; *)
(* DefRangeFun[tmax_,T_]:={Round[0.2*tmax/T,2],Round[tmax/T]-Round[0.3*tmax/T,2]}; *)
DefRangeFun[tmax_,T_]:={0,Round[tmax/T]};

Options[FLXColumnPlot] = Join[{ColorAssoc->None, PlotRangeFunction->DefRangeFun, Scale->1, LabelTransform->Identity}, Options[ListLinePlot]];
FLXColumnPlot[all_, opt :OptionsPattern[]]:=Block[{dat,tmax,T, labels},
    (* find all keys *)
    
    labels = DeleteDuplicates[Flatten[Keys /@ Values[all[[All, 1]] ] ]];
    GraphicsColumn[
    Table[
        LOG[magenta,  "Making joined plot for...", label];
        Show[KeyValueMap[(
            Check[dat = First[#2][label];
                tmax=Last[#2]["tmax"];
                T=Last[#2]["T"];
                (* Debug to check that field values are right *)
                (* If[label==="A",Print[2 \[Pi] Max[dat]/T] ]; *)
                ListLinePlot[
                    If[FLXFilter@label,OptionValue[FLXColumnPlot,Scale],1] dat,
                    Evaluate@FilterRules[{opt}, Options[ListLinePlot] ],
                    PlotStyle -> OptionValue[FLXColumnPlot,ColorAssoc][#1],
                    AspectRatio -> 1/2.5, 
                    GridLines -> {Range[IntegerPart[tmax/T] ],Automatic},
                    PlotLegends->None,
                    Frame->True,
                    PlotRange -> {OptionValue[FLXColumnPlot,PlotRangeFunction][tmax, T],Full},
                    FrameTicks->{Automatic,{If[Last@labels===label,Automatic,None],None}},
                    DataRange -> {0, tmax/T}
                ], 
                Nothing]) &, all], 
        (* PlotRange -> Full,  *)
        PlotRange -> {OptionValue[FLXColumnPlot,PlotRangeFunction][tmax, T],Full},
        FrameLabel->{{None,OptionValue[FLXColumnPlot,LabelTransform][label]},{None,None}},
        (* ImagePadding -> {{30, 20}, {30, 10}}, *)
        ImagePadding -> {{30, 20}, {30 , 10 }} OptionValue[FLXColumnPlot,ImageSize]/400,
        (* ImageSize -> 400, *)
        Frame -> True], 
    {label, labels}],
    (* Alignment -> Right, *)
    Spacings->{-50.0,-45.5 }OptionValue[FLXColumnPlot,ImageSize]/400,
    ImagePadding->{{0,0},{0,15}}
    (* ImagePadding->20 OptionValue[FLXColumnPlot,ImageSize]/400 *)
    ]
    
];
val = None;
Options[FLXJoin] = Join[{}, Options[FLXColumnPlot]];
FLXJoin[all_, opt :OptionsPattern[]]:=Block[
    {legend,size,colors},
    (* legend for data on the same plot and colors*)
    legend=Keys[all];
    size=Length[legend];
    colors=AssociationThread[legend -> ColorData[97, "ColorList"][[;;size]] ];
    
    (* Normalize the FLX data *)
    val = If[val===None,Flatten[Values[KeySelect[#,FLXFilter] ]& /@ Values[all[[All, 1]] ] ], val];
    Labeled[
        FLXColumnPlot[all,
            Evaluate@FilterRules[{opt}, Options[FLXColumnPlot] ],
            ColorAssoc->colors, Scale->10.0 /StandardDeviation[val] 
        ],
        If[size>1,LineLegend[
            Values[colors], 
            FontFix/@Keys[colors], 
            LegendLabel -> FontFix[""],
            LegendLayout -> "Row"],""], 
        Top]
];

    (* ListLinePlot[
            flx,
            (* Prepend[flx,data[[All,3]] ], *)
            PlotRange -> {Full, Full},
            DataRange -> {0, tmax/T},
            ImagePadding -> {{0, 10}, {0, 0}},
            PlotRangePadding -> {None, None},
            PlotLayout->{"Column", UpTo[1]},
            AspectRatio -> 1/3, 
            GridLines -> Automatic,
            ImageSize -> Large
            ,Method -> {"Spacings" -> {0, 0}, 
                PlotRange -> {Full, Automatic},
                "ColumnLabels" -> "time [cycles]" 
                (* ,"RowLabels" -> names *)
                ,"RowLabels" -> Rest[names]
                }
            
            ] *)