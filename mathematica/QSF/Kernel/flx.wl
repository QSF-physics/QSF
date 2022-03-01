
BeginPackage["QSF`flx`", {"cmdline`log`","cmdline`opt`"}];

OpenBin:=CHA[OpenRead[#,BinaryFormat->True], "File "<>#<>" cannot be opened."] &;
(* For extracting wavelength data *)
LambdaFromPath:=ToExpression[StringExtract[#,"nm_"->2,"/"->1] ]&;
FWHMFromPath:=ToExpression[StringExtract[#,"fwhm_cycles_"->2, "/"->1] ]&;

FLXDefMetadata[path_, labeled_]:=<|
    "tmax"->Last[labeled["time"] ],
    "dt"->First[Rest[labeled["time"] ] ], 
    "T"->2*\[Pi]/(45.563352529/LambdaFromPath[path]),
    "fwhm"->FWHMFromPath[path]
|>;

FLXFilter:=StringContainsQ[#,(StartOfString~~ (LetterCharacter|PunctuationCharacter)~~"2"~~ (LetterCharacter|"_")..~~EndOfString)]&;

FLXRegions[l_]:= (l=AssociateTo[l,"AA"->"LK"]; Print[Keys[l]]);

(* AssociateTo[l, { "*2S"-> (l["N2S"] - l["S2CAP"] - l["S2D_SYM"]- l["S2D_ASYM"]),  *)
    (* "*2D"-> (l["N2D_SYM"] + l["N2D_ASYM"]+ l["S2D_SYM"]+ l["S2D_ASYM"] - l["D2CAP"]) }]; *)

FAFLXFilter:=StringContainsQ[#,(StartOfString~~"A"~~EndOfString)|(StartOfString~~"F"~~EndOfString)|(StartOfString~~ (LetterCharacter|PunctuationCharacter)~~"2"~~ (LetterCharacter|"_")..~~EndOfString)]&;
Options[FLXLoad] = {"PreprocessTable"->Identity, 
                    "MetadataExtract"->FLXDefMetadata, 
                    "FilterLabeled"->FAFLXFilter,
                    "PostprocessLabeled"->FLXRegions};

FLXLoad[path_, OptionsPattern[{QSFcmdline,FLXLoad}]]:=Block[{st,meta, labels, data, labeled},
    st=OpenBin[path];
    (* TODO: placeholder - should be loaded from binary header *)
    labels = <|"step"->1,"time"->2,"A"->3,"norm"->4,"eta"->5,"N2S"->6,"N2D_SYM"->7,"N2D_ASYM"->8,"S2D_SYM"->9,"S2D_ASYM"->10,"S2CAP" -> 11, "D2CAP" -> 12|>;
    data = BinaryReadList[st, "Real64"];
    Close[st];
    data=Partition[data, 12];
    (* optional sparsification *)
    data=OptionValue[{QSFcmdline,FLXLoad},"PreprocessTable"][data];
    (* data time series is accessible by label index *)
    (* labels=KeySelect[labels, If[labelFilter===Identity, True &,labelFilter] ]; *)
    data=Transpose[data];
    labeled = Map[data[[#]]&, labels];
    meta=OptionValue[{QSFcmdline,FLXLoad},"MetadataExtract"][path, labeled];
    labeled = KeySelect[OptionValue[{QSFcmdline,FLXLoad},"FilterLabeled"] ][labeled];

    AssociateTo[labeled, { "*2S"-> (labeled["N2S"] - labeled["S2D_SYM"]- labeled["S2D_ASYM"] (* - labeled["S2CAP"]  *)),
    
    "*2D"-> (labeled["N2D_SYM"] + labeled["N2D_ASYM"]+ labeled["S2D_SYM"]+ labeled["S2D_ASYM"] (* - labeled["D2CAP"] *)) }];

    labeled=KeySelect[labeled,!StringContainsQ[#,"CAP"]&];
    (* OptionValue["PostprocessLabeled"][labeled]; *)
    {labeled, meta}
];

Options[FLXSmooth] = {"BlurRadius"->0.1};
FLXSmooth[data_, meta_, OptionsPattern[{QSFcmdline,FLXSmooth}]]:=GaussianFilter[data, OptionValue[{QSFcmdline,FLXSmooth},"BlurRadius"] meta["T"]/meta["dt"] ];

Options[ProcessAvgs] = {ShowPopulations->False, SmoothFunction->FLXSmooth};
ProcessAvgs[toAvg_, OptionsPattern[{QSFcmdline,ProcessAvgs}]]:=Block[{lambda,all,data,meta},
    (* Loading *)
    all=Map[FLXLoad, toAvg];
    (* TODO: check that the labels are identical *)
    (* TODO: check that the metadata are identical or at least similar*)
    data=Merge[all[[All,1]], Mean];
    meta=Last[Last[all]];
    (* Smooth the data before std dev gets calculated *)
    data=MapAt[If[OptionValue[{QSFcmdline,ProcessAvgs},ShowPopulations], meta["dt"] Accumulate[#], Identity[#] ] &, data, {#}&/@Select[Keys[data],FLXFilter] ];
    (* F= OptionValue[SmoothFunction]; *)
    (* Print[F]; *)
    data=MapAt[OptionValue[{QSFcmdline,ProcessAvgs},SmoothFunction][#, meta] &, data, {#}&/@Select[Keys[data],FLXFilter] ];
    
    (* Map[
        (data[#]=
        (If[integrateFlux, meta["dt"] Accumulate[#], Identity[#] ] &)@
        GaussianFilter[data[#],blur meta["T"]/meta["dt"] ];) &, 
        Select[Keys[data],FLXFilter]
    ]; *)
    {data, meta}
]

(* DefRangeFun[tmax_,T_]:={Round[0.4*tmax/T,2],Round[tmax/T]-Round[0.2*tmax/T,2]}; *)
(* FLXDefRangeFun[tmax_,T_]:={Round[0.2*tmax/T,2],Round[tmax/T]-4-Round[0.2*tmax/T,2]}; *)
FLXDefRangeFun[meta_]:=Block[{cycles=meta["tmax"]/meta["T"], fwhm=meta["fwhm"]},
    (* Print[cycles, " §FFFF ",fwhm];
    (* Print@{(cycles-4)/2-fwhm, (cycles-4)/2+fwhm}; *) *)
    (* Print@{(cycles-4)/2-fwhm, (cycles-4)/2+fwhm};  *)
    (* {(cycles-4)/2-fwhm, (cycles-4)/2+fwhm} *)
    {(cycles-4)/2-(fwhm/2), (cycles-4)/2+fwhm}
    ];
(* FLXDefRangeFun[meta_]:={0,10}; *)

Options[FLXColumnPlot] = Join[{ColorAssoc->None, PlotRangeFunction->FLXDefRangeFun, Scale->1, LabelTransform->Identity}, Options[ListLinePlot]];
FLXColumnPlot[all_, opt :OptionsPattern[{QSFcmdline,FLXColumnPlot}]]:=Block[{dat,tmax,T, labels},
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
                    If[FLXFilter@label,OptionValue[{QSFcmdline, FLXColumnPlot},Scale],1] dat,
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


FLXJoin[all_, opt :OptionsPattern[{QSFcmdline,FLXColumnPlot}] ]:=Block[{legend,size,colors},
    
    (* legend for data on the same plot and colors*)
    legend=Keys[all];
    size=Length[legend];
    colors=AssociationThread[legend -> ColorData[97, "ColorList"][[;;size]] ];
    
    (* Normalize the FLX data *)
    AddOpts["Scale"-> 10/StandardDeviation[Flatten[Values[KeySelect[#,FLXFilter] ]& /@ Values[all[[All, 1]] ] ] ] ];
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

EndPackage[];
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