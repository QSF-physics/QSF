
OpenBin:=Check[OpenRead[#,BinaryFormat->True],Abort[]] &;

FLXLoad[path_]:=Block[{st,labels, data},
    st=OpenBin[path];
    (* placeholder *)
    (* {3, A}, {5, NORM}, *)
    labels = {{3, A},   {6, N2S}, {7, N2D_SYM}, {8, N2D_ASYM}, {9, S2D_SYM}, {10, S2D_ASYM}};
    data = BinaryReadList[st, "Real64"];
    Close[st];
    data=Partition[data, 10];
    (* [[;; 5]]; *)
    {labels, data}
];
Options[WFPlot] = Join[{}, Options[ListLinePlot]];
FLXPlot[data_, labels_, lambda_:1, blur_:0.1, opt :OptionsPattern[]]:=Block[{last,tmax,dt, T,flx},
    tmax = data[[-1, 2]];
    dt = data[[2, 2]];
    T = 2*\[Pi]/(45.563352529/lambda);
    indices=First[Transpose[labels] ];
    names=Last[Transpose[labels] ];
    Print[tmax/T];
    flx=Transpose[data][[Rest@indices]];
    flx=GaussianFilter[#,blur T/dt] & /@flx;
    (* flx=flx/StandardDeviation[Flatten[flx]]; *)
    flx=10.0 flx/StandardDeviation[Flatten[flx]];
    (* /Max[Abs[flx]]; *)
    (* Map[If[Last[#]===A, data[[All,First[#]]], Standardize[data[[All,First[#] ]] ] ]&, labels] *)
    (* for currents *)
    (*dt Accumulate[GaussianFilter[#,0.05 T/dt]]& *)
    (* fun := If[#2 === A || #2 === NORM, Identity, GaussianFilter[#, blur T/dt] &]; *)
    (* grids[min_, max_] := Join[Range[Ceiling[min], Floor[max]], Table[{j + .5, Dashed}, {j, Round[min], Round[max - 1], 1}]]; *)
    GraphicsColumn[
        MapThread[ListLinePlot[#1,
            Evaluate@FilterRules[{opt}, Options[ListLinePlot] ],
            PlotRange -> {{tmax/T/2.5, 0.75 tmax/T},Full},
            DataRange -> {0, tmax/T},
            (* ImagePadding -> {{60, 0}, {If[Last@names===#2,20,0], 0}}, *)
            (* ImagePadding -> {{60, 20}, {20, 20}}, *)
            ImagePadding -> {{30, 20}, {30, 10}},
            PlotRangePadding -> {None, Scaled[0.1]},
            (* PlotLayout->{"Column", UpTo[1]}, *)
            Frame->True,
            AspectRatio -> 1/2.5, 
            FrameLabel->{{None,#2},{None,None}},
            GridLines -> {Range[IntegerPart[tmax/T]],Automatic},
            (* GridLines->Automatic, *)
            FrameTicks->{Automatic,{If[Last@names===#2,Automatic,None],None}},
            (* Background->White, *)
            (* Prolog-> {White, Rectangle[ImageScaled[{0, 0}], ImageScaled[{1, 1}]] }, *)
            ImageSize -> 400]&,{Prepend[flx,data[[All,3]] ],names}],
            Spacings->{0,-46},
            ImagePadding->20
    ]
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
];