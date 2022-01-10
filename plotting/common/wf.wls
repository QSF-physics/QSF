
BoolInv := # == 0 &;
Bool := # != 0 &;
FromEdge[ind_, ns_] := MapThread[If[#1 < #2/2, -#1, #1 - #2 + 1] &, {ind, ns}];
WFMask[hd_,ratio_:0.5] := Block[{ns, nCAP, eta, CAP, bd},
   ns = hd["ns"]/hd["downscale"];
   nCAP = IntegerPart[ns*ratio/2];
   bd = Bool /@ hd["bounded"];
   eta = Pick[Power[3.0/nCAP, 4], bd];
   CAP := Exp[-eta . Map[If[# > 0, #^4, 0] &, #]] &;
   Array[
    (s = Pick[(nCAP + FromEdge[List[##] - 1, ns]), bd];
      If[AllTrue[s, Negative], 1, CAP[s]]) &, ns]
   ];
WFInvMask[hd_,ratio_:0.5]:= 1 - WFMask[hd,ratio];
WFShiftTrick[hd_]:= Array[Mod[Total[List[##]] + 1, 2] 2 - 1 &, hd["ns"]/hd["downscale"] ];



OpenBin:=Check[OpenRead[#,BinaryFormat->True],Abort[]] &;

(* For internal use with open input string *)
WFHeader[st_InputStream]:=Block[{dim=BinaryRead[st, "UnsignedInteger8"]},
    <|"dim"-> dim,"rep"->BinaryRead[st, "UnsignedInteger8"],
    "isComplex"->(1===BinaryRead[st, "Byte"]),
    "downscale"->BinaryRead[st, "Integer64"],
    "ns"->BinaryReadList[st, "Integer64",dim],
    "dxs"->BinaryReadList[st, "Real64",dim],
    "bounded"->BinaryReadList[st, "Byte",dim]|>
];

WFHeader[path_String]:=Block[{st=OpenBin[path]},
    hd=WFHeader[st]; Close[path];
    hd
];

WFData[st_, hd_]:=Block[{data},
    data=BinaryReadList[st, If[hd["isComplex"],"Complex128","Real64"] ];
    ArrayReshape[data,hd["ns"]/hd["downscale"] ] 
];

WFDataStep[hd_]:=Block[{},
    dxs=hd["dxs"];
    dps=2\[Pi]/(ToExpression[hd["ns"]] ToExpression[hd["dxs"] ]);
    If[hd["rep"]==2,dps,dxs]*hd["downscale"]
]

WFNorm[hd_, data_]:=Total[data,\[Infinity] ] Times@@(WFDataStep[hd]); 

Options[WFLoad] = {Probability -> False};
WFLoad[path_,opt :OptionsPattern[]]:=Block[{st,hd,data,norm},
    st=OpenBin[path];
    hd=WFHeader[st];
    data=WFData[st,hd];
    Close[path];
    {hd,If[OptionValue[WFLoad, Probability] && hd["isComplex"], Abs[data]^2, data]}
]


WFDataRange[hd_]:=Block[{Xs,Ps,dps},
    dps=2\[Pi]/(ToExpression[hd["ns"]] ToExpression[hd["dxs"]]);
    Xs=0.5 hd["dxs"] (hd["ns"]-1);
    Ps=\[Pi]/ToExpression[hd["dxs"]];
    MapThread[{-#1,#2}&, If[hd["rep"]==2,{Ps, Ps-dps}, {Xs, Xs}] ]
]

Needs["ForScience`"];
(* BlurRadius and PlotRange are in a.u. *)
Options[WFPlot] = Join[{BlurRadius -> None, QuarterMerge->None}, Options[ListLinePlot],Options[ArrayPlot],Options[ListDensityPlot3D] ];
WFPlot[hd_, data_, opt :OptionsPattern[] ]:=Block[{drng,len,subsets,step,res},
    drng=WFDataRange[hd];
    step=WFDataStep[hd];
    res=If[OptionValue[WFPlot, BlurRadius]===None,data,
    GaussianFilter[data,{OptionValue[WFPlot, BlurRadius]/step}] ];
    
    res=If[OptionValue[WFPlot, QuarterMerge]===None,
        (* Print["--> Not merging quaters"]; *)
        res,
        len = Length[Dimensions[res] ];
        subsets=Subsets[Range[len ] ];
        subsets = If[OptionValue[WFPlot,QuarterMerge]===Correlated, 
            (* Print["--> Merging correlated quaters"]; *)
            Select[subsets, Length[#] == len || Length[#] == 0 &], 
            (* Print["--> Merging all quaters"];  *)
            subsets];
        Total[Reverse[res, #] & /@ subsets]/Length[subsets] ];
    tenPower=Table[10^(-i),{i,0,20}];
        
    {min, max} = {Min[res],Max[res]};
    exp = Round[Log[10, max],1]-2;
    Clear[f];
    f[x_] := x /. {NumberForm[y_, {w_, z_}] :> NumberForm[Round[y/(10^exp), 1], 0]};
    (* NumberForm[PaddedForm[y/(10^exp), {1, 1}], 2]}; *)
    Switch[hd["dim"],
        1, ListLinePlot[res, 
        Evaluate@FilterRules[{opt}, Options[ListLinePlot] ],
        Frame->True,
        FrameTicks->True,
        DataRange->drng],
        2, ArrayPlot[Reverse@res,
        Evaluate@FilterRules[{opt}, Options[ArrayPlot]],
        DataReversed->False,
        FrameTicks->{{True, None},{True, None}},
        
        PlotLegends->Placed[BarLegend[{Automatic, {min,max}},
            LegendFunction -> f, 
            LegendMarkerSize -> {136,10},
            LegendMargins -> {{0, 20}, {0, -5}},
            (* LegendMarkerSize -> {10, 15}, *)
            LegendLabel -> Placed[DisplayForm[SuperscriptBox[ToString["\[Times]10"], exp] ], Right]
            ],
            Below], 
        ColorFunction->ForScience`PlotUtils`Jet,
        DataRange->drng],
        3, ListDensityPlot3D[res,
        Evaluate@FilterRules[{opt}, Options[ListDensityPlot3D] ],
        DataReversed->{False, True},
        FrameTicks->True,
        PlotLegends->Placed[Automatic,Right], 
        ColorFunction->ForScience`PlotUtils`Jet,
        DataRange->drng,
        PlotRangePadding -> 0,
        ViewProjection -> "Orthographic",
        BoxStyle -> {LightGray}]  
    ]
   (* Epilog ->Inset[If[norm == 0.0, "-\[Infinity]", Round[Log10@norm, 0.01]], Scaled[{0.95, 0.95}], Scaled[{1, 1}]],
   ColorFunction -> (Color[norm] &), ColorFunctionScaling -> False,
   OpacityFunction -> Transparency, OpacityFunctionScaling -> True,  *)
];


