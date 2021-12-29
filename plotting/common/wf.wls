OpenBin:=Check[OpenRead[#,BinaryFormat->True],Abort[]] &;
WFHeader[st_]:=Block[{dim=BinaryRead[st, "UnsignedInteger8"]},
<|"dim"-> dim,
"rep"->BinaryRead[st, "UnsignedInteger8"],
"isComplex"->(1===BinaryRead[st, "Byte"]),
"downscale"->BinaryRead[st, "Integer64"],
"ns"->BinaryReadList[st, "Integer64",dim],
"dxs"->BinaryReadList[st, "Real64",dim],
"bounded"->BinaryReadList[st, "Byte",dim]|>];

WFData[st_, hd_]:=Block[{data},
    data=BinaryReadList[st, If[hd["isComplex"],"Complex128","Real64"] ];
    data=If[hd["isComplex"], Abs[data]^2, data];
    ArrayReshape[data,hd["ns"]/hd["downscale"] ] 
];

WFNorm[hd_, data_]:=Total[data,\[Infinity] ] If[hd["rep"]==2,
    (2\[Pi])^(hd["dim"])/(Times@@(hd["dxs"]))/(Times@@(hd["ns"])),
    Times@@(hd["dxs"]) 
];

WFLoad[path_]:=Block[{st,hd,data,norm},
    st=OpenBin[path];
    hd=WFHeader[st];
    data=WFData[st,hd];
    Close[path];
    {hd,data}
]

WFDataRange[hd_]:=Block[{Xs,Ps,dps},
    dps=2\[Pi]/(ToExpression[hd["ns"]] ToExpression[hd["dxs"]]);
    Xs=0.5 hd["dxs"] (hd["ns"]-1);
    Ps=\[Pi]/ToExpression[hd["dxs"]];
    MapThread[{-#1,#2}&, If[hd["rep"]==2,{Ps, Ps-dps}, {Xs, Xs}] ]
]

WFDataStep[hd_]:=Block[{},
    dxs=hd["dxs"];
    dps=2\[Pi]/(ToExpression[hd["ns"]] ToExpression[hd["dxs"]]);
    If[hd["rep"]==2,dps,dxs]*hd["downscale"]
]

(* BlurRadius and PlotRange are in a.u. *)
Options[WFPlot] = Join[{BlurRadius -> 1, QuarterMerge->None}, Options[ListLinePlot],Options[ArrayPlot],Options[ListDensityPlot3D] ];
(* Join[{BlurRadius -> 1, QuarterMerge->None},Options[ListLinePlot],Options[ArrayPlot],Options[ListDensityPlot3D] ]; *)
WFPlot[hd_, data_, opt :OptionsPattern[] ]:=Block[{drng,len,subsets,step,res},
    drng=WFDataRange[hd];
    step=WFDataStep[hd];
    res=GaussianFilter[data,{OptionValue[WFPlot, BlurRadius]/step}];
    
    res=If[OptionValue[WFPlot, QuarterMerge]===None,
        Print["--> Not merging quaters"]; res,
        len = Length[Dimensions[res] ];
        subsets=Subsets[Range[len ] ];
        subsets = If[OptionValue[WFPlot,QuarterMerge]===Correlated, 
            Print["--> Merging correlated quaters"];
            Select[subsets, Length[#] == len || Length[#] == 0 &], 
            Print["--> Merging all quaters"]; 
            subsets];
        Total[Reverse[res, #] & /@ subsets]/Length[subsets] ];

    Switch[hd["dim"],
        1, ListLinePlot[res, 
        Evaluate@FilterRules[{opt}, Options[ListLinePlot] ],
        Frame->True,
        FrameTicks->True,
        DataRange->drng],
        2, ArrayPlot[res,
        Evaluate@FilterRules[{opt}, Options[ArrayPlot]],
        DataReversed->{False, True},
        FrameTicks->True,
        PlotLegends->Placed[Automatic,Right], 
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

PrependName[path_, what_]:=Block[{fsp=FileNameSplit[path]}, FileNameJoin[Join[Most@fsp, {what<>fsp[[-1]]} ] ] ];

ParseInput[filterFun_:Identity]:=(
    (* Read first cmd line arg *)
    inp=$ScriptCommandLine[[2]];
    (* If it points directly to file or is not a directory, leave as is, otherwise add *.psi* pattern *)
    inp=If[FileExistsQ[inp]|| !DirectoryQ[inp], inp, inp<>"*.psi*"];
    (* Read second command line arg, with png format as default *)
    format=Quiet[Check[$ScriptCommandLine[[3]],"png"] ];
    format=If[StringStartsQ[format,"."],format,"."<>format];
    (* Look for patterns in the input *)
    sinp=FileNameSplit[inp];
    Print[sinp];
    pattern=FirstPosition[sinp, n_ /; StringContainsQ[Print[n]; n, "*"], Heads -> False];
    Print[pattern];
    (* Generate a list of files to process, excluding images *)
    inPaths=filterFun[FileNames[inp]];
    Print[inp <> Directory[] ];
    inPaths=Select[inPaths, ! (StringEndsQ[#,".png"] || StringEndsQ[#,".pdf"]) & ];
    Print["Will process following files: "];
    Print/@inPaths;

    relPaths=If[StringStartsQ[ExpandFileName[#],Directory[] ], StringDelete[ExpandFileName[#], Directory[]<>"/" ],#] & /@inPaths;
    outPaths=(If[ pattern === Missing["NotFound"], #, StringRiffle[FileNameSplit[#],"_"] ] <>format) &/@relPaths;
    outPathsMerge=FileNameSplit/@relPaths;
    outPathsMerge= (StringRiffle[If[ pattern === Missing["NotFound"], #,
    pos=First[pattern] - Length[sinp] - 1;
    replacement=If[pos== -1, StringDrop[#[[pos]],-1] <> "_joined",
    StringRiffle[Most[StringSplit[#[[pos]],"_"] ],"_"]<>"_AVG"];
    ReplacePart[#,pos->replacement] ],"_" ]<>format) & /@ outPathsMerge;
);


