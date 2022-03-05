BeginPackage["QSF`wf`", {"cmdline`log`","cmdline`opt`", "QSF`wf`styling`"}];

BoolInv := # == 0 &;
Bool := # != 0 &;
OpenBin:=Check[OpenRead[#,BinaryFormat->True],Abort[]] &;

FromEdge[ind_, ns_] := MapThread[If[#1 < #2/2, -#1, #1 - #2 + 1] &, {ind, ns}];
WFCAPMask[hd_,ratio_:0.5] := Block[{ns, nCAP, eta, CAP, bd},
   ns = hd["ns"]/hd["downscale"];
   nCAP = IntegerPart[ns*ratio/2];
   bd = Bool /@ hd["bounded"];
   eta = Pick[Power[3.0/nCAP, 4], bd];
   CAP := Exp[-eta . Map[If[# > 0, #^4, 0] &, #]] &;
   Array[
    (s = Pick[(nCAP + FromEdge[List[##] - 1, ns]), bd];
      If[AllTrue[s, Negative], 1, CAP[s]]) &, ns]
   ];
WFInvCAPMask[hd_,ratio_:0.5]:= 1 - WFCAPMask[hd,ratio];

(* Makes a fourier transform at dimension d *)
FourierAt[tens_List, d_Integer, inv_ : False] := Block[
    {dim = Length[Dimensions[tens]]},

    Transpose[
        Map[If[inv, InverseFourier[#], Fourier[#] ] &, 
            Transpose[tens, d <-> dim], 
            {dim - 1}], 
    d <-> dim]
];

decorator[LOGF]@
FourierAt[tens_List, dims_List, inv_ : False] := Block[
    {pos = Flatten@Position[dims, 1], dim = Length[Dimensions[tens]], ret = tens},
    If[Length[dims] != dim, LOGE[dims," ",dim, " ",Dimensions[tens]]; Abort[]];
   
    If[AllTrue[dims, 1 == # &], 
        If[inv, InverseFourier[ret], Fourier[ret] ],
        Fold[FourierAt[#1, #2, inv] &, ret, pos] 
    ]
];
cache=<||>;

Options[WFMaskEdgeAt] = {"WFMaskEdgeAtRatio"->0.5};
decorator[LOGF]@
WFMaskEdgeAt[tens_List, dims_List, opt :OptionsPattern[]] := 
 Block[{path,
        ns = Dimensions[tens],
        ratio=OptionValue[{QSFcmdline, WFMaskEdgeAt},"WFMaskEdgeAtRatio"],
        nCAP = IntegerPart[ns*ratio/2], 
        eta = Power[3.0/nCAP , 4], 
        CAPk := Exp[-Clip[(eta) . Map[If[# < 0, #^4, 0] &, #], {0, 100}] ] &},
    path="/tmp/ns"<>StringRiffle[ns,"_"]<>"_d"<>StringRiffle[dims,"_"]<>"_r"<>ToString[ratio]<>".mx";
    
    LOGV[nCAP, " ", ratio];
    If[MissingQ[cache[path] ],  
        
        cache[path]=If[FileExistsQ[path],
            LOG["Loading mask from ", path];
            Import[path],
            LOG["Making new ", path, " mask"];
            Array[N@CAPk[Clip[nCAP + FromEdge[{##} - 1, ns] dims, {-100000, 0}] ] &, ns]
        ];
        If[!FileExistsQ[path],Export[path,cache[path] ];];
    ];
    cache[path] tens 
];


Bounded[hd_Association]:=hd["bounded"];
BoundedQ[hd_Association]:=Total[Bounded[hd]] == 0;

(* WF Modifiers *)
SetAttributes[RemoveBoundedPart,{Listable}];
decorator[LOGF]@
RemoveBoundedPart[WF[hd_Association, data_List] ] :=If[
    BoundedQ[hd], LOG["Pass"]; WF[hd,data], 
    LOGV["Removing Bounded Part ", hd]; 
    WF[
        MapAt[#*0 &,hd,Key["bounded"] ], 
        FourierAt[WFMaskEdgeAt[FourierAt[data, Bounded[hd] ], Bounded[hd] ], Bounded[hd], True]
    ]
];


SetAttributes[Evenify,{Listable}];
Evenify[n_Integer]:=If[OddQ[n],n+1,n];

Options[TrimMargins]={"TrimMarginsPercent"->0.1};
SetAttributes[TrimMargins,{Listable}];
decorator[LOGF]@
TrimMargins[WF[hd_Association, data_List] ] := Block[{p, po},
    p=OptionValue[{QSFcmdline, TrimMargins},"TrimMarginsPercent"]; 
    po=Evenify@Round[hd["ns"] p];
    If[Not[BoundedQ[hd] ], WF[hd,data],
        LOG["Removing margins: ",p, "%, points: ", po]; 
        WF[
            ReplacePart[hd,{"ns"->hd["ns"] - 2 po}],
            (* MapAt[# - 2 Evenify@Round[# p] &, hd, Key["ns"] ], *)
            Apply[data[[##]] &, Map[(1 + Evenify@Round[# p] ;; # - Evenify@Round[# p]) &, hd["ns"] ] ]
        ]
    ]
];



Dim[hd_Association]:=hd["dim"];


cf = Compile[{{data, _Real, 2}, {n, _Integer},{h, _Integer}}, 
    Table[
        Sum[data[[j - i + h, j]], {j, i, n}], 
        {i, h, n}
    ], 
    CompilationTarget -> "C", RuntimeOptions -> "Speed"
];

SetAttributes[TransverseDiagSum,{Listable}];
decorator[LOGF]@
TransverseDiagSum[WF[hd_Association, data_List] ]:=Block[
    {n=First@hd["ns"], h=1+n/2},
    If[Dim[hd]!=2, LOGE["TransverseDiagSum currently supports only 2D wf's"] ];
    If[Apply[UnsameQ,hd["ns"] ], LOGE["TransverseDiagSum requires square wf's"] ];
    If[Apply[UnsameQ,hd["bounded"] ], LOGE["TransverseDiagSum requires wf's of uniform boundedness"] ];
    WF[ReplacePart[hd, 
        {"dim" -> 1, 
        "ns" -> Most[hd["ns"] ]/2,
        "bounded"->Most[Bounded[hd] ],
        "dxs"->Sqrt[2]*Most[hd["dxs"] ],
        "positiveAxisOnly"->True
    }], cf[data,n, h] ]
];

WFSameSignMask[hd_, ratio_ : 0.2] := Block[{ns, nCAP, CAP, s},
   ns = hd["ns"]/hd["downscale"];
   nCAP = IntegerPart[ns*ratio/2];
   CAP := Exp[-Power[3.0/nCAP, 4] . Map[If[# > 0, #^4, 0] &, #]] &;
   Array[(s = List[##] - (hd["ns"] + 1) 0.5;
      If[(Total[UnitStep[s] ] == hd["dim"]) || (Total[UnitStep[s] ] == 
          0), 1, 1 - CAP[nCAP - Abs[s] ] ]) &, ns] ];


WFShiftTrick[hd_]:= Array[Mod[Total[List[##]] + 1, 2] 2 - 1 &, hd["ns"]/hd["downscale"] ];


WFInFourierSpace[Op_,hdData_]:= Block[{hd=First@hdData,data=Last@hdData, trick=WFShiftTrick[First@hdData],pscale, fp},
    fp={0, If[hd["rep"]==2,1,-1]};
    {hd, trick InverseFourier[ Op[hd] Fourier[trick data, FourierParameters->fp], FourierParameters->fp]}
];
WFKeepCorrelated[hdData_]:= WFInFourierSpace[WFSameSignMask,hdData];
WFDropCorrelated[hdData_]:= WFInFourierSpace[(1-WFSameSignMask[#]&),hdData]





(* For internal use with open input string *)
WFHeader[st_InputStream]:=Block[{dim=BinaryRead[st, "UnsignedInteger8"]},
    <|"dim"-> dim,"rep"->BinaryRead[st, "UnsignedInteger8"],
    "isComplex"->(1===BinaryRead[st, "Byte"]),
    "downscale"->BinaryRead[st, "Integer64"],
    "ns"->BinaryReadList[st, "Integer64",dim],
    "dxs"->BinaryReadList[st, "Real64",dim],
    "bounded"->BinaryReadList[st, "Byte",dim],
    "positiveAxisOnly"->False
    |>
];
WFHeader[path_String]:=Block[{st=OpenBin[path]},
    hd=WFHeader[st]; Close[path];
    hd
];

WFData[st_, hd_]:=Block[{data},
    data=BinaryReadList[st, If[hd["isComplex"],"Complex128","Real64"] ];
    ArrayReshape[data,hd["ns"]/hd["downscale"] ]
];

WFDataStep[hd_Association]:=Block[{},
    dxs=hd["dxs"];
    dps=2\[Pi]/(ToExpression[hd["ns"]] ToExpression[hd["dxs"] ]);
    If[hd["rep"]==2,dps,dxs]*hd["downscale"]
]

WFNorm[hd_, data_]:=Total[data,\[Infinity] ] Times@@(WFDataStep[hd]); 


SetAttributes[WFLoad, {Listable}];
WFLoad[path_String,opt :OptionsPattern[]]:=Block[{st,hd,data,norm},
    st=OpenBin[path];
    hd=WFHeader[st];
    hd["path"]=path;
    data=WFData[st,hd];
    Close[path];
    WF[hd,data]
];


WFDataRange[hd_Association]:=Block[{Xs,Ps,dps},
    dps=2\[Pi]/(hd["ns"] hd["dxs"]);
    Xs=0.5 hd["dxs"] (hd["ns"]-1);
    Ps=\[Pi]/hd["dxs"];
    
    If[hd["dim"]==1,First[#],Identity[#] ]&@
    If[hd["positiveAxisOnly"],
        MapThread[{#1,#2}&, If[hd["rep"]==2,{Ps*0, Ps-dps}, {hd["dxs"]*0.5, Xs}] ],
        MapThread[{-#1,#2}&, If[hd["rep"]==2,{Ps, Ps-dps}, {Xs, Xs}] ]
    ]
];



(* Extractors *)
SetAttributes[Data, {Listable}];
Data[WF[hd_Association, data_List]] := data;
SetAttributes[Header, {Listable}];
Header[WF[hd_Association, data_List]] := hd;

SetAttributes[AbsSquare, {Listable}];
AbsSquare[WF[hd_Association, data_List]] := 
  Block[{HD = hd}, If[hd["isComplex"], HD["isComplex"] = False;
    WF[HD, Abs[data]^2], WF[hd, data]]];
AbsSquare[ass_Association]:=Map[AbsSquare,ass];

(* Actions on Associations *)
FirstHeader[ass_] := Merge[Level[Header[ass], {-2}], First];




decorator[LOGF]@
Average[inp_] := 
  Block[{r = AbsSquare[inp], fh, dp = ArrayDepth[inp, AllowedHeads -> Association]},
    fh = FirstHeader[r];
    (* LOGW["PURE", {dp}, {-1 - fh["dim"]}]; DP[r]; *)
    (* Mathematica works with value {-1 - fh["dim"]}
    while Wolfram Script seems to require {-1}, BUG? *)
    r=Level[Data[r], {-1}];
    (* LOGW["DATA AND LEVEL", {dp}, {-1 - fh["dim"]}]; DP[r]; *)
    WF[fh, Mean@r ]
];


(* "BlurRadius" and PlotRange are in a.u. *)

Options[MergeOrthants]={"Orthants"->All};
SetAttributes[MergeOrthants,{Listable}];
decorator[LOGF]@
MergeOrthants[WF[hd_Association, data_List],opt :OptionsPattern[{QSFcmdline,MergeOrthants}] ]:=
Block[{mo,res=data},mo=OptionValue[{QSFcmdline, MergeOrthants},"Orthants"];
    res=If[mo===None,LOG["Not merging orthants"]; res,
        subsets=Subsets[Range[hd["dim"] ] ];
        subsets=If[mo==="Correlated", 
            LOG["Merging correlated orthants"];
            Select[subsets, Length[#] == hd["dim"] || Length[#] == 0 &], 
            LOG["Merging all orthants"]; subsets];
        Total[Reverse[res, #] & /@ subsets]
        (* /Length[subsets]  *)
    ];
    WF[hd,data]
];

Needs["ForScience`"];
SetAttributes[GaussianBlur,{Listable}];
Options[GaussianBlur] = {"GaussianBlurRadius" -> None};
GaussianBlur[WF[hd_Association, data_List], opt :OptionsPattern[{QSFcmdline,GaussianBlur}] ]:=
If[
    OptionValue["GaussianBlurRadius"]===None, WF[hd,data],
    LOG["Applying GaussianBlur with radius ",OptionValue["GaussianBlurRadius"]];
    WF[hd,GaussianFilter[data,{OptionValue["GaussianBlurRadius"]/WFDataStep[hd]} ] ]
]; 

Options[WFPlot] = {"ColorIndex"->1};
SetAttributes[WFPlot,{Listable}];
decorator[LOGF]@
WFPlot[WF[hd_Association, data_List|data_Legended], opt :OptionsPattern[{QSFcmdline, WFPlot, ListLinePlot, ArrayPlot, ListDensityPlot3D}] ]:=Block[
    {drng,step, min, max, subsets, HD=hd, res=data, leg},
    drng=WFDataRange[hd];
    step=WFDataStep[hd];
    
    (* res= If[HD["isComplex"], HD["isComplex"]=False; Abs[res]^2, res]; *)    

    {min, max} = MinMax[data /. {Legended[x_, _] :> x}];
    exp = Round[Log[10, max],1]-2;
    Clear[f];
    f[x_] := x /. {NumberForm[y_, {w_, z_}] :> NumberForm[Round[y/(10^exp), 1], 0]};
    (* leg=OptionValue[{QSFcmdline, WFPlot},"Legend"];
    LOGI[leg]; *)
    LOG["Rules: ",ToString@FilterRules[{Options[QSFcmdline], opt}, Options[ListLinePlot] ] ];
    (* NumberForm[PaddedForm[y/(10^exp), {1, 1}], 2]}; *)
    Switch[HD["dim"],
        1,  
        ListLinePlot[
        res,
        FilterRules[{Options[QSFcmdline], opt}, Options[ListLinePlot] ],
        Frame->True,
        FrameTicks->True,
        DataRange->drng],
        2, 
        ArrayPlot[Reverse@res,
        FilterRules[{Options[QSFcmdline], opt}, Options[ArrayPlot] ],
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
            (* PlotLegends->Automatic, *)
        ColorFunction->ForScience`PlotUtils`Jet,
        DataRange->drng],
        3, ListDensityPlot3D[res,
        FilterRules[{Options[QSFcmdline], opt}, Options[ListDensityPlot3D] ],
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

(* Substitution for Show *)
Options[Combine]={"PlotRange"->Full};
Combine[ass_Association, keys_List : {}, opt :OptionsPattern[{QSFcmdline, Combine}] ] :=Show[KeyValueMap[Combine[#2, Flatten[Append[keys, #1] ] ] &, ass], PlotRange->OptionValue["PlotRange"] ];

decorator[LOGF]@
Combine[WF[hd_Association, data_List], keys_List:{}]:=WFPlot[WF[hd, Legended[data,keys] ],PlotStyle->cmdline`opt`GetColor[keys]];

GraphicsQ[x_Legended|x_Graphics]:=True;
GraphicsQ[x_]:=False;

(* SetAttributes[WFExport,{Listable}]; *)
decorator[LOGF]@
WFExport[gr_?GraphicsQ]:=LOG@Export["wf00/"<>options["groupOutput"]<>".png",gr];

Multicolumn[ass_Association, opt : OptionsPattern[{QSFcmdline, Multicolumn}]]:=
Multicolumn[KeyValueMap[Labeled[#2, #1] &, ass] ];

Options[WFMultiExport]={"OptionGroups"-><|""->{}|>,"OutputType"->".png"};
WFMultiExport[processed_]:=Block[{outputLocation,folder},
    KeyValueMap[
        (
            BackupOpts[];
            LOG["Loading opts [", ToString[#1],"]:\n|- ", #2];
            UpdateOpts[Normal@#2];
            folder="wf0/"<>ToString[#1]<>"/";
            Export[folder<>"opts.txt", ToString[options] ];
            outputLocation=folder<>options["groupOutput"];
            Export[
                outputLocation<>OptionValue[{QSFcmdline, WFMultiExport}, "OutputType" ], 
                Multicolumn[
                    KeyValueMap[(
                        (* LOGI[Keys[#2]]; *)
                        colors=AssociationThread[Keys[#2] -> ColorData[97, "ColorList"][[;;Length[Keys[#2] ] ]] ];
                        Labeled[
                            Show[KeyValueMap[
                                (
                                    (* LOGI[#1, " ", leg[#1]]; *)
                                    (* UpdateOpts[]; *)
                                    (* UpdateOpts[PlotStyle->leg[#1]]; *) 
                                    {hd,data}=#2;
                                    WFPlot[hd,data, "Legend"->#1,PlotStyle->colors[#1] ]
                                
                                )&
                                ,#2] 
                            ,PlotRange->Full]
                            ,FontFix@{#1},{Top} ]
                        )&, processed]//Transpose,
                    3, Appearance -> "Horizontal"]
            ];
            LOGW["|- Exporting WF plot: ", outputLocation];
            RestoreOpts[];
        )&
        ,OptionValue[{QSFcmdline, WFMultiExport}, "OptionGroups"]
    ];
];


(* ImportExport`RegisterImport["QSF-wf",WFLoad]; *)
(* ResourceFunction["RegisterFormat"]["QSF-wf", <|"Extensions" -> Table[StringTemplate["psi``"][i], {i, 1, 16}]|>] *)

EndPackage[];

