BeginPackage["QSF`wf`", {"cmdline`log`","cmdline`opt`", "QSF`styling`", "PlotGrid`","QSF`ColorFunctions`"}];

Begin["`Private`"]
Bool := # != 0 &;
End[];

ProbabilityQ[hd_Association] := Not[hd["isComplex"]];

OpenBin:=Check[OpenRead[#,BinaryFormat->True],Abort[]] &;
WFFileNameQ:=StringEndsQ[#, ".psi" ~~ DigitCharacter..] &;
FromEdge[ind_, ns_] := MapThread[If[#1 < #2/2, -#1, #1 - #2 + 1] &, {ind, ns}];
WFCAPMask[hd_,ratio_:0.5] := Module[{ns, nCAP, eta, CAP, bd},
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
FourierAt[tens_List, d_Integer, inv_ : False] := Module[
    {dim = Length[Dimensions[tens]]},

    Transpose[
        Map[If[inv, InverseFourier[#], Fourier[#] ] &, 
            Transpose[tens, d <-> dim], 
            {dim - 1}], 
    d <-> dim]
];

decorator[LOGF]@
FourierAt[tens_List, dims_List, inv_ : False] := Module[
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
 Module[{path,ns = Dimensions[tens],ratio,nCAP, eta,CAPk},
    ratio=OptionValue[{QSFcmdline, WFMaskEdgeAt},"WFMaskEdgeAtRatio"];
    nCAP = IntegerPart[ns*ratio/2];
    eta = Power[3.0/nCAP , 4];
    CAPk := Exp[-Clip[(eta) . Map[If[# < 0, #^4, 0] &, #], {0, 100}] ] &;
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
    LOG["Removing Bounded Part"]; 
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
TrimMargins[WF[hd_Association, data_List] ] := Module[{p, po},
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
TransverseDiagSum[WF[hd_Association, data_List] ]:=Module[
    {n=First@hd["ns"],h},
    h=1+n/2;
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

WFSameSignMask[hd_, ratio_ : 0.2] := Module[{ns, nCAP, CAP, s},
   ns = hd["ns"]/hd["downscale"];
   nCAP = IntegerPart[ns*ratio/2];
   CAP := Exp[-Power[3.0/nCAP, 4] . Map[If[# > 0, #^4, 0] &, #]] &;
   Array[(s = List[##] - (hd["ns"] + 1) 0.5;
      If[(Total[UnitStep[s] ] == hd["dim"]) || (Total[UnitStep[s] ] == 
          0), 1, 1 - CAP[nCAP - Abs[s] ] ]) &, ns] ];


WFShiftTrick[hd_]:= Array[Mod[Total[List[##]] + 1, 2] 2 - 1 &, hd["ns"]/hd["downscale"] ];


WFInFourierSpace[Op_,hdData_]:= Module[{hd=First@hdData,data=Last@hdData, trick=WFShiftTrick[First@hdData],pscale, fp},
    fp={0, If[hd["rep"]==2,1,-1]};
    {hd, trick InverseFourier[ Op[hd] Fourier[trick data, FourierParameters->fp], FourierParameters->fp]}
];
WFKeepCorrelated[hdData_]:= WFInFourierSpace[WFSameSignMask,hdData];
WFDropCorrelated[hdData_]:= WFInFourierSpace[(1-WFSameSignMask[#]&),hdData]





(* For internal use with open input string *)
WFHeader[st_InputStream]:=Module[{dim=BinaryRead[st, "UnsignedInteger8"]},
    <|"dim"-> dim,"rep"->BinaryRead[st, "UnsignedInteger8"],
    "isComplex"->(1===BinaryRead[st, "Byte"]),
    "downscale"->BinaryRead[st, "Integer64"],
    "ns"->BinaryReadList[st, "Integer64",dim],
    "dxs"->BinaryReadList[st, "Real64",dim],
    "bounded"->BinaryReadList[st, "Byte",dim],
    "positiveAxisOnly"->False
    |>
];
WFHeader[path_String]:=Module[{st=OpenBin[path]},
    hd=WFHeader[st]; Close[path];
    hd
];

WFData[st_, hd_]:=Module[{data},
    data=BinaryReadList[st, If[hd["isComplex"],"Complex128","Real64"] ];
    ArrayReshape[data,hd["ns"]/hd["downscale"] ]
];

WFDataStep[hd_Association]:=Module[{},
    dxs=hd["dxs"];
    dps=2\[Pi]/(ToExpression[hd["ns"]] ToExpression[hd["dxs"] ]);
    If[hd["rep"]==2,dps,dxs]*hd["downscale"]
]

WFNorm[hd_, data_]:=Total[data,\[Infinity] ] Times@@(WFDataStep[hd]); 


SetAttributes[WFLoad, {Listable}];
decorator[LOGF]@
WFLoad[path_String,opt :OptionsPattern[]]:=Module[{st,hd,data,norm},
    st=OpenBin[path];
    hd=WFHeader[st];
    hd["path"]=path;
    data=WFData[st,hd];
    Close[path];
    WF[hd,data]
];


WFDataRange[hd_Association]:=Module[{Xs,Ps,dps},
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
Data[WF[hd_Association, Legended[data_List]]] := data;
Data[WF[hd_Association, data_List]] := data;
SetAttributes[Header, {Listable}];
Header[WF[hd_Association, data_List]] := hd;

SetAttributes[AbsSquare, {Listable}];
AbsSquare[WF[hd_Association, data_List]] := 
  Module[{HD = hd}, If[hd["isComplex"], HD["isComplex"] = False;
    WF[HD, Abs[data]^2], WF[hd, data]]];
AbsSquare[ass_Association]:=Map[AbsSquare,ass];

(* Actions on Associations *)
(* FirstHeader[ass_] := Merge[Level[Header[ass], {-2}], First]; *)




decorator[LOGF]@
Average[inp_] := 
  Module[{r = AbsSquare[inp], fh, dp = ArrayDepth[inp, AllowedHeads -> Association]},
    WF[
        Merge[Cases[r, WF[hd_, data_] :> hd, \[Infinity] ], First], 
        (* Merge[Cases[r, WF[hd_, data_] :> hd, \[Infinity] ], First],  *)
        Mean@Cases[r, WF[_, data_] :> data, \[Infinity] ] 
    ]
];


(* "BlurRadius" and PlotRange are in a.u. *)

Options[MergeOrthants]={"Orthants"->All};
SetAttributes[MergeOrthants,{Listable}];
decorator[LOGF]@
MergeOrthants[WF[hd_Association, data_List],opt :OptionsPattern[{QSFcmdline,MergeOrthants}] ]:=
Module[{mo,res=data},mo=OptionValue[{QSFcmdline, MergeOrthants},"Orthants"];
    res=If[mo===None,LOG["Not merging orthants"]; res,
    =Subsets[Range[hd["dim"] ] ];
    =If[mo==="Correlated", 
            LOG["Merging correlated orthants"];
            Select, Length[#] == hd["dim"] || Length[#] == 0 &], 
            LOG["Merging all orthants"];];
        Total[Reverse[res, #] & /@]
        (* /Length]  *)
    ];
    WF[hd,res]
];


BestLegendPlacement[gr_] := Module[{cm, fo, fo2},
   cm = Counts[
     Join[Position[OptionsV[gr, FrameTicks], None], 
      Position[OptionsV[gr, FrameLabel], None]]];
   fo = Ordering@
      Flatten[AbsoluteOptionsV[gr, ImagePadding]] /. {1 -> {1, 1}, 
      2 -> {1, 2}, 3 -> {2, 1}, 4 -> {2, 2}};
   If[MemberQ[cm, First@fo] && cm[First@fo] > 1, First@fo, 
    First@Keys[TakeLargest[cm, 1]]]
   ];

AbsoluteOptionsV := #2 /. AbsoluteOptions[#1, #2] &;
OptionsV := #2 /. Options[#1, #2] &;
$SidePositions = <|Left -> {1, 1}, Right -> {1, 2}, Bottom -> {2, 1}, 
   Top -> {2, 2}|>;
$PlacePositions = <|{1, 1} -> {{None, Null}}, {1, 
     2} -> {{Null, None}}, {2, 1} -> {{Null}, {None}}, {2, 
     2} -> {{None}, {Null}}|>;

IntegerChop = With[{r = Round[#]}, r + Chop[# - r]] &;
BetterTicks[minmax_, aspectRatio_:1, backgroundColor_:White, expRescale_:False, minor_:{0.01,0.000,10},major_:{0.025,0.00,5}]:=
Module[{ticks, mD, exp=0, expN=1,min=First@minmax, max=Last@minmax, MC, mC},
    MC=ColorNegate@ Apply[RGBColor,ConstantArray[backgroundColor /.{GrayLevel[x_]->x, RGBColor[x___] -> Round[Total[{x}]/3]}, 3]];
    mC=Blend[{MC,Gray}];
    ticks=N@FindDivisions[{min, max}, Last@major, Method -> "ExtendRange"];
    If[expRescale,
        exp=Last@MantissaExponent[ticks[[2]]];
        expN=Power[10,-exp];
        ticks=ticks*expN;
    ];
    mD=0.1*(ticks[[2]]-ticks[[1]]);
    ticks=Flatten[
        Table[If[(i+j*mD>max*expN) || (i+j*mD<min*expN),
            Nothing, 
            If[j==0, {i,Style[IntegerChop[i],BaseStyle /. Options[ArrayPlot]], aspectRatio*Most[major], MC},{i+j*mD,"",aspectRatio*Most[minor], mC}]],
       {i,ticks},
       {j,If[i==First@ticks,-Last[minor]+1,0],Last[minor]-1,1}
    ],1];
    {exp,ticks}
];

AutoBarLegend[gr_Graphics, colorFn_, {min_, max_,n_:9}] := 
  Module[{bp=BestLegendPlacement[gr], fbpi, ip, is, isL, isG, ipG, ipL, ft, arL, ticks, exp, m, M, expText,
  auto={Automatic, Automatic}, none={None,None}},
    (* bp: {Horizontal: Left=1/Right=2, Vertical: Bottom=1, Top=2} *)
   fbpi = Mod[First@bp, 2] + 1;
   ip= AbsoluteOptionsV[gr, ImagePadding];
   is= AbsoluteOptionsV[gr, ImageSize];
   ipL = ReplacePart[{auto, auto}, fbpi->ip[[fbpi]]];
   ft = ReplacePart[{none, none}, bp->All];
   arL = Power[2, 4*(2*fbpi-3)];
    {exp,ticks}=BetterTicks[{min,max}, If[First[bp]==1,1,arL], White, True,{0.2,0.0,10},{0.5,0.0,5}];
    {m, M} = {min, max} Power[10,-exp];
    (* Legend label with common exponent: exp *)
    expText=Style[DisplayForm@SuperscriptBox["\[Times]10", exp], BaseStyle /. Options[ArrayPlot],Background -> Transparent];
    extraPadding=Most[Rasterize[Text[expText],"BoundingBox"]];
    expText=Text[expText, RotateLeft[{0.5*(M+m),M},First@bp-1], RotateLeft[{-0.5,-1},First@bp-1] ];
    extraPadding=Max[ip[[{fbpi,2}]], extraPadding[[fbpi]] ];
    ipG=ReplacePart[ip, {fbpi,2} -> extraPadding];
    isG=is+Map[Total,ipG];
    isL = ReplacePart[auto, fbpi -> isG[[fbpi]] ];
    ipL=ReplacePart[ipL, {fbpi,2} -> extraPadding]; 
    
    (* DeleteCases[ ____, Rule[ImageSize, _], \[Infinity]] *)
    Grid[$PlacePositions[bp] /. {
        Null -> Show[gr,ImagePadding ->ipG, ImageSize->isG],
        None -> DensityPlot[{y,x}[[First@bp]], {x, m, M}, {y, m, M}, 
            PlotRange -> {{m, M}, {m, M}}, 
            ColorFunction -> colorFn, 
            AspectRatio -> arL, 
            FrameTicks -> ReplacePart[{none, none}, bp -> ticks],
            ImagePadding -> ipL, 
            ImageSize -> isL, 
            PlotRangeClipping -> False, 
            FrameStyle->Black,
            Epilog -> expText, 
            BaseStyle -> ImageSizeMultipliers -> 1
            ] 
        }]
   ];

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



WFPlot[WF[hd_Association, leg_Legended]]:=Legended[WFPlot[hd,RemoveLegend[data]], leg /.{Legended[a_, b___] :> b}];

decorator[LOGF]@
WFPlot[WF[hd_Association, data_List], opt :OptionsPattern[{QSFcmdline, WFPlot, ListLinePlot, ArrayPlot, ListDensityPlot3D}] ]:=Module[
    {drng,step, min, max, HD=hd, res=data, leg,pr},
    drng=WFDataRange[hd];
    step=WFDataStep[hd];
    res= If[HD["isComplex"], HD["isComplex"]=False; Abs[res]^2, res];    

    {min, max} = MinMax[RemoveLegends[data]];
    exp = Round[Log[10, max],1]-2;
    MC=Apply[RGBColor, ConstantArray[Round[ColorNegate@ ColorData["Jet"][0] /. RGBColor[x__] -> Total[{x}]/3], 3]];
    mC=Blend[{MC,Gray}];
    pr=OptionValue[{Options[QSFcmdline], opt, {PlotRange->Full}}, PlotRange];
    pr=MapIndexed[If[First[#2]<HD["dim"] && #1===Full, drng[[First[#2]]], #1]&,pr];
    Switch[HD["dim"],
        1,  
        ListLinePlot[
        res,
        FilterRules[{Options[QSFcmdline], opt}, Options[ListLinePlot] ],
        Frame->True,
        FrameTicks->True,
        DataRange->drng],
        2, 
        AutoBarLegend[#,"Jet",{min,max}]&@ 
        ArrayPlot[Reverse@res,
        FilterRules[{Options[QSFcmdline], opt}, Options[ArrayPlot] ],
        DataReversed->False,
        FrameTicks->{{First@Rest@BetterTicks[pr[[1]],1, ColorData["Jet"][0]], None},
                     {First@Rest@BetterTicks[pr[[2]],1, ColorData["Jet"][0]], None}},
        ColorFunction->"Jet",
        ImageSize -> 300,
        DataRange->drng,
        PlotRangePadding->0,
        BaseStyle -> ImageSizeMultipliers -> 1],
        3, ListDensityPlot3D[res,
        FilterRules[{Options[QSFcmdline], opt}, Options[ListDensityPlot3D] ],
        DataReversed->{False, True},
        FrameTicks->True,
        PlotLegends->Placed[Automatic,Right], 
        ColorFunction->Jet,
        DataRange->drng,
        PlotRangePadding -> 0,
        ViewProjection -> "Orthographic",
        BoxStyle -> {LightGray}]  
    ]
   (* Epilog ->Inset[If[norm == 0.0, "-\[Infinity]", Round[Log10@norm, 0.01]], Scaled[{0.95, 0.95}], Scaled[{1, 1}]],
   ColorFunction -> (Color[norm] &), ColorFunctionScaling -> False,
   OpacityFunction -> Transparency, OpacityFunctionScaling -> True,  *)
];


fixNestedLegends = {Legended[Legended[k_, c___], dd___] :> Legended[k, Flatten[{c, dd}, 1] ]};
(* Substitution for Show *)
Options[WFCombine]={"PlotRange"->Full,"LegendPlacement"->Bottom};
WFCombine[ass_Association, keys_List : {}, opt :OptionsPattern[{QSFcmdline, WFCombine}] ] :=
Show[KeyValueMap[WFCombine[#2, Flatten[Append[keys, #1] ], opt ] &, ass], PlotRange->OptionValue["PlotRange"] ] //. fixNestedLegends;
(* Show[WFPlot[Cases[ass,WF[hd_Association, data_List]:>WF[hd,Legended[data,"elo"] ],\[Infinity] ], PlotStyle->GetColor[keys] ] ] //. fixNestedLegends; *)

decorator[LOGF]@
WFCombine[WF[hd_Association, data_List], keys_List:{}, opt :OptionsPattern[{QSFcmdline, WFCombine}] ]:=
WFPlot[WF[hd, Legended[data, Placed[keys, OptionValue[{QSFcmdline, WFCombine},"LegendPlacement"] ] ] ], PlotStyle->GetColor[keys] ];

GraphicsQ[x_Legended|x_Graphics|x_Grid]:=True;
GraphicsQ[x_]:=False;

NameRows[row_,lbl_]:=MapAt[Labeled[#, {lbl}, {Left}]&,row, 1];
GriddedLeaves[ass_Association] :=Module[{ret},
ret=KeyValueMap[If[AssociationQ[#2], Flatten[GriddedLeaves[#2],1], GriddedLeaves[#2] ] &, ass];
LOG["Composed a grid of dimensions", Dimensions[ret]];
ret
];
GriddedLeaves[any_] := {any};
(* Needs["ForScience`PlotUtils`"]; *)
(* PlotGrid[ass_Association,o:OptionsPattern[Join[Options[PlotGrid],Options[Graphics] ] ] ]:=ForScience`PlotUtils`PlotGrid[GriddedLeaves[ass],o]; *)

(* SetAttributes[WFExport,{Listable}]; *)
Options[WFExport]={"ExportPath"->"plots","FileFormat"->".png"};
decorator[LOGF]@
WFExport[gr_?GraphicsQ, keys_:{}, opt :OptionsPattern[{QSFcmdline, WFExport}] ]:=
LOG@Export[
    StringJoin[{
        FileNameJoin[Flatten[{
            StringJoin[Flatten[
                {
                    OptionValue[{QSFcmdline, WFExport},"ExportPath"], 
                    OptionValue[{QSFcmdline, WFExport},"path"]
                }]], 
                keys 
        }]],
        StringPadLeft[OptionValue[{QSFcmdline, WFExport},"FileFormat"],4, "."]
    }]
,gr];

WFExport[ass_Association, keys_:{}]:=MapIndexed[WFExport[#1,Flatten@Join[#2/.{Key[x_]->x}, keys]]&, ass];
(* MapIndexed[(AppendToKey["keys"->#2]; res=WFExport[#1]; PopKey["keys"];res)&,ass]; *)

Multicolumn[ass_Association, opt : OptionsPattern[{QSFcmdline, Multicolumn}]]:=
Multicolumn[KeyValueMap[Labeled[#2, #1] &, ass] ];
SameLegendQ:=SameQ[Cases[#1, LineLegend[_, x_, ___] :> x, \[Infinity] ], Cases[#2, LineLegend[_, y_, ___] :> y, \[Infinity] ] ] &;
(* WFGrid[x_List]:=Legended[Grid[(x//.fixNestedLegends)/.{Legended[a_, b___] :> a},BaseStyle->ImageSizeMultipliers->1], Flatten@Union[Flatten[(x//.fixNestedLegends) /.{Legended[a_, b___] :> b}] ] ]; *)

SubKeys[ass_Association]:=DeleteDuplicatesBy[DeleteCases[Rest/@Position[ass, x_Graphics, \[Infinity], Heads -> True], _Integer, \[Infinity] ],Last];
UnifyLegends[x_List]:=Flatten@Union[Flatten[Flatten[x] /.{Legended[a_, b___] :> b}],SameTest -> SameLegendQ];
RemoveLegends[x_]:=x/.{Legended[a_, b___] :> a};
Options[WFGrid]={"RowLabels"->{}, "ColumnLabels"->{}};
WFGrid[ass_Association]:=WFGrid[GriddedLeaves[ass],"RowLabels"->Keys[ass], "ColumnLabels"->SubKeys[ass] ];


LO[]:=PlotGrid0[{{Plot[x, {x, 0, 1}, Frame -> True], 
   Plot[x, {x, 0, 1}, Frame -> True]}}];

WFGrid[x_List,opt :OptionsPattern[{QSFcmdline, WFGrid}]]:=Legended[PlotGrid1[RemoveLegends[x]], UnifyLegends[x] ];
(* WFGrid[x_List,opt :OptionsPattern[{QSFcmdline, WFGrid}]]:=Legended[
    PlotGrid0[
        Prepend[
        MapThread[
            Prepend[#1,#2]&, 
            {
                RemoveLegends[x],
                OptionValue[{QSFcmdline, WFGrid},"RowLabels"]
            } 
        ], 
        Prepend[OptionValue[{QSFcmdline, WFGrid},"ColumnLabels"] ,""]
        ]
    ], 
    UnifyLegends[x]
]; *)



(* ImportExport`RegisterImport["QSF-wf",WFLoad]; *)
(* ResourceFunction["RegisterFormat"]["QSF-wf", <|"Extensions" -> Table[StringTemplate["psi``"][i], {i, 1, 16}]|>] *)

EndPackage[];

