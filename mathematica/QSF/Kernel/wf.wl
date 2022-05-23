BeginPackage["QSF`wf`", {"cmdline`log`","cmdline`opt`", "QSF`StyleUtils`", "PlotGrid`","QSF`ColorFunctions`","QSF`DataAnalysis`", "QSF`Physics`"}];

(* Whether given filename has a signature consistent with QSF package *)
WF;
WFFileNameQ;
WFDataStep;
WFLoad;
TrimMargins;
Average;
RemoveBoundedPart;
MergeOrthants;
PixelCircle;
PolarIntegral;
TransverseDiagSum;
WFCombine;
WFExport;
WFPlot;
WFPlotGrid;
GridKeys;

Begin["`Private`"];

Bool:=#!=0&;
ProbabilityQ[hd_Association] := Not[hd["isComplex"]];
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


cache=<||>;

Options[WFMaskEdgeAt] = {"WFMaskEdgeAtRatio"->0.5};
decorator[LOGF]@
WFMaskEdgeAt[tens_List, dims_List, opt:OptionsPattern[]] := 
 Module[{path,ns = Dimensions[tens],ratio,nCAP, eta,CAPk},
	ratio=COptionValue[WFMaskEdgeAt,"WFMaskEdgeAtRatio"];
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
		ReplacePart[hd,"bounded"->ConstantArray[0,hd["dim"]]], 
		FourierAt[WFMaskEdgeAt[FourierAt[data, Bounded[hd] ], Bounded[hd] ], Bounded[hd], True]
	]
];


SetAttributes[Evenify,{Listable}];
Evenify[n_Integer]:=If[OddQ[n],n+1,n];

Options[TrimMargins]={"TrimMarginsPercent"->0.1};
(* SetAttributes[TrimMargins,{Listable}]; *)
decorator[LOGF]@
TrimMargins[WF[hd_Association,data_List]]:=Module[{p,po,HD},
	p=COptionValue[TrimMargins,"TrimMarginsPercent"]; 
	po=Evenify@Round[hd["ns"] p];
	If[Not[BoundedQ[hd]],WF[hd,data],
    HD=ReplacePart[hd,{"ns"->hd["ns"] - 2 po}];
		LOG["Removing margins: ",p*100, "%, points: ", po, " (extracting ranges [", Map[(1+Evenify@Round[# p];;#-Evenify@Round[# p])&,hd["ns"]], "])"]; 
		WF[
			HD,
			Apply[data[[##]]&,Map[(1+Evenify@Round[# p];;#-Evenify@Round[# p])&,hd["ns"]]]
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
  Print[Most[hd["ns"]]/2];
  Print[Most[Bounded[hd]]];
  Print[Most[hd["dxs"]]];
	h=1+n/2;
	If[Dim[hd]!=2, LOGE["TransverseDiagSum currently supports only 2D wf's"] ];
	If[Apply[UnsameQ,hd["ns"] ], LOGE["TransverseDiagSum requires square wf's"] ];
	If[Apply[UnsameQ,hd["bounded"] ], LOGE["TransverseDiagSum requires wf's of uniform boundedness"] ];
	WF[ReplacePart[hd, 
		{"dim" -> 1, 
		"ns" -> Most[hd["ns"]]/2,
		"bounded"->Most[hd["bounded"]],
		"dxs"->Sqrt[2]*Most[hd["dxs"]],
		"positiveAxisOnly"->True
	}], cf[data,n, h] ]
];


(* PixelCircle60deg = Compile[{{data, _Real, 2}, {rad, _Real}, {n, _Integer}}, 
	Table[If[(i - n/2) (j - n/2) > 0 && 
      Max@Abs[{i - n/2, j - n/2}]/Min@Abs[{i - n/2, j - n/2}] < 
       2 + Sqrt[3], UnitStep[-Sqrt[(i - n/2)^2 + (j - n/2)^2] + rad], 
     0], {i, 1, n}, {j, 1, n}] data, 
	CompilationTarget -> "C", RuntimeOptions -> "Speed"
]; *)

PixelCircle60deg = Compile[{{data, _Real, 2}, {rad, _Real}, {n, _Integer}, {sign, _Integer}}, 
	Table[If[sign i j > 0 && 
      Max@Abs[{i, j}]/Min@Abs[{i, j}] < 2 + Sqrt[3], 
      data[[i+n+1,j+n+1]] UnitStep[-Sqrt[i^2+j^2]+rad], 
     0.0], {i,-n,n-1}, {j,-n,n-1}], 
	CompilationTarget -> "C", RuntimeOptions -> "Speed"
];

Options[PixelCircle]={"Orthants"->All};
SetAttributes[PixelCircle,{Listable}];
PixelCircle[WF[hd_,data_],opt:OptionsPattern[]]:=WF[hd,PixelCircle60deg[data, Length[data]/2, Length[data]/2,
  If[COptionValue[{opt,MergeOrthants},"Orthants"]==="Correlated",1,-1]] ];


SetAttributes[PolarIntegral,{Listable}];
Options[PolarIntegral]={"Orthants"->All};
decorator[LOGF]@
PolarIntegral[WF[hd_Association, data_List],opt:OptionsPattern[]]:=Module[
	{n=First@hd["ns"],h,pts,rs,tmp},
	h=1+n/2;
	If[Dim[hd]!=2, LOGE["TransverseDiagSum currently supports only 2D wf's"] ];
	If[Apply[UnsameQ,hd["ns"] ], LOGE["TransverseDiagSum requires square wf's"] ];
	If[Apply[UnsameQ,hd["bounded"] ], LOGE["TransverseDiagSum requires wf's of uniform boundedness"] ];
  pts=N@Table[k, {k, n/2 + 0.5, 0.5, -1}];
  rs=Most@pts+Differences[pts]/2;
	WF[ReplacePart[hd, 
		{"dim"->1, 
		"ns"->Most[hd["ns"]]/2,
		"bounded"->Most[hd["bounded"]],
		"dxs"->Most[hd["dxs"]],
		"positiveAxisOnly"->True
	}], 
  Reverse[Rest@Last@Last@Reap@Fold[
    (tmp=PixelCircle60deg[#1,#2,Length[#1]/2,If[COptionValue[{opt,MergeOrthants},"Orthants"]==="Correlated",1,-1]];Sow[Total[Total[#1-tmp]]];Print[#2];
    If[Length[#1]/2-#2 >2,tmp[[2;;-2,2;;-2]],tmp])&,data,pts]
    ]
  ]
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

WFDataStep[hd_Association]:=Module[{dxs,dps},
	dxs=hd["dxs"];
	dps=2\[Pi]/(ToExpression[hd["ns"]] ToExpression[hd["dxs"] ]);
	If[hd["rep"]==2,dps,dxs]*hd["downscale"]
]

WFNorm[hd_, data_]:=Total[data,\[Infinity] ] Times@@(WFDataStep[hd]); 


SetAttributes[WFLoad, Listable];
decorator[LOGF]@
WFLoad[path_,opt:OptionsPattern[]]:=Module[{st,hd,data},
	st=OpenBin[path];
	hd=WFHeader[st];
	hd["path"]=path;
	data=WFData[st,hd];
	Close[path];
	WF[hd,data]
];

Options[WFDataRange]={"RangeScale"->1};
WFDataRange[hd_Association,opt:OptionsPattern[]]:=Module[{Xs,Ps,dps},
	dps=2\[Pi]/(hd["ns"] hd["dxs"]);
	Xs=0.5 hd["dxs"] (hd["ns"]-1);
	Ps=\[Pi]/hd["dxs"];
	
	(* If[hd["dim"]==1,First[#],Identity[#] ]&@ *)
	If[hd["positiveAxisOnly"],
		MapThread[{#1,#2}&, If[hd["rep"]==2,{Ps*0, Ps-dps}, {hd["dxs"]*0.5, Xs}] ],
		MapThread[{-#1,#2}&, If[hd["rep"]==2,{Ps, Ps-dps}, {Xs, Xs}] ]
	]/ReleaseHold[(COptionValue[{opt,WFDataRange},"RangeScale"] /.hd)]
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



decorator[LOGF]@
Average[inp_]:=Module[
{r=AbsSquare[inp],dp=ArrayDepth[inp,AllowedHeads->Association]},
	WF[
		Merge[Cases[r,WF[hd_,data_]:>hd,\[Infinity]],First],
		Mean@Cases[r,WF[_,data_]:>data,\[Infinity]] 
	]
];


(* "BlurRadius" and PlotRange are in a.u. *)

Options[MergeOrthants]={"Orthants"->All};
SetAttributes[MergeOrthants,{Listable}];
decorator[LOGF]@
MergeOrthants[WF[hd_Association, data_List],opt:OptionsPattern[]]:=
Module[{mo,res=data},mo=COptionValue[{opt,MergeOrthants},"Orthants"];
	res=If[mo===None,LOG["Not merging orthants"]; res,
		subsets=Subsets[Range[hd["dim"] ] ];
		subsets=If[mo==="Correlated", 
			LOG["Merging correlated orthants"];
			Select[subsets, Length[#] == hd["dim"] || Length[#] == 0 &], 
			LOG["Merging all orthants"]; subsets];
		Total[Reverse[res, #] & /@ subsets]
		(* /Length[subsets]  *)
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




SetAttributes[WFPlot,{Listable}];
Options[WFPlot] = {"ColorIndex"->1};

decorator[LOGF]@
WFPlot[WF[hd_Association, data_List|data_Legended],opt:OptionsPattern[]]:=Module[
	{min, max, HD=hd, res=data, prng},
  With[{drng=WFDataRange[hd], step=WFDataStep[hd]},
  (* Default PlotRange will show the whole DataRange *)
	prng=COptionValue[{opt, PlotRange->Append[Full][drng]}, PlotRange];
	(* Attempt to keep PlotRange symbol free *)
  prng=MapThread[If[ListQ[#1],#1,#2]&,{prng,Append[Full][drng]}];

	res=If[HD["isComplex"], HD["isComplex"]=False; Abs[RemoveLegend[res]]^2, res];    
	{min, max} = MinMax[res];
	MC=Apply[RGBColor, ConstantArray[Round[ColorNegate@ ColorData["Jet"][0] /. RGBColor[x__] -> Total[{x}]/3], 3]];
	mC=Blend[{MC,Gray}];
  (* colorf=COptionValue[{opt,WFPlot},"ColorFunction"] *)
	Switch[HD["dim"],
		1,  
		ListLinePlot[
		res,
		FilterRules[{Options[QSFcmdline], opt}, Options[ListLinePlot] ],
    (* PlotRangePadding->Scaled[0.08], *)
		Frame->True,
		FrameTicks->True,
    PlotTheme -> "Detailed",
		DataRange->First@drng],
		2, 
		AutoBarLegend[#,"Jet",{min,max}]&@ 
		ArrayPlot[Reverse@res,
		FilterRules[{Options[QSFcmdline], opt}, Options[ArrayPlot] ],
		DataReversed->False,
		FrameTicks->{{First@Rest@BetterTicks[prng[[1]],1, ColorData["Jet"][0]], None},
					 {First@Rest@BetterTicks[prng[[2]],1, ColorData["Jet"][0]], None}},
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
]];



(* Substitution for Show *)
Options[WFCombine]={
  "PlotRange"->Full,"LegendLabels"->Identity
  ,"LegendPlacement"->Bottom,"LeafPath"->{}};
(* WFCombine[ass_Association,opt:OptionsPattern[]] :=
Show[KeyValueMap[WFCombine[#2, Flatten[Append[COptionValue[{opt,WFCombine},"LeafPath"], #1]],opt]&,ass],PlotRange->OptionValue["PlotRange"] ] //. FixNestedLegends; *)


(* [WFExport[#1,"LeafPath"->Flatten@Join[COptionValue[{opt,WFExport},"LeafPath"],#2/.{Key[x_]->x}]]&, ass]; *)


WFCombine[WF[hd_Association,data_List],opt:OptionsPattern[]]:=
WFPlot[
	WF[
		hd,
		Legended[
			data, 
			Placed[
				Style[
          COptionValue[{opt,WFCombine},"LegendLabels"]@
          Flatten[COptionValue[{opt,WFCombine},"LeafPath"]//.{Key[a__]:>a}]
        , defaultStyle]
				,COptionValue[{opt,WFCombine},"LegendPlacement"]
			]
		]
	]
	,PlotStyle->GetColor[COptionValue[{opt,WFCombine},"LeafPath"]]
];

decorator[LOGF]@
WFCombine[ass_Association,opt:OptionsPattern[]]:=
Show[
	Cases[
		MapIndexed[
			If[MatchQ[#1,_WF],WFCombine[#1,"LeafPath" -> #2],#1]&
			,ass
			,ArrayDepth[ass, AllowedHeads -> Association]
		], _Legended, ArrayDepth[ass, AllowedHeads -> Association]
	]
	,PlotRange->COptionValue[{opt,WFCombine},"PlotRange"]
]//.FixNestedLegends;

(* 
decorator[LOGF]@
WFCombine[x_Graphics,opt:OptionsPattern[]]:=
Legended[x, Placed[COptionValue[{opt,WFCombine},"LeafPath"],COptionValue[WFCombine,"LegendPlacement"]]];

WFCombine[ass_Association,opt:OptionsPattern[]]:=Show[Cases[
MapIndexed[If[MatchQ[#1,_Graphics],WFCombine[#1,"LeafPath" -> #2],#1]&,ass, 3], x_Legended, 3
],PlotRange->COptionValue[{opt,WFCombine},"PlotRange"]]; *)


(* Needs["ForScience`PlotUtils`"]; *)
(* PlotGrid[ass_Association,o:OptionsPattern[Join[Options[PlotGrid],Options[Graphics] ] ] ]:=ForScience`PlotUtils`PlotGrid[GriddedLeaves[ass],o]; *)

(* SetAttributes[WFExport,{Listable}]; *)
Options[WFExport]={"ExportPath"->"wf_plots/", "TreePath"->{}, "LeafPath"->{},"FileFormat"->".png"};
decorator[LOGF]@
WFExport[gr_?ExportableQ,opt:OptionsPattern[]]:=
(LOG@Export[
	StringJoin[{
		FileNameJoin[Flatten[{
			StringJoin[Flatten[
				{
					LOG["TreePath: ",COptionValue[WFExport,"TreePath"]];
					LOG["LeafPath: ",COptionValue[{opt,WFExport},"LeafPath"]];
					COptionValue[{opt,WFExport},"ExportPath"] 
					,COptionValue[{opt,WFExport},"TreePath"]
				}]] 
				,COptionValue[{opt,WFExport},"LeafPath"]
		}]],
		StringPadLeft[COptionValue[WFExport,"FileFormat"],4, "."]
	}]
,gr];gr)

WFExport[ass_Association,opt:OptionsPattern[]]:=MapIndexed[WFExport[#1,"LeafPath"->Flatten@Join[COptionValue[{opt,WFExport},"LeafPath"],#2/.{Key[x_]->x}]]&, ass];

Multicolumn[ass_Association,opt:OptionsPattern[]]:=
Multicolumn[KeyValueMap[Labeled[#2, #1] &, ass] ];


(* SameQ[Cases[#1, LineLegend[_, x_, ___] :> x, \[Infinity] ], Cases[#2, LineLegend[_, y_, ___] :> y, \[Infinity] ] ] &; *)
(* WFGrid[x_List]:=Legended[Grid[(x//.FixNestedLegends)/.{Legended[a_, b___] :> a},BaseStyle->ImageSizeMultipliers->1], Flatten@Union[Flatten[(x//.FixNestedLegends) /.{Legended[a_, b___] :> b}] ] ]; *)



Options[WFPlotGrid]={"GridLabels"->{},"GridTranspose"->False};
WFPlotGrid[ass_Association,opt:OptionsPattern[]]:=
Module[{g=Transpose@GriddedLeaves[ass]},
  If[Head[g[[1,1]]]===Grid,
    Grid[g,Spacings -> Scaled[-0.04]],
    WFPlotGrid[g,"GridLabels"->GridKeys[ass]]
  ]
];

WFPlotGrid[x_List,opt:OptionsPattern[]]:=If[MatrixQ[x]
  ,Legended[PlotGrid1[RemoveLegends[x],opt],UnifyLegends[x]]
  ,Grid[{x},Spacings -> Scaled[-0.04]]
];




(* WFGrid[x_List,opt :OptionsPattern[]]:=Legended[
	PlotGrid0[
		Prepend[
		MapThread[
			Prepend[#1,#2]&, 
			{
				RemoveLegends[x],
				COptionValue[WFGrid,"RowLabels"]
			} 
		], 
		Prepend[COptionValue[WFGrid,"ColumnLabels"] ,""]
		]
	], 
	UnifyLegends[x]
]; *)



(* ImportExport`RegisterImport["QSF-wf",WFLoad]; *)
(* ResourceFunction["RegisterFormat"]["QSF-wf", <|"Extensions" -> Table[StringTemplate["psi``"][i], {i, 1, 16}]|>] *)
End[];

EndPackage[];