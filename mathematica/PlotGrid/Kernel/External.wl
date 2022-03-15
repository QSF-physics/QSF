BeginPackage["PlotGrid`External`"];

PlotGrid;
GraphicsInformation;

Begin["`Private`"];

ListLookup[expr_, pos_, def_, h_ : Identity] :=
 h @@ Quiet[Check[Extract[expr, pos, HoldComplete], {def}]]
ListLookup[expr_, pos : {___List}, def_, h_ : Identity] :=
 
 ListLookup[expr, #, def, h] & /@ pos
ListLookup[pos_, def_][expr_] :=
 ListLookup[expr, pos, def]

ReverseY[{x_, y_}] :=
 {x, Reverse@y}

AccumulateShifts[vals_, shifts_] :=
 
 ReverseY[
  FoldList[Plus, 0, #] & /@ (Most /@ ReverseY@vals + ReverseY@shifts)]

XYLookup[{i_, j_}][{x_, y_}] :=
 {x[[j]], y[[i]]}

$SidePositions = <|Left -> {1, 1}, Right -> {1, 2}, Bottom -> {2, 1}, 
   Top -> {2, 2}|>;
$PositionSides = First /@ PositionIndex@$SidePositions;

MapAtHead[f_][h_[args___]] :=
 MapAtHead[f][h][args]
MapAtHead[f_][h_] :=
 f[h]

HeadPattern[p_, n_] :=
 Flatten[p | Through@HeadPattern[p, n - 1][___]]
HeadPattern[p_, 0] :=
 Alternatives[p]

Attributes[ThreadOver] = {HoldAll};

(f : ThreadOver[_, _, _, pat_])[struct_] :=
 
 f /@ struct /; MatchQ[struct, pat[___]]
ThreadOver[i_, _, src_, _][_] :=
 src[[++i]] /; i < Length@src
ThreadOver[_, miss_, _, _][_] :=
 (miss++; Missing["NotAvailable"])

ApplyStructure::tooFew = 
  "Missing `` elements while reshaping `` to the structure of ``.";
ApplyStructure::tooMany = 
  "`` elements discared while reshaping `` to the structure of ``.";

ApplyStructure[tmpl_][list_] :=
 ApplyStructure[list, tmpl]
ApplyStructure[list_, tmpl_] :=
 
 ApplyStructure[list, tmpl, Head@tmpl]
ApplyStructure[list_, tmpl_, head_] :=
 Module[
  {i = 0, miss = 0, res},
  res = ThreadOver[i, miss, list, head][tmpl];
  Which[
   miss > 0,
   ResourceFunction["ResourceFunctionMessage"][ApplyStructure::tooFew,
     miss, list, tmpl],
   i < Length@list,
   ResourceFunction["ResourceFunctionMessage"][
    ApplyStructure::tooMany, Length@list - i, list, tmpl]
   ];
  res
  ]

(* GraphicsInformation code adapted from Carl Woll's answer here: \
https://mathematica.stackexchange.com/a/138907/36508*)

ToNotebook[gr_] :=
 Notebook[
  {
   Cell[
    BoxData@ToBoxes@InstrumentGraphics[ExtractGraphics /@ gr],
    "Output",
    ImageSizeMultipliers -> {1, 1}
    ]
   },
  WindowSize -> AbsoluteCurrentValue[EvaluationNotebook[], WindowSize],
  Evaluator -> AbsoluteCurrentValue[EvaluationNotebook[], Evaluator]
  ]

SowRange[k_, dir_] :=
 (Sow[dir -> {##}, k]; None) &

InstrumentGraphics[gr : {__Graphics}] :=
 MapIndexed[
  Show[
    #1,
    GridLines -> {SowRange[#2, "x"], SowRange[#2, "y"]},
    Epilog -> {
      Annotation[Rectangle[Scaled@{0, 0}, Scaled@{1, 1}], #2, "pr"],
      Annotation[
       Rectangle[ImageScaled@{0, 0}, ImageScaled@{1, 1}], #2, "is"]
      }
    ] &,
  gr
  ]

GraphicsObj = _Graphics | _Legended;

$NullMarker;
If[$Notebooks,
  GraphicsInformation[gr : {GraphicsObj ..}] :=
  Query[Transpose][
    <|
      ImagePadding -> Abs[#is - #pr],
      ImageSize -> Abs[Subtract @@@ #is],
      "PlotRangeSize" -> Abs[Subtract @@@ #pr],
      PlotRange -> {#x, #y}
      |> & /@
    Values@KeySort@<|
        Last@Reap[
          Sow[#[[2]] -> #2, #[[1]]] & @@@
          Cases[
            "Regions" /. FrontEndExecute@ExportPacket[
              ToNotebook[gr],
              "BoundingBox",
              Verbose -> True
              ],
            {{_, "is" | "pr"}, _}
            ],
          _,
          # -> <|#2|> &
          ]
        |>
    ];
  GraphicsInformation[gr : GraphicsObj | Null] := First /@ GraphicsInformation[{gr}];

  GraphicsInformation[gr_List] :=
  With[
    {
    flat = DeleteCases[Null]@Flatten@gr,
    struct = gr /. {
        Null -> $NullMarker[],
        List -> List,
        Except[_List] -> 0
        }
    },
    (
        ApplyStructure[#, struct, 
          List | $NullMarker] /. $NullMarker[] -> Null
        ) & /@ GraphicsInformation[flat] /;
    
    MatchQ[flat, {GraphicsObj ..}]
    ];
,
  SetAttributes[GInformation, {Listable}];
  GInformation[gr : GraphicsObj] := 
    Module[{ass, padder}, 
    ass = 
      Association@
      AbsoluteOptions[
        gr, {ImagePadding, ImageSize, PlotRange, PlotRangePadding}];
    padder = Map[{-First[#], Last[#]} &, ass[PlotRangePadding], {-2}];
    ass[PlotRange] = ass[PlotRange] + padder;
    ass["PlotRangeSize"] = (ass[ImageSize] - 
        Total[(ass[ImagePadding] + 2 ass[PlotRangePadding]), {-1}]);
    ass
  ];
  GInformation[gr:Null]:=<|ImagePadding->Null,ImageSize->Null,PlotRange->Null|>;
  GraphicsInformation[gr_List]:=Merge[Map[GraphicsInformation ,gr],Identity];
  GraphicsInformation[gr:GraphicsObj]:=GInformation@gr;
  GraphicsInformation[gr:Null]:=GInformation[Null];
];
Options[ApplyToWrapped] = {Method -> Blank};

SyntaxInformation[
   ApplyToWrapped] = {"ArgumentsPattern" -> {_, _, _, _., 
     OptionsPattern[]}};

ApplyToWrapped::noMatch = 
  "Expression `` is not a wrapped expression matching ``.";

ApplyToWrapped[func_, expr_, target_, extract_ : None, 
  Longest[OptionsPattern[]]] :=
 ReleaseHold[
  Hold@IApplyToWrapped[expr, target, extract, {}] //.
    
    DownValues@IApplyToWrapped /. {
    IApplyToWrapped[e_, c_] :>
     With[
      {
       r = Which[
         extract === None,
         func@e,
         OptionValue@Method === Blank,
         func[e, c],
         True,
         With[
          {
           
           cf = 
            Unevaluated@c /. 
             w_[Verbatim@_, rest___] :> (w[#, rest] &)
           },
          func[e, cf]
          ]
         ]
       },
      r /; True
      ],
    _Hold?(MemberQ[#, 
         HoldPattern@
          IApplyToWrapped[_, _, _, _], {0, \[Infinity]}] &) :>
     Hold[
      ResourceFunction["ResourceFunctionMessage"][
       ApplyToWrapped::noMatch, HoldForm@expr, target];
      expr
      ]
    }
  ]

Attributes[IApplyToWrapped] = {HoldFirst};

IApplyToWrapped[expr_, target_, extract_, coll_] /;
  
  MatchQ[Unevaluated@expr, target] :=
 IApplyToWrapped[expr, coll]
IApplyToWrapped[expr : w_[wrapped_, args___], target_, 
   extract : Except[None], {coll___}] /;
  
  MatchQ[Unevaluated@expr, extract] :=
 
 IApplyToWrapped[wrapped, target, extract, {coll, w[_, args]}]
IApplyToWrapped[w_[wrapped_, args___], target_, extract_, coll_] :=
 
 w[IApplyToWrapped[wrapped, target, extract, coll], args]

ThreadWrappers[expr_, wrap_, nmax_ : Infinity] :=
 
 iThreadWrappers[expr, wrap, nmax]

iThreadWrappers[expr_, _, 0] :=
 expr
iThreadWrappers[expr_List, wrap_, nmax_] :=
 
 iThreadWrappers[#, wrap, nmax - 1] & /@ expr
iThreadWrappers[h_[expr_List, args___], wrap_, nmax_] /; 
  MatchQ[h, wrap] :=
 
 iThreadWrappers[h[#, args], wrap, nmax - 1] & /@ expr
iThreadWrappers[h_[expr_, args___], wrap_, nmax_] /; MatchQ[h, wrap] :=

  With[
  {proc = iThreadWrappers[expr, wrap, nmax]},
  iThreadWrappers[h[proc, args], wrap, nmax] /; proc =!= expr
  ]
iThreadWrappers[expr_, _, _] :=
 expr

ValidGraphicsQ[Legended[expr_, __]] :=
 ValidGraphicsQ[expr]
ValidGraphicsQ[_Graphics] :=
 True
ValidGraphicsQ[_] :=
 False

ExtractGraphics[gr_Graphics] :=
 gr
ExtractGraphics[Legended[expr_, __]] :=
 expr

GraphicsOpt[g_Graphics, opt_] :=
 GraphicsOpt[Options@g, opt]
GraphicsOpt[l_Legended, opt_] :=
 GraphicsOpt[ExtractGraphics@l, opt]
GraphicsOpt[opts_, opt_] :=
 OptionValue[Graphics, opts, opt]

NormalizeFrameSetting[val_] :=
 Replace[
  Replace[
   val,
   {None -> False, Automatic -> True},
   {2}
   ],
  Except@Table[True | False, 2, 2] -> Table[False, 2, 2]
  ]

NormalizeGraphicsOpt[_][s : {{_, _}, {_, _}}] :=
 s
NormalizeGraphicsOpt[FrameLabel][b_] :=
 {{None, None}, {b, None}}
NormalizeGraphicsOpt[FrameLabel][{b_, l_ : None, t_ : None, 
   r_ : None, ___}] :=
 {{l, r}, {b, t}}
NormalizeGraphicsOpt[FrameTicks][a_] :=
 {{a, a}, {a, a}}
NormalizeGraphicsOpt[FrameTicks][{h_, v_}] :=
 {{v, v}, {h, h}}
NormalizeGraphicsOpt[FrameStyle][
  a : Except[_List]] :=
 {{a, a}, {a, a}}
NormalizeGraphicsOpt[FrameStyle][{v_, h_}] :=
 {{v, v}, {h, h}}
NormalizeGraphicsOpt[FrameStyle][{b_, l_, t_, r_}] :=
 {{l, r}, {b, t}}
NormalizeGraphicsOpt[FrameStyle][_] :=
 
 NormalizeGraphicsOpt[FrameStyle][None]
NormalizeGraphicsOpt[FrameTicksStyle] :=
 
 NormalizeGraphicsOpt[FrameStyle]
NormalizeGraphicsOpt[Frame][a_] :=
 
 NormalizeFrameSetting@{{a, a}, {a, a}}
NormalizeGraphicsOpt[Frame][{a_}] :=
 
 NormalizeFrameSetting@NormalizeGraphicsOpt[FrameStyle][True]
NormalizeGraphicsOpt[Frame][{h_, v_}] :=
 
 NormalizeFrameSetting@{{v, v}, {h, h}}
NormalizeGraphicsOpt[Frame][{b_, l_, t_, r_ : True}] :=
 
 NormalizeFrameSetting@{{l, r}, {b, t}}
NormalizeGraphicsOpt[PlotRange][m_?NumericQ] :=
 {{-m, m}, {-m, m}}
NormalizeGraphicsOpt[PlotRange][s_] :=
 {{s, s}, {s, s}}
NormalizeGraphicsOpt[PlotRange][{min_, max_}] :=
 {{min, 
   max}, {Automatic, Automatic}}
NormalizeGraphicsOpt["Custom"][a_] :=
 {{a, a}, {a, a}}
NormalizeGraphicsOpt["Custom"][{h_, v_}] :=
 {{h, h}, {v, v}}

NormalizedOptionValue[g_, opt_List] :=
 
 NormalizedOptionValue[g, #] & /@ opt
NormalizedOptionValue[g_, opt_] :=
 
 NormalizeGraphicsOpt[opt][GraphicsOpt[g, opt]]

RotateGraphicsOpt[_, ang_][s : {{_, _}, {_, _}}] /; 
  Mod[ang, 2 \[Pi]] == 0 :=
 s
RotateGraphicsOpt[opt_, ang_][{{l_, r_}, {b_, t_}}] /; 
  Mod[ang, \[Pi]/2] == 0 :=
 
 RotateGraphicsOpt[opt, ang - \[Pi]/2][{{t, b}, {l, r}}]
RotateGraphicsOpt[opt : PlotRange, ang_][{{l_, r_}, {b_, t_}}] /; 
  Mod[ang, \[Pi]/2] == 0 :=
 RotateGraphicsOpt[opt, ang - \[Pi]/2][
  {Replace[{t, b}, lim_?NumericQ :> -lim, 1], {l, r}}
  ]

NormalizeTickSpec[tick : {_, _, {_, _}, _}] :=
 tick
NormalizeTickSpec[{x_, lbl_, len_, sty_}] :=
 {x, lbl, {len, len}/2, 
  sty}
NormalizeTickSpec[{x_, lbl_, len_ : 0.00625}] :=
 
 NormalizeTickSpec@{x, lbl, len, {}}
NormalizeTickSpec[x_] :=
 NormalizeTickSpec@{x, x}

$BuiltinTickSpecs = Automatic | None | False | True | All;

HideTickLabels[t : $BuiltinTickSpecs] :=
 t
HideTickLabels[ticks_List] :=
 {#, Spacer@0, ##3} & @@@ 
  NormalizeTickSpec /@ ticks
HideTickLabels[f_][lims__] :=
 HideTickLabels[f@lims]

ScaleFrameTicks[_][t : $BuiltinTickSpecs] :=
 t
ScaleFrameTicks[s_][ticks_List] :=
 {#, #2, s #3, #4} & @@@ 
  NormalizeTickSpec /@ ticks
ScaleFrameTicks[s_][fun_] :=
 ScaleFrameTicks[s, fun]
ScaleFrameTicks[s_, fun_][lims__] :=
 ScaleFrameTicks[s][fun@lims]

LoadScaledTicks[] := (
  Charting`ScaledTicks[{Identity, Identity}][0, 1];
  LoadScaledTicks[] = Null;
  )

CustomTicks[limits__] :=
 With[
  {
   rLimits = 
    Round[{limits}, 10.^(Round@Log10[-Subtract[limits]] - 4)]
   },
  Replace[
   LoadScaledTicks[];
   Internal`InheritedBlock[
    {Charting`SimplePadding, Charting`CommonDump`simplePadding},
    (*Fix for 0 tick label not having same number of decimal places \
as other labels in Charting`ScaledTicks*)
    
    (*Put any zeroes from the list of zeros to the list of medium \
numbers This causes zeros to be printed with the same number of \
decimal places as those numbers simplePadding is patched below to \
properly handle zeros Since simplePadding[
    "Medium",{0.},...] returns {0},
    this doesn't break the case where there are no medium numbers*)
  
      Unprotect@Charting`SimplePadding;
    DownValues[Charting`SimplePadding] = 
     DownValues[
       Charting`SimplePadding] /. {HoldPattern[
         z : Charting`CommonDump`zero = _] :> (z = {}), 
       HoldPattern[
         m : Charting`CommonDump`medium = Select[l_, c_ &]] :> (m = 
          Select[l, c || PossibleZeroQ@# &])};
    (*Chop 0. to 0 in MantissaExponent calls in simplePadding[
    "Medium",...] to prevent errors/broken formatting*)
    DownValues[Charting`CommonDump`simplePadding] = 
     DownValues[Charting`CommonDump`simplePadding] /. {
       HoldPattern[MantissaExponent[l_, 10]] :> 
        MantissaExponent[Chop@l, 10],
       HoldPattern@NumberForm[#, {Infinity, d_}] :> 
        NumberForm[N@#, {Infinity, d}]
       };
    Charting`ScaledTicks[{Identity, Identity}] @@ rLimits
    ],
   {
    {_?(Not@*Between[{limits}]), _Spacer, __} :> 
     Nothing, {x_, lbl : Except@_Spacer, rest__} :>
     {Clip[x, Sort@{limits}], ToString[lbl, TraditionalForm], rest}
    },
   1
   ]
  ]

CustomTicksNoLabels[limits__] := {#, Spacer@0, ##3} & @@@ 
  CustomTicks[limits]

(* replace Scaled with "Scaled" to prevent issues with Graphics \
typesetting *)
"Scaled";

HideOuterLabels[_, _, _, None | False] := None
HideOuterLabels[_, _, {-\[Infinity], -\[Infinity]}, t_] := t
HideOuterLabels[{_, 2}, size_, thres_, t : Automatic | True] :=
 t
HideOuterLabels[side_, size_, thres_, Automatic | True | All] :=
 
 HideOuterLabels[side, size, thres, CustomTicks]
HideOuterLabels[side_, size_, thres_, ticks_List][min_, max_] :=
 With[
  {t = (max - min) Replace[
      thres,
      {
       "Scaled"@thr_ :> thr ,
       thr_?NumericQ :> thr/size[[3 - side[[1]]]]
       },
      1
      ]
   },
  Replace[
   NormalizeTickSpec /@ 
    ticks, {x_, lbl_, rest__} /; ! 
      min + t[[1]] <= x <= max - t[[2]] :> {x, Spacer@0, rest},
   1
   ]
  ]
HideOuterLabels[side_, size_, thres_, f_][lims__] :=
 
 HideOuterLabels[side, size, thres, f@lims][lims]

EnsureCustomTicks[_, None | False] := None
EnsureCustomTicks[_, All] := CustomTicks
EnsureCustomTicks[{_, 1}, Automatic | True] := CustomTicks
EnsureCustomTicks[{_, 2}, Automatic | True] := CustomTicksNoLabels
EnsureCustomTicks[_, t_] := t

$ModifiedDefaultTickPattern = 
  Hold[EnsureCustomTicks][_, HideOuterLabels[_, _, _, CustomTicks]];

$CustomTicks = 
  HideTickLabels | ScaleFrameTicks | CustomTicks | 
   CustomTicksNoLabels | HideOuterLabels;

$CustomTicksNames = SymbolName /@ $CustomTicks;

$CustomTicksAssoc = 
  AssociationThread[List @@ $CustomTicksNames -> List @@ $CustomTicks];

TransformCustomTicksHeads[f_, p_, t_] /; MatchQ[t, p] :=
 
 MapAtHead[f]@Replace[
   t,
   h_[args___, t2 : p] :>
    h[args, TransformCustomTicksHeads[f, p, t2]]
   ]

(* force the resource function to be loaded before evaluating ticks \
by redirecting ticks functions over ResourceFunction["PlotGrid"] *)

InjectCustomTicksLoading[expr_] :=
 
 expr /.
  t : HeadPattern[$CustomTicks, 1] :>
   With[
    (* convert tick symbols to strings to ensure they work across \
different versions of PlotGrid *)
    {f = 
      TransformCustomTicksHeads[SymbolName, 
       HeadPattern[$CustomTicks, 1], t]},
    ResourceFunction["PlotGrid"][f][##] &
    ]

PlotGrid[t : HeadPattern[$CustomTicksNames, 1]][args___] :=
 
 TransformCustomTicksHeads[$CustomTicksAssoc, 
   HeadPattern[$CustomTicksNames, 1], t][args]

ClipFrameLabels[gr_?ValidGraphicsQ, sides_List] :=
 With[
  {
   drop = Complement[Keys@$SidePositions, sides] /. $SidePositions,
   frame = 
    Replace[NormalizedOptionValue[gr, Frame], None -> False, {2}],
   frameTicks = 
    Replace[NormalizedOptionValue[gr, FrameTicks], 
     None -> False, {2}]
   },
  Show[
   gr,
   FrameLabel -> 
    ReplacePart[NormalizedOptionValue[gr, FrameLabel], drop -> None],
   FrameTicks -> MapAt[HideTickLabels, frameTicks, drop],
   FrameTicksStyle -> MapIndexed[
     Which[
       (* only apply styling if the frame is shown, 
       to work around the bug that styles from hidden frame sides are \
inherited by the following frame sides *)
       
       TrueQ[! Extract[frame, #2]] || TrueQ[! Extract[frameTicks, #2]],
       None,
       MemberQ[drop, #2] && 
        MatchQ[Extract[frameTicks, #2], $BuiltinTickSpecs],
       Directive[FontSize -> 0, FontOpacity -> 0, #],
       True,
       #
       ] &,
     NormalizedOptionValue[gr, FrameTicksStyle],
     {2}
     ]
   ]
  ]
ClipFrameLabels[expr_, _] := expr

LookupNeighbors[grid_, {i_, j_}, {h_, w_}] :=
 Map[
  ListLookup[grid, {i, j} + # & /@ Thread@# /. (0 -> Null), Null] &,
  {
   {
    {Range[0, h], -1},
    {Range[0, h], 1 + w}
    }, {
    {1 + h, Range[0, w]},
    {-1, Range[0, w]}
    }
   },
  {2}
  ]

PlotsPresent[grid_, pos_, dim_] :=
 With[
  {neigh = LookupNeighbors[grid, pos, dim]},
  MemberQ[Extract[neigh, #], Except@Null] & /@ $SidePositions
  ]

LookupSpacings[spacings_, {i_, j_}, {h_, w_}] :=
 
 ListLookup[spacings, # /. (0 -> Null), Missing[]] & /@
  <|
   
   Left -> {1, j - 1}, Right -> {1, j + w}, Bottom -> {2, i + h}, 
   Top -> {2, i - 1}
   |>

SpacingsPresent[spacings_, pos_, dim_] :=
 # =!= 0 & /@ 
  LookupSpacings[spacings, pos, dim]

ApplyFrameSettings[pos : {i_, j_}, dims_, showFrameLabels_, 
  hideOuterLabels_, mergeAxes_, effSize_, grid_, spacings_] := 
 Module[
  {
   gr = grid[[i, j]],
   plotsPresent = PlotsPresent[grid, pos, dims],
   spacingsPresent = SpacingsPresent[spacings, pos, dims],
   autoValues,
   thresholds
   },
  (
    autoValues = <|
       Automatic -> ! plotsPresent[#],
       Full -> ! plotsPresent[#] || spacingsPresent[#]
       |> &;
    thresholds = MapIndexed[
      If[#[[1]] && (#[[2]] /. Not /@ autoValues@$PositionSides@#2),
        #[[3]],
        -Infinity
        ] &,
      hideOuterLabels,
      {2}
      ];
    gr = ClipFrameLabels[
      Show[
       grid[[i, j]],
       Frame -> MapThread[
         If[#2 =!= False, False, #] &,
         {
          NormalizedOptionValue[gr, Frame],
          mergeAxes
          },
         2
         ]
       ],
      If[
         #2 /. autoValues@#,
         #,
         Nothing
         ] & @@@ showFrameLabels
      ];
    Show[
     gr,
     FrameTicks -> MapIndexed[
       Hold[EnsureCustomTicks][
         #2,
         HideOuterLabels[#2, effSize, thresholds[[3 - #2[[1]]]], #]
         ] &,
       NormalizedOptionValue[gr, FrameTicks],
       {2}
       ]
     ]
    ) /; ValidGraphicsQ@gr
  ]
ApplyFrameSettings[{i_, j_}, _, _, _, _, _, grid_, _] := grid[[i, j]]

(* based on Legending`LegendDump`placedHelper - used to determine \
whether legends will be placed inside or outside the plot *)

iExplicitPosPat = _?NumericQ | (Scaled | ImageScaled)[_?NumericQ];
explicitPosPat = 
  iExplicitPosPat | Offset[_?NumericQ] | 
   Offset[_?NumericQ, iExplicitPosPat];

normalizeExplicitPos[(t : Scaled | ImageScaled)[p_?NumericQ]] :=
 {t, 
  0, p}
normalizeExplicitPos[p_?NumericQ] :=
 {Scaled, 0, p}
normalizeExplicitPos[Offset[o_?NumericQ]] :=
 {Scaled, o, 0}
normalizeExplicitPos[Offset[o_?NumericQ, p : iExplicitPosPat]] :=
 
 ReplacePart[normalizeExplicitPos[p], 2 -> o]

posPat = 
  Left | Right | Top | Bottom | Center | After | "After" | Before | 
   "Before" | Above | Below;
xposPat = 
  Before | Left | Center | Right | After | "Axis" | explicitPosPat;
yposPat = 
  Below | Bottom | Center | Top | Above | "Axis" | explicitPosPat;

eposPat = Left | Right | Top | Bottom | Center;
exposPat = Left | Center | Right | _?NumericQ;
eyposPat = Bottom | Center | Top | _?NumericQ;

placedHelper[pos : (Scaled | ImageScaled)[_List]] :=
  
  placedHelper@Thread@pos;
placedHelper[(Scaled | ImageScaled)[pos_]] :=
  placedHelper@pos;
placedHelper[
   pos : 
    Offset[_List, _List | (Scaled | ImageScaled)[_List] | 
      PatternSequence[]]] :=
  
  placedHelper@Thread@MapAt[Thread, 2]@pos;
placedHelper[{xpos : (Scaled | ImageScaled)[_List], (Scaled | 
        ImageScaled)[ypos_] | ypos_List}] :=
  
  placedHelper@{Thread@xpos, ypos};
placedHelper[{(t1 : Scaled | ImageScaled)[xpos_] | 
      xpos_List, (t2 : Scaled | ImageScaled)[ypos_] | 
      ypos_List}] /;
   {t1, t2} =!= {} :=
  placedHelper@{xpos, ypos};
placedHelper[pos : posPat] :=
  {posHelper@namedPlaceHelper@pos, 
   eposHelper@namedPlaceHelper@pos};
placedHelper[pos : {xposPat, yposPat}] :=
  {posHelper@pos, 
   automaticEposHelper@pos};
placedHelper[{pos : {xposPat, yposPat}}] :=
  {posHelper@pos, 
   automaticEposHelper@pos};
placedHelper[{pos : posPat, epos : eposPat}] :=
  {posHelper@
    namedPlaceHelper@pos, eposHelper@namedPlaceHelper@epos};
placedHelper[{pos : {xposPat, yposPat}, 
    epos : eposPat}] :=
  {posHelper@pos, 
   eposHelper@namedPlaceHelper@epos};
placedHelper[{pos : posPat, 
    epos : {exposPat, eyposPat}}] :=
  {posHelper@
    namedPlaceHelper@pos, eposHelper@epos};
placedHelper[{pos : {xposPat, yposPat}, 
    epos : {exposPat, eyposPat}}] :=
  {posHelper@pos, 
   eposHelper@epos};
placedHelper[pos_] :=
  {posHelper@namedPlaceHelper@After, 
   eposHelper@namedPlaceHelper@Left};

namedPlaceHelper[pos_] :=
  Switch[pos,
   Right | "After" | After, {After, Center},
   Left | "Before" | Before, {Before, Center},
   Center, {Center, Center},
   Bottom | Below, {Center, Below},
   Top | Above, {Center, Above},
   _, {After, Center}
   ];

$defOffset = 7;
posHelper[{xpos_, ypos_}] :=
  normalizeExplicitPos /@ {
    Switch[xpos,
     After, Offset[$defOffset, 1],
     Right, Offset[-$defOffset, 1],
     Before, Offset[-$defOffset, 0],
     Left, Offset[$defOffset, 0],
     Center, 0.5,
     explicitPosPat, xpos,
     "Axis", 1,
     _, Offset[$defOffset, 1]
     ],
    Switch[ypos,
     Below, Offset[-$defOffset, 0],
     Bottom, Offset[$defOffset, 0],
     Top, Offset[-$defOffset, 1],
     Above, Offset[$defOffset, 1],
     Center, 0.5,
     explicitPosPat, ypos,
     "Axis", 1,
     _, 0.5
     ]
    };

eposHelper[{expos_, eypos_}] :=
  {
   Switch[expos,
    After, 0,
    Right, 1,
    Before, 1,
    Left, 0,
    Center, 0.5,
    _?NumericQ, expos,
    _, 0.5
    ],
   Switch[eypos,
    Center, 0.5,
    Below, 1,
    Bottom, 0,
    Top, 1,
    Above, 0,
    _?NumericQ, eypos,
    _, 0.5
    ]
   };

coordSign[_ | PatternSequence[], o_, x_] /; 0 < x < 1 :=
 0
coordSign[_ | PatternSequence[], o_, x_] /; 1 < x || x < 0 :=
 Sign@x
coordSign[_ | PatternSequence[], o_, x_] /; x == 0 :=
 
 If[o > 0, -0.5, -1]
coordSign[_ | PatternSequence[], o_, x_] /; x == 1 :=
 
 If[o < 0, 0.5, 1]
coordSign[{_ | PatternSequence[], o_, x_}] :=
 coordSign[o, x]

coordInsideQ[_ | PatternSequence[], o_, x_] :=
 -0.5 <= 
  coordSign[o, x] <= 0.5
coordInsideQ[{_ | PatternSequence[], o_, x_}] :=
 coordInsideQ[o, x]

automaticEposHelper[{xpos_, ypos_}] :=
  {
   Switch[xpos,
    After, 0,
    Right, 1,
    Before, 1,
    Left, 0,
    Center, 0.5,
    explicitPosPat, Switch[coordSign@normalizeExplicitPos[xpos],
     0, 0.5,
     1 | -0.5, 0,
     -1 | 0.5, 1
     ],
    _, 0.5
    ],
   Switch[ypos,
    Center, 0.5,
    Below, 1,
    Bottom, 0,
    Top, 1,
    Above, 0,
    explicitPosPat, Switch[coordSign@normalizeExplicitPos[ypos],
     0, 0.5,
     1 | -0.5, 0,
     -1 | 0.5, 1
     ],
    _, 0.5
    ]
   };

LegendInsideQ[Placed[
   _,
   pos_,
   ___
   ]] :=
 AllTrue[
  First@placedHelper@pos,
  Apply[coordInsideQ]
  ]
LegendInsideQ[_] :=
 False

Options[ResolveCoordinatesTool] = {"CopiedValueFunction" -> Identity, 
   "DisplayFunction" -> Automatic};

ResolveCoordinatesTool[Automatic][_] :=
 Identity
ResolveCoordinatesTool[OptionsPattern[]][
  cvf : "CopiedValueFunction"] :=
 
 Replace[OptionValue@cvf, Automatic -> Identity]
ResolveCoordinatesTool[OptionsPattern[]][df : "DisplayFunction"] :=
 
 Replace[OptionValue@df, 
  Automatic -> OptionValue@"CopiedValueFunction"]

(*get CoordinateToolOptions from a Graphics object. Since the setting \
can be directly inside the Graphics or inside Method,we need to \
combine both (and handle Automatic settings for both \
Method,CoordinateToolOptions and DisplayFunction along the way*)

GetCoordinatesToolOptions[gr_] :=
 
 ResolveCoordinatesTool@GraphicsOpt[
   {
    FilterRules[
     Replace[GraphicsOpt[gr, Method], Automatic -> {}],
     CoordinatesToolOptions
     ],
    Options@gr
    },
   CoordinatesToolOptions
   ]

implicitAutomatic;
unspecified;

LastExplicit[specs_] := 
 FirstCase[Reverse@specs, Except@implicitAutomatic, implicitAutomatic]

ExpandSeqSpec[{spec_, rules : {(_Integer -> _) ..}}, n_, m_] :=
 
 MapThread[
  If[#2 === unspecified,
    #,
    m@{##}
    ] &,
  {
   ExpandSeqSpec[spec, n, m],
   ReplacePart[Table[unspecified, n], rules]
   }
  ]
ExpandSeqSpec[{spec_, rule : (_Integer -> _)}, n_, m_] :=
 
 ExpandSeqSpec[{spec, {rule}}, n, m]
ExpandSeqSpec[rules : {(_Integer -> _) ..}, n_, m_] :=
 
 ExpandSeqSpec[{{}, rules}, n, m]
ExpandSeqSpec[rule : (_Integer -> _), n_, m_] :=
 
 ExpandSeqSpec[{{}, {rule}}, n, m]
ExpandSeqSpec[{start___, cycle : {__}, end___}, n_, _] :=
 With[
  {
   rem = Max[n - Length@{start}, 0]
   },
  Join[
   Take[{start}, UpTo@n],
   PadRight[{}, Max[rem - Length@{end}, 0], cycle],
   Take[{end}, -Min[Length@{end}, rem]]
   ]
  ]
ExpandSeqSpec[{start___, {}, end___}, n_, m_] :=
 
 ExpandSeqSpec[{start, {implicitAutomatic}, end}, n, m]
ExpandSeqSpec[{start___}, n_, m_] :=
 
 ExpandSeqSpec[{start, {}}, n, m]
ExpandSeqSpec[spec_, n_, m_] :=
 ExpandSeqSpec[{{spec}}, n, m]

Expand2DSpec[{wspec_ : {}, hspec_ : {}, ___}, {m_, n_}, 
  merger_] :=
 {ExpandSeqSpec[wspec, m, merger], 
  ExpandSeqSpec[hspec, n, merger]}
Expand2DSpec[spec_, {m_, n_}, merger_] :=
 
 Expand2DSpec[{spec, spec}, {m, n}, merger]

ExpandGridSpec[spec_, {m_, n_}, merger_] :=
 Module[
  {wspec, hspec},
  {wspec, hspec} = Expand2DSpec[spec, {m, n}, merger];
  MapThread[
   merger@*List,
   {
    Transpose@Table[hspec, m],
    Table[wspec, n]
    },
   2
   ]
  ]
ExpandGridSpec[{wspec_, hspec_, posRule_Rule}, {m_, n_}, merger_] :=
 
 ExpandGridSpec[{wspec, hspec, {posRule}}, {m, n}, merger]
ExpandGridSpec[{wspec_, hspec_, posRules : {__Rule}}, {m_, n_}, 
  merger_] :=
 Module[
  {
   grid = ExpandGridSpec[{wspec, hspec}, {m, n}, merger]
   },
  (grid = MapAt[val |-> merger@{val, #2}, grid, Span @@@ #]) & @@@ 
   posRules;
  grid
  ]

ExpandFrameLabelSpec[Directive[spec___]] := 
 ExpandFrameLabelSpec[{spec}]
ExpandFrameLabelSpec[{def : Except[_Rule] : implicitAutomatic, 
   rules___Rule}] :=
 # -> (Lookup[<|rules|>, #, def] /. {All -> True,
        None -> False}) & /@ {Left, Right, Bottom, Top}
ExpandFrameLabelSpec[
  s : Automatic | Full | True | False | All | None | 
    implicitAutomatic] :=
 ExpandFrameLabelSpec[{s}]
ExpandFrameLabelSpec[{h_, v_}] :=
 {Left -> v, Right -> v, 
  Bottom -> h, Top -> h}
ExpandFrameLabelSpec[{{l_, r_}, {b_, t_}}] :=
 {Left -> l, Right -> r,
   Bottom -> b, Top -> t}

$FrameLabelMerger = 
  Map[ExpandFrameLabelSpec]/*Merge[LastExplicit]/*Normal;

SidedFrameOptionFromWrappers[wrappers_, opt_] :=
 Map[
  Cases[
   Frame[_, opts__] :>
    Replace[
     opt /. {opts},
     {
      opt -> Nothing,
      val_ :> NormalizeGraphicsOpt["Custom"]@val
      }
     ]
   ],
  wrappers,
  {2}
  ]

FillInStructure[structure_, fromWrap_] :=
 Module[
  {grid = fromWrap, h},
  (grid[[Sequence @@ Span @@@ Transpose@Accumulate@{##}]] = 
      h[Extract[fromWrap, #]]) & @@@ structure;
  grid /. h -> Identity
  ]

ListToPerSideLists[{}] :=
 Table[{}, 2, 2]
ListToPerSideLists[l_] :=
 Transpose[l, {3, 1, 2}]

ApplyStructureToOption::usage = 
  "ApplyStructureToOption[structure,opt,fromWrap,m] merges the \
relevant settings for a sided property based on the structure of the \
grid, where settings come from a Spacings-like option and per-side \
wrappers";
ApplyStructureToOption[structure_, opt_, fromWrap_, m_] :=
 Module[
  {
   h,
   filledFromWrap = FillInStructure[
     structure,
     Map[
      ListToPerSideLists,
      fromWrap,
      {2}
      ]
     ]
   },
  MapThread[
     m@*Prepend,
     {
      Transpose[
          {
            Extract[filledFromWrap, #],
            MapIndexed[
             Extract[#, Prepend[#2 {1, -1} + {0, 3}, All]] &,(* 
             get settings for opposite sides of neighbors *)
         
                 
             Replace[LookupNeighbors[filledFromWrap, ##], 
              Null -> Nothing, {3}],
             {2}
             ]
            } // MapAt[Map[h, #, {3}] &, {{1}, {2, All}}],
          {3, 1, 2}
          ] //
         MapAt[Reverse, {{1, 2}, {2, 1}}] //
        
        Map@Map@Flatten // ReplaceAll[h -> Identity],
      {opt[[1, #[[2]] + {0, #2[[2]] + 1}]], 
       opt[[2, #[[1]] + {#2[[1]] + 1, 0}]]}
      },
     2
     ] & @@@ structure
  ]

$OuterLabelModeSpec = Automatic | Full | False | True | All | None;
$OuterLabelSideSpec = 
  Automatic | Left | Right | Bottom | Top | First | Last | All | 
   None;
$OuterLabelMarginSpec = _?NumericQ | Scaled@_?NumericQ;

iExpandOuterLabelSpec[Directive[spec___]] :=
 
 iExpandOuterLabelSpec[{spec}]
iExpandOuterLabelSpec[implicitAutomatic] :=
 {implicitAutomatic, 
  implicitAutomatic, implicitAutomatic}
iExpandOuterLabelSpec[spec : Except@_List] :=
 
 iExpandOuterLabelSpec[{spec}]
iExpandOuterLabelSpec[spec : {_, _, _}] :=
 spec
iExpandOuterLabelSpec[{
   s : $OuterLabelSideSpec : implicitAutomatic,
   Longest[m : $OuterLabelModeSpec : implicitAutomatic],
   mar : $OuterLabelMarginSpec : implicitAutomatic
   }] :=
 {s, m, mar}

postProcessOuterLabelSpec[hor_, speed_][{s_, m_, mar_}] :=
 {
  Replace[
   s, {Automatic -> If[hor, First, Last], Left | Top -> First, 
    Right | Bottom -> Last, All -> True, None -> False}],
  Replace[m, {All -> True, None -> False}],
  mar /. {Automatic -> If[speed, -\[Infinity], 10], Scaled -> "Scaled"}
  }

ExpandOuterLabelSpec[structure_, spec_, fromWrap_, {m_, n_}, speed_] :=

  ApplyStructureToOption[
     structure,
     Expand2DSpec[
       spec,
       {m, n} - 1,
       Map[iExpandOuterLabelSpec]/*MapThread[LastExplicit@*List]
       ] //
      Map[ArrayPad[#, {{1, 1}}, {{False, False, 0}}] &](* 
     add settings for outer edges *),
     fromWrap,
     Map[
       iExpandOuterLabelSpec]/*(MapThread[LastExplicit@*List, 
         Prepend[{Automatic, Full, Automatic}]@#] &)
     ] //
    
    Query[(* resolve Automatic side spec & None/All *)
     All,
     {
      1 -> Map@postProcessOuterLabelSpec[True, speed],
      2 -> Map@postProcessOuterLabelSpec[False, speed]
      }
     ] //
   
   MapAt[
    Replace@{First -> False, Last -> True}, {{All, 1, 1, 1}, {All, 2, 
      2, 1}}] //(* resolve First/Last specs *)
  
  MapAt[Replace@{First -> True, Last -> False}, {{All, 1, 2, 1}, {All,
      2, 1, 1}}]

$NamedAxesMergeIcons = <|
   "ZigZag" -> (
     Graphics[
       {
        CapForm["Square"],
        Style[
         Line@{{0, 0}, Offset[{-#, 0}, {1, 0}]},
         Antialiasing -> False
         ],
        Line@{Offset[{-#, 0}, {1, 0}], 
          Offset[{-#/2, 0}, {1, -1}], {1, 0}}
        },
       AspectRatio -> Full,
       PlotRange -> {{0, 1}, 1.2 {-1, 1}}
       ] &),
   "Cut" -> (Graphics[
       {
        CapForm["Square"],
        Style[
         Line@{{0, 0}, Offset[{-#/2, 0}, {1, 0}]},
         Antialiasing -> False
         ],
        Line@{Offset[{-#, 0}, {1, -1}], {1, 1}}
        },
       AspectRatio -> Full,
       PlotRange -> {{0, 1}, 1.2 {-1, 1}}
       ] &),
   "Wave" -> (Graphics[
       {
        CapForm["Square"],
        Style[
         Line@{{0, 0}, Offset[{-#, 0}, {1, 0}]},
         Antialiasing -> False
         ],
        (*Circle[Offset[{-#/2,0},{1,0}],Offset[{#/2,0},{0,1}],{-\[Pi],
        0}]*)
        BezierCurve@{Offset[{-#, 0}, {1, 0}], 
          Offset[{-#/2, 0}, {1, 0}], Offset[{-# 3/5, 0}, {1, -1}], 
          Offset[{-# /3, 0}, {1, -1}], 
          Offset[{-# /6, 0}, {1, -1}], {1, 0}}
        },
       AspectRatio -> Full,
       PlotRange -> {{0, 1}, 1.2 {-1, 1}}
       ] &),
   "Dotted" -> (Graphics[
       {
        CapForm["Square"],
        Style[
         Line@{{0, 0}, Offset[{-#, 0}, {1, 0}]},
         Antialiasing -> False
         ],
        AbsoluteDashing@{1, 3},
        Line@{Offset[{-#, 0}, {1, 0}], {1, 0}}
        },
       AspectRatio -> Full,
       PlotRange -> {{0, 1}, 1.2 {-1, 1}}
       ] &),
   "Line" -> (Graphics[
       {
        CapForm["Square"],
        Style[
         Line@{{0, 0}, {1, 0}},
         Antialiasing -> False
         ]
        },
       AspectRatio -> Full,
       PlotRange -> {{0, 1}, 1.2 {-1, 1}}
       ] &)
   |>;
Show[#@40, ImageSize -> {80, 40}] & /@ $NamedAxesMergeIcons

NormalizeMergePrimitive[False | implicitAutomatic] :=
 False
NormalizeMergePrimitive[None] :=
 {{None, 0}, {None, 0}}
NormalizeMergePrimitive[Directive[spec__]] :=
 
 NormalizeMergePrimitive[{spec}]
NormalizeMergePrimitive[spec : Except[{_, _?NumericQ}]] :=
 
 NormalizeMergePrimitive[{spec, 10}]
NormalizeMergePrimitive[{name_String | (name_String -> size_), 
    width_?NumericQ}] /; KeyMemberQ[$NamedAxesMergeIcons, name] :=
 
 ConstantArray[{$NamedAxesMergeIcons[name][First@{size, 10}/2], 
   width}, 2]
NormalizeMergePrimitive[{{left_, right_}, 
   width_?NumericQ}] :=
 {{left, 
   width}, {Replace[right, Automatic -> left], width}}
NormalizeMergePrimitive[{prim_, width_?NumericQ}] :=
 {{prim, 
   width}, {implicitAutomatic, 0}}

ComputeCorners::usage = 
  "ComputeCorners[item] computers the indices of the four corners of \
an element of the structure specification";
ComputeCorners[{pos_, size_}] :=
 Table[
  pos + {i, j} (size + 1) - 1,
  {i, 0, 1},
  {j, 0, 1}
  ]

ComputeCorners::usage = 
  "ComputeCorners[structure,{m,n}] finds all points in the where four \
corners join, and returns the indices of the participating items and \
the relevant corner";
GetCornerMembers[structure_, {m_, n_}] :=
 With[
  {cor = ComputeCorners /@ structure},
  Catenate@Table[
    With[
     {pos = Position[cor, {i, j}]},
     If[Length@pos != 4,
      Nothing,
      pos
      ]
     ],
    {i, m - 1},
    {j, n - 1}
    ]
  ]

FindMergeLoops::usage = 
  "FindMergeLoops[structure,spec,{m,n}] finds points where four \
corners meet and all relevant axes are merged, which leads to a tiny \
floating ring of merge indicators";
FindMergeLoops[structure_, spec_, {m_, n_}] :=
 
 Map[First] /@ Select[Cases[Rest /@ #, False, {2}] == {} &]@Apply[
    {{##}, spec[[#, 1, #3]], spec[[#, 2, 3 - #2]]} &,
    GetCornerMembers[structure, {m, n}],
    {2}
    ]

ExpandMergeAxesSpec[grid_, structure_, spec_, fromWrap_, {m_, n_}] :=

  With[(* remove floating loops of merge indicators *)
    {loops = 
      GroupBy[First -> Rest]@
       Catenate@FindMergeLoops[structure, #, {m, n}]},
    MapThread[
     ReplacePart[
       Replace[#, s : Except[False] :> {s, s}, {2}],
       Catenate[
        {
           {1, #[[2]], 3 - #[[1]]} -> {{None, 0}, {None, 0}},
           {2, 3 - #[[1]], #[[2]]} -> {{None, 0}, {None, 0}}
           } & /@ #2
        ]
       ] &,
     {
      #,
      Lookup[loops, Range@Length@structure, {}]
      }
     ]
    ] &@MapThread[
   MapThread[
     If[##, False] &,
     {
      SparseArray@
       Normal@KeyMap[$SidePositions]@PlotsPresent[grid, Sequence @@ #],
      #2
      },
     2
     ] &,
   {
    structure,
    ApplyStructureToOption[
      structure,
      Expand2DSpec[
        spec,
        {m, n} - 1,
        LastExplicit
        ] //
       Map[ArrayPad[#, {1, 1}, False] &],(* 
      add settings for outer edges *)
      fromWrap,
      LastExplicit
      ] // Map[NormalizeMergePrimitive, #, {3}] &
    }
   ]

ApplyStyles[prim_, style_, ang_] :=
 
 iApplyStyles[
  prim, {"Graphics", "GraphicsFrame", 
   Replace[style, None -> Nothing]}, ang]
iApplyStyles[gr : Graphics[prim_, opts___], style_, ang_] :=
 Graphics[
  Rotate[prim, ang, {0, 0}],
  PlotRange -> RotateGraphicsOpt[PlotRange, ang]@
    NormalizedOptionValue[gr, PlotRange],
  If[Mod[ang, \[Pi]] == 0,
   {},
   {
    AspectRatio -> 
     Replace[GraphicsOpt[gr, AspectRatio], o_?NumericQ :> 1/o],
    ImageSize -> 
     Replace[
      GraphicsOpt[gr, ImageSize], {{h_, v_} :> {v, h}, 
       s : Except@_UpTo | Automatic :> {Automatic, s}}]
    }
   ],
  opts,
  BaseStyle -> style
  ]
iApplyStyles[prim_, style_, ang_] :=
 Style[prim, style]

$SideAngles = <|Left -> \[Pi], Right -> 0, Top -> \[Pi]/2, 
   Bottom -> 3 \[Pi]/2|>;

ConstructMergeInset[_, _, _, _, _, False, _, _] :=
 {}
ConstructMergeInset[{_, _, pos_, posOff_, size_, sizeOff_}, {spacing_,
    spacingOff_}, frame_, padding_, dir_, spec_, styles_, fullSpec_] :=

  With[
  {
   off = {# - 1, 2 - #} &@$SidePositions[dir][[1]],
   posDelta = (#2 - 1) {2 - #, # - 1} & @@ $SidePositions[dir],
   effSizeOff = sizeOff,
   scale = If[#[[2, 1]] === implicitAutomatic, 1, 1/2] & /@ spec
   },
  MapThread[
   With[
     {effSpec = spec[[#2 + 1, 3 - $SidePositions[dir][[2]]]]},
     If[
      effSpec[[1]] === None,
      {},
      If[
       frame[[3 - $SidePositions[dir][[1]], #2 + 1]] || 
        fullSpec[[3 - $SidePositions[dir][[1]], #2 + 1]] =!= False,
       Inset[
        ApplyStyles[effSpec[[1]], #, $SideAngles[dir]],
        Offset[posOff + (off #2 + posDelta) effSizeOff, 
         Scaled[pos + (off #2 + posDelta) size]],
        If[Head@effSpec[[1]] === Graphics, {0, 0}, 
         ImageScaled[{0.5, 0.5} - 0.5 AngleVector@$SideAngles[dir]]],
        Offset[
         off effSpec[[2]] + (1 - off) spacingOff scale, (1 - 
            off) spacing scale]
        ],
       {}
       ]
      ]
     ] &,
   {
    styles[[3 - $SidePositions[dir][[1]]]],
    {0, 1}
    }
   ]
  ]

ComputeRangeMinMax[ranges_] :=
 {Reverse@MinMax@#, MinMax@#2} & @@ 
  Transpose@ranges
ComputeRangeMinMax[{}] :=
 {{\[Infinity], -\[Infinity]}, \
{-\[Infinity], \[Infinity]}}

ComputePlotRangeMinMax[plotRanges_] :=
 Map[
  ComputeRangeMinMax,
  ReplaceAll[{None, None} -> Nothing]@{
    Transpose@plotRanges[[All, All, 1]],
    plotRanges[[All, All, 2]]
    },
  {2}
  ]

ComputeEffectivePlotRanges[plotRanges_, dims_, spec_] :=
 
 ReplacePart[
  iComputeEffectivePlotRanges[
   plotRanges /. Null -> {{None, None}, {None, None}},
   dims,
   spec
   ],
  Position[plotRanges, Null, {2}] -> Null
  ]

ExpandXYPlotRanges[{xranges_, yranges_}] :=
 Outer[
  {#2, #} &,
  yranges,
  xranges,
  1
  ]

iComputeEffectivePlotRanges[plotRanges_, {m_, n_}, spec_] :=
 
 MapThread[
  Switch[#,
    Inherited,
    #2,
    Min,
    #3[[1]],
    Max,
    #3[[2]],
    _,
    #
    ] &,
  {
   ExpandXYPlotRanges@Replace[
     Expand2DSpec[spec, {m, n}, LastExplicit] /. Directive -> List,
     {
      a_?NumericQ :> {-a, a},
      s : Except@_List :> {s, s}
      },
     {2}
     ],
   plotRanges,
   ExpandXYPlotRanges@ComputePlotRangeMinMax[plotRanges]
   },
  4
  ]

$SpanMarkers = {SpanFromLeft, SpanFromAbove, SpanFromBoth};
ProcessGridStructure[grid_] := 
 DeleteDuplicates[(* get the first (i.e. 
    biggest & highest priority) sub-grid containing each cell *)
    #,
    And @@ MapThread[LessEqual, {#[[1]], #2[[1]], #[[1]] + #[[2]]}] &
    ] &@
  Cases[(* extract all valid spanned sub-grids. SpanFrom* 
    are only allowed as 1x1 *)
    {pos_, 
       dim_, {{p_, 
         SpanFromLeft ...}, {SpanFromAbove, SpanFromBoth ...} ...}} /;
       dim == {0, 0} || 
       MatchQ[p, Except[Alternatives @@ $SpanMarkers]] :>
     {pos, dim}
    ]@Flatten[
    Table[(* generate all possible sub-
       grids (those with higher priority first) *)
       {{x, y}, {h,
          w}, grid[[x ;; x + h, y ;; y + w]]},
       {x, #},
       {y, #2},
       {w, #2 - y, 0, -1},
       {h, # - x, 0, -1}
       ] & @@ Dimensions@grid,
    3
    ]

ComputeEffectivePaddings[{nx_, ny_}, struc_, paddings_] :=
 MapIndexed[
  Table[
    Max[
     Cases[
      struc,
      {{i_, j_}, {h_, w_}} /; #[j, i, w, h] == idx :> 
       ListLookup[paddings, {i, j, Sequence @@ #2}, 0]
      ],
     0
     ],
    {idx, {nx, ny}[[#2[[1]]]]}
    ] &,
  {
   {# &, # + #3 &},
   {#2 + #4 &, #2 &}
   },
  {2}
  ]

NormalizeLabelSpec[Labeled[Verbatim[_], lbl_]] :=
 {lbl}
NormalizeLabelSpec[
  Labeled[Verbatim[_], lbl_, pos_]] :=
 {Placed[lbl, pos]}
NormalizeLabelSpec[Labeled[Verbatim[_], lbl_List, pos_List]] /; 
  Length@lbl == Length@pos :=
 
 Join @@ MapThread[NormalizeLabelSpec[Labeled[_, ##]] &, {lbl, pos}]
NormalizeLabelSpec[Labeled[Verbatim[_], lbl_List, All]] /; 
  Length@lbl == 4 :=
 
 NormalizeLabelSpec[Labeled[_, lbl, {Bottom, Left, Top, Right}]]
NormalizeLabelSpec[
  Labeled[Verbatim[_], 
   Placed[lbl_, pos_, f_ : Identity]]] :=
 {Placed[f@lbl, pos]}

NormalizeLabelOption[n_, opt_] :=
 Map[
  ApplyToWrapped[
    Last[#2, # &]@# &,
    #,
    Except@_Placed,
    _Placed,
    Method -> Function
    ] &,
  iNormalizeLabelOption[n, ThreadWrappers[opt, Placed, 3]],
  {2}
  ]

iNormalizeLabelOption[n_, opt : {___, {___, {___}, ___}, ___}] :=
 
 MapThread[Join, iiNormalizeLabelOption[n, #] & /@ opt]
iNormalizeLabelOption[n_, opt_] :=
 iiNormalizeLabelOption[n, opt]

iiNormalizeLabelOption[n_, None] :=
 Table[{}, n]
iiNormalizeLabelOption[n_, lbl_] :=
 Table[{lbl}, n]
iiNormalizeLabelOption[n_, lbls_List] :=
 
 PadRight[Replace[lbls, lbl : Except[_List] :> {lbl}, 1], n, {{}}]

ComputeLabelIndices[plots_] :=
 Transpose[
    Ordering@*DeleteCases[Null] /@ {Join @@ #, Join @@ Transpose@#}
    ] &@MapIndexed[
   If[# === Null, Null, #2] &,
   plots,
   {2}
   ]

CollectLabels[plots_, fromWrap_, fromOpt_, defPos_] :=
 
 Association@MapThread[
   # -> GroupBy[Last -> First]@
      Map[Replace[lbl : Except[_Placed] :> Placed[lbl, defPos]]]@
       Join[##2] &,
   {ComputeLabelIndices@plots, fromOpt, 
    Pick[fromWrap, Join @@ Map[# =!= Null &, plots, {2}]]}
   ]

ProcessAutoLabel[idx : {_, _}, lbls_] :=
 
 Map[iProcessAutoLabel[#, idx] &, lbls, {2}]

iProcessAutoLabel[lbl_, _] :=
 lbl
iProcessAutoLabel[
  t : Automatic | Row | Column -> spec_, {rowIdx_, colIdx_}] :=
 
 iiProcessAutoLabel[spec, If[t === Row, rowIdx, colIdx]]
iiProcessAutoLabel[spec_, n_] :=
 Quiet[
  Check[
   ApplyToWrapped[
    StringReplace[
      #,
      {
       "i" :> ToLowerCase@RomanNumeral[n],
       "I" :> RomanNumeral[n],
       "a" :> StringJoin@FromLetterNumber@IntegerDigits[n, 26],
       "A" :> 
        ToUpperCase@StringJoin@FromLetterNumber@IntegerDigits[n, 26],
       "1" :> ToString[n]
       },
      1
      ] &,
    spec,
    _String
    ],
   spec,
   {ApplyToWrapped::noMatch}
   ],
  {ApplyToWrapped::noMatch}
  ]
iiProcessAutoLabel[spec_Function, n_] :=
 spec[n]

ConstructLabelInset[{_, _, pos_, posOff_, size_, sizeOff_}, 
  pad_, {{lblPosType : {_, _}, lblPosOff : {xo_, yo_}, 
    lblPos : {x_, y_}}, anchor_}, lbls_, sty_, comb_] :=
 Inset[
  Pane[
   Replace[
     Replace[
      comb,
      Automatic -> (Automatic -> Identity)
      ],
     (Automatic -> func_) :> func@*Which[
        Length@lbls == 1,
        First,
        Abs[x - 0.5] == Abs[y - 0.5],(* if on the diagonals, 
        break tie with offsets*)
        
        If[xo Sign[x - 0.5] > yo Sign[y - 0.5],
         Column,
         Row[#, " "] &
         ],
        Abs[x - 0.5] > Abs[y - 0.5],(* if not on diagonal, 
        use column for left/right quadrant & row otherwise*)
        
        Column,
        True,
        Row[#, " "] &
        ]
     ]@lbls,
   BaseStyle -> {"GraphicsLabel", sty},
   ContentPadding -> False
   ],
  Offset[#, Scaled@#2] & @@ Transpose@MapThread[
     {ipos, iposOff, isize, isizeOff, ipad, ilblPosType, ilblPos, 
       ilblPosOff, isign, iotherInsideQ} |-> {
       iposOff + isizeOff*ilblPos +
        If[ilblPosType === ImageScaled,
         Total@ipad*ilblPos - ipad[[1]],
         If[isign == -1, -ipad[[1]], 0] +
          If[isign == 1, ipad[[2]], 0] +
          
          If[-0.5 <= isign <= 0.5 && ! iotherInsideQ, 
           Total@ipad*ilblPos - ipad[[1]], 0]
         ] +
        ilblPosOff,
       ipos + isize*ilblPos
       },
     {
      pos, posOff,
      size, sizeOff,
      pad,
      lblPosType, lblPos, lblPosOff,
      MapThread[coordSign, {lblPosOff, lblPos}],
      Reverse@MapThread[coordInsideQ, {lblPosOff, lblPos}]
      }
     ],
  ImageScaled@anchor
  ]

$ValidPlaceholders = {Null, None}~Join~$SpanMarkers;

$ValidWrappers = Frame | Labeled;

ValidGridItemQ[gr_?ValidGraphicsQ] :=
 True
ValidGridItemQ[$ValidWrappers[expr_, ___]] :=
 ValidGridItemQ[expr]
ValidGridItemQ[_] :=
 False

PlotGrid::noScaled = 
  "Invalid item size spec ``. At least one column/row dimension must \
be relative.";
PlotGrid::cantScale = 
  "Cannot use Scaled/ImageScaled size for `` ``, as no plots are \
available to estimate the size.";

Options[PlotGrid] = {FrameStyle -> Automatic, ItemSize -> Automatic, 
   Spacings -> None, "ShowFrameLabels" -> Full, 
   "HideOuterTickLabels" -> Full, "MergeAxes" -> False, 
   Method -> Automatic, PlotLabels -> None, 
   LabelingFunction -> Automatic, PerformanceGoal -> "Quality", 
   PlotRange -> Inherited};

PlotGrid[
  l_?(MatrixQ[#, 
      ValidGridItemQ@# || MemberQ[$ValidPlaceholders, #] &] &),
  o : OptionsPattern[Join[Options[PlotGrid], Options[Graphics]]
    ]
  ] := Module[
  {
   nx, ny,
   plots,
   structure,
   wrappers,
   paddingData,
   padding,
   method,
   favorSpeed = OptionValue[PerformanceGoal] === "Speed",
   fixTicks,
   allCustomTicks,
   frameStyle,
   rotateLabel = 
    If[OptionValue[RotateLabel], Rotate[#, \[Pi]/2] &, Identity],
   frameInsets,
   effectivePlotRanges,
   labelData,
   labelInsets,
   gi,
   giNoLabels,
   grid,
   legends,
   computeSizeMeans,
   sizes,
   rangeSizes,
   imageSizes,
   sizeTotals,
   rawSpacings,
   spacings,
   positions,
   positionOffsets,
   sizeOffsets,
   layoutData,
   showFrameLabels,
   hideOuterLabels,
   mergeAxes,
   mergeAxesInsets,
   effectiveSizes,
   options,
   effectiveGridSize
   },
  {ny, nx} = Dimensions@l;
  rawSpacings = 
   Expand2DSpec[OptionValue[Spacings], {nx, ny} - 1, 
     Last] /. {implicitAutomatic | None | Automatic -> 0};
  {plots, wrappers} = Transpose[
    Map[
     ApplyToWrapped[List, #, 
       Except[$ValidWrappers[___]], $ValidWrappers[___]] &,
     l,
     {2}
     ],
    {2, 3, 1}
    ];
  structure = ProcessGridStructure[l];
  showFrameLabels = MapThread[
    ExpandFrameLabelSpec@$FrameLabelMerger@
       Prepend[#]@
        Cases[#2, 
         Frame[_, opts__] :> 
          Replace["ShowFrameLabels" /. {opts}, 
           "ShowFrameLabels" -> Nothing]]
     &,
    {
     ExpandGridSpec[
       OptionValue["ShowFrameLabels"],
       {nx, ny},
       $FrameLabelMerger
       ] /. implicitAutomatic -> Full,
     wrappers
     },
    2
    ];
  hideOuterLabels = ExpandOuterLabelSpec[
    structure,
    OptionValue["HideOuterTickLabels"],
    SidedFrameOptionFromWrappers[wrappers, "HideOuterTickLabels"],
    {nx, ny},
    favorSpeed
    ];
  mergeAxes = ExpandMergeAxesSpec[
    plots,
    structure,
    OptionValue["MergeAxes"],
    SidedFrameOptionFromWrappers[wrappers, "MergeAxes"],
    {nx, ny}
    ];
  method = 
   DeleteDuplicatesBy[First]@
    Flatten@{Replace[OptionValue[Method], Automatic -> {}], 
      "FixFrameTicks" -> True, "AllCustomTicks" -> Automatic};
  fixTicks = "FixFrameTicks" /. method;
  allCustomTicks = "AllCustomTicks" /. method;
  plots = Table[
    With[
     {idx = First@FirstPosition[structure, {{i, j}, _}]},
     ApplyFrameSettings[
      {i, j},
      ListLookup[structure, {idx, 2}, {0, 0}],
      showFrameLabels[[i, j]],
      ListLookup[hideOuterLabels, {idx}, None],
      ListLookup[mergeAxes, {idx}, None],
      Indexed[effectiveSizes, {i, j}],
      plots,
      rawSpacings
      ]
     ],
    {i, ny}, {j, nx}
    ];
  plots = 
   plots /. 
    Hold[EnsureCustomTicks] -> 
     If[TrueQ[
       allCustomTicks /. 
        Automatic :> ! favorSpeed && 
          MemberQ[plots, $ModifiedDefaultTickPattern, All]],
      EnsureCustomTicks,
      #2 &
      ];
  plots = 
   Replace[plots, Alternatives @@ $ValidPlaceholders -> Null, {2}];
  gi = GraphicsInformation@plots;
  effectivePlotRanges = ComputeEffectivePlotRanges[
    gi[PlotRange],
    {nx, ny},
    OptionValue[PlotRange]
    ];
  sizes = 
   Expand2DSpec[OptionValue[ItemSize], {nx, ny}, Last] /. 
    implicitAutomatic -> Automatic;
  computeSizeMeans = Map[Mean@*Replace[{} -> {Indeterminate}]] /@ (
      MapAt[Transpose, 1]@Transpose[
         ReplacePart[
           #,
           
           Position[SparseArray[Rule @@@ structure], Except@0, {3}] ->
             Null
           ] /.
          Null -> {Null, Null},
         {2, 3, 1}
         ] /. Null -> Nothing
      ) &;
  rangeSizes = 
   computeSizeMeans[Apply[Abs@*Subtract, effectivePlotRanges, {3}]];
  imageSizes = computeSizeMeans[gi[ImageSize]];
  sizeOffsets = 
   Replace[sizes, {Offset[off_, _ : 0] :> off, _ -> 0}, {2}];
  sizes = Replace[sizes, Offset[_, sz_ : 0] :> sz, {2}];
  If[MemberQ[Total /@ sizes, 0],
   ResourceFunction["ResourceFunctionMessage"][PlotGrid::noScaled, 
    OptionValue[ItemSize]];
   Return@$Failed
   ];
  sizeTotals = 
   Total@
      Replace[#, {_[s_] :> s, Scaled | ImageScaled | Automatic -> 1}, 
       1] & /@ sizes;
  sizes = MapThread[
      Replace[
        #,
        {
         Scaled[s_] :> s #3,
         ImageScaled[s_] :> s #4,
         Automatic -> 1/#2,
         Scaled -> #3,
         ImageScaled -> #4,
         s_ :> s/#2
         }
        ] &,
      #
      ] & /@ Transpose@{
      sizes,
      sizeTotals + 0 sizes,
      Normalize[#, Total@*Replace[Indeterminate -> 1]] & /@ rangeSizes,
      Normalize[#, Total@*Replace[Indeterminate -> 1]] & /@ 
       imageSizes
      };
  If[MemberQ[Total /@ sizes, 0],
   ResourceFunction["ResourceFunctionMessage"][PlotGrid::noScaled, 
    OptionValue[ItemSize]];
   Return@$Failed
   ];
  If[MemberQ[sizes, Indeterminate, {2}],
   ResourceFunction["ResourceFunctionMessage"][PlotGrid::cantScale, 
      If[# == 1, "row", "column"], #2] & @@ 
    FirstPosition[sizes, Indeterminate];
   Return@$Failed
   ];
  sizes = Normalize[#, Total] & /@ sizes;
  positionOffsets = Replace[rawSpacings, _Scaled -> 0, {2}];
  spacings = Replace[rawSpacings, {Scaled@s_ :> s, _ -> 0}, {2}];
  sizeOffsets += -(Total /@ positionOffsets + Total /@ sizeOffsets)*
    sizes;
  sizes *= (1 - Total /@ spacings);
  positionOffsets = AccumulateShifts[sizeOffsets, positionOffsets];
  positions = AccumulateShifts[sizes, spacings];
  layoutData = With[
      {
       xyLookupBL = XYLookup[# + #2 {1, 0}],
       xyLookupTR = XYLookup[# + #2 {0, 1}]
       },
      {
       ##,
       xyLookupBL@positions,
       xyLookupBL@positionOffsets,
       xyLookupTR@sizes + xyLookupTR@positions - 
        xyLookupBL@positions,
       xyLookupTR@sizeOffsets + xyLookupTR@positionOffsets - 
        xyLookupBL@positionOffsets
       }
      ] & @@@ structure;
  options = {
    PlotRange -> Automatic,
    PlotRangePadding -> None,
    AspectRatio -> Replace[
      OptionValue[AspectRatio],
      {
       m : Mean | GeometricMean :> 
        ny/nx/
         m[
          Divide @@@ 
           DeleteCases[Null]@Flatten[gi["PlotRangeSize"], 1]],
       Automatic -> 
        1/Divide @@ sizeTotals/
         GeometricMean[
          Divide @@@ DeleteCases[Null]@Flatten[gi["PlotRangeSize"], 1]]
       }
      ],
    FilterRules[FilterRules[{o}, Options@Graphics], 
     Except[FrameLabel | Method | PlotRange]],
    Method -> 
     FilterRules[method, Except["FixFrameTicks" | "AllCustomTicks"]],
    ImagePadding -> Automatic
    };
  {plots, legends} = Reap@Map[
     If[# =!= Null,
       ApplyToWrapped[
        (Sow@#2; #) &,
        #,
        _Graphics,
        Legended[_, Except@_?LegendInsideQ],
        Method -> Function
        ]
       ] &,
     plots,
     {2}
     ];
  Block[{effectiveSizes = Table[{\[Infinity], \[Infinity]}, ny, nx]},
   gi = GraphicsInformation@
     {
      Graphics[{}, options],
      With[
       {pl = MapIndexed[
          If[! MatchQ[#, Alternatives @@ $SpanMarkers | Null | None],
            Show[
             #,
             PlotRangePadding -> None,
             PlotRange -> Extract[effectivePlotRanges, #2]
             ]
            ] &,
          plots,
          {2}
          ]},
       {
        pl,
        Replace[(* 
         we also need GraphicsInformation for the plots without \
PlotLabel, 
         since PlotLabel is placed outside explicit ImagePadding \
settings *)
         pl,
         {
          
          p : Except[Null] /; GraphicsOpt[p, PlotLabel] =!= None :> 
           Show[p, PlotLabel -> None],
          _ -> Null
          },
         {2}
         ]
        }
       ]
      }
   ];
  effectiveGridSize = gi[["PlotRangeSize", 1]];
  effectiveSizes = ReplacePart[
    Table[{0, 0}, ny, nx],
    # -> #6 + #5*effectiveGridSize & @@@ layoutData
    ];
  (* giNoLabels = MapThread[If[#2 === Null, ##] &, #, 2] & /@ gi[[All, 2]]; *)
  (* RANZA FIX *)
  giNoLabels = MapThread[If[# === Null, ##] &, #, 2] & /@ gi[[All, 2]];
  gi = gi[[All, 2, 1]];
  
  paddingData = 
   ComputeEffectivePaddings[{nx, ny}, structure, gi[ImagePadding]];
  padding = 
   MapThread[Construct, {{{First, Last}, {Last, First}}, paddingData},
     2];
  frameStyle = NormalizeGraphicsOpt[FrameStyle]@Replace[
     OptionValue[FrameStyle],
     Automatic -> 
      GraphicsOpt[FirstCase[plots, _Graphics, {}, All], FrameStyle]
     ];
  frameInsets = MapThread[
    If[
      # =!= None,
      Inset[
       #3@Pane[#, BaseStyle -> {"GraphicsLabel", "GraphicsFrame", #2}],
       Offset[(5 + #5) #4, Scaled[{0.5, 0.5} + 0.5 #4]],
       ImageScaled[{0.5, 0.5} - 0.5 #4]
       ],
      Nothing
      ] &,
    {
     NormalizeGraphicsOpt[FrameLabel]@OptionValue[FrameLabel],
     frameStyle,
     {{rotateLabel, rotateLabel}, {Identity, Identity}},
     {{{-1, 0}, {1, 0}}, {{0, -1}, {0, 1}}},
     padding
     },
    2
    ];
  labelData = CollectLabels[
    plots,
    Catenate@*Map[NormalizeLabelSpec]@*Cases[_Labeled] /@ 
     Join @@ wrappers,
    NormalizeLabelOption[Count[plots, Except[Null], {2}], 
     OptionValue[PlotLabels]],
    {Left, Top}
    ];
  labelData = KeyValueMap[ProcessAutoLabel, labelData];
  labelInsets = MapThread[
    With[
      {
       layout = #,
       pad = 
        MapThread[
         Extract, {paddingData, 
          ReverseY@Reverse@Transpose@{#[[1]], #[[1]] + #[[2]]}}, 2]
       },
      KeyValueMap[
       ConstructLabelInset[
         layout,
         pad,
         MapAt[Transpose, 1]@placedHelper@#,
         #2,
         OptionValue[LabelStyle],
         OptionValue[LabelingFunction]
         ] &,
       #2
       ]
      ] &,
    {
     Select[layoutData, Extract[plots, #[[1]]] =!= Null &],
     labelData
     }
    ];
  mergeAxesInsets = MapThread[
    With[
      {
       gr = Extract[plots, #[[1]]],
       spac = LookupSpacings[rawSpacings, #[[1]], #[[2]]]
       },
      If[gr =!= Null,
       Table[
        ConstructMergeInset[
         #,
         
         Replace[
          spac[dir], {Scaled@s_ :> {s, 0}, Missing[] -> {0, 0}, 
           s_ :> {0, s}}],
         NormalizedOptionValue[gr, Frame],
         Extract[gi[ImagePadding], #[[1]]],
         dir,
         Extract[#2, $SidePositions[dir]],
         NormalizedOptionValue[gr, FrameStyle],
         #2
         ],
        {dir, Keys@$SidePositions}
        ],
       Nothing
       ]
      ] &,
    {
     layoutData,
     mergeAxes
     }
    ];
  
  grid = Graphics[
    {
     Inset[
      Graphics[
       If[
          ! MatchQ[Extract[l, #], Null | None],
          Annotation[
           Inset[
            If[! MatchQ[Extract[l, #], Alternatives @@ $SpanMarkers],
             Show[
              Extract[plots, #],
              PlotRangePadding -> None,
              PlotRange -> Extract[effectivePlotRanges, #],
              ImagePadding -> Extract[giNoLabels[ImagePadding], #],
              AspectRatio -> Full
              ],
             Graphics[
              Table[
               Inset[
                Graphics@Disk[],
                Scaled[
                 0.5 + 0.2 i*Switch[
                    Extract[l, #],
                    SpanFromLeft, {1, 0},
                    SpanFromAbove, {0, 1},
                    _, {1, -1}
                    ]
                 ],
                {0, 0},
                Scaled@0.2
                ],
               {i, -1, 1, 1}
               ],
              PlotRange -> 1,
              AspectRatio -> Full
              ]
             ],
            Offset[#4, #3],
            Scaled[{0, 0}],
            
            Offset[
             Total /@ (Extract[gi[ImagePadding], #] /. 
                 Null -> 0) + #6, #5]
            ],
           {"PlotGridItem", ##}
           ]
          ] & @@@ layoutData,
       PlotRange -> {{0, 1}, {0, 1}},
       ImagePadding -> padding,
       AspectRatio -> Full
       ],
      Offset[-padding[[All, 1]], Scaled@{0, 0}],
      ImageScaled@{0, 0},
      Offset[Total /@ padding, Scaled@{1, 1}]
      ],
     frameInsets,
     labelInsets,
     mergeAxesInsets
     },
    options
    ];
  effectiveSizes = ReplacePart[
    Table[{0, 0}, ny, nx],
    # -> #6 + #5*effectiveGridSize & @@@ layoutData
    ];
  If[fixTicks,
   grid = grid /.
     Annotation[
       Inset[
        gr_,
        rest___
        ],
       a : {"PlotGridItem", idx_, _, _, _, size_, sizeOff_}
       ] :>
      Annotation[
       Inset[
        Show[
         gr,
         FrameTicks -> MapAt[
           
           ScaleFrameTicks[
            Divide @@ (sizeOff + size*effectiveGridSize)/
             Divide @@ Subtract @@@ Extract[effectivePlotRanges, idx]],
           NormalizedOptionValue[gr, FrameTicks],
           {2, All}
           ]
         ],
        rest
        ],
       a
       ]
   ];
  grid = Show[
    grid,
    CoordinatesToolOptions -> With[
      {
       layoutData = layoutData,
       imagePadding = First /@ padding,
       plotRanges = effectivePlotRanges,
       df = Map[
         If[# =!= Null,
           GetCoordinatesToolOptions[#]@"DisplayFunction"
           ] &,
         plots,
         {2}
         ],
       cvf = Map[
         If[# =!= Null,
           GetCoordinatesToolOptions[#]@"CopiedValueFunction"
           ] &,
         plots,
         {2}
         ]
       },
      Function[
        {type, updateScale, funcs, valFunc, def},
        type -> (
          (
            If[updateScale || ! ListQ@plotScaleData,
             (* 
             Try to compute plot size in absolute units needed for \
coordinate transformations, 
             using image corner or plot range corner as reference. 
             Since "CopiedValueFunction" is only called once the \
values are copied, 
             we need to cache this value for the case that the mouse \
is no longer on top of the plot *)
             
             plotScaleData = {MousePosition["GraphicsImageScaled"], 
               MousePosition["GraphicsAbsolute"]}
             ];
            With[
             {
              plotData = FirstCase[
                layoutData,
                {idx_, _, pos_, posOff_, size_, sizeOff_} :> With[
                  {
                   cursor = #*plotScaleData[[2]],
                   effPos = Total[{posOff, pos}*plotScaleData],
                   effSize = Total[{sizeOff, size}*plotScaleData]
                   },
                  {idx, cursor, effPos, effPos + effSize} /;
         
                             
                   And @@ Thread[effPos < cursor < effPos + effSize]
                  ],
                {0, 0}
                ]
              },
             With[
              {
               id = First@plotData
               },
              If[
               FreeQ[id, 0] && Extract[plotRanges, id] =!= Null,
               valFunc[
                id,
                
                If[# =!= None, #, Nothing &] &[Extract[funcs, id]]@
                 MapThread[
                  Rescale,
                  {
                   plotData[[2]],
                   Transpose@plotData[[3 ;;]],
                   Extract[plotRanges, id]
                   }
                  ]
                ],
               def
               ]
              ]
             ]
            ) &)
        ] @@@ {
        {"DisplayFunction", True, df, 
         Column@{Row@{"Plot: ", #}, #2} &, "Not inside a plot"},
        {"CopiedValueFunction", False, cvf, List, 
         Missing["OutsidePlot"]}
        }
      ]
    ];
  (RightComposition @@ Flatten@legends)@InjectCustomTicksLoading@grid
  ];
End[];
EndPackage[];