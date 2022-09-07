#!/usr/bin/env wolframscript
BeginPackage["PlotGrid`",{"PlotGrid`External`"}];

PlotGrid0;
PlotGrid1;
PlotGrid2;
TrimTicksAndLabels;

Begin["`Private`"];


$SidePositions = <|System`Left->{1,1},System`Right->{1,2}
                  ,System`Bottom->{2, 1},System`Top ->{2,2}|>;
$SidePlacement = AssociationMap[Reverse,$SidePositions];

ClearTickLabel := ReplacePart[#, 2 -> ""] &;
AbsoluteOptionsV := #2 /. AbsoluteOptions[#1, #2] &;

ClearLastTickLabels[ticks_List, n_ : 2] :=  Map[MapAt[ClearTickLabel, #, -n ;;] &, ticks, {2}];

TrimTicksAndLabels[g_Grid, sides_List,opt:OptionsPattern[]]:=g;
(* Options[TrimTicksAndLabels]={"PlotGridPadding"->{{30,20},{10,20}}}; *)
TrimTicksAndLabels[gr_Graphics, sides_List,opt:OptionsPattern[]] := 
  Show[gr, 
    FrameLabel -> 
        ToExpression@
        ReplacePart[AbsoluteOptionsV[gr, FrameLabel], sides -> None] 
    ,FrameTicks -> 
        ClearLastTickLabels@
        MapAt[
            ReplacePart[#, {_, 2, ___} -> ""] &, 
            AbsoluteOptionsV[gr, FrameTicks]
        , sides]
    ,FrameStyle->{Black}
    ,ImagePadding->ReplacePart[
      AbsoluteOptionsV[gr,ImagePadding]
      ,{sides->0.3,Complement[Values@$SidePositions,sides]->OptionValue[{opt,PlotGrid1},"PlotGridPadding"]}]
];

GridEdge[pos_List,dims_List]:=
Position[MapAt[Reverse,Reverse@Transpose@{pos-1,dims-pos},{2}],0];
GridEdgeComplement[pos_List,dims_List]:=Complement[Values@$SidePositions,GridEdge[pos,dims]];


RemovedImagePadding[gr_Graphics,loc_List,dims_List]:=
ReplacePart[AbsoluteOptionsV[gr,ImagePadding],{GridEdgeComplement[loc,dims]->0.3,
GridEdge[loc,dims]->40}];



KeepVisible[l_,edges_]:=MapIndexed[If[MemberQ[edges,{First@#2,1}|{First@#2,2}],#1,Nothing]&, l];

PPair[x_]:=If[(x===System`Top)||(x===System`Bottom),{System`Center,x},{x,System`Center}];
PPairRel[x_]:=If[(x===System`Top)||(x===System`Bottom),{System`Center,0},{0,System`Center}];
PRot[x_]:=If[x===System`Left || x===System`Right,Rotate[#,90 Degree]&,Identity];

Options[PlotGrid1]={"GridLabels"->{},"GridTranspose"->False,"LabelPlacement"->{System`Right,System`Top},"PlotGridPadding"->55};

PlotGrid1[pl_?MatrixQ,opt:OptionsPattern[]]:=
Module[{vis, ims, pref,gLab, dims,g},
pref=OptionValue[{opt,PlotGrid1},"LabelPlacement"];
gLab=OptionValue[{opt,PlotGrid1},"GridLabels"];
pgp=OptionValue[{opt,PlotGrid1},"PlotGridPadding"];
dims=Dimensions@pl;
Grid[
  MapIndexed[
    (
    vis=Intersection[pref/.$SidePositions,GridEdge[#2,dims]];
    ims=If[Last@#2==1||Last@#2==Last@dims,{300+pgp,Automatic},{300,Automatic}];
    (* Overlay[{ *)
      Show[Show[
        TrimTicksAndLabels[#1,GridEdgeComplement[#2,dims],opt]
        (* ,ImagePadding->RemovedImagePadding[#1,#2,dims] *)
        ,ImageSize->ims
        ],Epilog->Inset[Framed[Style[
            "("<>ToString[Part[Alphabet[],(dims[[2]] (#2[[1]] - 1)) + #2[[2]]]]<>")",  FontFamily->"LatinModernMath",FontSize->14,FontColor->Black]
            ,Background -> LightGray,FrameStyle->Directive[Thin, LightGray]
          ], Scaled[{1.01, 1.03}], {Right, Top}]]
      (* ,Graphics[{MapThread[
        Inset[
          Framed[
            PRot[#2]@Style[First@#1,FontSize->12]
            ,Background->LightYellow]
        ,PPair[#2],PPair[#2]
        ]&
        ,{KeepVisible[Apply[Part[gLab,#1,{1,#2+1}]&,#2],vis]
        ,KeepVisible[pref,vis]/.$SidePlacement}]}
        ,AbsoluteOptions[g,{ImageSize}]
        ,AspectRatio->Apply[Divide, Reverse@AbsoluteOptionsV[g, ImageSize]]
        ,PlotRangeClipping->False] 
    }] *)
    )& 
    ,pl
    ,{2}]
,BaseStyle->ImageSizeMultipliers->1
,Spacings -> {0, 0.1}
]];


PlotGrid0[pl__]:=PlotGrid[pl];
PlotGrid2[pl_?MatrixQ]:=GraphicsGrid[MapIndexed[
    TrimTicksAndLabels[#1, GridEdgeComplement[#2, Dimensions@pl]] &, pl, {2}], 
    Spacings -> Scaled[-0.1]];



EmptyQ:=Length[#]==0&;
RandomPlotGrid[width_:2,height_:2]:=
  Table[ListLinePlot[RandomReal[RandomInteger[10],10],Frame->True],{width},{height}];

RandomAssoc[size_:4,depth_:3,balanced_:True]:=If[depth>0,AssociationThread[RandomWord[#,Language->"English"],
  Table[RandomAssoc[size,depth-1,balanced],#]]&@If[balanced||(depth!=1),size,RandomInteger[size]],
  Style[FileNameJoin[RandomWord[2,Language->"English"]],Background->LightBlue]];

  
End[];


EndPackage[];