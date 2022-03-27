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

ClearLastTickLabels[ticks_List, n_ : 3] :=  Map[MapAt[ClearTickLabel, #, -n ;;] &, ticks, {2}];

TrimTicksAndLabels[gr_Graphics, sides_List] := 
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
    (* ,ImagePadding ->  *)
      (* ReplacePart[AbsoluteOptionsV[gr, ImagePadding], sides->10]  *)
];

GridEdge[pos_List, dims_List] := 
  Position[MapAt[Reverse, Reverse@Transpose@{pos - 1, dims - pos}, {2}], 0];
GridEdgeComplement[pos_List, dims_List] := Complement[Values@$SidePositions, GridEdge[pos,dims] ];


RemovedImagePadding[gr_Graphics,loc_List,dims_List]:=
ReplacePart[AbsoluteOptionsV[gr,ImagePadding],{GridEdgeComplement[loc,dims]->0.3,
GridEdge[loc,dims]->40}];



PlotGrid0[pl__]:=PlotGrid[pl];

(* PlotGrid1[pl_?MatrixQ]:=Grid[MapIndexed[
   Show[#1, ImagePadding -> RemovedImagePadding[#1, #2, Dimensions@pl] ] &, pl, {2}], 
  Spacings -> 0, BaseStyle -> ImageSizeMultipliers -> 1]; *)
  (* //If[COptionValue[{opt,WFGrid},"GridTranspose"],Transpose,Identity], *)
KeepVisible[l_,edges_]:=MapIndexed[If[MemberQ[edges,{First@#2,1}|{First@#2,2}],#1,Nothing]&, l];

PPair[x_]:=If[Print[x];(x===System`Top)||(x===System`Bottom),{System`Center,x},{x,System`Center}];
PRot[x_]:=If[x===System`Left || x===System`Right,Rotate[#,90 Degree]&,Identity];
Options[PlotGrid1]={"GridLabels"->{},"GridTranspose"->False,"LabelPlacement"->{System`Right,System`Top}};
PlotGrid1[pl_?MatrixQ,opt:OptionsPattern[]]:=
Module[{vis, ims, pref,gLab, dims,xx,xxx},
pref=OptionValue[{opt,PlotGrid1},"LabelPlacement"];
gLab=OptionValue[{opt,PlotGrid1},"GridLabels"];
dims=Dimensions@pl;
Grid[
  MapIndexed[
    (
    vis=Intersection[pref/.$SidePositions,GridEdge[#2,dims]];
    ims=If[Last@#2==1||Last@#2==Last@dims,{340,Automatic},{300,Automatic}];
    g=Overlay[{
      g=Show[
        TrimTicksAndLabels[#1,GridEdgeComplement[#2,dims]]
        ,ImagePadding->RemovedImagePadding[#1,#2,dims]
        ,ImageSize->ims
        ]
      ,Graphics[ {MapThread[
        Inset[
          Framed[
            PRot[#2]@Style[First@#1,FontSize->12]
            ,Background->LightYellow]
        ,PPair[#2],PPair[#2]
        ]&
        ,{KeepVisible[Apply[Part[gLab,#1,{1,#2+1}]&,#2],vis]
        ,KeepVisible[pref,vis]/.$SidePlacement}]}
        ,AbsoluteOptions[g,{ImageSize}]
        ,AspectRatio->Apply[Divide, AbsoluteOptionsV[g, ImageSize]]
        ,PlotRangeClipping->False] 
    }];
    xx="/tmp/"<>RandomWord[];
    Export[xx<>"_g.txt",FullForm@g];
    Export[xx<>"_g.png",g];
    g    
    )& 
    ,pl
    ,{2}]
(* ,BaseStyle->ImageSizeMultipliers->1 *)
,Spacings -> {0, 0.1}
]];
(* Grid[MapIndexed[
  Show[
    TrimTicksAndLabels[#1, GridEdgeComplement[#2, dims]]
    , ImagePadding -> RemovedImagePadding[#1, #2, dims]
    , ImageSize -> 
     If[Last@#2 == 1 || Last@#2 == Last@dims, {340, 
       Automatic}, {300, Automatic}]
    ] &, mod, {2}] *)
  

(* Show[#1,ImagePadding -> RemovedImagePadding[#1,#2,Dimensions@pl]],
    Labeled[
      
      (* Apply Labels *)
     ,KeepVisible[Apply[Part[OptionValue[{opt,PlotGrid1},"GridLabels"],#1,{1,#2+1}]&,#2],visible]
      (* Only one vertical and horizontal, the one that's free *)
     ,KeepVisible[pref,visible]/.$SidePlacement
     ,RotateLabel->True] *)
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