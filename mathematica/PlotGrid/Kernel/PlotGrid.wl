#!/usr/bin/env wolframscript
BeginPackage["PlotGrid`",{"PlotGrid`External`"}];

PlotGrid0;
PlotGrid1;
PlotGrid2;

Begin["`Private`"];


$SidePositions = <|Left->{1,1},Right->{1,2},Bottom->{2, 1},Top ->{2,2}|>;

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
    ,ImagePadding -> 
        MapAt[(0 #) + 1 &, AbsoluteOptionsV[gr, ImagePadding], sides ] 
    
];

GridEdge[pos_List, dims_List] := 
  Position[MapAt[Reverse, Reverse@Transpose@{pos - 1, dims - pos}, {2}], 0];
GridEdgeComplement[pos_List, dims_List] := Complement[Values@$SidePositions, GridEdge[pos,dims] ];

RemovedImagePadding[gr_Graphics, loc_List, dims_List] := 
  MapAt[0 &, AbsoluteOptionsV[gr, ImagePadding], GridEdgeComplement[loc,dims]];


PlotGrid0[pl__]:=PlotGrid[pl];

PlotGrid1[pl_?MatrixQ]:=Grid[MapIndexed[
   Show[#1, ImagePadding -> RemovedImagePadding[#1, #2, Dimensions@pl] ] &, pl, {2}], 
  Spacings -> 0, BaseStyle -> ImageSizeMultipliers -> 1];

PlotGrid2[pl_?MatrixQ]:=GraphicsGrid[MapIndexed[
    TrimTicksAndLabels[#1, GridEdgeComplement[#2, Dimensions@pl]] &, pl, {2}], 
    Spacings -> Scaled[-0.1]];




End[];


EndPackage[];