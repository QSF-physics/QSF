BeginPackage["filestruct`",{"cmdline`log`","cmdline`opt`"}];

StructKeyMap;
StructPrint;
StructPeak;
StructMap;
StructProcess;
Options[ParsePattern] = {"ParseFilter"->Identity};
ParsePattern;

Begin["`Private`"];

NumberSort[str_String]:=Map[If[StringMatchQ[#, NumberString],ToExpression[#],#] &, StringSplit[str, {"_"," "}]];
NumberSort[list_List]:=NumberSort/@Flatten[list];
KeyCleanup[key_String]:=StringReplace[If[StringContainsQ[key,"."],StringDelete[key, "0" .. ~~ EndOfLine],key], "_"->" "];
KeyCleanup[keys_List]:=If[Length[keys]==0, "", StringRiffle[KeyCleanup/@Flatten[keys],", "] ];
(* KeyCleanup[assoc_Association]:=KeyMap[KeyCleanup,assoc]; *)
PrependName[path_, what_]:=Block[{fsp=FileNameSplit[path]}, FileNameJoin[Join[Most@fsp, {what<>fsp[[-1]]} ] ] ];


(* Helper functions for creating relative paths and flatterning *)
Relativizer:=If[StringStartsQ[ExpandFileName[#],Directory[] ], StringDelete[ExpandFileName[#], Directory[]<>"/" ],#] &;

Options[StructKeyMap]={"StructKeyMapFn"->Identity};
StructKeyMap[assoc_Association, op: OptionsPattern[{StructKeyMap}] ]:=Map[
    StructKeyMap[#,op]&, OptionValue[{StructKeyMap},"StructKeyMapFn"][assoc], -1];
StructKeyMap[else_,op: OptionsPattern[{StructKeyMap}] ]:=else;

Options[StructPrint]= {"StructPrintFn"->LOG,"Heads"->False};
StructPrint[assoc_Association, op: OptionsPattern[{QSFcmdline,StructPrint}]] := KeyValueMap[
    (OptionValue[{QSFcmdline,StructPrint},"StructPrintFn"][
        If[OptionValue["Heads"],ToString[#1,InputForm]<>" -> "<>ToString[Head@#2],#1]
    ]; 
    LL[]++;
    StructPrint[#2,op]; 
    LL[]--;
    ) &, assoc];

StructPeak[ass_Association]:=StructPrint[ass,"Heads"->True];

(* StructPrint[lis_] := OptionValue[{QSFcmdline,StructPrint},"StructPrintFn"]["[", lis, "]"]; *)


(* StructMap[f_, ass_, lvl_] := Block[{dpth = StructDepth[ass]},
  Map[f, ass, {If[# < 0, # + dpth, # - 1] &@lvl}]
]; *)

StructMap[ass_, rule_Rule] := Module[
    {dpth = ArrayDepth[ass, AllowedHeads -> Association],lvl,fn,res,tmp, TreePath},
    (* lvl=If[(First[rule]===All||First[rule]===0), {0}, {dpth-First[rule]+1}]; *)
    (* lvl=DeleteDuplicates[lvl]; *)
    lvl=First[rule];
    fn=Last[rule];
    
    res=MapIndexed[
        (
            LL[]++;
            TreePath=Flatten[{#2}//.{Key[k_]->k}];
            LOG["Applying [", fn,"] at TreePath [", If[TreePath=={},"root", ToString@TreePath],"]"];
            UpdateOpts["TreePath"->TreePath];
            tmp=fn[#1];
            LL[]--;
            tmp
        )&, ass, {lvl}];
    res
];
currentOperations={};

LookForCache[ops_,chi_]:=Module[{rops},
(* Ignore functions, suboperations and output paths *)
rops=DeleteCases[ops,_List|_String|_Symbol];
AssociationThread[Flatten[Position[ops,Rule[_Integer, _],1]]->Map[chi+Hash[Part[rops,;;#]]&,Flatten[Position[rops,Rule[_Integer, _],1]]]]
];
OpName[x_]:=If[Head@x===Composition, StringRiffle[Map[ToString[#]&,Level[x,1]],"_"], ToString[x]];
Options[StructProcess] = {"Operations" -> {}, "CacheHashIndex"->0};
StructProcess[ass_Association, op: OptionsPattern[{StructProcess}]] := 
Module[{ops,opN,chi,ch,H=0,res=ass,cn=1,tmpl=StringTemplate["[``/``] "]},
  Block[{cmdline`opt`options=cmdline`opt`options},
  ops=DeleteCases[OptionValue["Operations"],Null];
  ch=LookForCache[ops,If[OptionValue["CacheHashIndex"]==0,Hash[ass],OptionValue["CacheHashIndex"]]];
  LL[]++; 
  Do[
    step=tmpl[cn++,Length[ops] ];
    Switch[o,
      _Rule,Switch[First[o]
        ,_Integer,H=ch[cn-1];
          opN=StringTemplate["/tmp/mmac/``_``.mx"][OpName[Last[o]],H];
          If[FileExistsQ[opN],
            LOG[step,"Operation: ", ToString[Last@o,InputForm], " cached from ", opN]; 
            res=Import[opN];,
            LOG[step,"Operation: ", ToString[Last@o,InputForm]]; 
            res=StructMap[res,o]; LOG["Caching results to: ",opN];
            Export[opN,res];];
        ,_,LOG[step,"Option: ", ToString[o,InputForm]]; 
          UpdateOpts[o];
        ]; 
      ,_String, AppendToKey["ExportPath"->o]; LOG[step,"Current ExportPath set to: ","ExportPath"/.Options[QSFcmdline]];
      ,_Symbol, LOG[step, "Function ", o]; Catch[o[res],o[]]; 
      ,_List, LOG[step, "Inner List"]; StructProcess[res, "Operations"->o, "CacheHashIndex"->H];
      ,_ , LOGE[step, "Unrecognized command: ", o];
    ];
  ,
  {o, ops}];
  LL[]--; 
  ];
];
LookIn[inp_]:=If[$Notebooks, Quiet@FileNameJoin[{Check[NotebookDirectory[],Directory[]],inp}], inp];


ParsePattern[inp_, op: OptionsPattern[{QSFcmdline,ParsePattern}]]:=Block[
    {inps, ptrn, poses, grFn},
    LOGJ["Parsing Input Patern: ", inp];
        
    inps = FileNameSplit[inp];
    (* Generate a list of files to process, excluding images *)
    
    
    inputFiles=Select[FileNames[LookIn[inp]], 
        ! (StringEndsQ[#,".png"] || StringEndsQ[#,".pdf"]) & ];
    inputFiles=OptionValue[{QSFcmdline, ParsePattern}, "ParseFilter"][inputFiles];
    (* Abort if no files found *)
    CHA[First[inputFiles],"No valid files found in " <> Directory[]];
    UpdateOpts["inputFiles"->inputFiles]; 
    (* Find GROUP pattern positions up to level 10 *)
    ptrn = AssociationMap[Position[inps, el_ /; StringCount[el, "*"] == #, Heads -> False] &, Range[10]];
    ptrn = KeyDrop[ptrn, TakeWhile[Reverse@Range[10], ptrn[#] == {} &]]; (* Remove empty tail *)
    
    UpdateOpts["groupOutput"->Relativizer@inp];

    poses = Flatten /@ Values@Most@ptrn;
    grFn= Table[f[i], {i, poses}] /. {f[a_] -> (Part[FileNameSplit[#], a] &)};
    groupInput = GroupBy[inputFiles,grFn];
    (* Transform List at last level to Association *)
    groupInput = Map[AssociationThread[#->#]&,groupInput,{-2}];
    
    (* Create simpler keys *)
    (* groupInput=StructKeyMap[groupInput, "StructKeyMapFn"->(KeyMap[KeyCleanup])]; *)
    groupInput=StructKeyMap[groupInput, "StructKeyMapFn"->(KeySortBy[NumberSort])];
    (* groupInput=Map[groupInput,0->If[ListQ[#] && (Length[#]==1),First@#,#] ]; *)
    
    UpdateOpts["groupInput"->groupInput];
    UpdateOpts["groupInput_"->groupInput];
    UpdateOpts["groupCount"->Length[Flatten[Keys[groupInput] ] ] ];
    UpdateOpts["maxJoinCount"->Max[Length/@Keys[Values[groupInput] ] ] ];
    (* LOGW["Identified structure sizes:", Dimensions[groupInput, AllowedHeads -> All]]; *)
    LOGW["Identified structure:"];
    StructPrint[groupInput];
    groupInput
];
End[];
EndPackage[];