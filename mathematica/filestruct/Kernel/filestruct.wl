BeginPackage["filestruct`",{"cmdline`log`","cmdline`opt`"}];

StructKeyMap;
StructPrint;
StructPeak;
StructMap;
StructProcess;
Options[ParsePattern] = {"FilterFunction"->Identity};
ParsePattern;

Begin["`Private`"];

NumberSort[str_String]:=Map[If[StringMatchQ[#, NumberString],ToExpression[#],#] &, StringSplit[str, {"_"," "}]];
NumberSort[list_List]:=NumberSort/@Flatten[list];
KeyCleanup[key_String]:=StringReplace[If[StringContainsQ[key,"."],StringDelete[key, "0" .. ~~ EndOfLine],key], "_"->" "];
KeyCleanup[keys_List]:=If[Length[keys]==0, "", StringRiffle[KeyCleanup/@Flatten[keys],", "] ];
(* KeyCleanup[assoc_Association]:=KeyMap[KeyCleanup,assoc]; *)
PrependName[path_, what_]:=Block[{fsp=FileNameSplit[path]}, FileNameJoin[Join[Most@fsp, {what<>fsp[[-1]]} ] ] ];


Options[ParsePattern] = {"FilterFunction"->Identity};

(* Helper functions for creating relative paths and flatterning *)
Relativizer:=If[StringStartsQ[ExpandFileName[#],Directory[] ], StringDelete[ExpandFileName[#], Directory[]<>"/" ],#] &;

Options[StructKeyMap]={"StructKeyMapFn"->Identity};
StructKeyMap[assoc_Association, op: OptionsPattern[{StructKeyMap}] ]:=Map[
    StructKeyMap[#,op]&, OptionValue[{StructKeyMap},"StructKeyMapFn"][assoc], -1];
StructKeyMap[else_,op: OptionsPattern[{StructKeyMap}] ]:=else;

Options[StructPrint]= {"StructPrintFn"->LOG,"Heads"->False};
StructPrint[assoc_Association, op: OptionsPattern[{QSFcmdline,StructPrint}]] := KeyValueMap[
    (OptionValue[{QSFcmdline,StructPrint},"StructPrintFn"][
        If[OptionValue["Heads"],ToString[#1]<>" -> "<>ToString[Head@#2],#1]
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

StructMap[ass_, rule_Rule] := Block[
    {dpth = ArrayDepth[ass, AllowedHeads -> Association],lvl,fn,res, keys,tmp},
    (* lvl=If[(First[rule]===All||First[rule]===0), {0}, {dpth-First[rule]+1}]; *)
    (* lvl=DeleteDuplicates[lvl]; *)
    lvl=First[rule]-1;
    fn=Last[rule];
    
    res=MapIndexed[
        (
            keys=Flatten[{#2}//.{Key[k_]->k}];
            LOG["Applying to group of [",Head@#1,"] at key path [", If[keys=={},"root", keys],"]"];
            LL[]++;
            tmp=fn[#1];
            LL[]--;
            tmp
        )&, ass, {lvl}];
    res
];
Options[StructProcess] = {"Operations" -> {}};
StructProcess[ass_Association, op: OptionsPattern[{QSFcmdline,StructProcess}]] := 
Module[{res = ass, cn=1, tmpl=StringTemplate["[``/``]"]},   
    Do[
        step=tmpl[cn++,Length[OptionValue["Operations"] ] ];
        Switch[o,
            _Rule, 
                
                Switch[First[o],
                        
                    _Integer, 
                        LL[]++; 
                        LOG[step," Operation ", ToString[o]]; 
                        res=StructMap[res,o];,
                        LL[]--;
                    _, 
                        LL[]++; 
                        LOG[step," Option ", ToString[o]]; 
                        cmdline`opt`UpdateOpts[o];
                        LL[]--;
                ]; 
                ,
            _Symbol, LOG[step, " Function ", o]; 
                LL[]++; Catch[o[res],o[]]; LL[]--;,
            _List, LOG[step, " Inner process (to avoid saving results)"]; 
                LL[]++;
                StructProcess[res, "Operations"->o];,
                LL[]--;
            _ , LOGE[step, " Unrecognized command: ", o];
        ];
    ,
    {o, OptionValue["Operations"]}];
];

ParsePattern[inp_, op: OptionsPattern[{QSFcmdline,ParsePattern}]]:=Block[
    {inps, ptrn, poses, grFn},
    LOGII["Parsing Input Patern: ", inp];
        
    inps = FileNameSplit[inp];
    (* Generate a list of files to process, excluding images *)
    inputFiles=Select[FileNames[
        If[$Notebooks, FileNameJoin[{NotebookDirectory[],inp}], inp]], ! (StringEndsQ[#,".png"] || StringEndsQ[#,".pdf"]) & ];
    inputFiles=OptionValue[{QSFcmdline, ParsePattern}, "FilterFunction"][inputFiles];
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
    LOGW["Key count at particular levels: ", Dimensions[groupInput, AllowedHeads->Association],"\n"];
    groupInput
];
End[];
EndPackage[];