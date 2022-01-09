SetOptions[Plot, BaseStyle -> {FontFamily -> "Latin Modern Math", FontSize -> 12}];
SetOptions[ListLinePlot, BaseStyle -> {FontFamily -> "Latin Modern Math", FontSize -> 12}];
SetOptions[Labeled, BaseStyle -> {FontFamily -> "Latin Modern Math", FontSize -> 12}];
FontFix:=Style[#,FontFamily -> "Latin Modern Math"]&;
NumberSort[str_String]:=ToExpression/@Select[StringSplit[str, {"_"," "}], NumberQ[ToExpression[#]] &];
NumberSort[list_List]:=NumberSort/@Flatten[list];
KeyCleanup[key_String]:=StringReplace[If[StringContainsQ[key,"."],StringDelete[key, "0" .. ~~ EndOfLine],key], "_"->" "];
KeyCleanup[keys_List]:=If[Length[keys]==0, "", StringRiffle[KeyCleanup/@Flatten[keys],", "] ];
PrependName[path_, what_]:=Block[{fsp=FileNameSplit[path]}, FileNameJoin[Join[Most@fsp, {what<>fsp[[-1]]} ] ] ];

esc = Association[reset -> "\033[1;0m",
   black -> "\033[1;30m", red -> "\033[1;31m",
   green -> "\033[1;32m", yellow -> "\033[1;33m",
   blue -> "\033[1;34m", magenta -> "\033[1;35m"];

LOG[c_,test___]:=Print[esc[c],test, esc[reset]];

ParseInput[defaultPattern_, filterFun_:Identity]:=(
    (* Helper functions for creating relative paths and flatterning *)
    Relativizer:=If[StringStartsQ[ExpandFileName[#],Directory[] ], StringDelete[ExpandFileName[#], Directory[]<>"/" ],#] &;
    PathFlattener:=(StringReplace[#,"/"->"_"]<>outFormat ) &;
    (* Read first cmd line arg *)
    inp=$ScriptCommandLine[[2]];
    (* If it points directly to file or is not a directory, leave as is, otherwise add *.psi* pattern *)
    inp=If[FileExistsQ[inp]|| !DirectoryQ[inp], inp, inp<>defaultPattern];
    LOG[green, "Parser will look files of kind:",defaultPattern]
    LOG[green, "Search pattern (* - AVG, @ - group, @@ - join):", inp];
    
    (* Read second command line arg, with png format as default *)
    outFormat=Quiet[Check[$ScriptCommandLine[[3]],"png"] ];
    outFormat=If[StringStartsQ[outFormat,"."],outFormat,"."<>outFormat];
    LOG[blue, "Output outFormat:",outFormat];
    (* Look for blur parameter *)
    blur=Quiet[ToExpression@Check[$ScriptCommandLine[[4]],0.1] ];
    LOG[blue, "Blur radius [au]:",blur];
    
    inps = FileNameSplit[inp];
    (* Generate a list of files to process, excluding images *)
    inputFiles=Select[FileNames[inp], ! (StringEndsQ[#,".png"] || StringEndsQ[#,".pdf"]) & ];
    inputFiles=filterFun[inputFiles];
    (* Abort if no files found *)
    If[Length[inputFiles]==0,LOG[red,"No valid files found. Aborting"]; Abort[]];

    (* Find JOIN patterns (putting on the same plot) *)
    joinIDs = Flatten[Position[inps, n_ /; StringContainsQ[n, "@@"], Heads -> False]];
    (* Find GROUP patterns (putting plots next to each other)*)
    groupIDs = Flatten[Position[inps, n_ /; StringCount[n, "@"]==1, Heads -> False]];
    (* Find AVG patterns (averaging the data) *)
    avgIDs = Flatten[Position[inps, n_ /; StringContainsQ[n, "*"], Heads -> False] ];
    (* Output file name *)
    groupOutput = StringReplace[inp,{"*"->"AVG", "@@"->"JOIN", "@"->"GROUP","/"->"_"}]<>outFormat;
    (* Create groups with embedded joins *)
    groupInput = GroupBy[inputFiles, {Table[FileNameSplit[#][[i]], {i, groupIDs}] &, Table[FileNameSplit[#][[i]], {i, joinIDs}] &} ];
    (* Create simpler keys *)
    groupInput=KeyMap[KeyCleanup,groupInput];
    groupInput=Map[KeyMap[KeyCleanup,#]&, groupInput];
    groupInput=KeySortBy[NumberSort][groupInput];
    groupInput=Map[KeySortBy[NumberSort][#]&, groupInput];


    groupCount=Length[Flatten[Keys[groupInput]]];
    maxJoinCount=Max[Length/@Keys[Values[groupInput] ] ];
   
    LOG[blue,"Identified structure:\n"];
    KeyValueMap[
        (LOG[green,"GROUP:",#1]; 
        KeyValueMap[(
            LOG[yellow,"  JOIN:", #1]; 
            LOG[magenta,"    AVG:", #]&/@#2; 
        )&,#2] 
    ) &,groupInput];
);

