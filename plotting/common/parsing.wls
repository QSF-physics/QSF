SetOptions[Plot, BaseStyle -> {FontFamily -> "Latin Modern Math", FontSize -> 12}];
SetOptions[ListLinePlot, BaseStyle -> {FontFamily -> "Latin Modern Math", FontSize -> 12}];
SetOptions[Labeled, BaseStyle -> {FontFamily -> "Latin Modern Math", FontSize -> 12}];

PrependName[path_, what_]:=Block[{fsp=FileNameSplit[path]}, FileNameJoin[Join[Most@fsp, {what<>fsp[[-1]]} ] ] ];

ParseInput[defaultPattern_, filterFun_:Identity]:=(
    (* Read first cmd line arg *)
    inp=$ScriptCommandLine[[2]];
    (* If it points directly to file or is not a directory, leave as is, otherwise add *.psi* pattern *)
    inp=If[FileExistsQ[inp]|| !DirectoryQ[inp], inp, Print["Looking for wf files (*.psi*)"]; inp<>defaultPattern];
    Echo[inp, "Using search pattern (* - AVG, @ - group):"];
    
    (* Read second command line arg, with png format as default *)
    outFormat=Quiet[Check[$ScriptCommandLine[[3]],"png"] ];
    outFormat=If[StringStartsQ[outFormat,"."],outFormat,"."<>outFormat];
    Echo[outFormat, "Output outFormat:"];
    (* Look for blur parameter *)
    blur=Quiet[ToExpression@Check[$ScriptCommandLine[[4]],0.1] ];
    Echo[blur, "Blur radius [au]:"];
    
    inps = FileNameSplit[inp];
    (* Find JOIN patterns (putting on the same plot) *)
    joinIDs = Flatten[Position[inps, n_ /; StringContainsQ[n, "@@"], Heads -> False]];
    (* Find GROUP patterns (putting plots next to each other)*)
    groupIDs = Flatten[Position[inps, n_ /; StringContainsQ[n, "@"], Heads -> False]];
    (* Find AVG patterns (averaging the data) *)
    avgIDs = Flatten[Position[inps, n_ /; StringContainsQ[n, "*"], Heads -> False] ];
    (* Generate a list of files to process, excluding images *)
    inputFiles=Select[FileNames[inp], ! (StringEndsQ[#,".png"] || StringEndsQ[#,".pdf"]) & ];
    inputFiles=filterFun[inputFiles];
    (* Create groups, filter trash *)
    groups = GroupBy[inputFiles, Table[FileNameSplit[#][[i]], {i, groupIDs}] & ];
    (* Create simpler keys *)
    keys = StringDelete[#, "0" .. ~~ EndOfLine] & /@ (StringTrim[#, StringSplit[#, "@"] & /@ inps[[groupIDs]]]) & /@ Keys[groups];
    (* Function for creating relative paths *)
    Relativizer:=If[StringStartsQ[ExpandFileName[#],Directory[] ], StringDelete[ExpandFileName[#], Directory[]<>"/" ],#] &;
    PathFlattener:=(StringRiffle[FileNameSplit[#],"_"]<>outFormat ) &;
    outs = Map[
        FileNameJoin@ReplacePart[
            FileNameSplit[Relativizer[#]], 
            AssociationMap[StringReplace[inps[[#]], "*" -> "AVG"] &, avgIDs]
            ] &, Values[groups], {-1}];
    (* groupOutput = Map[
        FileNameJoin@ReplacePart[
            FileNameSplit[#], 
            AssociationMap[StringReplace[inps[[#]], "@" -> ""] &, groupIDs]
            ] &, outs, {-1}];
    groupOutput = PathFlattener@First[DeleteDuplicates[Flatten[groupOutput]]]; *)
    groupOutput = StringReplace[inp,{"*"->"AVG", "@"->"GROUP","/"->"_"}]<>outFormat;

    outs = Map[PathFlattener[First[#]]&, outs, {-2}];
    groupInput = AssociationThread[Transpose[Join[Transpose[keys], {outs}] ], Values[groups] ];
    groupCount=Length[Transpose[Keys[groupInput] ] ]-1;

    Print["Will process following files from: "];
    KeyValueMap[Echo[#2,StringRiffle[Most@#1,"::"]<>":\n"] &,groupInput];
);

