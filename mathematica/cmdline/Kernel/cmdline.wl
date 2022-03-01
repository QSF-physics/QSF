(* If[$Notebooks, Abort[]]; *)
BeginPackage["cmdline`",{"cmdline`opt`", "cmdline`log`"}];
If[!$Notebooks, 
    (* << "log.wl"; *)
    (* << "opt.wl"; *)
    LOGI["Using CommandLine tools"];
];
EndPackage[];
(* If[!$Notebooks, *)
(* $ContextAliases["log`"] = "cmdline`log`"; *)
(* $ContextAliases["opt`"] = "cmdline`opt`"; *)
(* ]; *)