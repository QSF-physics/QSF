#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Build*)


(* ::Subsection:: *)
(*Unlink live package*)


(* ::Input:: *)
(*(*move 'live' package outside*)*)
(*CurrentValue[$FrontEnd,"DefaultStyleDefinitions"]=Inherited;*)
(*RenameFile[$UserBaseDirectory<>"/Autoload/StrongFieldUtils","~/Downloads/StrongFieldUtils"];*)
(*PacletDataRebuild[];*)


(* ::Subsection:: *)
(*Create docs*)


(* ::Input:: *)
(*Needs["ForScience`PacletUtils`"]*)


(* ::Input:: *)
(*printAbort[str_]:=(Print["ABORTING: ",Style[str,Red,Bold]];Quit[])If[$VersionNumber<10.0,printAbort["Mathematica 10.0 or later required."]]*)


(* ::Input:: *)
(*$appName="StrongFieldUtils";*)
(*$docDir=FileNameJoin[{ParentDirectory[NotebookDirectory[]],$appName,"Documentation"}];*)
(*If[Not@DirectoryQ[$docDir],Print["Documentation directory not found.  Aborting."];*)
(*Abort[]]*)


(* ::Input:: *)
(*path=FileNameJoin[{$docDir,"English"}];*)
(*ref=FileNameJoin[{path,"ReferencePages","Symbols"}];*)


(* ::Input:: *)
(*DocumentationHeader[EfromIAtomic]=*)
(*{*)
(*"StrongFieldUtils",*)
(*Red,*)
(*"StrongFieldUtils documentation page"*)
(*};*)
(*Usage[EfromIAtomic]=*)
(*"EfromIAtomic[I] returns the value of electric field magnitude in atomic units from intensity ```I``` given in W/cm^2 unit.";*)


(* ::Input:: *)
(*doc1=DocumentationBuilder[EfromIAtomic];*)
(*NotebookSave[doc1,FileNameJoin[{ref,"EfromIAtomic.nb"}]];*)


(* ::Input:: *)
(*myGuide=Guide["StrongFieldUtils guide"];*)
(*DocumentationHeader[myGuide]={"GUIDE",Orange};*)
(*(*Abstract[myGuide]=*)
(*"StrongFieldUtils gathers some frequently used routines used in strong field physics.";*)*)
(*GuideSections[myGuide]={*)
(*{SectionTitle["Unit conversion functions"],*)
(*{EfromIAtomic,IfromEAtomic,\[Omega]from\[Lambda]Atomic,Tfrom\[Omega]Atomic,Text["Usually SI\[LeftRightArrow]Atomic"]}},*)
(*{SectionTitle["Calculation of often used quantities"],*)
(*{QuiverAtomic,PonderoAtomic,CorkumCutoff(*,Text["Calculating often used quantities"]*)}},*)
(*{SectionTitle["HHG Spectrum"],*)
(*{CalcSpectrum,PresentSpectrum,GaborSpectrum}*)
(*}  *)
(*};*)
(*guide=DocumentationBuilder[myGuide];*)
(*NotebookSave[guide,FileNameJoin[{path,"Guide.nb"}]];*)


(* ::Input:: *)
(*Needs["DocumentationSearch`"]*)
(*SetDirectory[path];*)
(*CreateDocumentationIndex[path]*)


(* ::Input:: *)
(*indexDir=CreateDirectory["Index"]*)
(*ind=NewDocumentationNotebookIndexer[indexDir]*)
(*AddDocumentationDirectory[ind,"ReferencePages/Symbols"]*)
(*CloseDocumentationNotebookIndexer[ind]*)
(*spellIndexDir=CreateDirectory["SpellIndex"]*)
(*CreateSpellIndex[indexDir,spellIndexDir]*)


(* ::Subsection:: *)
(*Update version*)


(* ::Input:: *)
(*(*load*)*)
(*PacletDirectoryLoad@ParentDirectory[NotebookDirectory[]];*)


(* ::Input:: *)
(*(*update version*)*)
(*exportPath=(First@PacletFind["StrongFieldUtils"])@"Location"<>"/PacletInfo.m";*)
(*pacletE=Get[exportPath] /.(Version->version_):>(Version->(StringRiffle[Date[][[;;4]],"."]/.r:Except[_String]:>Throw[r,tag]));*)
(*Block[{Internal`$ContextMarks=False},Export[exportPath,pacletE,"Package",PageWidth->80]];*)
(*Version/.List@@pacletE*)


(* ::Input:: *)
(*(*check status*)*)
(*PacletDataRebuild[];*)
(*PacletFind["StrongFieldUtils"]*)


(* ::Subsection:: *)
(*Pack*)


(* ::Input:: *)
(*(*view junk*)*)
(*(*PacletDirectoryLoad@ParentDirectory[NotebookDirectory[]];*)*)
(*FileNames[".*",ParentDirectory[NotebookDirectory[]],Infinity]*)
(*(*delete junk*)*)
(*DeleteFile/@FileNames[".*",ParentDirectory[NotebookDirectory[]],Infinity];*)
(*CreatePacletArchive[ParentDirectory[NotebookDirectory[]]<>"/StrongFieldUtils",NotebookDirectory[]]*)


(* ::Subsection:: *)
(*Unload temporary package used for packing*)


(* ::Input:: *)
(*(*unload temporary package*)*)
(*PacletDirectoryUnload@ParentDirectory[NotebookDirectory[]];*)
(*(*this should be empty now*)*)
(*PacletFind["StrongFieldUtils"]*)


(* ::Subsection:: *)
(*Install from paclet file*)


(* ::Input:: *)
(*(*Preview the installed version*)*)
(*SetDirectory@NotebookDirectory[];*)
(*Last@FileNames["*.paclet"]*)


(* ::Input:: *)
(*(*Install it*)*)
(*PacletInstall@Last@FileNames["*.paclet"];*)
(*PacletFind["StrongFieldUtils"]*)


(* ::Subsection:: *)
(*Remove it*)


(* ::Input:: *)
(*PacletUninstall/@PacletFind["StrongFieldUtils"];*)
(*(*This should be empty now*)*)
(*PacletFind["StrongFieldUtils"]*)


(* ::Subsection:: *)
(*Restore live package*)


(* ::Input:: *)
(*(*move 'live' package inside*)*)
(*RenameFile["~/Downloads/StrongFieldUtils",$UserBaseDirectory<>"/Autoload/StrongFieldUtils"];*)
(*(*Look for loaded paclets*)*)
(*PacletDataRebuild[];*)
(*PacletFind["StrongFieldUtils"]*)
(*PacletDataRebuild[];*)
