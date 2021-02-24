#This Tcl file is written for demostrating tcl commands for Heat Transfer module in OpenSees


wipe

set sectionType protected
# set sectionType exposed

set tf 0.0108
set tw 0.0076
set h 0.45
set dp 0.016
set b 0.15
set Centrex 0.0
set centrey [expr $h/2.0]
set midBotFlange [expr 0.5*$tf]
set midTopFlange [expr $h-0.5*$tf]
set midWeb [expr 0.5*$h]

set elemFx 10; # Flange x-elements
set elemFy 8
set elemWx 8
set elemWy 20
set elemdp 4

HeatTransfer 2D;       #HeatTransfer activates the HTModule. 2D ,or 2d, or 3D or 3d indicate the model dimension. 

#Defining HeatTransfer Material with Material tag 1.
HTMaterial CarbonSteelEC3 1;
HTMaterial SFRM 2 1;


#HTEntity Isection $tag $centreX $centreY $Bf $Hb  $Tw  $Tf;
# HTEntity $EntityTYpe $EntityTag $Centrelocx $Centrelocy 		bf  	hb 		 tw		tf    	coat
if {$sectionType == "protected"} {
  HTEntity ProtectedIsection  	1		0.0				0.0		$b		$h		$tw	    $tf     $dp
  # HTMesh $meshTag $EntityTag $MaterialTag  <-SecondMat $secMatTag> -phaseChange $PCTag <$PCTag2> -MeshCtrls $FlangeEleSizeX    $FlangeEleSizeY 		$WebEleSizeX 		$WebEleSizeY 					$CoatLayerThick
  HTMesh	1		1			1			-SecondMat  2 		   -phaseChange	1 		0	     -MeshCtrls	[expr ($b-$tw-2*$dp)/$elemFx]  [expr $tf/$elemFy]	[expr $tw/$elemWx]	[expr ($h-2*$tf-2*$dp)/$elemWy]		[expr $dp/$elemdp]
  puts "Using the protected section."
} elseif {$sectionType == "exposed"} {
  HTEntity Isection  	1		0.0				0.0		$b		$h		$tw	    $tf  
  HTMesh	1		1			1			-phaseChange	1 	   -MeshCtrls	[expr ($b-$tw)/$elemFx]  [expr $tf/$elemFy]	[expr $tw/$elemWx]	[expr ($h-2*$tf)/$elemWy]		[expr $dp/$elemdp]
  puts "Using the UNprotected section."
} else {
puts "no section type selected. Aborting."
return
}
HTMeshAll



SetInitialT 293.15;
HTNodeSet 1 -Locx 0.0 -Locy $midBotFlange
HTNodeSet 2 -Locx 0.0 -Locy $midWeb
HTNodeSet 3 -Locx 0.0 -Locy $midTopFlange
HTConstants 1 35.0 293.15 0.85 0.85

# thermal load assignment 
set fileName AST1.dat
FireModel	UserDefined	1	-file	$fileName -type 1
# FireModel standard 1
# FireModel hydroCarbon 1



HTPattern fire 1 model 1 {
	HeatFluxBC -HTEntity 1 -face 1 4 5 6 7 8 9 -type -ConvecAndRad -HTConstants 1
}

HTRecorder -file BotFlange.dat -NodeSet 1
HTRecorder -file theWeb.dat -NodeSet 2
HTRecorder -file TopFlange.dat -NodeSet 3
HTAnalysis HeatTransfer TempIncr 0.1 1000 2 Newton
HTAnalyze 10 10
wipeHT


