#This Tcl file is written for demostrating tcl commands for Heat Transfer module in OpenSees


wipe

set tf 17.4e-3
set tw 10.16e-3
set h 215.9e-3
set b 205.994e-3
set dp 25e-3
set midBotFlange [expr -0.5*$h+0.5*$tf]
set midTopFlange [expr +0.5*$h-0.5*$tf]
set midWeb 0

set elemFx 20
set elemFy 12
set elemWx 10
set elemWy 30
set elemdp 5
# set midWeb [expr -($h-2*$tf)/$elemWy]

HeatTransfer 2D;       #HeatTransfer activates the HTModule. 2D ,or 2d, or 3D or 3d indicate the model dimension. 

# HTMaterial CarbonSteelEC3 1;       #Defining HeatTransfer Material with Material tag 1.

# HTMaterial ConcreteEC2 2 0.0; 
HTMaterial CarbonSteelEC3 1
# HTMaterial CarbonSteelEC3 2
HTMaterial SFRM 2


#HTEntity Isection $tag $centreX $centreY $Bf $Hb  $Tw  $Tf;
# HTEntity $EntityTYpe $EntityTag $Centrelocx $Centrelocy 		bf  	hb 		 tw		tf    	coat
  HTEntity ProtectedIsection  	1		0				0		$b		$h		$tw	    $tf     $dp
	# HTEntity Isection  	1		0				0		$b		$h		$tw	    $tf  

# HTMesh $meshTag $EntityTag $MaterialTag  <-SecondMat $secMatTag> -phaseChange $PCTag <$PCTag2> -MeshCtrls $FlangeEleSizeX    $FlangeEleSizeY 		$WebEleSizeX 		$WebEleSizeY 					$CoatLayerThick
  HTMesh	1		1			1			-SecondMat  2 		   -phaseChange	1 		2	     -MeshCtrls	[expr $b/$elemFx]  [expr $tf/$elemFy]	[expr $tw/$elemWx]	[expr ($h-2*$tf)/$elemWy]		[expr $dp/$elemdp]
  # HTMesh	1		1			1			-phaseChange	1 	   -MeshCtrls	[expr $b/$elemFx]  [expr $tf/$elemFy]	[expr $tw/$elemWx]	[expr ($h-2*$tf)/$elemWy]		[expr $dp/$elemdp]

#HTMesh $MeshTag $EntityTag  $MaterialTag -SecondMat 2



HTMeshAll

#RefineMesh 1 -seed 1 ;
SetInitialT 293.15
# HTNodeSet 	1	-Entity	1	 -face	7
# HTNodeSet $NodeSetTag -Locx $locxMin <$locxMax> <-Locy $locyMin $locyMax>

HTNodeSet 1 -Locx 0 0 -Locy $midBotFlange $midBotFlange
HTNodeSet 2 -Locx 0.0 0.0 -Locy 0.0 0.0
HTNodeSet 3 -Locx 0 0 -Locy $midTopFlange $midTopFlange
HTConstants 1 25.0 293.15 0.9 0.9


#FireModel standard 1; 
set fileName AST1.dat
# FireModel	UserDefined	1	-file	$fileName -type 1

FireModel hydroCarbon 1



HTPattern fire 1 model 1 {
	HeatFluxBC -HTEntity 1 -face 1 -type ConvecAndRad -HTConstants 1
	HeatFluxBC -HTEntity 1 -face 2 -type ConvecAndRad -HTConstants 1
	HeatFluxBC -HTEntity 1 -face 3 -type ConvecAndRad -HTConstants 1
	HeatFluxBC -HTEntity 1 -face 4 -type ConvecAndRad -HTConstants 1
	HeatFluxBC -HTEntity 1 -face 5 -type ConvecAndRad -HTConstants 1
	HeatFluxBC -HTEntity 1 -face 6 -type ConvecAndRad -HTConstants 1
	HeatFluxBC -HTEntity 1 -face 7 -type ConvecAndRad -HTConstants 1
	HeatFluxBC -HTEntity 1 -face 8 -type ConvecAndRad -HTConstants 1
	HeatFluxBC -HTEntity 1 -face 9 -type ConvecAndRad -HTConstants 1
	HeatFluxBC -HTEntity 1 -face 10 -type ConvecAndRad -HTConstants 1
	HeatFluxBC -HTEntity 1 -face 11 -type ConvecAndRad -HTConstants 1
}

HTRecorder -file BotFlange.dat -NodeSet 1;
HTRecorder -file theWeb.dat -NodeSet 2;
HTRecorder -file TopFlange.dat -NodeSet 3;

HTAnalysis HeatTransfer Residual 0.1 1000 1 Newton
HTAnalyze 1800 10

wipeHT;


