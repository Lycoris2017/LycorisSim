// $Id: DetectorConstruction.cc 94 2010-01-26 13:18:30Z adotti $
/**
 * @file
 * @brief Implements mandatory user class DetectorConstruction.
 */

#include "DetectorConstruction.hh"
//#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"


// New include files - used for magnetic field
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnitsTable.hh"
//#include "SensitiveDetector.hh"
#include "G4SDManager.hh"

DetectorConstruction::DetectorConstruction()
{
	//--------- Material definition ---------
	DefineMaterials();

	//--------- Sizes of the principal geometrical components (solids)  ---------
	ComputeParameters();
}
 
DetectorConstruction::~DetectorConstruction()
{
}
 
void DetectorConstruction::DefineMaterials() 
{
	//Get Materials from NIST database
	G4NistManager* man = G4NistManager::Instance();
	man->SetVerbose(0);

	// define NIST materials
	vacuum  = man->FindOrBuildMaterial("G4_Galactic");

	//Tracker
	air     = man->FindOrBuildMaterial("G4_AIR");
	ar = man->FindOrBuildMaterial("G4_Ar");
	cu = man->FindOrBuildMaterial("G4_Cu");
	steel = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
	kapton = man->FindOrBuildMaterial("G4_KAPTON");
	si = man->FindOrBuildMaterial("G4_Si");
}
 
void DetectorConstruction::ComputeParameters() 
{

	//This function defines the defaults
	//of the geometry construction

	// ** world **
	halfWorldLength = 5.3* m;
	halfFieldStripsLength = 305*mm;
	innerKaptonRadius = 382.5*mm;
	outerKaptonRadius = (innerKaptonRadius+0.011*kapton->GetRadlen());
	innerInnerFieldStripsRadius = 382.55*mm;
	outerInnerFieldStripsRadius = 382.65*mm;
	innerOuterFieldStripsRadius = 382.75*mm;
	outerOuterFieldStripsRadius = 382.85*mm;
	halfLPLength = 305*mm;
	innerMagnetRadius = 425*mm;
	outerMagnetRadius = (425*mm+0.2*cu->GetRadlen());
	gasRadius=innerKaptonRadius;
	noOfFieldStrips=210;  //number of field strips per layer //default:210
	halfFieldStripWidth=1.15*mm;
	electrodeRadius = 379.5*mm;
	overlapdistance = 0.*mm; //0.001*mm;
	//	siSensorThickness= 0.25*mm;  // original
	//	siSensorThickness= 0.100*mm;
	siSensorThickness= 0.3*mm;
	minimumForTest=0.025*mm;
	moveInnerSiLayer= 2.1*mm;//0.1*mm;
	moveOuterSiLayer= 6.9*mm;//15*mm;

	moveInnerSiLayerBack= 2.1*mm;//0.1*mm;
	moveOuterSiLayerBack= 6.9*mm;//15*mm;
	//	overlapdistanceBack = 28.*mm; //0.001*mm;
	overlapdistanceBack = 2.*mm; //0.001*mm;

	//	ExtraFirstPosition = (innerMagnetRadius - outerKaptonRadius - moveInnerSiLayer - moveOuterSiLayer)/4.; // for 3 additional layers (ie 2+3=5)
	//	ExtraFirstPosition = (innerMagnetRadius - outerKaptonRadius - moveInnerSiLayer - moveOuterSiLayer)/3.; // for 2 additional layers
	ExtraFirstPosition = (innerMagnetRadius - outerKaptonRadius - moveInnerSiLayer - moveOuterSiLayer)/2.; // for 1 additional layer
	Layers=false; // change also flag in Analysis.cc
	Nlayers=false; // change also flag in Analysis.cc
	NNlayers=false; // change also flag in Analysis.cc
}

void DetectorConstruction::ConstructField() 
{
  static G4TransportationManager* trMgr= 
    G4TransportationManager::GetTransportationManager();

  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager* pZeroFieldMgr= trMgr->GetFieldManager();
  G4MagneticField* pNullField= (G4MagneticField*) 0;
  pZeroFieldMgr->SetDetectorField(pNullField); 


  G4FieldManager*  magnetFieldMgr= new G4FieldManager(); 
  G4MagneticField* myField;
  G4ThreeVector  fieldVector( 0 , 0 , 1*tesla );
//  G4ThreeVector  fieldVector( 0 , 0 , 0 );
  myField = new G4UniformMagField( fieldVector ); 
  magnetFieldMgr->SetDetectorField(myField);
  magnetFieldMgr->CreateChordFinder(myField);  
//  tpcGasLogic->SetFieldManager( magnetFieldMgr, true );
  airInMagnetLogic->SetFieldManager( magnetFieldMgr, true );
}
 
G4VPhysicalVolume* DetectorConstruction::Construct()
{

	//This function is called by G4 when the detector has to be created
	//--------- Definitions of Solids, Logical Volumes, Physical Volumes ---------

  
	//------------------------------
	// World
	//------------------------------
	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.*halfWorldLength);
	G4cout << "Computed tolerance = "
			<< G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
			<< " mm" << G4endl;


	G4Box * solidWorld= new G4Box("world",halfWorldLength,halfWorldLength,halfWorldLength);
	logicWorld= new G4LogicalVolume( solidWorld, air, "World", 0, 0, 0);
  
	//  Must place the World Physical volume unrotated at (0,0,0).
	//

	G4VPhysicalVolume * physiWorld = new G4PVPlacement(0,               // no rotation
			G4ThreeVector(), // at (0,0,0)
			logicWorld,      // its logical volume
			"World",         // its name
			0,               // its mother  volume
			false,           // no boolean operations
			0);              // copy number
				 


	
	G4Color
		green(0.0,1.0,0.0),
		red(1.0,0.0,0.0),
		blue(0.0,0.0,1.0),
		brown(0.4,0.4,0.1),
		white(1.0,1.0,1.0);
	std::cout<<"test"<<std::endl;
	ConstructGeometry();
	ConstructField(); 

	//--------- Visualization attributes -------------------------------

	//logicWorld -> SetVisAttributes(new G4VisAttributes(white));
	logicWorld -> SetVisAttributes(G4VisAttributes::Invisible);
    
	//always return the physical World
	return physiWorld;
}

G4VPhysicalVolume* DetectorConstruction::ConstructGeometry()
{
	/*copy numbers
	 * 0: world with air
	 * 1: air between magnet and tpc
	 * 2: kapton of field cage
	 * 3: magnet solid
	 * 4: tpc gas
	 * 5: anode
	 * 6: cathode
	 * 7: first si-Sensor
	 * 8: second si-Sensor
	 * 9: third si-Sensor
	 * 10: fourth si-Sensor
	 * 11+: field strips of the field cage
	 */


	//creating air in magnet
	G4Tubs* airInMagnetSolid= new G4Tubs( "AirInMagnetSolid",//its name
			0, //inner radius
			innerMagnetRadius,//outer radius
			halfLPLength, //half length in z: number of layers
			0, //Starting angle in phi
			CLHEP::twopi);//Ending angle in phi, pi defined in

	airInMagnetLogic = new G4LogicalVolume( airInMagnetSolid,//its solid
			                         		     air,//its material
								     "AirInMagnetLogic");//its name

	airInMagnetPart = new G4PVPlacement( 0, //no rotation
					     G4ThreeVector(), //translation
					     airInMagnetLogic,//its logical volume
					     "AirInMagnetPart",//its name
					     logicWorld,//its mother volume
					     false,
				  	     1);//copy number

	//craeting Kapton layer
	G4Tubs* kaptonSolid = new G4Tubs( "KaptonSolid",//its name
 					innerKaptonRadius, //inner radius
					outerKaptonRadius,//outer radius
					halfLPLength, //half length in z: number of layers
					0, //Starting angle in phi
					CLHEP::twopi);//Ending angle in phi, pi defined in
	kaptonLogic = new G4LogicalVolume( kaptonSolid,//its solid
		                         		     kapton,//its material
							     "KaptonLogic");//its name
	
	kaptonPart = new G4PVPlacement( 0, //no rotation
				     G4ThreeVector(), //translation
				     kaptonLogic,//its logical volume
				     "KaptonPart",//its name
				     airInMagnetLogic,//its mother volume
				     false,
			  	     2);//copy number

	//creating Magnet	
	G4Tubs* magnetSolid = new G4Tubs( "MagnetSolid",//its name
 					innerMagnetRadius, //inner radius
					outerMagnetRadius,//outer radius
					halfLPLength, //half length in z: number of layers
					0, //Starting angle in phi
					CLHEP::twopi);//Ending angle in phi, pi defined in
	magnetLogic = new G4LogicalVolume( magnetSolid,//its solid
		                         		     cu,//its material
							     "MagnetLogic");//its name
	
	magnetPart = new G4PVPlacement( 0, //no rotation
				     G4ThreeVector(), //translation
				     magnetLogic,//its logical volume
				     "MagnetPart",//its name
				     logicWorld,//its mother volume
				     false,
			  	     3);//copy number

	//creating drift gas	
	G4Tubs* tpcGas = new G4Tubs( "TPCGas",//its name
 					0, //inner radius
					gasRadius,//outer radius
					halfLPLength, //half length in z: number of layers
					0, //Starting angle in phi
					CLHEP::twopi);//Ending angle in phi, pi defined in
	tpcGasLogic = new G4LogicalVolume( tpcGas,//its solid
		                         		     ar,//its material
							     "TPCGasLogic");//its name
	G4double maxStep = 6*mm; //small step length to create many hits
	stepLimit = new G4UserLimits(maxStep);
	tpcGasLogic ->SetUserLimits(stepLimit);

	tpcGasPart = new G4PVPlacement( 0, //no rotation
				     G4ThreeVector(), //translation
				     tpcGasLogic,//its logical volume
				     "TPCGasPart",//its name
				     airInMagnetLogic,//its mother volume
				     false,
			  	     4);//copy number

	//creating anode and kathode plane
	G4Tubs* anodeSolid = new G4Tubs( "AnodeSolid",//its name
 					0, //inner radius
					electrodeRadius,//outer radius
					1*mm, //half length in z: number of layers
					0, //Starting angle in phi
					CLHEP::twopi);//Ending angle in phi, pi defined in
	anodeLogic = new G4LogicalVolume( anodeSolid,//its solid
		                         		     cu,//its material
							     "AnodeLogic");//its name
	
	anodePart = new G4PVPlacement( 0, //no rotation
				     G4ThreeVector(0,0,-halfLPLength), //translation
				     anodeLogic,//its logical volume
				     "AnodePart",//its name
				     logicWorld,//its mother volume
				     false,
			  	     5);//copy number

	G4Tubs* cathodeSolid = new G4Tubs( "CathodeSolid",//its name
 					0, //inner radius
					electrodeRadius,//outer radius
					1*mm, //half length in z: number of layers
					0, //Starting angle in phi
					CLHEP::twopi);//Ending angle in phi, pi defined in
	cathodeLogic = new G4LogicalVolume( cathodeSolid,//its solid
		                         		     cu,//its material
							     "CathodeLogic");//its name
	
	cathodePart = new G4PVPlacement( 0, //no rotation
				     G4ThreeVector(0,0,halfLPLength), //translation
				     cathodeLogic,//its logical volume
				     "CathodePart",//its name
				     logicWorld,//its mother volume
				     false,
			  	     6);//copy number

	//creating first layer of silicon sensor
	G4Tubs* firstSiSensorSolid= new G4Tubs( "FirstSiSensorSolid",//its name
			    innerMagnetRadius-overlapdistance-siSensorThickness-moveOuterSiLayer, //inner radius
				innerMagnetRadius-overlapdistance-moveOuterSiLayer,//outer radius
				halfLPLength, //half length in z: number of layers
				0*deg, //Starting angle in phi
				90*deg);//Ending angle in phi, pi defined in

	firstSiSensorLogic = new G4LogicalVolume( firstSiSensorSolid,//its solid
				                         		     si,//its material
									     "FirstSiSensorLogic");//its name

	G4RotationMatrix * rot1  = new G4RotationMatrix();
	rot1->rotateZ(45*deg);
	firstSiSensorPart = new G4PVPlacement( rot1, //no rotation
						     G4ThreeVector(), //translation
						     firstSiSensorLogic,//its logical volume
						     "FirstSiSensorPart",//its name
						     airInMagnetLogic,//its mother volume
						     false,
					  	     7);//copy number
	if ( Layers==true ) {
	//creating first extra layer of silicon sensor
	G4Tubs* firstExtraSiSensorSolid= new G4Tubs( "FirstExtraSiSensorSolid",//its name
			        innerMagnetRadius-siSensorThickness-ExtraFirstPosition-moveOuterSiLayer, //inner radius
				innerMagnetRadius-ExtraFirstPosition-moveOuterSiLayer,//outer radius
				halfLPLength, //half length in z: number of layers
				0*deg, //Starting angle in phi
				90*deg);//Ending angle in phi, pi defined in

	firstExtraSiSensorLogic = new G4LogicalVolume( firstExtraSiSensorSolid,//its solid
				                         		     si,//its material
									     "FirstExtraSiSensorLogic");//its name

	firstExtraSiSensorPart = new G4PVPlacement( rot1, //no rotation
						     G4ThreeVector(), //translation
						     firstExtraSiSensorLogic,//its logical volume
						     "FirstExtraSiSensorPart",//its name
						     airInMagnetLogic,//its mother volume
						     false,
					  	     1000);//copy number
	}

	//creating first extra N layer of silicon sensor
	if ( Nlayers== true ) {
	G4Tubs* firstExtraNSiSensorSolid= new G4Tubs( "FirstExtraNSiSensorSolid",//its name
				innerMagnetRadius-siSensorThickness-2.*ExtraFirstPosition-moveOuterSiLayer, //inner radius
				innerMagnetRadius-2.*ExtraFirstPosition-moveOuterSiLayer,//outer radius
				halfLPLength, //half length in z: number of layers
				0*deg, //Starting angle in phi
				90*deg);//Ending angle in phi, pi defined in

	firstExtraNSiSensorLogic = new G4LogicalVolume( firstExtraNSiSensorSolid,//its solid
				                         		     si,//its material
									     "FirstExtraNSiSensorLogic");//its name

	firstExtraNSiSensorPart = new G4PVPlacement( rot1, //no rotation
						     G4ThreeVector(), //translation
						     firstExtraNSiSensorLogic,//its logical volume
						     "FirstExtraNSiSensorPart",//its name
						     airInMagnetLogic,//its mother volume
						     false,
					  	     1001);//copy number
	}

	if ( NNlayers == true ) {

	//creating first extra NN layer of silicon sensor
	G4Tubs* firstExtraNNSiSensorSolid= new G4Tubs( "FirstExtraNNSiSensorSolid",//its name
    		                innerMagnetRadius-siSensorThickness-3.*ExtraFirstPosition-moveOuterSiLayer, //inner radius
				innerMagnetRadius-3.*ExtraFirstPosition-moveOuterSiLayer,//outer radius
				halfLPLength, //half length in z: number of layers
				0*deg, //Starting angle in phi
				90*deg);//Ending angle in phi, pi defined in

	firstExtraNNSiSensorLogic = new G4LogicalVolume( firstExtraNNSiSensorSolid,//its solid
				                         		     si,//its material
									     "FirstExtraNNSiSensorLogic");//its name

	firstExtraNNSiSensorPart = new G4PVPlacement( rot1, //no rotation
						     G4ThreeVector(), //translation
						     firstExtraNNSiSensorLogic,//its logical volume
						     "FirstExtraNNSiSensorPart",//its name
						     airInMagnetLogic,//its mother volume
						     false,
					  	     1002);//copy number
	}
	//creating second layer of silicon sensor
	G4Tubs* secondSiSensorSolid= new G4Tubs( "SecondSiSensorSolid",//its name
						 outerKaptonRadius+moveInnerSiLayer, //inner radius
						 outerKaptonRadius+siSensorThickness+moveInnerSiLayer,//outer radius
						 halfLPLength, //half length in z: number of layers
						 0*deg, //Starting angle in phi
						 90*deg);//Ending angle in phi, pi defined in

	secondSiSensorLogic = new G4LogicalVolume( secondSiSensorSolid,//its solid
				                         		si,//its material
				                         		"SecondSiSensorLogic");//its name

	secondSiSensorPart = new G4PVPlacement( rot1, //no rotation
						     G4ThreeVector(), //translation
						     secondSiSensorLogic,//its logical volume
						     "SecondSiSensorPart",//its name
						     airInMagnetLogic,//its mother volume
						     false,
					  	     8);//copy number

	//creating third layer of silicon sensor
	G4Tubs* thirdSiSensorSolid= new G4Tubs( "ThirdSiSensorSolid",//its name
						outerKaptonRadius+moveInnerSiLayerBack, //inner radius
						outerKaptonRadius+siSensorThickness+moveInnerSiLayerBack,//outer radius
						halfLPLength, //half length in z: number of layers
						0*deg, //Starting angle in phi
						90*deg);//Ending angle in phi, pi defined in

	thirdSiSensorLogic = new G4LogicalVolume( thirdSiSensorSolid,//its solid
					                         si,//its material
										     "ThirdSiSensorLogic");//its name
	G4RotationMatrix * rot2  = new G4RotationMatrix();
	rot2->rotateZ(225*deg);
	thirdSiSensorPart = new G4PVPlacement( rot2, //no rotation
							     G4ThreeVector(), //translation
							     thirdSiSensorLogic,//its logical volume
							     "ThirdSiSensorPart",//its name
							     airInMagnetLogic,//its mother volume
							     false,
						  	     9);//copy number

	if (Layers==true ) {
	//creating third extra layer of silicon sensor
	G4Tubs* fourthExtraSiSensorSolid= new G4Tubs( "FourthExtraSiSensorSolid",//its name
				 innerMagnetRadius-siSensorThickness-ExtraFirstPosition-moveOuterSiLayer, //inner radius
				 innerMagnetRadius-ExtraFirstPosition-moveOuterSiLayer,//outer radius
				 halfLPLength, //half length in z: number of layers
				 0*deg, //Starting angle in phi
				 90*deg);//Ending angle in phi, pi defined in

	fourthExtraSiSensorLogic = new G4LogicalVolume( fourthExtraSiSensorSolid,//its solid
					                         si,//its material
										     "FourthExtraSiSensorLogic");//its name

	fourthExtraSiSensorPart = new G4PVPlacement(rot2, //no rotation
							     G4ThreeVector(), //translation
							     fourthExtraSiSensorLogic,//its logical volume
							     "FourthExtraSiSensorPart",//its name
							     airInMagnetLogic,//its mother volume
							     false,
						  	     1010);//copy number
	}

	if ( Nlayers==true ) {
	//creating third extra layer of silicon sensor
	G4Tubs* fourthExtraNSiSensorSolid= new G4Tubs( "FourthExtraNSiSensorSolid",//its name
				 innerMagnetRadius-siSensorThickness-2.*ExtraFirstPosition-moveOuterSiLayer, //inner radius
				 innerMagnetRadius-2.*ExtraFirstPosition-moveOuterSiLayer,//outer radius
				 halfLPLength, //half length in z: number of layers
				 0*deg, //Starting angle in phi
				 90*deg);//Ending angle in phi, pi defined in

	fourthExtraNSiSensorLogic = new G4LogicalVolume( fourthExtraNSiSensorSolid,//its solid
							 si,//its material
							 "FourthExtraNSiSensorLogic");//its name

	fourthExtraNSiSensorPart = new G4PVPlacement(rot2, //no rotation
							     G4ThreeVector(), //translation
							     fourthExtraNSiSensorLogic,//its logical volume
							     "FourthExtraNSiSensorPart",//its name
							     airInMagnetLogic,//its mother volume
							     false,
						  	     1011);//copy number
	}
	if (NNlayers==true ) {
	//creating third extra layer of silicon sensor
	G4Tubs* fourthExtraNNSiSensorSolid= new G4Tubs( "FourthExtraNNSiSensorSolid",//its name
				 innerMagnetRadius-siSensorThickness-3.*ExtraFirstPosition-moveOuterSiLayer, //inner radius
				 innerMagnetRadius-3.*ExtraFirstPosition-moveOuterSiLayer,//outer radius
				 halfLPLength, //half length in z: number of layers
				 0*deg, //Starting angle in phi
				 90*deg);//Ending angle in phi, pi defined in

	fourthExtraNNSiSensorLogic = new G4LogicalVolume( fourthExtraNNSiSensorSolid,//its solid
							si,//its material
							"FourthExtraNNSiSensorLogic");//its name

	fourthExtraNNSiSensorPart = new G4PVPlacement(rot2, //no rotation
							     G4ThreeVector(), //translation
							     fourthExtraNNSiSensorLogic,//its logical volume
							     "FourthExtraNNSiSensorPart",//its name
							     airInMagnetLogic,//its mother volume
							     false,
						  	     1012);//copy number

	}
	//creating fourth layer of silicon sensor
	G4Tubs* fourthSiSensorSolid= new G4Tubs( "FourthSiSensorSolid",//its name
						 innerMagnetRadius-overlapdistanceBack-siSensorThickness-moveOuterSiLayerBack, //inner radius
						 innerMagnetRadius-overlapdistanceBack-moveOuterSiLayerBack,//outer radius
				 halfLPLength, //half length in z: number of layers
				 0*deg, //Starting angle in phi
				 90*deg);//Ending angle in phi, pi defined in

	fourthSiSensorLogic = new G4LogicalVolume( fourthSiSensorSolid,//its solid
					                         si,//its material
										     "FourthSiSensorLogic");//its name

	fourthSiSensorPart = new G4PVPlacement(rot2, //no rotation
							     G4ThreeVector(), //translation
							     fourthSiSensorLogic,//its logical volume
							     "FourthSiSensorPart",//its name
							     airInMagnetLogic,//its mother volume
							     false,
						  	     10);//copy number



	std::cout<<" 1st inner radius "<<innerMagnetRadius-overlapdistance-siSensorThickness-moveOuterSiLayer
		 <<" 1st outer radius "<<innerMagnetRadius-overlapdistance-moveOuterSiLayer<<std::endl;
	std::cout<<" 2nd inner radius "<<outerKaptonRadius+moveInnerSiLayer
		 <<" 2nd outer radius "<<outerKaptonRadius+siSensorThickness+moveInnerSiLayer<<std::endl;
	std::cout<<" 3rd inner radius "<<outerKaptonRadius+moveInnerSiLayerBack
		 <<" 3rd outer radius "<<outerKaptonRadius+siSensorThickness+moveInnerSiLayerBack<<std::endl;
	std::cout<<" 4th inner radius "<<innerMagnetRadius-overlapdistanceBack-siSensorThickness-moveOuterSiLayerBack
		 <<" 4th outer radius "<<innerMagnetRadius-overlapdistanceBack-moveOuterSiLayerBack<<std::endl;
	std::cout<<" innerMagnetRadius "<<innerMagnetRadius<<" , outerMagnetRadius "<<outerMagnetRadius<<std::endl;
	std::cout<<" innerKaptonRadius "<<innerKaptonRadius<<" , outerKaptonRadius "<<outerKaptonRadius<<std::endl;

	//creating field strips
	G4Tubs* innerFieldStripsSolid = new G4Tubs( "InnerFieldStripsSolid", innerInnerFieldStripsRadius , outerInnerFieldStripsRadius , halfFieldStripWidth , 0 ,  CLHEP::twopi);
	G4Tubs* outerFieldStripsSolid = new G4Tubs( "OuterFieldStripsSolid", innerOuterFieldStripsRadius , outerOuterFieldStripsRadius , halfFieldStripWidth , 0 , CLHEP::twopi);
	innerFieldStripsLogic = new G4LogicalVolume(innerFieldStripsSolid , cu , "InnerFieldStripsLogic");
	outerFieldStripsLogic = new G4LogicalVolume(outerFieldStripsSolid , cu , "OuterFieldStripsLogic");
	//Translation of one Layer with respect previous Layer
	G4double position;
	G4int stripCopyNum = 10;
	for ( int layerIdx = 0 ; layerIdx < (noOfFieldStrips) ; layerIdx++ )
	{
		position = (layerIdx+0.5)*halfFieldStripsLength*2/noOfFieldStrips-halfLPLength;
		new G4PVPlacement(0,G4ThreeVector(0,0,position),innerFieldStripsLogic,"InnerFieldStrips",kaptonLogic,false,++stripCopyNum);
		position = (layerIdx+1)*halfFieldStripsLength*2/noOfFieldStrips-halfLPLength;
		new G4PVPlacement(0,G4ThreeVector(0,0,position),outerFieldStripsLogic,"OuterFieldStrips",kaptonLogic,false,++stripCopyNum);
	}
	
	G4Color
		green(0.0,1.0,0.0),
		red(1.0,0.0,0.0),
		blue(0.0,0.0,1.0),
		brown(0.4,0.4,0.1),
		white(1.0,1.0,1.0);

	G4VisAttributes* innerFieldStripsAttributes = new G4VisAttributes(brown);
	innerFieldStripsAttributes->SetForceSolid(true);
	innerFieldStripsLogic->SetVisAttributes(innerFieldStripsAttributes);
	G4VisAttributes* outerFieldStripsAttributes = new G4VisAttributes(red);
	outerFieldStripsAttributes->SetForceSolid(true);
	outerFieldStripsLogic->SetVisAttributes(outerFieldStripsAttributes);
	G4VisAttributes* kaptonAttributes = new G4VisAttributes(green);
	kaptonLogic->SetVisAttributes(kaptonAttributes);
	G4VisAttributes* tpcGasAttributes = new G4VisAttributes(blue);
	tpcGasLogic->SetVisAttributes(tpcGasAttributes);
	G4VisAttributes* magnetAttributes = new G4VisAttributes(blue);
	magnetLogic->SetVisAttributes(magnetAttributes);
	G4VisAttributes* anodeAttributes = new G4VisAttributes(red);
	anodeLogic->SetVisAttributes(anodeAttributes);
	G4VisAttributes* cathodeAttributes = new G4VisAttributes(red);
	cathodeLogic->SetVisAttributes(cathodeAttributes);
	G4VisAttributes* airInMagnetAttributes = new G4VisAttributes(red);
	airInMagnetLogic->SetVisAttributes(airInMagnetAttributes);
	G4VisAttributes* firstSiSensorAttributes = new G4VisAttributes(white);//red);
	firstSiSensorLogic->SetVisAttributes(firstSiSensorAttributes);
	G4VisAttributes* secondSiSensorAttributes = new G4VisAttributes(white);//blue);
	secondSiSensorLogic->SetVisAttributes(secondSiSensorAttributes);
	G4VisAttributes* thirdSiSensorAttributes = new G4VisAttributes(white);//brown);
	thirdSiSensorLogic->SetVisAttributes(thirdSiSensorAttributes);
	G4VisAttributes* fourthSiSensorAttributes = new G4VisAttributes(white);//green);
	fourthSiSensorLogic->SetVisAttributes(fourthSiSensorAttributes);

	if ( Layers==true ) {
	  G4VisAttributes* firstExtraSiSensorAttributes = new G4VisAttributes(white);
	  firstExtraSiSensorLogic->SetVisAttributes(firstExtraSiSensorAttributes);
	  G4VisAttributes* fourthExtraSiSensorAttributes = new G4VisAttributes(white);
	  fourthExtraSiSensorLogic->SetVisAttributes(fourthExtraSiSensorAttributes);
	}
	if (Nlayers==true ) {
	  G4VisAttributes* firstExtraNSiSensorAttributes = new G4VisAttributes(white);
	  firstExtraNSiSensorLogic->SetVisAttributes(firstExtraNSiSensorAttributes);
	  G4VisAttributes* fourthExtraNSiSensorAttributes = new G4VisAttributes(white);
	  fourthExtraNSiSensorLogic->SetVisAttributes(fourthExtraNSiSensorAttributes);
	}
	if (NNlayers==true ) {
	  G4VisAttributes* firstExtraNNSiSensorAttributes = new G4VisAttributes(white);
	  firstExtraNNSiSensorLogic->SetVisAttributes(firstExtraNNSiSensorAttributes);
	  G4VisAttributes* fourthExtraNNSiSensorAttributes = new G4VisAttributes(white);
	  fourthExtraNNSiSensorLogic->SetVisAttributes(fourthExtraNNSiSensorAttributes);
	}
	return kaptonPart;
}

