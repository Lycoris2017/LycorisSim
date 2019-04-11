// $Id: DetectorConstruction.cc 94 2010-01-26 13:18:30Z adotti $
/**
 * @file
 * @brief Implements mandatory user class DetectorConstruction.
 *
 * Origin @Dimitra
 * Modify @Mengqing, 2019-01-28
 * Dev    @Mengqing, 2019-04-11
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

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

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
	mNumOfLayers = 6;
	mTRKSiliconXY = 92*mm;
	mSiSensorThickness = 0.320*mm;

	mInnerMagnetRadius = 425*mm;
	mOuterMagnetRadius = (425*mm+0.2*cu->GetRadlen()); // magnet wall is 20% X0 in copper thick

	//--------- Material definition ---------
	DefineMaterials(); // mengqing

	//--------- Sizes of the principal geometrical components (solids)  ---------
	ComputeParameters(); 

	//--------- Extra: Mengqing  ---------
	fCheckOverlaps=true;
}
 
DetectorConstruction::~DetectorConstruction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

void DetectorConstruction::ComputeParameters(){ 
	/*
	 * This function defines the defaults of the geometry parameters
	 * TODO@too dirty, to modify!
	 */
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
	
	gasRadius=innerKaptonRadius;
	noOfFieldStrips=210;  //number of field strips per layer //default:210
	halfFieldStripWidth=1.15*mm;
	electrodeRadius = 379.5*mm;
	overlapdistance = 0.*mm; //0.001*mm;
	//	mSiSensorThickness= 0.25*mm;  // original
	//	mSiSensorThickness= 0.100*mm;
	//  mSiSensorThickness= 0.3*mm;
	minimumForTest=0.025*mm;
	moveInnerSiLayer= 2.1*mm;//0.1*mm;
	moveOuterSiLayer= 6.9*mm;//15*mm;

	moveInnerSiLayerBack= 2.1*mm;//0.1*mm;
	moveOuterSiLayerBack= 6.9*mm;//15*mm;
	//	overlapdistanceBack = 28.*mm; //0.001*mm;
	overlapdistanceBack = 2.*mm; //0.001*mm;

	//	ExtraFirstPosition = (mInnerMagnetRadius - outerKaptonRadius - moveInnerSiLayer - moveOuterSiLayer)/4.; // for 3 additional layers (ie 2+3=5)
	//	ExtraFirstPosition = (mInnerMagnetRadius - outerKaptonRadius - moveInnerSiLayer - moveOuterSiLayer)/3.; // for 2 additional layers
	ExtraFirstPosition = (mInnerMagnetRadius - outerKaptonRadius - moveInnerSiLayer - moveOuterSiLayer)/2.; // for 1 additional layer
	Layers=false; // change also flag in Analysis.cc
	Nlayers=false; // change also flag in Analysis.cc
	NNlayers=false; // change also flag in Analysis.cc
	
}

G4VPhysicalVolume* DetectorConstruction::Construct(){
	/*
	 * Main function to construct the geometry	
	 * This function is called by G4 when the detector has to be created
	 * Definitions of Solids, Logical Volumes, Physical Volumes
	 */
	
	/*--------- Material definition ---------*/
    DefineMaterials();

    /*--------- Volumes definition ---------*/
    return DefineVolumes();
	
}


 
G4VPhysicalVolume* DetectorConstruction::DefineVolumes(){
	/*
	 * Define all the boxes:
	 * Version 1: 
	 *   World Box <- upstream Box + downstream Box
	 * Version 2:
	 *   World Box <- Magnet World <- upstream Box + downstream Box
	 */
	
	/*--------- World Box ---------*/
	G4Box * solidWorld= new G4Box("world", halfWorldLength, halfWorldLength, halfWorldLength);
	// World Logical Volume (Associate with the world box)
	flogicWorld= new G4LogicalVolume( solidWorld, air, "World", 0, 0, 0);
	// Placement of the World Logical Volume
	fphysiWorld = new G4PVPlacement(0,                // no rotation
	                                G4ThreeVector(),  // at (0,0,0)
	                                flogicWorld,    // its logical volume
	                                "World",          // its name
	                                0,                // its mother volume 
	                                false,			  // no boolean operations
	                                0, 				  // copy number
	                                fCheckOverlaps);  // flag on/off check overlaps

	/*--------- Upstream Box ---------*/
	// creating 1st layer of silicon sensor
	G4Tubs* firstSiSensorSolid= new G4Tubs( "FirstSiSensorSolid",//its name
	                                        mInnerMagnetRadius-overlapdistance-mSiSensorThickness-moveOuterSiLayer, //inner radius
	                                        mInnerMagnetRadius-overlapdistance-moveOuterSiLayer,//outer radius
	                                        halfLPLength, //half length in z: number of layers
	                                        0*deg, //Starting angle in phi
	                                        90*deg);//Ending angle in phi, pi defined in
	
	firstSiSensorLogic = new G4LogicalVolume( firstSiSensorSolid,//its solid
	                                          si,//its material
	                                          "FirstSiSensorLogic");//its name
	
	G4RotationMatrix * rot1  = new G4RotationMatrix();
	rot1->rotateZ(45*deg); // very large degree!
	firstSiSensorPart = new G4PVPlacement( rot1, //no rotation
	                                       G4ThreeVector(), //translation
	                                       firstSiSensorLogic,//its logical volume
	                                       "FirstSiSensorPart",//its name
	                                       airInMagnetLogic,//its mother volume
	                                       false,
	                                       7, fCheckOverlaps);//copy number

	// creating 2nd layer of silicon sensor
	G4Tubs* secondSiSensorSolid= new G4Tubs( "SecondSiSensorSolid",//its name
	                                         outerKaptonRadius+moveInnerSiLayer, //inner radius
	                                         outerKaptonRadius+mSiSensorThickness+moveInnerSiLayer,//outer radius
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
	                                        8, fCheckOverlaps);//copy number
	
	/*--------- Downstream Box ---------*/
	// creating 3rd layer of silicon sensor
	G4Tubs* thirdSiSensorSolid= new G4Tubs( "ThirdSiSensorSolid",//its name
	                                        outerKaptonRadius+moveInnerSiLayerBack, //inner radius
	                                        outerKaptonRadius+mSiSensorThickness+moveInnerSiLayerBack,//outer radius
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
	                                       9, fCheckOverlaps);//copy number
	
	// creating 4th layer of silicon sensor
	G4Tubs* fourthSiSensorSolid= new G4Tubs( "FourthSiSensorSolid",//its name
	                                         mInnerMagnetRadius-overlapdistanceBack-mSiSensorThickness-moveOuterSiLayerBack, //inner radius
	                                         mInnerMagnetRadius-overlapdistanceBack-moveOuterSiLayerBack,//outer radius
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
	                                       10, fCheckOverlaps);//copy number
	


	//-- For visulization:
	
	G4Color
		green(0.0,1.0,0.0),
		red(1.0,0.0,0.0),
		blue(0.0,0.0,1.0),
		brown(0.4,0.4,0.1),
		white(1.0,1.0,1.0);

	/*
	std::cout<<"test"<<std::endl;
	ConstructGeometry();
	ConstructField();
	*/
	return fphysiWorld;

}


void DetectorConstruction::ConstructField() {
	/* 
	 * 1 Tesla B field.
	 */
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
	                                      mInnerMagnetRadius,//outer radius
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
	                                     flogicWorld,//its mother volume
	                                     false,
	                                     1, fCheckOverlaps);//copy number
	
	//craeting Kapton layer -- TPC

	//creating Magnet	

	//creating drift gas  -- TPC
	//creating anode and kathode plane -- TPC

	//creating first layer of silicon sensor

	// tracker box
	auto trkBox = new G4Box("", mTRKSiliconXY/2, mTRKSiliconXY/2, mSiSensorThickness);
	// tracker logical volume
	mpTrkLogic = new  G4LogicalVolume( trkBox, si, "Tracker");//its name
	// tracker physical volume
	mpTrkPhysical = new G4PVPlacement( 0,                  //no rotation
	                                   G4ThreeVector(),    //translation
	                                   mpTrkLogic,         //its logical volume
	                                   "FirstSiSensorPart",//its name
	                                   airInMagnetLogic,   //its mother volume
	                                   false,
	                                   7, fCheckOverlaps); //copy number
	
	//creating first extra layer of silicon sensor
	//creating first extra N layer of silicon sensor
	//creating first extra NN layer of silicon sensor

	//creating second layer of silicon sensor

	//creating third layer of silicon sensor
	//creating third extra layer of silicon sensor
	//creating third extra layer of silicon sensor
	//creating third extra layer of silicon sensor
	//creating fourth layer of silicon sensor

	//creating field strips -- TPC field cage copper strips

    //Translation of one Layer with respect previous Layer

	return airInMagnetPart;
}

void DetectorConstruction::PrintParameters()
{
	// was part of Construct() from Dimitra's code
	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.*halfWorldLength);
	G4cout << "Computed tolerance = "
	       << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
	       << " mm" << G4endl;

}
