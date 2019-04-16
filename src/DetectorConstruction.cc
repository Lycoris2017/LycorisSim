// $Id: DetectorConstruction.cc 94 2010-01-26 13:18:30Z adotti $
/**
 * @file
 * @brief Implements mandatory user class DetectorConstruction.
 *
 * Origin @Dimitra
 * Modify @Mengqing, 2019-01-28
 * Dev    @Mengqing, 2019-04-11
 * Note by Mengqing:
 * -- Magnetic field is hard coded in it, but one can set it to 0 tesla.
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
    //--------- Material definition ---------
    DefineMaterials(); // mengqing

	mNumOfLayers = 6;
	mSiSensorYZ = 92*mm;
	mSiSensorThickness = 0.320*mm;

	mCassetteX=33*mm;
	mCassetteY=121*mm;
	mCassetteZ=321*mm;

	mInnerMagnetRadius = 425*mm;
	mOuterMagnetRadius = (425*mm+0.2*cu->GetRadlen()); // magnet wall is 20% X0 in copper thick

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
	air= man->FindOrBuildMaterial("G4_AIR");
	ar = man->FindOrBuildMaterial("G4_Ar");
	cu = man->FindOrBuildMaterial("G4_Cu");
	//kapton = man->FindOrBuildMaterial("G4_KAPTON");
	si = man->FindOrBuildMaterial("G4_Si");
}

void DetectorConstruction::ComputeParameters(){
	/*
	 * This function defines the defaults of the geometry parameters
	 * TODO@too dirty, to modify!
	 */
	// ** world **
	halfWorldLength = 5.3* m;
	halfLPLength = 305*mm;

	//minimumForTest=0.025*mm;

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

	/* copy numbers (unique index)
	 * 0: world with air
	 * 1: air in magnet
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

	/*--------- World Box ---------*/
	G4Box * solidWorld= new G4Box("world", halfWorldLength, halfWorldLength, halfWorldLength);
	// World Logical Volume (Associate with the world box)
	pLogicWorld= new G4LogicalVolume( solidWorld, air, "World", 0, 0, 0);
	// Placement of the World Logical Volume
	pPhysiWorld = new G4PVPlacement(0,                // no rotation
	                                G4ThreeVector(),  // at (0,0,0)
	                                pLogicWorld,      // its logical volume
	                                "World",          // its name
	                                0,                // its mother volume
	                                false,			  // pMany, not used, Set it to false...
	                                0, 				  // copy number, unique arbitrary index
	                                fCheckOverlaps);  // optional overlap check

	/*--------- Magnet : 1st daughter of mother World  ---------*/
	G4Tubs* magnetSolid = new G4Tubs( "MagnetSolid",      //its name
	                                  mInnerMagnetRadius, //inner radius
	                                  mOuterMagnetRadius, //outer radius
	                                  halfLPLength,      //half length in z: number of layers
	                                  0,                 //Starting angle in phi
	                                  CLHEP::twopi);     //Ending angle in phi, pi defined in

	pMagnetLogic = new G4LogicalVolume( magnetSolid, cu, "MagnetLogic");//its name
	pMagnetPhysi= new G4PVPlacement( 0, G4ThreeVector(), pMagnetLogic, "MagnetPart",
	                                pLogicWorld, false, 1, fCheckOverlaps);

	//create a cylinder magnetic field region
	G4Tubs* airInMagnetSolid= new G4Tubs( "AirInMagnetSolid",// name
	                                      0,                 // inner radius
	                                      mInnerMagnetRadius,// outer radius
	                                      halfLPLength,      // half height (z axis)
	                                      0,                 // Starting angle in phi
	                                      CLHEP::twopi);     // Ending angle in phi, pi defined in

	pAirInMagnetLogic = new G4LogicalVolume( airInMagnetSolid, air, "AirInMagnetLogic");

	pAirInMagnetPhysi = new G4PVPlacement(0, G4ThreeVector(), pAirInMagnetLogic,"AirInMagnetPart",
	                                      pLogicWorld,  false, 2, fCheckOverlaps);

	/*--------- Upstream Cassette : daughter of AirInMagnetLogic ---------*/

	auto CassetteSolid = new G4Box("Cassette Solid", mCassetteX/2, mCassetteY/2, mCassetteZ/2); //half-length
	pCassetteLogic = new G4LogicalVolume( CassetteSolid, air, "Cassette Logic");
	pCassettePhysi = new G4PVPlacement( 0, G4ThreeVector(), pCassetteLogic, "Cassette Physical",
	                                    pAirInMagnetLogic,  false, 3, fCheckOverlaps);


	// Create sensor layers
	auto SensorSolid = new G4Box("Sensor solid", mSiSensorThickness/2, mSiSensorYZ/2, mSiSensorYZ/2 );
	pSensorLogic = new G4LogicalVolume( SensorSolid, si, "Cassette Logic");
	// TODO: orientation and placement!

	/*--------- Downstream Box ---------*/



	//** Construct an uniform local B field
	ConstructField();

	//** Visualization attributes
	G4Color
		green(0.0,1.0,0.0),
		red(1.0,0.0,0.0),
		blue(0.0,0.0,1.0),
		brown(0.4,0.4,0.1),
		white(1.0,1.0,1.0);

	pLogicWorld -> SetVisAttributes(G4VisAttributes::Invisible);
	auto MagnetVisAtt = new G4VisAttributes(blue);
	auto CassetteVisAtt = new G4VisAttributes(white);
	auto SensorVisAtt = new G4VisAttributes(red);

	pMagnetLogic  -> SetVisAttributes(MagnetVisAtt);
	pCassetteLogic-> SetVisAttributes(CassetteVisAtt);
	pSensorLogic  -> SetVisAttributes(SensorVisAtt);

	std::cout<<"test"<<std::endl;
	// PrintParameters();

	//always return the physical World
	return pPhysiWorld;

}


void DetectorConstruction::ConstructField() {
	/*
	 * 1 Tesla B field.
	 * using an uniform LOCAL B field to represent the inner-solenoid.
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

	// create a uniform magnetic field along Z axis
	myField = new G4UniformMagField( fieldVector );
	// set this field as the global field
	magnetFieldMgr->SetDetectorField(myField);
	// prepare the propagation with dfault parameters and other choices
	magnetFieldMgr->CreateChordFinder(myField);

	pAirInMagnetLogic->SetFieldManager( magnetFieldMgr, true );

}



// G4VPhysicalVolume* DetectorConstruction::ConstructGeometry()
// {

// }

void DetectorConstruction::PrintParameters()
{
	// was part of Construct() from Dimitra's code
	G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.*halfWorldLength);
	G4cout << "Computed tolerance = "
	       << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()/mm
	       << " mm" << G4endl;

}
