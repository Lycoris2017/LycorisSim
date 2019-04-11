// $Id: DetectorConstruction.hh 33 2010-01-14 17:08:18Z adotti $
/*
 * Origin @Dimitra
 * Modify @Mengqing, 2019-01-28
 */

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

/**
 * @file
 * @brief Defines mandatory user class DetectorConstruction.
 */

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4StepLimiter.hh"
#include "G4Colour.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UserLimits;
//class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*!
  \brief This mandatory user class defines the geometry.
  
  It is responsible for
  - Definition of material, and
  - Construction of geometry
  
  \sa Construct()
*/
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
	//! Constructor
	DetectorConstruction();
	//! Destructor
	~DetectorConstruction();
public:
	//! Construct geometry of the setup
	G4VPhysicalVolume* Construct();
	
private:

	
	//! define needed materials
	void DefineMaterials();
	//! initialize geometry parameters
	void ComputeParameters();
	//! Construct geometry
	G4VPhysicalVolume* ConstructGeometry();
	G4VPhysicalVolume* DefineVolumes();
	
	void ConstructField();
	void PrintParameters();
private:
	// Extra Mengqing
	bool fCheckOverlaps;
	
	//! \name Materials
	//@{
	G4Material* air;
	G4Material* vacuum;
	G4Material* steel;
	G4Material* kapton;
	G4Material* cu;
	G4Material* ar;
	G4Material* si;
	//@}
	
	//! \name global mother volume
	//@{
	G4LogicalVolume * flogicWorld;
	G4VPhysicalVolume * fphysiWorld;
	G4double halfWorldLength;
	//@}
	
	//! \name Parameters for LP
	//@{
	bool Layers;
	bool NNlayers;
	bool Nlayers;
	
	G4double moveInnerSiLayerBack;
	G4double moveOuterSiLayerBack;
	G4double overlapdistanceBack;
	
	G4int noOfFieldStrips;
	G4double halfFieldStripWidth;
	G4double halfFieldStripsLength;
	G4double overlapdistance;

	G4double mSiSensorThickness;
	G4double mTRKSiliconXY;
	G4int    mNumOfLayers;
	
	G4double minimumForTest;
	G4double moveInnerSiLayer;
	G4double moveOuterSiLayer;
	G4double makeVisibleGraphics;
	G4double ExtraFirstPosition;
	G4double innerKaptonRadius;
	G4double outerKaptonRadius;
	G4double gasRadius;
	G4double innerInnerFieldStripsRadius;
	G4double outerInnerFieldStripsRadius;
	G4double innerOuterFieldStripsRadius;
	G4double outerOuterFieldStripsRadius;
	G4double electrodeRadius;
	G4double halfLPLength;
	/*	G4VPhysicalVolume* kaptonPart;
	G4VPhysicalVolume* magnetPart;
	G4VPhysicalVolume* tpcGasPart;
	G4VPhysicalVolume* anodePart;
	G4VPhysicalVolume* cathodePart;*/
	G4VPhysicalVolume* airInMagnetPart;

	G4VPhysicalVolume* mpTrkPhysical;

	G4VPhysicalVolume* firstSiSensorPart;
	G4VPhysicalVolume* secondSiSensorPart;
	G4VPhysicalVolume* thirdSiSensorPart;
	G4VPhysicalVolume* fourthSiSensorPart;
	G4VPhysicalVolume* firstExtraSiSensorPart;
	G4VPhysicalVolume* fourthExtraSiSensorPart;
	G4VPhysicalVolume* firstExtraNSiSensorPart;
	G4VPhysicalVolume* fourthExtraNSiSensorPart;
	G4VPhysicalVolume* firstExtraNNSiSensorPart;
	G4VPhysicalVolume* fourthExtraNNSiSensorPart;
	
	
	//@}
	G4LogicalVolume* kaptonLogic;
	G4LogicalVolume* magnetLogic;
	G4LogicalVolume* tpcGasLogic;
	G4LogicalVolume* anodeLogic;
	G4LogicalVolume* cathodeLogic;
	G4LogicalVolume* innerFieldStripsLogic;
	G4LogicalVolume* outerFieldStripsLogic;
	G4LogicalVolume* airInMagnetLogic;

	G4LogicalVolume* mpTrkLogic;
	
	G4LogicalVolume* firstSiSensorLogic;
	G4LogicalVolume* secondSiSensorLogic;
	G4LogicalVolume* thirdSiSensorLogic;
	G4LogicalVolume* fourthSiSensorLogic;
	
	G4LogicalVolume* firstExtraSiSensorLogic;
	G4LogicalVolume* fourthExtraSiSensorLogic;
	G4LogicalVolume* firstExtraNSiSensorLogic;
	G4LogicalVolume* fourthExtraNSiSensorLogic;
	G4LogicalVolume* firstExtraNNSiSensorLogic;
	G4LogicalVolume* fourthExtraNNSiSensorLogic;
	
	
	
	//! \name Parameters for Magnet
	//@{
	G4double mInnerMagnetRadius;
	G4double mOuterMagnetRadius;
	//@}
	G4UserLimits* stepLimit;             // pointer to user step limits
};


#endif
