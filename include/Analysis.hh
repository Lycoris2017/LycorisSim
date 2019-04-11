/*
 * @file Analysis.hh
 *
 * @date: 9 Feb 2010
 * @author: adotti
 *
 * \brief Analysis class
 */

#ifndef ANALYSIS_HH_
#define ANALYSIS_HH_

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//ROOT
#include "TH2F.h"
#include "TH1F.h"
#include "TTree.h"
#include "TVector3.h"

/*!
 * \brief Analysis class
 * This class contains the code to collect information from
 * the different UserActions.
 * The class is designed as a singleton.
 * To access it you need to use:
 * Analysis* analysis = Analysis::GetInstance()
 */
class Analysis {
public:
	//! Singleton pattern
	static Analysis* GetInstance() {
		if ( Analysis::singleton == NULL ) Analysis::singleton = new Analysis();
		return Analysis::singleton;
	}
	//! destructor
	virtual ~Analysis() {};
	//! Should be called at the beginning of an event
	void PrepareNewEvent(const G4Event* anEvent);
	//! Should be called at the end of an event
	void EndOfEvent(const G4Event* anEvent);
	//! Should be called at the beginning of a run
	void PrepareNewRun(const G4Run* aRun);
	//! Should be called at the end of a run
	void EndOfRun(const G4Run* aRun);
	//! Increase number of secondaries
	void AddSecondary( G4int num ) { thisEventSecondaries += num; };
	void AddChargedSecondary( G4int num ) { noOfChargedPart += num; };
	void AddChargedSecondaryEnetringVolume( G4int num ) { noOfChargedPartEnteringVolume += num; };
	void AddEDepInGas( G4double eDep) { eDepTPCGas+=eDep/MeV; };
	void AddTrackLength( G4double length) { trackLength+=length/mm; };
	void AddEnergyOfChargedPart( G4double energyOfSec );
	void InsertEntryTPC( G4ThreeVector pos, G4double eDep);
	void InsertExitTPC( G4ThreeVector pos, G4double eDep);
	void InsertEntryGas( G4ThreeVector pos, G4double eDep);
	void InsertExitGas( G4ThreeVector pos, G4double eDep);
	void InsertEntryMagnet( G4ThreeVector pos, G4double eDep);
	void InsertExitMagnet( G4ThreeVector pos, G4double eDep);
	void InsertExitWorld( G4ThreeVector pos, G4double eDep);

        void InsertExtraFirstSi( G4ThreeVector pos, G4double eDep);
  	void InsertExtraThirdSi( G4ThreeVector pos, G4double eDep);
        void InsertExtraNFirstSi( G4ThreeVector pos, G4double eDep);
  	void InsertExtraNThirdSi( G4ThreeVector pos, G4double eDep);
  void InsertExtraNNFirstSi( G4ThreeVector pos, G4double eDep);
  void InsertExtraNNThirdSi( G4ThreeVector pos, G4double eDep);

	//add a simulated hit to the output
	void AddSimTrackHit(G4ThreeVector pos, G4double  eDep, G4int trackID);

private:
	//! Private construtor: part of singleton pattern
	Analysis();
	//! Singleton static instance
	static Analysis* singleton;
	G4ThreeVector entryPosition;
	TVector3* entryPositionT;
	G4ThreeVector entryPositionGas;
	TVector3* entryPositionGasT;
	G4ThreeVector entryPositionMagnet;
	TVector3* entryPositionMagnetT;
	G4ThreeVector exitPosition;
	TVector3* exitPositionT;
	G4ThreeVector exitPositionGas;
	TVector3* exitPositionGasT;
	G4ThreeVector exitPositionMagnet;
	TVector3* exitPositionMagnetT;
	G4ThreeVector exitPositionWorld;
	TVector3* exitPositionWorldT;

	G4ThreeVector entryPositionExtraFirstSi;
	TVector3* entryPositionExtraFirstSiT;
	G4ThreeVector entryPositionExtraThirdSi;
	TVector3* entryPositionExtraThirdSiT;

	G4ThreeVector entryPositionExtraNFirstSi;
	TVector3* entryPositionExtraNFirstSiT;
	G4ThreeVector entryPositionExtraNThirdSi;
	TVector3* entryPositionExtraNThirdSiT;

	G4ThreeVector entryPositionExtraNNFirstSi;
	TVector3* entryPositionExtraNNFirstSiT;
	G4ThreeVector entryPositionExtraNNThirdSi;
	TVector3* entryPositionExtraNNThirdSiT;

	G4double eDepTPCGas;
	G4double trackLength;
	G4double entryEnergyMagnet;
	G4double exitEnergyMagnet;
	G4double entryEnergyGas;
	G4double exitEnergyGas;
	G4double entryEnergy;
	G4double exitEnergy;
        G4double exitEnergyWorld;

	G4double entryEnergyExtraFirstSi;
	G4double entryEnergyExtraThirdSi;

	G4double entryEnergyExtraNFirstSi;
	G4double entryEnergyExtraNThirdSi;
  G4double entryEnergyExtraNNFirstSi;
  G4double entryEnergyExtraNNThirdSi;

	G4int noOfChargedPartEnteringVolume;
	G4int noOfChargedPart;
	G4int thisEventSecondaries;
	G4int thisRunTotSecondaries;
	std::vector<G4double>* energyOfSecChargedPart;
	G4bool enterMagnet;
	TTree* tree;
	std::vector<G4double>* hitXPositions;
	std::vector<G4double>* hitYPositions;
	std::vector<G4double>* hitZPositions;
	std::vector<G4double>* hitEnergy;
	G4int eventNumber;


};

#endif /* ANALYSIS_HH_ */
