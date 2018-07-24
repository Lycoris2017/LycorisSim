/*
 * Analysis.cc
 *
 *  Created on: 9 Feb 2010
 *      Author: adotti
 */

#include "Analysis.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
//ROOT Stuff
#include "TProfile.h"
#include "TFile.h"

Analysis* Analysis::singleton = 0;

Analysis::Analysis() :
	thisEventSecondaries(0),
	thisRunTotSecondaries(0)
{
}

void Analysis::PrepareNewEvent(const G4Event* /*anEvent*/)
{
	//Reset variables relative to this event
	thisEventSecondaries = 0;
	eDepTPCGas = 0;
	trackLength = 0;
	noOfChargedPartEnteringVolume = 0;
	noOfChargedPart = 0;
	enterMagnet = true;
	eventNumber++;
}


// The TPC stuff are relevant to the first 2 Si sensors the beam sees

void Analysis::InsertEntryTPC( G4ThreeVector pos, G4double eDep)
{
		entryPosition = pos/mm;
		entryEnergy = eDep/MeV;
}

void Analysis::InsertExitTPC( G4ThreeVector pos, G4double eDep)
{
		exitPosition = pos/mm;
		exitEnergy = eDep/MeV;
}


// The Gas stuff are relevant to the last 2 Si sensors the beam particles see

void Analysis::InsertEntryGas( G4ThreeVector pos, G4double eDep)
{
		entryPositionGas = pos/mm;
		entryEnergyGas = eDep/MeV;
}

void Analysis::InsertExitGas( G4ThreeVector pos, G4double eDep)
{
		exitPositionGas = pos/mm;
		exitEnergyGas = eDep/MeV;
}

// The first Extra sensor stuff are relevant to an extra layer added in the front
// and thirs Extra relevant to an extra layer added in the back
void Analysis::InsertExtraFirstSi( G4ThreeVector pos, G4double eDep)
{
		entryPositionExtraFirstSi = pos/mm;
		entryEnergyExtraFirstSi = eDep/MeV;
}

void Analysis::InsertExtraThirdSi( G4ThreeVector pos, G4double eDep)
{
		entryPositionExtraThirdSi = pos/mm;
		entryEnergyExtraThirdSi = eDep/MeV;
}
// N Extra
void Analysis::InsertExtraNFirstSi( G4ThreeVector pos, G4double eDep)
{
		entryPositionExtraNFirstSi = pos/mm;
		entryEnergyExtraNFirstSi = eDep/MeV;
}
void Analysis::InsertExtraNThirdSi( G4ThreeVector pos, G4double eDep)
{
		entryPositionExtraNThirdSi = pos/mm;
		entryEnergyExtraNThirdSi = eDep/MeV;
}

// NN Extra
void Analysis::InsertExtraNNFirstSi( G4ThreeVector pos, G4double eDep)
{
		entryPositionExtraNNFirstSi = pos/mm;
		entryEnergyExtraNNFirstSi = eDep/MeV;
}
void Analysis::InsertExtraNNThirdSi( G4ThreeVector pos, G4double eDep)
{
		entryPositionExtraNNThirdSi = pos/mm;
		entryEnergyExtraNNThirdSi = eDep/MeV;
}


void Analysis::InsertEntryMagnet( G4ThreeVector pos, G4double eDep)
{
	//check whether it is the frontside or backside of the magnet
	if (enterMagnet == true)
	{
		entryPositionMagnet = pos/mm;
		entryEnergyMagnet = eDep/MeV;
	}
	enterMagnet = false;
}

void Analysis::InsertExitMagnet( G4ThreeVector pos, G4double eDep)
{

		exitPositionMagnet = pos/mm;
		exitEnergyMagnet = eDep/MeV;

}

void Analysis::InsertExitWorld( G4ThreeVector pos, G4double eDep)
{

		exitPositionWorld = pos/mm;
		exitEnergyWorld = eDep/MeV;

}


// what is this? The SimTrackHit is in the TPC region. Where is this position? In the entrance?
void Analysis::AddSimTrackHit( G4ThreeVector pos, G4double eDep , G4int trackID )
{
	hitEnergy->push_back(eDep);
	hitXPositions->push_back(pos.getX()/mm);
	hitYPositions->push_back(pos.getY()/mm);
	hitZPositions->push_back(pos.getZ()/mm);
}

void Analysis::AddEnergyOfChargedPart( G4double energyOfSec)

{
	energyOfSecChargedPart->push_back( energyOfSec/MeV );
}

void Analysis::PrepareNewRun(const G4Run* /*aRun*/ )
{
	//Reset variables relative to the run
	thisRunTotSecondaries = 0;
	eventNumber=0;
	//Create an empty histogram
	tree = new TTree("tree", "TPC Simulation");
	tree->Branch("eventNumber", &eventNumber, "eventNumber/I");
	tree->Branch("entryEnergy", &entryEnergy, "entryEnergy/D");
	tree->Branch("exitEnergy", &exitEnergy, "exitEnergy/D");
	tree->Branch("entryEnergyGas", &entryEnergyGas, "entryEnergyGas/D");
	tree->Branch("exitEnergyGas", &exitEnergyGas, "exitEnergyGas/D");

	tree->Branch("entryEnergyExtraFirstSi", &entryEnergyExtraFirstSi, "entryEnergyExtraFirstSi/D");
	tree->Branch("entryEnergyExtraThirdSi", &entryEnergyExtraThirdSi, "entryEnergyExtraThirdSi/D");
	tree->Branch("entryEnergyExtraNFirstSi", &entryEnergyExtraNFirstSi, "entryEnergyExtraNFirstSi/D");
	tree->Branch("entryEnergyExtraNThirdSi", &entryEnergyExtraNThirdSi, "entryEnergyExtraNThirdSi/D");
	tree->Branch("entryEnergyExtraNNFirstSi", &entryEnergyExtraNNFirstSi, "entryEnergyExtraNNFirstSi/D");
	tree->Branch("entryEnergyExtraNNThirdSi", &entryEnergyExtraNNThirdSi, "entryEnergyExtraNNThirdSi/D");

	tree->Branch("eDepTPCGas", &eDepTPCGas, "depositedEnergy/D");
	tree->Branch("exitEnergyWorld", &exitEnergyWorld, "exitEnergyWorld/D");
	tree->Branch("track Length", &trackLength, "trackLength/D");
	tree->Branch("entryEnergyMagnet", &entryEnergyMagnet, "entryEnergyMagnet/D");
	tree->Branch("exitEnergyMagnet", &exitEnergyMagnet, "exitEnergyMagnet/D");
	tree->Branch("noOfChargedParticleEnteringDriftVolume", &noOfChargedPartEnteringVolume, "noOfChargedPartEnteringVolume/I");
	tree->Branch("noOfSecondaryChargedParticle", &noOfChargedPart, "noOfChargedPart/I");
	tree->Branch("noOfSecondaryParticle", &thisEventSecondaries, "thisEventSecondaries/I");
	tree->Branch("energyOfSecondaryChargedPart", "std::vector<double>", &energyOfSecChargedPart);
	entryPositionMagnetT=new TVector3;
	exitPositionMagnetT=new TVector3;
	entryPositionT=new TVector3;
	exitPositionT=new TVector3;
	entryPositionGasT=new TVector3;
	exitPositionGasT=new TVector3;
	exitPositionWorldT=new TVector3;

	entryPositionExtraFirstSiT=new TVector3;
	entryPositionExtraThirdSiT=new TVector3;
	entryPositionExtraNFirstSiT=new TVector3;
	entryPositionExtraNThirdSiT=new TVector3;
	entryPositionExtraNNFirstSiT=new TVector3;
	entryPositionExtraNNThirdSiT=new TVector3;

	hitEnergy=new std::vector<G4double>;
	hitXPositions=new std::vector<G4double>;
	hitYPositions=new std::vector<G4double>;
	hitZPositions=new std::vector<G4double>;

	energyOfSecChargedPart=new std::vector<G4double>;
	tree->Branch("entryPosition", "TVector3", &entryPositionT );
	tree->Branch("entryPositionGas", "TVector3", &entryPositionGasT );
	tree->Branch("entryPositionMagnet", "TVector3", &entryPositionMagnetT );
	tree->Branch("exitPosition", "TVector3", &exitPositionT );
	tree->Branch("exitPositionGas", "TVector3", &exitPositionGasT );
	tree->Branch("exitPositionMagnet", "TVector3", &exitPositionMagnetT );
	tree->Branch("exitPositionWorld.", "TVector3", &exitPositionWorldT );

	tree->Branch("entryPositionExtraFirstSi", "TVector3", &entryPositionExtraFirstSiT );
	tree->Branch("entryPositionExtraThirdSi", "TVector3", &entryPositionExtraThirdSiT );
	tree->Branch("entryPositionExtraNFirstSi", "TVector3", &entryPositionExtraNFirstSiT );
	tree->Branch("entryPositionExtraNThirdSi", "TVector3", &entryPositionExtraNThirdSiT );
	tree->Branch("entryPositionExtraNNFirstSi", "TVector3", &entryPositionExtraNNFirstSiT );
	tree->Branch("entryPositionExtraNNThirdSi", "TVector3", &entryPositionExtraNNThirdSiT );

	tree->Branch("hitXPositions","std::vector<double>", &hitXPositions);
	tree->Branch("hitYPositions","std::vector<double>", &hitYPositions);
	tree->Branch("hitZPositions","std::vector<double>", &hitZPositions);

}

void Analysis::EndOfEvent(const G4Event* anEvent)
{
	//Accumulate over the run
	thisRunTotSecondaries += thisEventSecondaries;
	entryPositionT->SetXYZ(entryPosition.getX(),entryPosition.getY(),entryPosition.getZ());
	entryPositionGasT->SetXYZ(entryPositionGas.getX(),entryPositionGas.getY(),entryPositionGas.getZ());
	entryPositionMagnetT->SetXYZ(entryPositionMagnet.getX(),entryPositionMagnet.getY(),entryPositionMagnet.getZ());
	exitPositionT->SetXYZ(exitPosition.getX(),exitPosition.getY(),exitPosition.getZ());
	exitPositionGasT->SetXYZ(exitPositionGas.getX(),exitPositionGas.getY(),exitPositionGas.getZ());
	exitPositionMagnetT->SetXYZ(exitPositionMagnet.getX(),exitPositionMagnet.getY(),exitPositionMagnet.getZ());
	exitPositionWorldT->SetXYZ(exitPositionWorld.getX(),exitPositionWorld.getY(),exitPositionWorld.getZ());

	bool Layers2 = false;
	if ( Layers2 == true ) {
	  entryPositionExtraFirstSiT->SetXYZ(entryPositionExtraFirstSi.getX(),entryPositionExtraFirstSi.getY(),entryPositionExtraFirstSi.getZ());
	  entryPositionExtraThirdSiT->SetXYZ(entryPositionExtraThirdSi.getX(),entryPositionExtraThirdSi.getY(),entryPositionExtraThirdSi.getZ());
	}
	bool Nlayers2 = false;
	if ( Nlayers2 == true ) {
	  entryPositionExtraNFirstSiT->SetXYZ(entryPositionExtraNFirstSi.getX(),entryPositionExtraNFirstSi.getY(),entryPositionExtraNFirstSi.getZ());
	  entryPositionExtraNThirdSiT->SetXYZ(entryPositionExtraNThirdSi.getX(),entryPositionExtraNThirdSi.getY(),entryPositionExtraNThirdSi.getZ());
	}
	bool NNlayers2 = false;
	if ( NNlayers2 == true ) {
	  entryPositionExtraNNFirstSiT->SetXYZ(entryPositionExtraNNFirstSi.getX(),entryPositionExtraNNFirstSi.getY(),entryPositionExtraNNFirstSi.getZ());
	  entryPositionExtraNNThirdSiT->SetXYZ(entryPositionExtraNNThirdSi.getX(),entryPositionExtraNNThirdSi.getY(),entryPositionExtraNNThirdSi.getZ());
	}
	else {
	  entryPositionExtraNNFirstSiT->SetXYZ(-999., -999., -999.);
	  entryPositionExtraNNThirdSiT->SetXYZ(-999., -999., -999.);
	}

	tree->Fill();
	energyOfSecChargedPart->clear();
	hitXPositions->clear();
	hitYPositions->clear();
	hitZPositions->clear();
	hitEnergy->clear();

}

void Analysis::EndOfRun(const G4Run* aRun)
{
	//Some print outs

	G4int numEvents = aRun->GetNumberOfEvent();

	
	G4cout<<"================="<<G4endl;
	G4cout<<"Summary for run: "<<aRun->GetRunID()<<G4endl;
	G4cout<<"\t Event processed: "<<numEvents<<G4endl;
	G4cout<<"\t Average number of secondaries: "<<thisRunTotSecondaries/numEvents<<G4endl;
	//for ( int layer = 0 ; layer < NUMLAYERS ; ++layer)
	//{
	//	G4cout<<"\t\t Average energy in Layer "<<layer<<": "<<G4BestUnit(thisRunTotHad[layer],"Energy")<<G4endl;
	//}

	//At the end of the run we can now save a ROOT file containing the histogram
	// From the aRun variable we get the run number and use it to create a unique name for the output file
	char filename[256];
	sprintf(filename,"run_%d.root",aRun->GetRunID() );
	TFile* outfile = TFile::Open(filename,"recreate");
	tree->Write();
	outfile->Close();

	delete tree;
	tree = 0;


}
