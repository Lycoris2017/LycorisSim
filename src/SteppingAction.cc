#include "SteppingAction.hh"
//#include "G4Track.hh"
#include "G4Step.hh"
//#include "G4ParticleDefinition.hh"
//#include "G4ParticleTypes.hh"
//#include "G4StepPoint.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh"
//#include "G4TouchableHistory.hh"
#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"
#include "Analysis.hh"

SteppingAction::SteppingAction()
{
}

SteppingAction::~SteppingAction()
{
}

void SteppingAction::UserSteppingAction( const G4Step * theStep ) {
	//Check if this is the first step of this event,

	//G4cout<<fpSteppingManager->GetFirstStep()<<G4endl;
	//if ( fpSteppingManager->GetFirstStep() )
	//{
	//	G4cout<<"Previous Event had : "<<G4BestUnit(totalEdepEM,"Energy")<<G4endl;
	//	PrepareNewEvent();
	//}

	//Get the ;
	G4int preVolCopyNum = theStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetCopyNo();
	G4VPhysicalVolume* postPhysicalVolume = theStep->GetPostStepPoint()->GetPhysicalVolume();
	G4int postVolCopyNum;
	//Get the trackID of the Step
	G4int trackID = theStep->GetTrack()->GetTrackID();
	// catches a problem if it is the last step of the track
	if (postPhysicalVolume == 0)
	{
		if ( trackID ==1 )
			{
			Analysis::GetInstance()->InsertExitWorld( theStep->GetPreStepPoint()->GetPosition() , theStep->GetPreStepPoint()->GetTotalEnergy() );
			}
		postVolCopyNum=preVolCopyNum;

//		return;

	}
	else
	{
		postVolCopyNum = theStep->GetPostStepPoint()->GetTouchable()->GetVolume()->GetCopyNo();
	}
	//Get the trackID of the Step

	//Check if the step is the transition between air and second sensor
	if ( preVolCopyNum == 1 && postVolCopyNum == 8 )
	{
		//Check if the track is the primary particle
		if (trackID ==1)
		{
			Analysis::GetInstance()->InsertEntryTPC( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
		}
		//Check if the secondary particle is charged
		else if ( theStep->GetTrack()->GetParticleDefinition()->GetPDGCharge() > 0 )
		{
			Analysis::GetInstance()->AddChargedSecondaryEnetringVolume( 1 );
		}
	}

	//Check if the step the primary particle in the transition between third sensor and air
	else if ( preVolCopyNum == 9 && postVolCopyNum == 1 && trackID ==1)
	{
		Analysis::GetInstance()->InsertExitTPC( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}
	//Check if the step the primary particle in the transition air and magnet (incoming magnet)
	else if ( preVolCopyNum == 0 && postVolCopyNum == 3 && trackID ==1 )
	{
		Analysis::GetInstance()->InsertEntryMagnet( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}
	//Check if the step the primary particle in the tpc
	else if ( preVolCopyNum == 4 && trackID ==1)
	{
		Analysis::GetInstance()->AddTrackLength( theStep->GetStepLength() );
		Analysis::GetInstance()->AddEDepInGas( theStep->GetPreStepPoint()->GetKineticEnergy()-theStep->GetPostStepPoint()->GetKineticEnergy() );
		Analysis::GetInstance()->AddSimTrackHit(theStep->GetPreStepPoint()->GetPosition(), theStep->GetPreStepPoint()->GetKineticEnergy()-theStep->GetPostStepPoint()->GetKineticEnergy(), trackID);
		//G4cout << "EnergyDeposit=" <<theStep->GetDeltaEnergy()<<G4endl;
		//G4cout << "EnergyDeposit2=" <<theStep->GetPreStepPoint()->GetKineticEnergy()-theStep->GetPostStepPoint()->GetKineticEnergy()<<G4endl;
	}
	//Check if the step the primary particle in the transition between magnet and air (outgoing magnet)
	else if (preVolCopyNum == 3 && postVolCopyNum == 0 && trackID ==1)
	{
		Analysis::GetInstance()->InsertExitMagnet( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}
	//Check if the step the primary particle in the transition between first sensor and air
	else if ( preVolCopyNum == 7 && postVolCopyNum == 1 &&trackID ==1)
	{
		Analysis::GetInstance()->InsertEntryGas( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}
	//Check if the step the primary particle in the transition between air and fourth sensor
	if ( preVolCopyNum == 1 && postVolCopyNum == 10 &&trackID ==1)
	{
	  Analysis::GetInstance()->InsertExitGas( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}
	



	//Check if the step the primary particle in the transition between first EXTRA sensor and air
	else if ( preVolCopyNum == 1 && postVolCopyNum == 1000 && trackID ==1) {
	  Analysis::GetInstance()->InsertExtraFirstSi( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}
	//Check if the step the primary particle in the transition between third EXTRA sensor and air
	else if ( preVolCopyNum == 1 && postVolCopyNum == 1010 && trackID ==1) {
	  Analysis::GetInstance()->InsertExtraThirdSi( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}
	//N Check if the step the primary particle in the transition between first EXTRA sensor and air
	else if ( preVolCopyNum == 1 && postVolCopyNum == 1001 && trackID ==1) {
	  Analysis::GetInstance()->InsertExtraNFirstSi( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}
	//Check if the step the primary particle in the transition between third EXTRA sensor and air
	else if ( preVolCopyNum == 1 && postVolCopyNum == 1011 && trackID ==1) {
	  Analysis::GetInstance()->InsertExtraNThirdSi( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}

	//NN Check if the step the primary particle in the transition between first EXTRA sensor and air
	else if ( preVolCopyNum == 1 && postVolCopyNum == 1002 && trackID ==1) {
	  Analysis::GetInstance()->InsertExtraNNFirstSi( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}
	//Check if the step the primary particle in the transition between third EXTRA sensor and air
	else if ( preVolCopyNum == 1 && postVolCopyNum == 1012 && trackID ==1) {
	  Analysis::GetInstance()->InsertExtraNNThirdSi( theStep->GetPostStepPoint()->GetPosition() , theStep->GetPostStepPoint()->GetTotalEnergy() );
	}



}
