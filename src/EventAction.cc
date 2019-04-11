// $Id: EventAction.cc 94 2010-01-26 13:18:30Z adotti $
/**
 * @file   EventAction.cc
 *
 * @date   17 Dec 2009
 * @author adotti
 * 
 * @brief  Implements user class EventAction.
 */

#include "EventAction.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4DigiManager.hh"
#include "G4Event.hh"
#include "G4UnitsTable.hh"
#include "Analysis.hh"

EventAction::EventAction()
{
	//We build the digitization module
	//SiDigitizer* digitizer = new SiDigitizer("SiDigitizer");
	//G4DigiManager * digiManager = G4DigiManager::GetDMpointer();
	//digiManager->AddNewModule( digitizer );
}


void EventAction::BeginOfEventAction(const G4Event* anEvent )
{
  if ( anEvent->GetEventID() % 5000 == 0 )
    {
      //		G4cout<<"Starting Event: "<<anEvent->GetEventID()<<G4endl;
      std::cout<<"Starting Event: "<<anEvent->GetEventID()<<std::endl;
    }
  Analysis::GetInstance()->PrepareNewEvent(anEvent);
	//Retrieve the ID for the hit collection
	//if ( hitsCollID == -1 )
	//{
	//	G4SDManager * SDman = G4SDManager::GetSDMpointer();
	//	hitsCollID = SDman->GetCollectionID(hitsCollName);
	//}
}

void EventAction::EndOfEventAction(const G4Event* anEvent)
{
	Analysis::GetInstance()->EndOfEvent(anEvent);
}

