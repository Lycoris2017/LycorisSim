// $Id: RunAction.cc 102 2010-01-26 16:45:56Z adotti $
/**
 * @file   RunAction.cc
 *
 * @date   17 Dec 2009
 * @author adotti
 *
 * @brief  Implements user class RunAction.
 */

#include "RunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "Analysis.hh"

#include <ctime>

RunAction::RunAction()
{
}

void RunAction::BeginOfRunAction(const G4Run* aRun )
{

	// start -- added Mengqing
	//to be sure to generate different events!
	long seeds[2];
	long systime = time(NULL);
	seeds[0] = systime;
	seeds[1] = systime * G4UniformRand();
	G4Random::setTheSeeds(seeds);
	/* long mylong = 1024;
	 const long* myseeds=&mylong;
	 G4Random::setTheSeeds(myseeds);*/
    // end -- added Mengqing
	
	G4cout<<"Starting Run: "<<aRun->GetRunID()<<G4endl;
	Analysis::GetInstance()->PrepareNewRun(aRun);

	
}

void RunAction::EndOfRunAction( const G4Run* aRun )
{
	Analysis::GetInstance()->EndOfRun(aRun);
}
