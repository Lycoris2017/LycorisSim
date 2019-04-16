// $Id: task2.cc 94 2010-01-26 13:18:30Z adotti $
/**
* @file
* @brief Main program.
*/

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "RunAction.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"
#include "PrimaryGeneratorAction.hh"

#include "PhysicsList.hh"
#include "QGSP_BERT.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"

#include "Randomize.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4PhysListFactory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {
    void PrintUsage() {
        G4cerr << " Usage: " << G4endl;
        G4cerr << " tpcsim [-m macro] [-u UIsession] [-h help]" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*!
\brief Main program

\callgraph

*/
int main(int argc,char** argv)
{
    // Evaluate arguments
    if ( argc > 7 ) {
        PrintUsage();
        return 1;
    }

    G4String macro;
    G4String session;

    for ( G4int i=1; i<argc; i=i+2 ) {
        if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
        else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
        else if ( G4String(argv[i]) == "-h" ) {
            PrintUsage();
            return 1;
        }
        else {
            PrintUsage();
            return 1;
        }
    }

    // Detect interactive mode (if no macro provided) and define UI session
    G4UIExecutive* ui = 0;

    if ( ! macro.size() ) {
        ui = new G4UIExecutive(argc, argv, session);
    }

    // Choose the Random engine
    G4Random::setTheEngine(new CLHEP::RanecuEngine);

    // Construct the default run manager
    G4RunManager * runManager = new G4RunManager;

    // Set mandatory initialization classes
    // Construction of the detector
    auto detConstruction = new DetectorConstruction();
    runManager->SetUserInitialization(detConstruction);

    // Choice of the physics List
    G4VUserPhysicsList* physics = new PhysicsList();
    runManager->SetUserInitialization(physics);

    // mandatory User Action classes
    G4VUserPrimaryGeneratorAction* gen_action = new PrimaryGeneratorAction();
    runManager->SetUserAction(gen_action);

    //Optional User Action classes
    //Stacking Action
    StackingAction* aStackingAction = new StackingAction();
    runManager->SetUserAction(aStackingAction);
    //Stepping Action
    SteppingAction* aSteppingAction = new SteppingAction();
    runManager->SetUserAction(aSteppingAction);
    //Event action (handles for beginning / end of event)
    EventAction* anEventAction = new EventAction();
    runManager->SetUserAction( anEventAction );
    //Run action (handles for beginning / end of event)
    RunAction* aRunAction = new RunAction();
    runManager->SetUserAction( aRunAction );

    // Initialize visualization
    auto visManager = new G4VisExecutive;
    visManager->Initialize();

    // Get the pointer to the User Interface manager
    auto UImanager = G4UImanager::GetUIpointer();

    // Process macro or start UI session
    if ( macro.size() ) {
        // batch mode
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command+macro);
    }
    else  {
        // interactive mode : define UI session
        UImanager->ApplyCommand("/control/execute init_vis.mac");
        if (ui->IsGUI()) {
            UImanager->ApplyCommand("/control/execute gui.mac");
        }
        ui->SessionStart();
        delete ui;
    }

    // Job termination
    // Free the store: user actions, physics_list and detector_description are
    // owned and deleted by the run manager, so they should not be deleted
    // in the main() program !
    delete visManager;
    delete runManager;
}
