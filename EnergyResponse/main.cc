
//// classes to include from the Geant4 framework
#include "G4RunManager.hh"
#include "Shielding.hh"
#include "G4ScoringManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//// user defined classes to include
#include "UserDefinedDetectorConstruction.hh"
#include "UserDefinedActionInitialization.hh"

//// stl libraries
// needed for getpid()
#include <unistd.h>

int main(int argc, char** argv){

  // the default randomness is not this great, so let's set somthing more sophisticated
  G4Random::setTheEngine(new CLHEP::RanecuEngine());
  G4Random::setTheSeed(time(NULL)*getpid()<<16);

  // create instance that sets up the simulation and guides the run
  G4RunManager* runManager = new G4RunManager();

  // world and detector setup for the simulation
  runManager->SetUserInitialization(new UserDefinedDetectorConstruction);

  // which physics processes should be used for the simulation
  runManager->SetUserInitialization(new Shielding);

  // actions during the simulation including the generation of particles
  runManager->SetUserInitialization(new UserDefinedActionInitialization);

  // set score writer to scoring manager
  G4ScoringManager::GetScoringManager();

  // run a macro file
  G4String command = "/control/execute ";
  G4String fileName;

  //Henning:
  std::ofstream ofs;
  ofs.open("SomeFile.txt", std::ofstream::out | std::ofstream::trunc);
  ofs << "Energy of Particle in MeV" << G4endl;

  // if no macro (argv[0] is the executable itself), run in GUI mode
  if (argc == 1)
  {
    // create an instance for the visualization and init it
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    // get interactive session
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);

    // initialize the view itself via macro
    fileName = "macros/init_vis.mac";
    G4UImanager::GetUIpointer()->ApplyCommand(command+fileName);

    ui->SessionStart();
    delete ui;
    delete visManager;
  }
  else // execute macro and close
  {
    // read in macro file name and run it
    fileName = argv[1];
    G4UImanager::GetUIpointer()->ApplyCommand(command+fileName);
  }

  delete runManager;
  return 0;
}
