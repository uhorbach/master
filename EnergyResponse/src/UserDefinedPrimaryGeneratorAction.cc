
// classes to include from the Geant4 framework
#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
// user defined classes to include
#include "UserDefinedPrimaryGeneratorAction.hh"

// call parent constructor and empty initialize the particle source
// using the initializer list
UserDefinedPrimaryGeneratorAction::UserDefinedPrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    gps(new G4GeneralParticleSource())
{}

// free the allocated memory
UserDefinedPrimaryGeneratorAction::~UserDefinedPrimaryGeneratorAction()
{
  delete gps;
}

// this is what happens when "BeamOn" is called
// simply tell the general particle source to generate the GeneratePrimaryVertex
void UserDefinedPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  gps->GeneratePrimaryVertex(anEvent);
  G4double ke = anEvent->GetPrimaryVertex()->GetPrimary()->GetKineticEnergy();
  std::ofstream ofs;
  ofs.open("SomeFile.txt", std::ofstream::app); // append instead of overwrite
  ofs << ke << G4endl;
  ofs.close();
}
