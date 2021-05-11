
#ifndef UserDefinedPrimaryGeneratorAction_h
#define UserDefinedPrimaryGeneratorAction_h 1

// classes to include from the Geant4 framework
#include "G4VUserPrimaryGeneratorAction.hh"

// use classes without importing it yet - forward definition
class G4Event;
class G4GeneralParticleSource;

class UserDefinedPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    UserDefinedPrimaryGeneratorAction();
    virtual ~UserDefinedPrimaryGeneratorAction();

    // this method is called at the beginning of each event
    virtual void GeneratePrimaries(G4Event*);

  private:
    // our particle source
    G4GeneralParticleSource* gps;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
