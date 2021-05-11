
#ifndef UserDefinedDetectorConstruction_h
#define UserDefinedDetectorConstruction_h 1

// classes to include from the Geant4 framework
#include "G4VUserDetectorConstruction.hh"

// use classes without importing it yet - forward definition
class G4VPhysicalVolume;

class UserDefinedDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    UserDefinedDetectorConstruction();
    virtual ~UserDefinedDetectorConstruction();

    // this is called, when run is initialized
    // it should return the "simulation world"
    virtual G4VPhysicalVolume* Construct();
};


#endif
