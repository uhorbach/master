
// local headers to include
#include "UserDefinedDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

// a constructor that does nothing except call its parent's constructor
UserDefinedDetectorConstruction::UserDefinedDetectorConstruction()
: G4VUserDetectorConstruction()
{}

// empty destructor
UserDefinedDetectorConstruction::~UserDefinedDetectorConstruction()
{}


  G4VPhysicalVolume* UserDefinedDetectorConstruction::Construct()
  {
    // Get nist material manager
    G4NistManager* nist = G4NistManager::Instance();

    // Envelope parameters
    //
    G4double env_sizeXY = 3*m, env_sizeZ = 30*cm;

    // Option to switch on/off checking of volumes overlaps
    //
    G4bool checkOverlaps = true;

    //
    // World
    //
    G4double world_sizeXY = 1.2*env_sizeXY;
    G4double world_sizeZ  = 1.2*env_sizeZ;
    G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

    G4Box* solidWorld =
      new G4Box("World",                       //its name
         0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

    G4LogicalVolume* logicWorld =
      new G4LogicalVolume(solidWorld,          //its solid
                          world_mat,           //its material
                          "World");            //its name

    G4VPhysicalVolume* physWorld =
      new G4PVPlacement(0,                     //no rotation
                        G4ThreeVector(),       //at (0,0,0)
                        logicWorld,            //its logical volume
                        "World",               //its name
                        0,                     //its mother  volume
                        false,                 //no boolean operation
                        0,                     //copy number
                        checkOverlaps);        //overlaps checking

    //
    // Envelope
    //
    G4Box* solidEnv =
      new G4Box("Envelope",                    //its name
          0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

    G4LogicalVolume* logicEnv =
      new G4LogicalVolume(solidEnv,            //its solid
                          world_mat,             //its material
                          "Envelope");         //its name

    new G4PVPlacement(0,                       //no rotation
                      G4ThreeVector(),         //at (0,0,0)
                      logicEnv,                //its logical volume
                      "Envelope",              //its name
                      logicWorld,              //its mother  volume
                      false,                   //no boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking










  return physWorld;
}
