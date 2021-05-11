//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file UserDefinedDetectorConstruction.cc
/// \brief Implementation of the UserDefinedDetectorConstruction class

#include "UserDefinedDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh" 
#include "G4Region.hh"
#include "G4UserLimits.hh"
#include "G4ProductionCuts.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UserDefinedDetectorConstruction::UserDefinedDetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

UserDefinedDetectorConstruction::~UserDefinedDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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



    //Detektoren
    // Vordefinitionen Material!
        // Aluminium-Legierung
        	G4double Z, a, density;
          G4String name, symbol;
          G4int ncomponents, natoms;
          G4double fractionmass;
          //Silizium
          a=28.09*g/mole;
          G4Element* elSi = new G4Element(name="Silizium", symbol="Si", Z=14.,a);
          //Eisen
          a=55.85*g/mole;
          G4Element* elFe = new G4Element(name="Eisen", symbol="Fe", Z=26.,a);
          //Kupfer
          a=63.55*g/mole;
          G4Element* elCu = new G4Element(name="Kupfer", symbol="Cu", Z=29.,a);
          //Mangan
          a=54.94*g/mole;
          G4Element* elMn = new G4Element(name="Mangan", symbol="Mn", Z=25.,a);
          //Magnesium
          a=24.31*g/mole;
          G4Element* elMg = new G4Element(name="Magnesium", symbol="Mg", Z=12.,a);
          //Chrom
          a=51.99*g/mole;
          G4Element* elCr = new G4Element(name="Chrom", symbol="Cr", Z=24.,a);
          //Zink
          a=65.38*g/mole;
          G4Element* elZn = new G4Element(name="Zink", symbol="Zn", Z=30.,a);
          //Titan
          a=47.88*g/mole;
          G4Element* elTi = new G4Element(name="Titan", symbol="Ti", Z=22.,a);
          //Aluminium
          a=26.98*g/mole;
          G4Element* elAl = new G4Element(name="Aluminium", symbol="Al", Z=13.,a);
          density = 2.7*g/cm3;
          G4Material* AlSi1MgMn = new G4Material(name="AlSi1MgMn", density, ncomponents=9);
          AlSi1MgMn->AddElement(elSi, fractionmass=1.0*perCent);
          AlSi1MgMn->AddElement(elFe, fractionmass=0.5*perCent);
          AlSi1MgMn->AddElement(elCu, fractionmass=0.1*perCent);
          AlSi1MgMn->AddElement(elMn, fractionmass=0.7*perCent);
          AlSi1MgMn->AddElement(elMg, fractionmass=0.9*perCent);
          AlSi1MgMn->AddElement(elCr, fractionmass=0.25*perCent);
          AlSi1MgMn->AddElement(elZn, fractionmass=0.2*perCent);
          AlSi1MgMn->AddElement(elTi, fractionmass=0.1*perCent);
          AlSi1MgMn->AddElement(elAl, fractionmass=96.25*perCent);


          //Keramikschicht-> Gemisch aus Al2O3 und TiO2
          //Sauerstoff
          a=16.*g/mole;
          G4Element* elO = new G4Element(name="Sauerstoff", symbol="O", Z=16.,a);
          //Aluminiumoxid
          density=3.94*g/cm3;
          G4Material* Al2O3 = new G4Material(name="Al2O3", density, ncomponents=2);
          Al2O3->AddElement(elAl, natoms=2);
          Al2O3->AddElement(elO, natoms=3);
          //Titanoxid
          density=4.24*g/cm3;
          G4Material* TiO2 = new G4Material(name="TiO2", density, ncomponents=2);
          TiO2->AddElement(elTi, natoms=1);
          TiO2->AddElement(elO, natoms=2);
          //Keramik
          density=4*g/cm3;
          G4Material* Keramik = new G4Material(name="Keramik", density, ncomponents=2);
          Keramik->AddMaterial(Al2O3, fractionmass=87*perCent);
          Keramik->AddMaterial(TiO2, fractionmass=13*perCent);



          //LiF
          G4Material* Lithiumfluorid = nist->FindOrBuildMaterial("G4_LITHIUM_FLUORIDE");

    //CODERING

    //Position
    G4double Px = 0*mm;
    G4double Py = 0*mm;
    G4double Pz = 0*mm;
    //face direction
    G4double Dx = 1;
    G4double Dy = 0;
    G4double Dz = 0;

    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateX(acos(Dz)+asin(Dx));  //um die x-Achse drehen
    rotm.rotateY(acos(Dy)+asin(Dz));   //um die y Achse drehen
    rotm.rotateZ(acos(Dx)+asin(Dy));  //um die z-Achse drehen, macht bei Zylinder keinen Sinn, da symmetrisch
    //Aluring oben

    	G4ThreeVector posa = G4ThreeVector(Px, Py, Pz);
    	G4Transform3D transforma = G4Transform3D(rotm,posa);

      G4double rmina = 5*mm;
      G4double rmaxa = 8.5*mm;
      G4double dza = 0.5*mm;
      G4double phimina = 0.*deg;
      G4double phimaxa = 360.*deg;
      G4Tubs* ringo =
       new G4Tubs("ringo", rmina, rmaxa, dza, phimina, phimaxa);

     G4LogicalVolume* logicalo = new G4LogicalVolume(ringo, AlSi1MgMn, "ringo");

      new G4PVPlacement(transforma,		//rot&pos
                        logicalo,             //its logical volume
                        "ringo",                //its name
                        logicEnv,             //its mother  volume
                        false,         	     //boolean operation
                        0,                       //copy number
                        checkOverlaps);          //overlaps checking

    //   Aluring unten

    G4ThreeVector posb = G4ThreeVector(Px+Dx*(-0.75)*mm,Py+Dy*(-0.75)*mm,Pz+Dz*(-0.75)*mm);
    G4Transform3D transformb = G4Transform3D(rotm,posb);

      G4double rminb = 3.25*mm;
      G4double rmaxb = 8.5*mm;
      G4double dzb = 0.25*mm;
      G4double phiminb = 0.*deg;
      G4double phimaxb = 360.*deg;
      G4Tubs* ringu =
       new G4Tubs("ringu", rminb, rmaxb, dzb, phiminb, phimaxb);

     G4LogicalVolume* logicalu = new G4LogicalVolume(ringu, AlSi1MgMn, "ringo");

      new G4PVPlacement(transformb,		//rot&pos
                        logicalu,             //its logical volume
                        "ringu",                //its name
                        logicEnv,             //its mother  volume
                        false,         	     //boolean operation
                        0,                       //copy number
                        checkOverlaps);          //overlaps checking




    // AluminiumTräger

    //ring
    G4double rminT1 = 3.04*mm;
    G4double rmaxT1 = 5*mm;
    G4double dzT1 = 0.13*mm;
    G4double phiminT1 = 0.*deg;
    G4double phimaxT1 = 360.*deg;
    G4Tubs* T1 =
     new G4Tubs("AluTräger1", rminT1, rmaxT1, dzT1, phiminT1, phimaxT1);

    G4LogicalVolume* logicT1 = new G4LogicalVolume(T1, AlSi1MgMn, "AluTräger1");

    new G4PVPlacement(transforma,		//rot&pos
                      logicT1,             //its logical volume
                      "AluTräger1",                //its name
                      logicEnv,             //its mother  volume
                      false,         	     //boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

    //platte

    G4ThreeVector posc = G4ThreeVector(Px+Dx*(-0.315)*mm,Py+Dy*(-0.315)*mm,Pz+Dz*(-0.315)*mm);
    G4Transform3D transformc = G4Transform3D(rotm,posc);

    G4double rminT2 = 0*mm;
    G4double rmaxT2 = 5*mm;
    G4double dzT2 = 0.185*mm;
    G4double phiminT2 = 0.*deg;
    G4double phimaxT2 = 360.*deg;
    G4Tubs* T2 =
     new G4Tubs("AluTräger2", rminT2, rmaxT2, dzT2, phiminT2, phimaxT2);

    G4LogicalVolume* logicT2 = new G4LogicalVolume(T2, AlSi1MgMn, "AluTräger2");

    new G4PVPlacement(transformc,		//rot&pos
                      logicT2,             //its logical volume
                      "AluTräger2",                //its name
                      logicEnv,             //its mother  volume
                      false,         	     //boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking




    // Keramikschicht (Flammspritzen)

    //ring
    G4double rminfs1 = 3.01*mm;
    G4double rmaxfs1 = 3.04*mm;
    G4double dzfs1 = 0.1*mm;
    G4double phiminfs1 = 0.*deg;
    G4double phimaxfs1 = 360.*deg;
    G4Tubs* fs1 =
     new G4Tubs("Flammspritzen1", rminfs1, rmaxfs1, dzfs1, phiminfs1, phimaxfs1);

    G4LogicalVolume* logicfs1 = new G4LogicalVolume(fs1, Keramik, "Flammspritzen1");

    new G4PVPlacement(transforma,		//rot&pos
                      logicfs1,             //its logical volume
                      "Flammspritzen1",                //its name
                      logicEnv,             //its mother  volume
                      false,         	     //boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking

    //platte

    G4ThreeVector posd = G4ThreeVector(Px+Dx*(-0.115)*mm, Py+Dy*(-0.115)*mm, Pz+Dz*(-0.115)*mm);
    G4Transform3D transformd = G4Transform3D(rotm,posd);

    G4double rminfs2 = 0*mm;
    G4double rmaxfs2 = 3.04*mm;
    G4double dzfs2 = 0.015*mm;
    G4double phiminfs2 = 0.*deg;
    G4double phimaxfs2 = 360.*deg;
    G4Tubs* fs2 =
     new G4Tubs("Flammspritzen2", rminfs2, rmaxfs2, dzfs2, phiminfs2, phimaxfs2);

    G4LogicalVolume* logicfs2 = new G4LogicalVolume(fs2, Keramik, "Flammspritzen2");

    new G4PVPlacement(transformd,		//rot&pos
                      logicfs2,             //its logical volume
                      "Flammspritzen2",                //its name
                      logicEnv,             //its mother  volume
                      false,         	     //boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking



    //LiF

    G4double rminLiF = 0*mm;
    G4double rmaxLiF = 3.01*mm;
    G4double dzLiF = 0.1*mm;
    G4double phiminLiF = 0.*deg;
    G4double phimaxLiF = 360.*deg;
    G4Tubs* Lif =
     new G4Tubs("LiF", rminLiF, rmaxLiF, dzLiF, phiminLiF, phimaxLiF);

    G4LogicalVolume* logicLif = new G4LogicalVolume(Lif, Lithiumfluorid, "LiF");

    new G4PVPlacement(transforma,		//rot&pos
                      logicLif,             //its logical volume
                      "LiF",                //its name
                      logicEnv,             //its mother  volume
                      false,         	     //boolean operation
                      0,                       //copy number
                      checkOverlaps);          //overlaps checking


  G4Region* Lithium = new G4Region("Li-Detector");
  Lithium->AddRootLogicalVolume(logicLif);
  logicLif->SetRegion(Lithium);
  G4ProductionCuts* myCuts = new G4ProductionCuts();
  myCuts->SetProductionCut(40*micrometer, "e-");
  myCuts->SetProductionCut(40*micrometer, "gamma");
  Lithium->SetProductionCuts(myCuts); 
  logicLif->SetUserLimits(new G4UserLimits(40*micrometer));






 //scoringvolume definieren

      fScoringVolume = logicLif;
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
