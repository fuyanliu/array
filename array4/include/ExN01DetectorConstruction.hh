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
// $Id$
//

#ifndef ExN01DetectorConstruction_H
#define ExN01DetectorConstruction_H 1
#include "globals.hh"
#include "ExN01MagneticField.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4FieldManager.hh"
class G4Box;
class G4Orb;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VPVParameterisation;
class G4UserLimits;
class ExN01DetectorMessenger;
class ExN01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    ExN01DetectorConstruction();
    virtual ~ExN01DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

  public:
    
    // Logical volumes
    //
    G4LogicalVolume* experimentalHall_log;
    G4LogicalVolume* pilotA_log;
    G4LogicalVolume* pilotB_log;
    G4LogicalVolume* pilotC_log;
    G4LogicalVolume* pilotD_log;
    G4LogicalVolume* kapton_log;
    G4LogicalVolume* sample_log;
    G4LogicalVolume* air_log;
    G4LogicalVolume* detector_log;
    // Physical volumes
    //
    G4VPhysicalVolume* experimentalHall_phys;
    G4VPhysicalVolume* pilotA_phys;
    G4VPhysicalVolume* pilotB_phys;
    G4VPhysicalVolume* pilotC_phys;
    G4VPhysicalVolume* pilotD_phys;
    G4VPhysicalVolume* kapton_phys;
    G4VPhysicalVolume* sample_phys;
    G4VPhysicalVolume* air_phys;
    G4VPhysicalVolume* detector_phys;
   // G4VPhysicalVolume* tracker_phys;
void SetMaxStep (G4double);     
void SetMagField (G4double);
     
   private:
G4UserLimits* stepLimit; 
ExN01MagneticField* fpMagField;
ExN01DetectorMessenger* detectorMessenger;


};

#endif

