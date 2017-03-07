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

#include "ExN01PrimaryGeneratorAction.hh"
#include "ExN01DetectorConstruction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
//#include "Randomize.hh"


ExN01PrimaryGeneratorAction::ExN01PrimaryGeneratorAction()
{
  ParticleGun = new G4GeneralParticleSource();
   

}

ExN01PrimaryGeneratorAction::~ExN01PrimaryGeneratorAction()
{
  delete ParticleGun;
}

void ExN01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 G4ParticleDefinition* particle = ParticleGun->GetParticleDefinition();
  
G4String n = particle->GetParticleName();

  ParticleGun->GeneratePrimaryVertex(anEvent);
}

















