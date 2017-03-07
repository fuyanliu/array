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

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "iostream"

#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
ExN01PrimaryGeneratorAction::ExN01PrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0),
   fpilotABox(0)
{

  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

 
//default kinematic
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();  
G4ParticleDefinition* particle
                    = particleTable->FindParticle("gamma");
  fParticleGun->SetParticleDefinition(particle);
  //fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(511*keV);
  /*G4double Radius = 1*mm;
  G4double x0 = Radius*2;
  G4double y0 = Radius*2;
  G4double z0 = Radius*2;
  
  while(((x0*x0)+(y0*y0)+(z0*z0)) > (Radius*Radius)) {
      x0 = G4UniformRand();
      y0 = G4UniformRand();
      z0 = G4UniformRand();
      
      x0 = (x0*2.*Radius) - Radius;
      y0 = (y0*2.*Radius) - Radius;
      z0 = (z0*2.*Radius) - Radius;
    }*/

                                                             
}



ExN01PrimaryGeneratorAction::~ExN01PrimaryGeneratorAction()
{
  delete fParticleGun;
}

void ExN01PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  /*G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
  if (particle == G4ChargedGeantino::ChargedGeantino()) {
    //fluorine 
    G4int Z = 11, A = 22;
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;
    
    G4ParticleDefinition* ion
       = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
   
   
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
  }*/
  G4double plASizeXY = 3*mm;
  G4double plASizeZ = 70*mm;


  /*if (!fpilotABox)
  {
    G4LogicalVolume* plALV
      = G4LogicalVolumeStore::GetInstance()->GetVolume("PilotA");
    if ( plALV ) fpilotABox = dynamic_cast<G4Box*>(plALV->GetSolid());
  }

  if ( fpilotABox ) {
    plASizeXY = fpilotABox->GetXHalfLength()*2.;
    plASizeZ = fpilotABox->GetZHalfLength()*2.;

  }  
  else  {
    G4ExceptionDescription msg;
    msg << "Envelope volume of box shape not found.\n"; 
    msg << "Perhaps you have changed geometry.\n";
    msg << "The gun will be place at the center.";
    G4Exception("ExN01PrimaryGeneratorAction::GeneratePrimaries()",
     "MyCode0002",JustWarning,msg);
  }*/
  //G4cout<<"fpilotABox is "<<fpilotABox<<"  plASizeXY is "<<plASizeXY<<"  plASizeZ is "<<plASizeZ<<G4endl;
  G4double x0 = plASizeXY * (G4UniformRand()-0.5);
  G4double y0 = plASizeXY * (G4UniformRand()-0.5);
  G4double z0 = 200 ;
  /*G4double Radius = 1*um;
  G4double x0 = Radius*2;
  G4double y0 = Radius*2;
  G4double z0 = Radius*2;
  
  while(((x0*x0)+(y0*y0)+(z0*z0)) > (Radius*Radius)) {
      x0 = G4UniformRand();
      y0 = G4UniformRand();
      z0 = G4UniformRand();
      
      x0 = (x0*2.*Radius) - Radius;
      y0 = (y0*2.*Radius) - Radius;
      z0 = (z0*2.*Radius) - Radius;
    }*/
  //G4cout<<"position is x:"<<x0/mm<<"  y:"<<y0/mm<<"  z:"<<z0/mm<<G4endl;
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(x0,y0,z0));
  fParticleGun->SetParticlePosition(G4ThreeVector(0,0,0.*cm));
  fParticleGun->GeneratePrimaryVertex(anEvent);
}



