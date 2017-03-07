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
// $Id: ExN01SteppingAction.cc 69899 1013-05-17 10:05:33Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN01SteppingAction.hh"
#include "ExN01EventAction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"
#include "G4SystemOfUnits.hh"
#include "iostream"
#include "G4EventManager.hh"
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01SteppingAction::ExN01SteppingAction(ExN01EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01SteppingAction::~ExN01SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01SteppingAction::UserSteppingAction(const G4Step* step)
{ 
  /*G4double reset = 142.0*ns; 
  G4double sameTime = reset; 
  G4EventManager::GetEventManager()->GetUserStackingAction()->SetSameTime(sameTime); */

  G4double parentID = step->GetTrack()->GetParentID();
  G4Track* aTrack = step->GetTrack();
 
  G4StepPoint * thePrePoint = step->GetPreStepPoint();
  G4VPhysicalVolume * thePrePV = thePrePoint->GetPhysicalVolume();
  G4String thePrePVname = thePrePV->GetName();
  G4double time_pre = thePrePoint->GetLocalTime(); 
  G4double time_post = step->GetPostStepPoint()->GetLocalTime();
  G4double time_glo = thePrePoint->GetGlobalTime();
  G4double time_track = step->GetTrack()->GetLocalTime();
  G4double ProperTime=step->GetTrack()->GetProperTime();
  G4String ProcessName = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  G4String ParticleName = step->GetTrack()->GetDynamicParticle()->GetDefinition()->GetParticleName();
  G4String volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume()->GetName();
  G4double kinEnergy = step->GetTrack()->GetKineticEnergy();
  G4double id = step->GetTrack()->GetTrackID();
  G4double x = step->GetTrack()->GetPosition().x();
  G4double y = step->GetTrack()->GetPosition().y();
  G4double z = step->GetTrack()->GetPosition().z();
  G4double edepStep = step->GetTotalEnergyDeposit();
  G4double annihillength = sqrt(x*x + y*y +z*z);
                //std::ofstream afile("volume.dat",std::ios::app);
                //afile <<" "<<thePrePVname<<" "<<ProcessName<<" "<<ParticleName<<" "<< kinEnergy/keV <<" "<< id << std::endl;
  
 /* if (thePrePVname == "sample" && ProcessName == "annihil" && ParticleName == "e+")
  {
                std::ofstream afile("anti_pre.dat",std::ios::app);
                afile << time_pre/ns<< std::endl;
                std::ofstream cfile("propertime.dat",std::ios::app);
                cfile << ProperTime<< std::endl;
  }
  if (thePrePVname == "sample" && ProcessName == "annihil" && ParticleName == "e+")
  {
                std::ofstream bfile("globaltime.dat",std::ios::app);
                bfile << time_glo<< std::endl;
  }*/
  /*if (thePrePVname == "air1" && ProcessName == "annihil" && ParticleName == "e+")
  {
                
                std::ofstream cfile("air.dat",std::ios::app);
                cfile << annihillength/mm<< std::endl;
  }
  if (thePrePVname == "air2" && ProcessName == "annihil" && ParticleName == "e+")
  {
                
                std::ofstream c1file("air.dat",std::ios::app);
                c1file << annihillength/mm<< std::endl;
  }
  
  if (thePrePVname == "kapton2" && ProcessName == "annihil" && ParticleName == "e+")
  {
                std::ofstream efile("kapton.dat",std::ios::app);
                efile << annihillength/mm<< std::endl;
  }
  if (thePrePVname == "base" && ProcessName == "annihil" && ParticleName == "e+")
  {
                std::ofstream ffile("base.dat",std::ios::app);
                ffile << annihillength/mm<< std::endl;
  }
  */
  
  /*if (thePrePVname == "detector" && ParticleName == "gamma")
  {
                //std::ofstream afile("volume.dat",std::ios::app);
                //afile <<" "<<thePrePVname<<" "<<ProcessName<<" "<<ParticleName<< std::endl;
                std::ofstream bfile("detector.dat",std::ios::app);
                bfile << kinEnergy/MeV << std::endl;
  }*/
  /*
  if (ProcessName == "phot" && ParticleName == "gamma")
  {
    std::ofstream bfile("phot.dat",std::ios::app);
                bfile << "a" << std::endl;
  }*/
  /*if (ParticleName == "gamma")
  {
                std::ofstream afile("volume.dat",std::ios::app);
                afile <<" "<<thePrePVname<<" "<<ProcessName<<" "<<ParticleName<<" "<< kinEnergy/keV << std::endl;
  }*/

  if (thePrePVname == "PilotA" && ParticleName == "e-" )
  {
                //fEventAction->AddEdep(edepStep);
                std::ofstream dafile("ao.dat",std::ios::app);
                dafile <<edepStep/keV << std::endl;
  }
  if (thePrePVname == "PilotB" && ParticleName == "e-" )
  {
                //fEventAction->AddEdep(edepStep);
                std::ofstream dbfile("bo.dat",std::ios::app);
                dbfile <<edepStep/keV << std::endl;
  }
  if (thePrePVname == "PilotC" && ParticleName == "e-" )
  {
                //fEventAction->AddEdep(edepStep);
                std::ofstream dcfile("co.dat",std::ios::app);
                dcfile <<edepStep/keV << std::endl;
  }
  if (thePrePVname == "PilotD" && ParticleName == "e-" )
  {
                //fEventAction->AddEdep(edepStep);
                std::ofstream ddfile("do.dat",std::ios::app);
                ddfile <<edepStep/keV << std::endl;
  }
  /*if (thePrePVname == "PilotB" && ParticleName == "e-" )
  {
                fEventAction->AddEdep(edepStep);
                std::ofstream efile("b.dat",std::ios::app);
                efile <<edepStep/keV << std::endl;
  }*/
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

