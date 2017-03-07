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
// $Id: ExN01EventAction.cc 69899 1013-05-17 10:05:33Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#include "ExN01EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "ExN01CrystalSD.hh"
#include "ExN01CrystalHit.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN01EventAction::ExN01EventAction()
: G4UserEventAction(),
  fEdepAHCID(-1),
  fEdepBHCID(-1),
  fEdepCHCID(-1),
  fEdepDHCID(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN01EventAction::~ExN01EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01CrystalHitsCollection*
ExN01EventAction::GetHitsCollection(G4int hcID,
                                    const G4Event* event) const
{
  ExN01CrystalHitsCollection* hitsCollection
    = static_cast<ExN01CrystalHitsCollection*>(
        event->GetHCofThisEvent()->GetHC(hcID));

  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID "<<hcID;
    G4Exception("ExN01EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }

  return hitsCollection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01EventAction::PrintEventStatistics(
                            G4double aEdep, G4double bEdep,
                            G4double cEdep, G4double dEdep) const
{
  // Print event statistics
  //
  G4cout
     << "   PilotA: total energy: " 
     << std::setw(7) << G4BestUnit(aEdep, "Energy")
     << "       PilotB total energy: " 
     << std::setw(7) << G4BestUnit(bEdep, "Energy")
     << G4endl
     << "   PilotC: total energy: " 
     << std::setw(7) << G4BestUnit(cEdep, "Energy")
     << "       PilotD total energy: " 
     << std::setw(7) << G4BestUnit(dEdep, "Energy")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExN01EventAction::BeginOfEventAction(const G4Event* event)
{
  /*fEdepA = 0.;
  fEdepB = 0.;
  fEdepC = 0.;
  fEdepD = 0.;*/
  /*G4int event_id = event->GetEventID();
  std::ofstream dfile("step.dat",std::ios::app);
               dfile << "****** "<<event_id<<"*********" << std::endl;*/
  // std::ofstream efile("bo.dat",std::ios::app);
  //               efile << event_id << std::endl;
  // std::ofstream ffile("co.dat",std::ios::app);
  //               ffile << event_id << std::endl;
  // std::ofstream gfile("do.dat",std::ios::app);
  //               gfile << event_id << std::endl;
  //std::ofstream afile("volume.dat",std::ios::app);
  //              afile <<event_id << std::endl;
  //std::ofstream bfile("detector.dat",std::ios::app);
  //              bfile << event_id <<std::endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExN01EventAction::EndOfEventAction(const G4Event* event)
{
  //Get hits collections IDs
  if ( fEdepAHCID == -1 ) {
    fEdepAHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("aPilotHitsCollection");
    fEdepBHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("bPilotHitsCollection");
    fEdepCHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("cPilotHitsCollection");
    fEdepDHCID
      = G4SDManager::GetSDMpointer()->GetCollectionID("dPilotHitsCollection");
  }

  ExN01CrystalHitsCollection* aHC = GetHitsCollection(fEdepAHCID, event);
  ExN01CrystalHitsCollection* bHC = GetHitsCollection(fEdepBHCID, event);
  ExN01CrystalHitsCollection* cHC = GetHitsCollection(fEdepCHCID, event);
  ExN01CrystalHitsCollection* dHC = GetHitsCollection(fEdepDHCID, event);

  ExN01CrystalHit* aHit = (*aHC)[aHC->entries()-1];
   //G4cout << "********aHC " << aHC->entries() << "*********" << std::endl;
  ExN01CrystalHit* bHit = (*bHC)[bHC->entries()-1];
  ExN01CrystalHit* cHit = (*cHC)[cHC->entries()-1];
  ExN01CrystalHit* dHit = (*dHC)[dHC->entries()-1];

  //Get sum values from hits collections
  G4double aEdep = aHit->GetEdep();
  G4double bEdep = bHit->GetEdep();
  G4double cEdep = cHit->GetEdep();
  G4double dEdep = dHit->GetEdep();

  G4int event_id = event->GetEventID();
  //G4cout<< fEdep/CLHEP::keV<<" ***** "<<event_id<<G4endl;
  std::ofstream dfile("ao.dat",std::ios::app);
              dfile << aEdep/CLHEP::keV << std::endl;
              //dfile << event_id << std::endl;
  std::ofstream efile("bo.dat",std::ios::app);
                //efile << bEdep/CLHEP::keV << std::endl;
                efile << bEdep/CLHEP::keV << std::endl;
  std::ofstream ffile("co.dat",std::ios::app);
                ffile << cEdep/CLHEP::keV << std::endl;
                //ffile << event_id << std::endl;
  std::ofstream gfile("do.dat",std::ios::app);
                gfile << dEdep/CLHEP::keV << std::endl;
                //gfile << event_id << std::endl;
  //std::ofstream afile("volume.dat",std::ios::app);
  //              afile <<event_id << std::endl;
  //std::ofstream bfile("detector.dat",std::ios::app);
  //              bfile <<event_id <<std::endl;  
  // get number of stored trajectories
  //
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  // periodic printing
  //
  if (event_id < 100 || event_id%100 == 0) {
    G4cout << ">>> Event " << event->GetEventID() << G4endl;
    /*PrintEventStatistics(
      aHit->GetEdep(), bHit->GetEdep(),
      cHit->GetEdep(), dHit->GetEdep());*/
    //G4cout << "    " << n_trajectories 
    //       << " trajectories stored in this event." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
