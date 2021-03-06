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
#include "ExN01Run.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN01EventAction::ExN01EventAction()
: G4UserEventAction(),
  fEdep(0.)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
ExN01EventAction::~ExN01EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExN01EventAction::BeginOfEventAction(const G4Event* evt)
{
  fEdep = 0.;
  G4int event_id = evt->GetEventID();
  std::ofstream dfile("ao.dat",std::ios::app);
                dfile << event_id << std::endl;
  std::ofstream efile("bo.dat",std::ios::app);
                efile << event_id << std::endl;
  std::ofstream ffile("co.dat",std::ios::app);
                ffile << event_id << std::endl;
  std::ofstream gfile("do.dat",std::ios::app);
                gfile << event_id << std::endl;
  //std::ofstream afile("volume.dat",std::ios::app);
  //              afile <<event_id << std::endl;
  //std::ofstream bfile("detector.dat",std::ios::app);
  //              bfile << event_id <<std::endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void ExN01EventAction::EndOfEventAction(const G4Event* evt)
{
  ExN01Run* run
    = static_cast<ExN01Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    run->AddEdep(fEdep);


  G4int event_id = evt->GetEventID();
  //G4cout<< fEdep/CLHEP::keV<<" ***** "<<event_id<<G4endl;
  std::ofstream dfile("ao.dat",std::ios::app);
                //dfile << fEdep/CLHEP::keV << std::endl;
              dfile << event_id << std::endl;
  std::ofstream efile("bo.dat",std::ios::app);
                efile << event_id << std::endl;
  std::ofstream ffile("co.dat",std::ios::app);
                ffile << event_id << std::endl;
  std::ofstream gfile("do.dat",std::ios::app);
                gfile << event_id << std::endl;
  //std::ofstream afile("volume.dat",std::ios::app);
  //              afile <<event_id << std::endl;
  //std::ofstream bfile("detector.dat",std::ios::app);
  //              bfile <<event_id <<std::endl;  
  // get number of stored trajectories
  //
  G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  
  // periodic printing
  //
  if (event_id < 100 || event_id%100 == 0) {
    G4cout << ">>> Event " << evt->GetEventID() << G4endl;
    //G4cout << "    " << n_trajectories 
    //       << " trajectories stored in this event." << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
