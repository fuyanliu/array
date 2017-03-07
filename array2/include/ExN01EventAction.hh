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
// $Id: ExN01EventAction.hh 69899 1013-05-17 10:05:33Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
#ifndef ExN01EventAction_h
#define ExN01EventAction_h 1
//#include "globals.hh"
//#include "G4VUserDetectorConstruction.hh"
#include "G4UserEventAction.hh"
#include "ExN01MagneticField.hh"
class G4Event;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExN01EventAction : public G4UserEventAction
{
  public:
    ExN01EventAction();
    virtual ~ExN01EventAction();

    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);
    void AddEdepA(G4double edep);
    void AddEdepB(G4double edep);
    void AddEdepC(G4double edep);
    void AddEdepD(G4double edep);
   
   private:
   	 G4double fEdepA, fEdepB, fEdepC, fEdepD;
};

inline void ExN01EventAction::AddEdepA(G4double edep) {
  fEdepA += edep; 
}

inline void ExN01EventAction::AddEdepB(G4double edep) {
  fEdepB += edep;
}

inline void ExN01EventAction::AddEdepC(G4double edep) {
  fEdepC += edep;
}

inline void ExN01EventAction::AddEdepD(G4double edep) {
  fEdepD += edep;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
