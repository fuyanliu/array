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
// $Id: ExN01CrystalHit.hh  2017-03-06 17:18:29  $
//
/// \file ExN01CrystalHit.hh
/// \brief Implementation of the ExN01CrystalHit class

#ifndef ExN01CrystalHit_h
#define ExN01CrystalHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"
#include "G4SystemOfUnits.hh"

/// Calorimeter hit class
///
/// It defines data members to store the the energy deposit
/// of charged particles in a selected volume:
/// - fEdep

class ExN01CrystalHit : public G4VHit
{
  public:
    ExN01CrystalHit();
    ExN01CrystalHit(const ExN01CrystalHit&);
    virtual ~ExN01CrystalHit();

    // operators
    const ExN01CrystalHit& operator=(const ExN01CrystalHit&);
    G4int operator==(const ExN01CrystalHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // methods to handle data
    void SetEdep(G4double de);
    void Add(G4double de);

    // get methods
    G4double GetEdep() const;
      
  private:
    G4double fEdep;        ///< Energy deposit in the sensitive volume
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<ExN01CrystalHit> ExN01CrystalHitsCollection;

extern G4ThreadLocal G4Allocator<ExN01CrystalHit>* ExN01CrystalHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* ExN01CrystalHit::operator new(size_t)
{
  if(!ExN01CrystalHitAllocator)
    ExN01CrystalHitAllocator = new G4Allocator<ExN01CrystalHit>;
  void *hit;
  hit = (void *) ExN01CrystalHitAllocator->MallocSingle();
  return hit;
}

inline void ExN01CrystalHit::operator delete(void *hit)
{
  if(!ExN01CrystalHitAllocator)
    ExN01CrystalHitAllocator = new G4Allocator<ExN01CrystalHit>;
  ExN01CrystalHitAllocator->FreeSingle((ExN01CrystalHit*) hit);
}
inline void ExN01CrystalHit::SetEdep(G4double de){
  fEdep = de;
}
inline void ExN01CrystalHit::Add(G4double de) {
  fEdep += de; 
  //G4cout << "*********@@@@@fEdep is "<<fEdep/CLHEP::keV << std::endl;
}

inline G4double ExN01CrystalHit::GetEdep() const { 
  return fEdep; 
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

