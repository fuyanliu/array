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
// $Id: ExN01CrystalSD.cc  2017-03-06 17:18:29  $
//
/// \file ExN01CrystalSD.cc
/// \brief Implementation of the ExN01CrystalSD class

#include "ExN01CrystalSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN01CrystalSD::ExN01CrystalSD(
						const G4String& name,
						const G4String& hitCollectionName)
 : G4VSensitiveDetector(name),
   fHitsCollection(0)
{	//hitCollectionName.pushback("aPilotHitsCollection");
	collectionName.insert(hitCollectionName);
	G4cout<< "******hitCollectionName is "<< hitCollectionName<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
ExN01CrystalSD::~ExN01CrystalSD(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN01CrystalSD::Initialize(G4HCofThisEvent* hce)
{
	// Create hits collection
	fHitsCollection
	  = new ExN01CrystalHitsCollection(SensitiveDetectorName, collectionName[0]);

	// Add this collection in hce
	  G4int hcID
	  	= G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
	  hce->AddHitsCollection(hcID, fHitsCollection);
	  //ExN01CrystalHit* hit = new ExN01CrystalHit();
	  fHitsCollection->insert(new ExN01CrystalHit());

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool ExN01CrystalSD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
	// energy deposit
    G4String ParticleName = step->GetTrack()->GetDynamicParticle()->GetDefinition()
    	->GetParticleName();
    
	G4double edep = step->GetTotalEnergyDeposit();
	/*if (ParticleName == "e-") {
	std::ofstream dfile("step.dat",std::ios::app);
              dfile << "Energy is "<<edep/keV <<" Particle is  "<< ParticleName<< std::endl;
    }*/
	//if (edep==0.) return false;

	/*G4TouchableHistory* touchable 
		= (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());*/
	ExN01CrystalHit* hit = new ExN01CrystalHit();
	ExN01CrystalHit* hittotal = (*fHitsCollection)[fHitsCollection->entries()-1];
	// Add values
	//hit->SetEdep(edep);
	if (ParticleName == "e-" & (edep != 0)) {
	//std::ofstream dfile("step.dat",std::ios::app);
    //          dfile << "Energy is "<<edep/keV <<" Particle is  "<< ParticleName<< std::endl;
	hit->Add(edep);
	hittotal->Add(edep);
	//fHitsCollection->insert(hittotal);
	}

	return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void ExN01CrystalSD::EndOfEvent(G4HCofThisEvent*)
{
	
}