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

#include "ExN01DetectorConstruction.hh"
#include "ExN01MagneticField.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
//#include "G4Tubs.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"

#include "G4UserLimits.hh"
#include "ExN01DetectorMessenger.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"
//#include "G4SDChargedFilter.hh"
#include "ExN01CrystalSD.hh"
//#include "G4VPrimitiveScorer.hh"
//#include "G4PSEnergyDeposit.hh"


ExN01DetectorConstruction::ExN01DetectorConstruction()
 :  experimentalHall_log(0), 
    pilotA_log(0),pilotB_log(0),pilotC_log(0),pilotD_log(0),
    kapton_log(0),
    sample_log(0), 
    air_log(0),
    detector_log(0),
    experimentalHall_phys(0), 
    pilotA_phys(0),pilotB_phys(0),pilotC_phys(0),pilotD_phys(0),
    kapton_phys(0), 
    sample_phys(0),
    air_phys(0),
    detector_phys(0),
    stepLimit(0), fpMagField(0)
{  
  fpMagField = new ExN01MagneticField();
  detectorMessenger = new ExN01DetectorMessenger(this);
}

ExN01DetectorConstruction::~ExN01DetectorConstruction()
{
  delete fpMagField;
  delete stepLimit;
  delete detectorMessenger;
  //delete chamberParam;
  
}

G4VPhysicalVolume* ExN01DetectorConstruction::Construct()
{

  //------------------------------------------------------ materials

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;
  G4double natoms;
  G4double ncomponents;
  G4double fractionmass;
  G4String symbol;

  G4NistManager* nist = G4NistManager::Instance();

 /*G4Material* Ar = 
  new G4Material("ArgonGas", z= 18., a= 39.95*g/mole, density= 1.781*mg/cm3);*/
  G4bool isotopes = false;
  G4Element*  O = nist->FindOrBuildElement("O" , isotopes); 
  G4Element* Si = nist->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = nist->FindOrBuildElement("Lu", isotopes);

  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* N  = new G4Element("Nitrogen",symbol="N" , z= 7., a= 14.01*g/mole);
  //G4Element* Tl = new G4Element("Thallium",symbol="Tl", z= 81., a= 204.3833*g/mole);
  G4Element* Gd = nist->FindOrBuildElement("Gd", isotopes);
  G4Element* F = nist->FindOrBuildElement("F", isotopes);
  G4Element* Ba = nist->FindOrBuildElement("Ba", isotopes);
  G4Element* Ce = nist->FindOrBuildElement("Ce", isotopes);
  G4Element* Cs = nist->FindOrBuildElement("Cs", isotopes);
  G4Element* I = nist->FindOrBuildElement("I", isotopes);
  G4Material* Tl = nist->FindOrBuildMaterial("G4_Tl");

  G4Material* Vacuum = 
  new G4Material("Vacuum", z=1., a=1.01*g/mole,density= CLHEP::universe_mean_density, kStateGas, 2.73*kelvin, 3.e-18*pascal);
  G4bool checkOverlaps = true;    
  
  G4double expHall_x = 50*cm;
  G4double expHall_y = 50*cm;
  G4double expHall_z = 50*cm;
  G4Material* airNist = nist->FindOrBuildMaterial("G4_AIR", isotopes);
  G4Material* film_mat = nist->FindOrBuildMaterial("G4_KAPTON");
  G4Material* VINYLTOLUENE = nist->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  G4Material* water = nist->FindOrBuildMaterial("G4_WATER");


  //--------------------------------------------------
  // CsI(Tl)
  //--------------------------------------------------
  G4Material* CsI = nist->FindOrBuildMaterial("G4_CESIUM_IODIDE"); 
  //G4Material* CsI = G4Material::GetMaterial("G4_CESIUM_IODIDE");
  /*G4Material* CsI = new G4Material("CsI", 4.51*g/cm3, 2);
  CsI->AddElement(Cs, 1);
  CsI->AddElement(I, 1);*/
  G4Material* CsI_Tl = new G4Material("CsI_Tl", density = 4.51*g/cm3, ncomponents = 2);
  CsI_Tl->AddMaterial(CsI,fractionmass = 99.6*perCent);
  CsI_Tl->AddMaterial(Tl,fractionmass = 0.4*perCent);
  //--------------------------------------------------
  // LSO（pure)
  //--------------------------------------------------
  G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
  LSO->AddElement(Lu, 2);
  LSO->AddElement(Si, 1);
  LSO->AddElement(O , 5); 
  /*G4Material* LSO = new G4Material("LSO_Ce", 7.4*g/m3, 2);
  LSO->AddMaterial(LSO_pure,99.9*perCent);
  LSO->AddElement(Tl,0.1*perCent);*/

  //--------------------------------------------------
  // GSO（pure)
  //--------------------------------------------------
  G4Material* GSO = new G4Material("GSO", 6.71*g/cm3, 3);
  GSO->AddElement(Gd, 2);
  GSO->AddElement(Si, 1);
  GSO->AddElement(O,  5);

  //-------------------------------------------------
  G4double density1;
  std::vector<G4int> natoms1;
  std::vector<G4String> elements;
  //--------------------------------------------------
  // polyethylene
  //--------------------------------------------------
 
  elements.push_back("C");     natoms1.push_back(2);
  elements.push_back("H");     natoms1.push_back(4);

  density1 = 1.200*g/cm3;

  G4Material* fPethylene = nist->
          ConstructNewMaterial("Pethylene", elements, natoms1, density1);

  elements.clear();
  natoms1.clear();

// optical properties

  //--------------------------------------------------
  //optical properties of pilot-u
  //--------------------------------------------------

  const G4int wlsnum = 5;
  G4double wls_Energy[] = {2.55*eV,3.01*eV, 3.18*eV, 3.29*eV,3.55*eV};
 
  G4double rIndexPstyrene[wlsnum]={ 1.58, 1.58, 1.58, 1.58, 1.58};
  G4double absorption1[wlsnum]={210.*cm, 210.*cm, 210.*cm, 210.*cm, 210.*cm};
  G4double scintilFast[wlsnum]={0.00, 0.60, 1.00, 0.60, 0.00};
  G4MaterialPropertiesTable* fMPTPStyrene = new G4MaterialPropertiesTable();
  fMPTPStyrene->AddProperty("RINDEX",wls_Energy,rIndexPstyrene,wlsnum);
  fMPTPStyrene->AddProperty("ABSLENGTH",wls_Energy,absorption1,wlsnum);
  fMPTPStyrene->AddProperty("FASTCOMPONENT",wls_Energy, scintilFast,wlsnum);
  fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10200./MeV);
  fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE",1);
  fMPTPStyrene->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
  VINYLTOLUENE->SetMaterialPropertiesTable(fMPTPStyrene);
  // Set the Birks Constant for the Polystyrene scintillator
  VINYLTOLUENE->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  //--------------------------------------------------
  //optical properties of CsI(Tl)
  //--------------------------------------------------
  G4double photonEnergies[wlsnum] = {2.761*eV,2.485*eV,2.258*eV,2.071*eV,1.911*eV}; 
  G4double refractiveIndexCesiumIodide[wlsnum] = {1.824,1.807,1.795,1.786,1.779};
  G4double absorptionLengthCesiumIodide[wlsnum] = {39.3*cm,39.3*cm,39.3*cm,39.3*cm,39.3*cm};
  G4double Scnt_SLOW[wlsnum] = { 0.29,0.75,1.0,0.72,0.39};

  // defining an object of the G4MaterialPropertiesTable specific to CsI
  G4MaterialPropertiesTable* mptCesiumIodide = new G4MaterialPropertiesTable();
  // add the optical properties of Cesium Iodide to this object
  //  mptCesiumIodide->AddProperty("FASTCOMPONENT", photonEnergies, Scnt_FAST, wlsnum);
  mptCesiumIodide->AddProperty("FASTCOMPONENT", photonEnergies, Scnt_SLOW, wlsnum);
  mptCesiumIodide->AddProperty("RINDEX",photonEnergies,refractiveIndexCesiumIodide,wlsnum);
  mptCesiumIodide->AddProperty("ABSLENGTH",photonEnergies,absorptionLengthCesiumIodide,wlsnum);
  mptCesiumIodide->AddConstProperty("SCINTILLATIONYIELD", 54000./MeV);
  mptCesiumIodide->AddConstProperty("YIELDRATIO", 0.0);
  mptCesiumIodide->AddConstProperty("FASTTIMECONSTANT", 1000.*ns);
  //mptCesiumIodide->AddConstProperty("SLOWTIMECONSTANT", 1000.*ns);
  mptCesiumIodide->AddConstProperty("RESOLUTIONSCALE", 1.0);
  // set these to Cesium Iodide
  CsI_Tl->SetMaterialPropertiesTable(mptCesiumIodide);
  CsI_Tl->GetIonisation()->SetBirksConstant(0.00152*mm/MeV);

  //--------------------------------------------------
  //optical properties of LSO
  //--------------------------------------------------

  const G4int LSO_NUMENTRIES = 5;
  G4double LSO_Energies[LSO_NUMENTRIES] = {3.544*eV, 3.101*eV, 2.757*eV, 2.481*eV, 2.255*eV};
  G4double LSO_RINDEX[LSO_NUMENTRIES] = { 1.82, 1.82, 1.82, 1.82, 1.82};
  G4double LSO_ABS_LENGTH[LSO_NUMENTRIES] = { 50.*m, 50.*m, 50.0*m, 50.0*m, 50.0*m};
  G4double LSO_scintFAST[LSO_NUMENTRIES] = {0.00, 0.65, 0.63, 0.14, 0.00};
  G4MaterialPropertiesTable* LSO_mpt = new G4MaterialPropertiesTable();
  LSO_mpt->AddProperty("FASTCOMPONENT", LSO_Energies, LSO_scintFAST, LSO_NUMENTRIES);
  LSO_mpt->AddProperty("RINDEX",        LSO_Energies, LSO_RINDEX, LSO_NUMENTRIES);
  LSO_mpt->AddProperty("ABSLENGTH",     LSO_Energies, LSO_ABS_LENGTH, LSO_NUMENTRIES);
  LSO_mpt->AddConstProperty("SCINTILLATIONYIELD",26000./MeV); 
  LSO_mpt->AddConstProperty("RESOLUTIONSCALE",4.41);
  LSO_mpt->AddConstProperty("FASTTIMECONSTANT",40.*ns);
  LSO_mpt->AddConstProperty("YIELDRATIO",1.0);
  LSO->SetMaterialPropertiesTable(LSO_mpt);
  LSO->GetIonisation()->SetBirksConstant(0.00700*mm/MeV);

  //--------------------------------------------------
  //optical properties of GSO
  //--------------------------------------------------

  const G4int GSO_NUMENTRIES = 5;
  G4double GSO_Energies[LSO_NUMENTRIES] = {3.101*eV, 2.885*eV, 2.697*eV, 2.532*eV, 2.386*eV};
  G4double GSO_RINDEX[LSO_NUMENTRIES] = { 1.87, 1.87, 1.87, 1.87, 1.87};
  G4double GSO_ABS_LENGTH[LSO_NUMENTRIES] = { 340.*cm, 340.*cm, 340.0*cm, 340.0*cm, 340.0*cm};
  G4double GSO_scintFAST[LSO_NUMENTRIES] = {0.29, 1.00, 0.93, 0.24, 0.18};
  G4MaterialPropertiesTable* GSO_mpt = new G4MaterialPropertiesTable();
  GSO_mpt->AddProperty("FASTCOMPONENT", GSO_Energies, GSO_scintFAST, GSO_NUMENTRIES);
  GSO_mpt->AddProperty("RINDEX",        GSO_Energies, GSO_RINDEX, GSO_NUMENTRIES);
  GSO_mpt->AddProperty("ABSLENGTH",     GSO_Energies, GSO_ABS_LENGTH, GSO_NUMENTRIES);
  GSO_mpt->AddConstProperty("SCINTILLATIONYIELD",8000./MeV); 
  GSO_mpt->AddConstProperty("RESOLUTIONSCALE",1.0);
  GSO_mpt->AddConstProperty("FASTTIMECONSTANT",60.*ns);
  GSO_mpt->AddConstProperty("YIELDRATIO",0.0);
  GSO->SetMaterialPropertiesTable(GSO_mpt);
  GSO->GetIonisation()->SetBirksConstant(0.00525*mm/MeV);

//------------------------------World----------------------------
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       expHall_x, expHall_y, expHall_z);     //its size
      
  experimentalHall_log =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        Vacuum,           //its material
                        "World");            //its name
                                   
  experimentalHall_phys = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      experimentalHall_log,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking


/*
//----------------------------------air gap
  G4double air_x = 0.5*cm;
  G4double air_y = 0.5*cm;
  G4double air_z = 3.5*um;
  G4Box* air_box =    
    new G4Box("air",                      //its name
             air_x,air_y, air_z); //its size
                
  air_log =                         
    new G4LogicalVolume(air_box,         //its solid
                        film_mat,          //its material
                        "air");           //its name
               
  air_phys =
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,7*um),                    //at position
                    air_log,             //its logical volume
                    "air",                //its name
                    experimentalHall_log,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                
G4VisAttributes* air_logVisAtt
    = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  air_log->SetVisAttributes(air_logVisAtt);
*/
//--------------------------------cysatal detector----------------
  G4double block_x = 5*mm;
  G4double block_y = block_x;
  G4double block_z = 25*mm;
  G4Box* pilotA_box =    
    new G4Box("PilotA",                    //its name
        block_x, block_y, block_z); //its size
      
  pilotA_log =                         
    new G4LogicalVolume(pilotA_box,            //its solid
                        VINYLTOLUENE,             //its material
                        "PilotALV",0,0,0);         //its name
               
  pilotA_phys =
    new G4PVPlacement(0 ,                      //no rotation
                    G4ThreeVector(0.*cm,0.*cm,(200.+block_z/mm)*mm),         //at (-1mm,0,0)
                    pilotA_log,                //its logical volume
                    "PilotA",              //its name
                    experimentalHall_log,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  G4Box* pilotB_box =    
    new G4Box("PilotB",                    //its name
        block_x, block_z, block_y); //its size
      
  pilotB_log =                         
    new G4LogicalVolume(pilotB_box,            //its solid
                        VINYLTOLUENE,             //its material
                        "PilotBLV",0,0,0);         //its name
               
  pilotB_phys =
    new G4PVPlacement(0 ,                      //no rotation
                    G4ThreeVector(0.*cm,(200.+block_z/mm)*mm,0.*cm),         //at (-1mm,0,0)
                    pilotB_log,                //its logical volume
                    "PilotB",              //its name
                    experimentalHall_log,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

G4Box* pilotC_box =    
    new G4Box("PilotC",                    //its name
        block_x, block_y, block_z); //its size
      
  pilotC_log =                         
    new G4LogicalVolume(pilotC_box,            //its solid
                        VINYLTOLUENE,             //its material
                        "PilotCLV",0,0,0);         //its name
               
  pilotC_phys =
    new G4PVPlacement(0 ,                      //no rotation
                    G4ThreeVector(0.*cm,0.*cm,-(200.+block_z/mm)*mm),         //at (-1mm,0,0)
                    pilotC_log,                //its logical volume
                    "PilotC",              //its name
                    experimentalHall_log,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  G4Box* pilotD_box =    
    new G4Box("PilotD",                    //its name
        block_x, block_z, block_y); //its size
      
  pilotD_log =                         
    new G4LogicalVolume(pilotD_box,            //its solid
                        VINYLTOLUENE,             //its material
                        "PilotDLV",0,0,0);         //its name
               
  pilotD_phys =
    new G4PVPlacement(0 ,                      //no rotation
                    G4ThreeVector(0.*cm,-(200.+block_z/mm)*mm,0.*cm),         //at (-1mm,0,0)
                    pilotD_log,                //its logical volume
                    "PilotD",              //its name
                    experimentalHall_log,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
/*//----------------------------- gap --------------------------------
  //
    
  //G4Material* sample_mat = nist->FindOrBuildMaterial("G4_Fe");

  //定义样品和kapton膜的间隙gap

  G4double sample_x = block_x;
  G4double sample_y = 100*um;
  G4double sample_z = 3.5*cm;
  G4Box* sample_box =    
    new G4Box("gap",                      //its name
             sample_x,sample_y, sample_z); //its size
                
  sample_log =                         
    new G4LogicalVolume(sample_box,         //its solid
                        Vacuum,          //its material
                        "gap");           //its name
               
  sample_phys =
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0.,(2*block_y/um + sample_y/um)*um,(sample_z/cm)*cm),                    //at position
                    sample_log,             //its logical volume
                    "gap",                //its name
                    experimentalHall_log,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
                
G4VisAttributes* sample_logVisAtt
    = new G4VisAttributes(G4Colour(0.88,1.0,0.0));
  sample_log->SetVisAttributes(sample_logVisAtt);
*/
/*//********************************************** PilotB ***********************************
  //G4Material* film_mat = nist->FindOrBuildMaterial("G4_KAPTON");
        
  // Conical section shape       
  G4double film_x = block_x;
  G4double film_y = film_x;
  G4double film_z = 3.5*cm;
  G4Box* kapton_box =    
    new G4Box("PilotB", 
    film_x, film_y, film_z);
                      
  kapton_log=                         
    new G4LogicalVolume(kapton_box,         //its solid
                        VINYLTOLUENE,          //its material
                        "PilotB");           //its name
  
  kapton_phys =             
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0.*um,(2*block_y/cm + 2*sample_y/cm + film_y/cm)*cm,(film_z/cm)*cm),                    //at position
                    kapton_log,             //its logical volume
                    "PilotB",                //its name
                    experimentalHall_log,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  G4VisAttributes* kaptonVisAtt
    = new G4VisAttributes(G4Colour(0.5,1.0,0.0));
  kapton_log->SetVisAttributes(kaptonVisAtt);*/

//------------------------define a detector for the gamma that don't deposit its all energy----------
  G4Sphere* detector_sph = new G4Sphere("detector_sph",
                            300.*mm,
                            310*mm,
                            0.*deg,
                            360.*deg,
                            0.*deg,
                            360.*deg);
  detector_log = new G4LogicalVolume(detector_sph,
                                      Vacuum,
                                      "detector");
  detector_phys = new G4PVPlacement(0,
                                      G4ThreeVector(0.,0.,0.),
                                      detector_log,
                                      "detector",
                                      experimentalHall_log,
                                      false,
                                      0, 
                                      checkOverlaps);
  G4VisAttributes* detectorVisAtt
    = new G4VisAttributes(G4Colour(0.5,0.8,0.8));
  detectorVisAtt->SetForceWireframe(checkOverlaps); 
  detectorVisAtt->SetForceAuxEdgeVisible(checkOverlaps);
  detector_log->SetVisAttributes(detectorVisAtt);

  //------------------------------------------------------------------
G4double maxStep = 1000*nm;
  stepLimit = new G4UserLimits(maxStep);
  experimentalHall_log->SetUserLimits(stepLimit);
  return experimentalHall_phys;
}

void ExN01DetectorConstruction::ConstructSDandField()
{
  ExN01CrystalSD* aDetectorSD
    = new ExN01CrystalSD("aPilotSD", "aPilotHitsCollection");
  SetSensitiveDetector("PilotALV", aDetectorSD);

  ExN01CrystalSD* bDetectorSD
    = new ExN01CrystalSD("bPilotSD", "bPilotHitsCollection");
  SetSensitiveDetector("PilotBLV", bDetectorSD);

  ExN01CrystalSD* cDetectorSD
    = new ExN01CrystalSD("cPilotSD", "cPilotHitsCollection");
  SetSensitiveDetector("PilotCLV", cDetectorSD);

  ExN01CrystalSD* dDetectorSD
    = new ExN01CrystalSD("dPilotSD", "dPilotHitsCollection");
  SetSensitiveDetector("PilotDLV", dDetectorSD);

}

void ExN01DetectorConstruction::SetMagField(G4double fieldValue)
{
  fpMagField->SetMagFieldValue(fieldValue);  
G4cout<<"Detector check fieldValue : "<<fieldValue<<G4endl;
}

void ExN01DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit)&&(maxStep>0.)) stepLimit->SetMaxAllowedStep(maxStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


