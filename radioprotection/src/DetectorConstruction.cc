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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

#include "DetectorConstruction.hh"
#include "DetectorConstructionMessenger.hh"

#include "G4Material.hh"
#include "G4Isotope.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4PhysicalConstants.hh"
#include "G4RunManager.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4AssemblyVolume.hh"
#include "G4Region.hh"

using namespace std;

DetectorConstruction::DetectorConstruction(AnalysisManager* analysis_manager):
	logicTreatmentRoom(0), physicalTreatmentRoom(0), logicAbsorber(0), logicTumor(0),
	aRegion(0)
{
	analysis = analysis_manager;
}

DetectorConstruction::~DetectorConstruction(){

}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	//messenger to move the TEPC
	dMessenger = new DetectorConstructionMessenger(this);

	// Sets default geometry and materials
	SetDefaultDimensions();

	// -----------------------------
	// Treatment room - World volume
	//------------------------------
	// Treatment room sizes
	const G4double worldX = 50.0 *cm;
	const G4double worldY = 50.0 *cm;
	const G4double worldZ = 50.0 *cm;
	G4bool isotopes = false;

	//G4Material* airNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR", isotopes);
	G4Material* vacuumNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic", isotopes);
	//G4Material* waterNist =  G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", isotopes);

	G4Material* worldMaterial = vacuumNist;

	G4Box* treatmentRoom = new G4Box("World",
			worldX,
			worldY,
			worldZ);

	logicTreatmentRoom = new G4LogicalVolume(treatmentRoom,
			worldMaterial,
			"World",
			0,0,0);

	physicalTreatmentRoom = new G4PVPlacement(0,
			G4ThreeVector(),
			"World",
			logicTreatmentRoom,
			0,
			false,
			0);


	// The treatment room is invisible in the Visualisation
	logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::GetInvisible());



	G4VSolid* solidAbsorber = new G4Tubs("Absorber",
			0,
			AbsR,
			0.5*AbsH,
			0.*deg,
			360.*deg);


	logicAbsorber = new G4LogicalVolume(solidAbsorber,
			AbsMaterial,
			"Absorber");

	new G4PVPlacement(0,                     //no rotation
			G4ThreeVector(0,0,0.5*AbsH),       //at (0,0,0)
			logicAbsorber,            //its logical volume
			"Absorber",               //its name
			logicTreatmentRoom,                     //its mother  volume
			false,                 //no boolean operation
			0,                     //copy number
			checkOverlaps);        //overlaps checking



	G4VSolid* solidTumor = new G4Tubs("Tumor",
                        0,
                        TumorR,
                        0.5*TumorH,
                        0.*deg,
                        360.*deg);

	logicTumor = new G4LogicalVolume(solidTumor,
			TumorMaterial,
			"Tumor");

	G4ThreeVector Tumor_pos = G4ThreeVector(0., 0., 0.);

	phyisicalTumor = new G4PVPlacement(0,
			Tumor_pos,
			logicTumor,
			"Tumor",
			logicAbsorber,
			false,
			0,
			checkOverlaps);


	return physicalTreatmentRoom; 

}
////////////////////////////////////////////////////////////////////////
void DetectorConstruction::SetDefaultDimensions()
{
	// Set of coulors that can be used
	white = new G4VisAttributes( G4Colour());
	white -> SetVisibility(true);
	white -> SetForceSolid(true);

	blue = new G4VisAttributes(G4Colour(0. ,0. ,1.));
	blue -> SetVisibility(true);
	blue -> SetForceSolid(true);

	gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5 ));
	gray-> SetVisibility(true);
	gray-> SetForceSolid(true);

	red = new G4VisAttributes(G4Colour(1. ,0. ,0.));
	red-> SetVisibility(true);
	red-> SetForceSolid(true);

	yellow = new G4VisAttributes(G4Colour(1., 1., 0. ));
	yellow-> SetVisibility(true);
	yellow-> SetForceSolid(true);

	green = new G4VisAttributes( G4Colour(25/255. , 255/255. ,  25/255. ));
	green -> SetVisibility(true);
	green -> SetForceSolid(true);

	darkGreen = new G4VisAttributes( G4Colour(0/255. , 100/255. ,  0/255. ));
	darkGreen -> SetVisibility(true);
	darkGreen -> SetForceSolid(true);

	darkOrange3 = new G4VisAttributes( G4Colour(205/255. , 102/255. ,  000/255. ));
	darkOrange3 -> SetVisibility(true);
	darkOrange3 -> SetForceSolid(false);

	skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255. ));
	skyBlue -> SetVisibility(true);
	skyBlue -> SetForceSolid(true);


	AbsR = 10*cm;
	AbsH = 5*cm;
	G4Material* PMMANist = G4NistManager::Instance()->FindOrBuildMaterial("G4_PLEXIGLASS", 0);
	//G4Material* waterNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", 0);

	//AbsMat = PMMANist;
	AbsMaterial = PMMANist;

	TumorR = 2.5*cm;
        TumorH = 5*cm;

	BrainNIST = G4NistManager::Instance()->FindOrBuildMaterial("G4_BRAIN_ICRP", false); 
	Cu63 = MaterialWithSingleIsotope("Copper63", "Cu63", 10.781791*g/cm3, 29, 63);
	Sc45 = MaterialWithSingleIsotope("Scandium45", "Sc45", 7.7023730*g/cm3, 21, 45);
	Y89 = MaterialWithSingleIsotope("Yttrium89", "Y89", 15.232262*g/cm3, 39, 89);
	F19 = MaterialWithSingleIsotope("Fluorine89", "F19", 3.25508*g/cm3, 9, 19);

	fFraction = 50;
	BrainMix = new G4Material("BrainMix", 1.04*g/cm3, 2);
	BrainMix -> AddMaterial(BrainNIST, (100-fFraction)*perCent);
	BrainMix -> AddMaterial(F19, (fFraction)*perCent);

	//TumorMaterial = PMMAF;
	//TumorMaterial = F19;
	TumorMaterial = BrainMix;
	//TumorMaterial = P31;
	//TumorMaterial = PMMANist;
	//TumorMaterial = waterNist;

	checkOverlaps = true;  //<---------------------
}



/////////////////////////////////////////////////////////////////////////////
/////////////////////////// MESSENGER ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void DetectorConstruction::SetAbsMaterial(G4String materialChoice)
{
	if (G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice, false) )
	{
		if (pttoMaterial)
		{
			AbsMaterial  = pttoMaterial;
			logicAbsorber -> SetMaterial(pttoMaterial);
			G4cout << "The material of the Absorber has been changed to " << materialChoice << G4endl;
		}
	}
	else
	{
		G4cout << "WARNING: material \"" << materialChoice << "\" doesn't exist in NIST elements/materials"
			" table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
		G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl;
	}

}

void DetectorConstruction::SetTumorMaterial(G4String materialChoice)
{
	if (G4Material* pttoMaterial = G4NistManager::Instance()->FindOrBuildMaterial(materialChoice, false) )
	{
		if (pttoMaterial)
		{
			TumorMaterial  = pttoMaterial;
			logicTumor -> SetMaterial(pttoMaterial);
			G4cout << "The material of the Tumor has been changed to " << materialChoice <<"\n NUMBER OF ATOMS PER VOLUME [cm-3] " <<TumorMaterial -> GetTotNbOfAtomsPerVolume()/(1/cm3) << G4endl;
		}
	}
	else if (materialChoice == "Cu63")
	{
		TumorMaterial  = Cu63;
                logicTumor -> SetMaterial(TumorMaterial);
                G4cout << "The material of the Tumor has been changed to " << materialChoice <<"\n NUMBER OF ATOMS PER VOLUME [cm-3] " <<TumorMaterial -> GetTotNbOfAtomsPerVolume()/(1/cm3) << G4endl;
	}
	else if (materialChoice == "Sc45")
        {
                TumorMaterial  = Sc45;
                logicTumor -> SetMaterial(TumorMaterial);
                G4cout << "The material of the Tumor has been changed to " << materialChoice <<"\n NUMBER OF ATOMS PER VOLUME [cm-3] " <<TumorMaterial -> GetTotNbOfAtomsPerVolume()/(1/cm3) << G4endl;
        }
	else if (materialChoice == "Y89")
        {
                TumorMaterial  = Y89;
                logicTumor -> SetMaterial(TumorMaterial);
                G4cout << "The material of the Tumor has been changed to " << materialChoice <<"\n NUMBER OF ATOMS PER VOLUME [cm-3] " <<TumorMaterial -> GetTotNbOfAtomsPerVolume()/(1/cm3) << G4endl;
        }
	else if (materialChoice == "F19")
        {
                TumorMaterial  = F19;
                logicTumor -> SetMaterial(TumorMaterial);
                G4cout << "The material of the Tumor has been changed to " << materialChoice <<"\n NUMBER OF ATOMS PER VOLUME [cm-3] " <<TumorMaterial -> GetTotNbOfAtomsPerVolume()/(1/cm3) << G4endl;
        }
	else if (materialChoice == "BrainMix50")
        {
		fFraction = 50;
        	BrainMix = new G4Material("BrainMix", 1.04*g/cm3, 2);
        	BrainMix -> AddMaterial(BrainNIST, (100-fFraction)*perCent);
        	BrainMix -> AddMaterial(Sc45, (fFraction)*perCent);
                
		TumorMaterial  = BrainMix;
                logicTumor -> SetMaterial(TumorMaterial);
                G4cout << "The material of the Tumor has been changed to " << materialChoice <<"\n NUMBER OF ATOMS PER VOLUME [cm-3] " <<TumorMaterial -> GetTotNbOfAtomsPerVolume()/(1/cm3) << G4endl;
        }
	else if (materialChoice == "BrainMix1")
        {
                fFraction = 1;
                BrainMix = new G4Material("BrainMix", 1.04*g/cm3, 2);
                BrainMix -> AddMaterial(BrainNIST, (100-fFraction)*perCent);
                BrainMix -> AddMaterial(Sc45, (fFraction)*perCent);
                
                TumorMaterial  = BrainMix;
                logicTumor -> SetMaterial(TumorMaterial);
                G4cout << "The material of the Tumor has been changed to " << materialChoice <<"\n NUMBER OF ATOMS PER VOLUME [cm-3] " <<TumorMaterial -> GetTotNbOfAtomsPerVolume()/(1/cm3) << G4endl;
        }
	else if (materialChoice == "BrainMix01")
        {
                fFraction = 0.1;
                BrainMix = new G4Material("BrainMix", 1.04*g/cm3, 2);
                BrainMix -> AddMaterial(BrainNIST, (100-fFraction)*perCent);
                BrainMix -> AddMaterial(Sc45, (fFraction)*perCent);
                
                TumorMaterial  = BrainMix;
                logicTumor -> SetMaterial(TumorMaterial);
                G4cout << "The material of the Tumor has been changed to " << materialChoice <<"\n NUMBER OF ATOMS PER VOLUME [cm-3] " <<TumorMaterial -> GetTotNbOfAtomsPerVolume()/(1/cm3) << G4endl;
        }


	else
	{
		G4cout << "WARNING: material \"" << materialChoice << "\" doesn't exist in NIST elements/materials"
			" table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
		G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl;
	}

}

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);

 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}


/*void DetectorConstruction::SetIsotopeFraction(G4double fraction)
  {
  if (fraction < 0)
  fFraction = 100;
  else
  {
  G4Material* 
  TumorMaterial  = pttoMaterial;
  logicTumor -> SetMaterial(pttoMaterial);
  G4cout << "The material of the Tumor has been changed to " << materialChoice << G4endl;
  }
  else
  {
  G4cout << "WARNING: material \"" << materialChoice << "\" doesn't exist in NIST elements/materials"
  " table [located in $G4INSTALL/source/materials/src/G4NistMaterialBuilder.cc]" << G4endl;
  G4cout << "Use command \"/parameter/nist\" to see full materials list!" << G4endl;
  }

  }
 */
