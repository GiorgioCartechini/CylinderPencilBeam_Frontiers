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

#ifndef DetectorConstruction_H 
#define DetectorConstruction_H 1
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4UnionSolid.hh"
#include "G4Trd.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4RotationMatrix.hh"
#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"

#include "G4Region.hh"

class G4VPhysicalVolume;
class DetectorConstructionMessenger;
class G4LogicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(AnalysisManager* analysis);
    ~DetectorConstruction();

    G4VPhysicalVolume* Construct();


    //static TrentoPassiveProtonBeamLine* GetInstance();
    
    
    void SetAbsMaterial(G4String materialChoice);
    // Allows to change the TEPC position with respect to the WP entrance windows

    void SetTumorMaterial(G4String materialChoice);
    // Allows to change the Tumor Material

    G4Material* MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A);

private:
  AnalysisManager* analysis;
  
  
  void SetDefaultDimensions();
  void ConstructTrentoPassiveProtonBeamLine();

  
  G4LogicalVolume* logicTreatmentRoom;
  G4VPhysicalVolume* physicalTreatmentRoom;
  
  G4LogicalVolume* logicAbsorber;

  G4LogicalVolume* logicTumor;
  G4VPhysicalVolume* phyisicalTumor;
	
	G4double AbsR;
	G4double AbsH;
	G4double TumorR;
	G4double TumorH;
	
	G4Material* AbsMaterial;
	G4Material* TumorMaterial;
	G4Material* Cu63;
	G4Material* Y89;
	G4Material* Sc45;
	G4Material* F19;
	G4double fFraction;
	G4Material* BrainMix;
	G4Material* BrainNIST;

    G4VisAttributes* blue;
    G4VisAttributes* gray;
    G4VisAttributes* white;
    G4VisAttributes* red;
    G4VisAttributes* yellow;
    G4VisAttributes* green;
    G4VisAttributes* darkGreen;
    G4VisAttributes* darkOrange3;
    G4VisAttributes* skyBlue;
    
    
    G4bool checkOverlaps;
    
    G4Region* aRegion;
    
    
    DetectorConstructionMessenger* dMessenger;
};
#endif
