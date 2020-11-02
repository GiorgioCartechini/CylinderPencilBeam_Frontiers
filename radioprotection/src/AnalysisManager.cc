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

#include <stdlib.h>
#include "AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

AnalysisManager::AnalysisManager(G4String fileName)
{
  fFileName = fileName;
}

AnalysisManager::~AnalysisManager() 
{
}

void AnalysisManager::CreateNtuples() 
{ 
  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  G4cout << "Using " << manager->GetType() << G4endl;
  
  manager->SetVerboseLevel(1);
  manager->SetNtupleMerging(true);
  
  
  //creating a ntuple, containing the information about secondary particles
  manager -> CreateNtuple("GammaExitTumor", "GammaExitTumor");
  manager -> CreateNtupleDColumn("Ekin");
  manager -> CreateNtupleDColumn("PosX");
  manager -> CreateNtupleDColumn("PosY");
  manager -> CreateNtupleDColumn("PosZ");
  manager -> CreateNtupleDColumn("DirX");
  manager -> CreateNtupleDColumn("DirY");
  manager -> CreateNtupleDColumn("DirZ");
  manager -> CreateNtupleDColumn("AngleCrossBeam");
  manager -> CreateNtupleDColumn("AngleAlongBeam");
  manager -> FinishNtuple();

  manager -> CreateNtuple("GammaExitAbsorber", "GammaExitAbsorber");
  manager -> CreateNtupleDColumn("Ekin");
  manager -> CreateNtupleDColumn("PosX");
  manager -> CreateNtupleDColumn("PosY");
  manager -> CreateNtupleDColumn("PosZ");
  manager -> CreateNtupleDColumn("DirX");
  manager -> CreateNtupleDColumn("DirY");
  manager -> CreateNtupleDColumn("DirZ");
  manager -> CreateNtupleDColumn("AngleCrossBeam");
  manager -> CreateNtupleDColumn("AngleAlongBeam");
  manager -> FinishNtuple();

  
  manager -> CreateNtuple("GammaVertex", "GammaVertex");
  manager -> CreateNtupleDColumn("Ekin");
  manager -> CreateNtupleDColumn("PosX");
  manager -> CreateNtupleDColumn("PosY");
  manager -> CreateNtupleDColumn("PosZ");
  manager -> CreateNtupleDColumn("DirX");
  manager -> CreateNtupleDColumn("DirY");
  manager -> CreateNtupleDColumn("DirZ");
  manager -> CreateNtupleDColumn("AngleCrossBeam");
  manager -> CreateNtupleDColumn("AngleAlongBeam");
  manager -> FinishNtuple();
  
  manager -> CreateNtuple("GammaVertexFrag", "GammaVertexFrag");
  manager -> CreateNtupleDColumn("Ekin");
  manager -> CreateNtupleDColumn("FinalA");
  manager -> CreateNtupleDColumn("FinalZ");
  manager -> FinishNtuple();
  
  //manager -> CreateNtuple("ProtonEkin", "ProtonEkin");
  //manager -> CreateNtupleDColumn("Ekin");
  //manager -> CreateNtupleDColumn("PosZ");
  //manager -> FinishNtuple();
 manager -> CreateH2("0","Ekin Vs PosZ", 100, 0.*mm, 100*mm, 100, 0*MeV, 100*MeV);


  //manager->SetFirstNtupleId(1);
   
}

void AnalysisManager::OpenFile()
{
	G4AnalysisManager* manager = G4AnalysisManager::Instance();
	
	// Create a root file
  //G4String fileName = "radioprotection";

  // Create directories  
  //manager->SetNtupleDirectoryName("radioprotection_ntuple");
  

  G4bool fileOpen = manager->OpenFile(fFileName);
  if (!fileOpen) {
    G4cout << "\n---> HistoManager::book(): cannot open " 
           << fFileName << G4endl;
    return;
  }
  
}



void AnalysisManager::FillExitTumor(G4double Ekin, G4double PositionX, G4double PositionY, G4double PositionZ, G4double DirectionX, G4double DirectionY, G4double DirectionZ, G4double AngleCrossBeam, G4double AngleAlongBeam)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleDColumn(0, 0, Ekin);
  manager -> FillNtupleDColumn(0, 1, PositionX);
  manager -> FillNtupleDColumn(0, 2, PositionY);
  manager -> FillNtupleDColumn(0, 3, PositionZ);
  manager -> FillNtupleDColumn(0, 4, DirectionX);
  manager -> FillNtupleDColumn(0, 5, DirectionY);
  manager -> FillNtupleDColumn(0, 6, DirectionZ);
  manager -> FillNtupleDColumn(0, 7, AngleCrossBeam);
  manager -> FillNtupleDColumn(0, 8, AngleAlongBeam);
  manager -> AddNtupleRow(0);
}


void AnalysisManager::FillExitAbsorber(G4double Ekin, G4double PositionX, G4double PositionY, G4double PositionZ, G4double DirectionX, G4double DirectionY, G4double DirectionZ, G4double AngleCrossBeam, G4double AngleAlongBeam)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleDColumn(1, 0, Ekin);
  manager -> FillNtupleDColumn(1, 1, PositionX);
  manager -> FillNtupleDColumn(1, 2, PositionY);
  manager -> FillNtupleDColumn(1, 3, PositionZ);
  manager -> FillNtupleDColumn(1, 4, DirectionX);
  manager -> FillNtupleDColumn(1, 5, DirectionY);
  manager -> FillNtupleDColumn(1, 6, DirectionZ);
  manager -> FillNtupleDColumn(1, 7, AngleCrossBeam);
  manager -> FillNtupleDColumn(1, 8, AngleAlongBeam);
  manager -> AddNtupleRow(1);
}


void AnalysisManager::FillVertex(G4double Ekin, G4double PositionX, G4double PositionY, G4double PositionZ, G4double DirectionX, G4double DirectionY, G4double DirectionZ, G4double AngleCrossBeam, G4double AngleAlongBeam)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleDColumn(2, 0, Ekin);
  manager -> FillNtupleDColumn(2, 1, PositionX);
  manager -> FillNtupleDColumn(2, 2, PositionY);
  manager -> FillNtupleDColumn(2, 3, PositionZ);
  manager -> FillNtupleDColumn(2, 4, DirectionX);
  manager -> FillNtupleDColumn(2, 5, DirectionY);
  manager -> FillNtupleDColumn(2, 6, DirectionZ);
  manager -> FillNtupleDColumn(2, 7, AngleCrossBeam);
  manager -> FillNtupleDColumn(2, 8, AngleAlongBeam);
  manager -> AddNtupleRow(2);
}

void AnalysisManager::FillVertexFrag(G4double Ekin, G4int FinalA, G4int FinalZ)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  manager -> FillNtupleDColumn(3, 0, Ekin);
  manager -> FillNtupleDColumn(3, 1, FinalA);
  manager -> FillNtupleDColumn(3, 2, FinalZ);
  manager -> AddNtupleRow(3);
}

void AnalysisManager::FillProton(G4double Ekin, G4double PositionZ)
{

  G4AnalysisManager* manager = G4AnalysisManager::Instance();
  //manager -> FillNtupleDColumn(3, 0, Ekin);
  //manager -> FillNtupleDColumn(3, 1, PositionZ);
  //manager -> AddNtupleRow(3);
  manager -> FillH2(0,PositionZ, Ekin);
}
 
void AnalysisManager::finish() 
{   
    G4AnalysisManager* manager = G4AnalysisManager::Instance();    
    manager -> Write();
    manager -> CloseFile();
}












