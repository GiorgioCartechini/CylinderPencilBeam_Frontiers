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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy


#include "DetectorConstruction.hh"
#include "DetectorConstructionMessenger.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4SystemOfUnits.hh"

   DetectorConstructionMessenger::DetectorConstructionMessenger(DetectorConstruction* detector)
:Abs_Tumor(detector)

{

    GeomDir = new G4UIdirectory("/changeAbsorberTumor/");
    GeomDir -> SetGuidance("Command to Change geometry");

    Abs_material = new G4UIcmdWithAString("/changeAbsorberTumor/AbsMaterial", this);
    Abs_material -> SetGuidance("Change the Absorber material");
    Abs_material -> SetParameterName("AbsorberMaterial", false);
    Abs_material -> SetDefaultValue("G4_A-150_TISSUE");
    Abs_material -> AvailableForStates(G4State_Idle);

    Tumor_material = new G4UIcmdWithAString("/changeAbsorberTumor/TumorMaterial", this);
    Tumor_material -> SetGuidance("Change the Tumor material");
    Tumor_material -> SetParameterName("TumorMaterial", false);
    Tumor_material -> SetDefaultValue("G4_A-150_TISSUE");
    Tumor_material -> AvailableForStates(G4State_Idle);

}

DetectorConstructionMessenger::~DetectorConstructionMessenger()
{ 
    delete Abs_material;
    delete Tumor_material;
    delete GeomDir; 
}




void DetectorConstructionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
    if( command == Abs_material )
    {
		Abs_Tumor -> SetAbsMaterial(newValue);
    }
    else if( command == Tumor_material)
    {
        Abs_Tumor -> SetTumorMaterial(newValue);
    }
}

