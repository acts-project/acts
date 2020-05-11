// This file is part of the Acts project.
//
// Copyright (C) 2017-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "SteppingAction.hpp"

#include <G4Material.hh>
#include <G4Step.hh>
#include <stdexcept>

#include "Acts/Utilities/Units.hpp"

using namespace ActsExamples;

SteppingAction* SteppingAction::s_instance = nullptr;

SteppingAction* SteppingAction::instance() {
  return s_instance;
}

SteppingAction::SteppingAction() : G4UserSteppingAction() {
  if (s_instance) {
    throw std::logic_error(
        "Attempted to duplicate the SteppingAction singleton");
  } else {
    s_instance = this;
  }
}

SteppingAction::~SteppingAction() {
  s_instance = nullptr;
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
  // get the material
  G4Material* material = step->GetPreStepPoint()->GetMaterial();

  if (material && material->GetName() != "Vacuum" &&
      material->GetName() != "Air") {
    // go through the elements of the material & weigh it with its fraction
    const G4ElementVector* elements = material->GetElementVector();
    const G4double* fraction = material->GetFractionVector();
    size_t nElements = material->GetNumberOfElements();
    double A = 0.;
    double Z = 0.;
    double X0 = material->GetRadlen();
    double L0 = material->GetNuclearInterLength();
    double rho = material->GetDensity() * CLHEP::mm3 / CLHEP::gram;
    double steplength = step->GetStepLength() / CLHEP::mm;
    if (nElements == 1) {
      A = material->GetA() * CLHEP::mole / CLHEP::gram;
      Z = material->GetZ();
    } else {
      for (size_t i = 0; i < nElements; i++) {
        A += elements->at(i)->GetA() * fraction[i] * CLHEP::mole / CLHEP::gram;
        Z += elements->at(i)->GetZ() * fraction[i];
      }
    }

    /*   G4cout << *material << G4endl;
       G4cout << "----G4StepMaterial----" << G4endl;
       /// @TODO remove output after testing
       G4cout << "Material: " << material->GetName() << G4endl;
       G4cout << "X0: " << X0 << G4endl;
       G4cout << "L0: " << L0 << G4endl;
       G4cout << "A: " << A << G4endl;
       G4cout << "Z: " << Z << G4endl;
       G4cout << "rho: " << rho << G4endl;
       G4cout << "steplength: " << steplength << G4endl;*/

    // create the RecordedMaterialProperties
    const auto& rawPos = step->GetPreStepPoint()->GetPosition();
    const auto& rawDir = step->GetPreStepPoint()->GetMomentum();
    Acts::MaterialInteraction mInteraction;
    mInteraction.position = Acts::Vector3D(rawPos.x(), rawPos.y(), rawPos.z());
    mInteraction.direction = Acts::Vector3D(rawDir.x(), rawDir.y(), rawDir.z());
    mInteraction.direction.normalized();
    mInteraction.materialProperties =
        Acts::MaterialProperties(X0, L0, A, Z, rho, steplength);
    m_materialSteps.push_back(mInteraction);

    //   // Get the track associated to the step
    //   G4Track* track = step->GetTrack();
    //   const auto& trkPos = track->GetPosition();
    //   const auto& trkTime = track->GetGlobalTime();
    //   // Get the particle from the track
    //   const auto& par = track->GetDynamicParticle();
    //   std::string parName = track->GetDefinition()->GetParticleName();
    //   const auto& par4Mom = par->Get4Momentum();
    //
    //   if (track->GetTrackID() == 1) {  // only consider primary tracks
    //     /*G4cout <<
    // "px: " << par4Mom.px() << "\t" <<
    // "py: " << par4Mom.py() << "\t" <<
    // "pz: " << par4Mom.pz() << "\t" <<
    // "e: " << par4Mom.e() << G4endl;
    //
    //     G4cout <<
    // "x: " << trkPos.x() << "\t" <<
    // "y: " << trkPos.y() << "\t" <<
    // "z: " << trkPos.z() << "\t" <<
    // "t: " << trkTime << G4endl;
    //     */
    //
    //     m_tracksteps.emplace_back(
    //         0,
    //         0,  // set Acts::GeometryID = 0 and Barcode = 0
    //         Acts::ActsVectorD<4>(trkPos.x() * Acts::UnitConstants::mm,
    //                              trkPos.y() * Acts::UnitConstants::mm,
    //                              trkPos.z() * Acts::UnitConstants::mm,
    //                              trkTime * Acts::UnitConstants::ns),  // pos4
    //         Acts::ActsVectorD<4>(
    //             par4Mom.px() * Acts::UnitConstants::MeV,
    //             par4Mom.py() * Acts::UnitConstants::MeV,
    //             par4Mom.pz() * Acts::UnitConstants::MeV,
    //             par4Mom.e() * Acts::UnitConstants::MeV),  // before4
    //         Acts::ActsVectorD<4>(0, 0, 0, 0));            // after4
    //   }
  }
}

void SteppingAction::clear() {
  m_materialSteps.clear();
  m_trackSteps.clear();
}
