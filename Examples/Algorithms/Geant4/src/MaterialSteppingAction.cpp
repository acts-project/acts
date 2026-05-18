// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/MaterialSteppingAction.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "ActsExamples/Geant4/AlgebraConverters.hpp"
#include "ActsExamples/Geant4/EventStore.hpp"
#include "ActsExamples/Geant4/UnitConversion.hpp"

#include <cstddef>
#include <ostream>
#include <unordered_map>
#include <utility>

#include <G4Material.hh>
#include <G4RunManager.hh>
#include <G4Step.hh>

namespace ActsExamples::Geant4 {

MaterialSteppingAction::MaterialSteppingAction(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : G4UserSteppingAction(), m_cfg(cfg), m_logger(std::move(logger)) {}

void MaterialSteppingAction::UserSteppingAction(const G4Step* stepPtr) {
  assert(stepPtr != nullptr);
  const G4Step& step = *stepPtr;

  // Get the material & check if it is present
  const G4Material* materialPtr = step.GetPreStepPoint()->GetMaterial();
  if (materialPtr == nullptr) {
    return;
  }
  const G4Material& material = *materialPtr;

  // First check for exclusion
  const std::string materialName = material.GetName();
  for (const std::string& emat : m_cfg.excludeMaterials) {
    if (emat == materialName) {
      ACTS_VERBOSE("Exclude step in material '" << materialName << "'.");
      return;
    }
  }

  ACTS_VERBOSE("Performing a step with step size = "
               << convertLengthToActs * step.GetStepLength());

  // Quantities valid for elemental materials and mixtures
  const double X0 = convertLengthToActs * material.GetRadlen();
  const double L0 = convertLengthToActs * material.GetNuclearInterLength();
  const double rho = convertDensityToActs * material.GetDensity();

  // Get{A,Z} is only meaningful for single-element materials (according to
  // the Geant4 docs). Need to compute average manually.
  const G4ElementVector* elements = material.GetElementVector();
  const G4double* fraction = material.GetFractionVector();
  const std::size_t nElements = material.GetNumberOfElements();
  double Ar = 0.;
  double Z = 0.;
  if (nElements == 1) {
    Ar = material.GetA() * (CLHEP::mole / CLHEP::gram);
    Z = material.GetZ();
  } else {
    for (std::size_t i = 0; i < nElements; i++) {
      Ar += elements->at(i)->GetA() * fraction[i] * (CLHEP::mole / CLHEP::gram);
      Z += elements->at(i)->GetZ() * fraction[i];
    }
  }

  // Construct passed material slab for the step
  const Acts::MaterialSlab slab(
      Acts::Material::fromMassDensity(X0, L0, Ar, Z, rho),
      convertLengthToActs * step.GetStepLength());

  // Create the RecordedMaterialSlab
  Acts::MaterialInteraction mInteraction;
  mInteraction.position =
      convertPosition(step.GetPreStepPoint()->GetPosition());
  mInteraction.direction =
      convertDirection(step.GetPreStepPoint()->GetMomentum()).normalized();
  mInteraction.materialSlab = slab;
  mInteraction.pathCorrection = step.GetStepLength() * convertLengthToActs;

  if (m_cfg.recordElementFractions) {
    if (nElements == 1) {
      mInteraction.elementZ.push_back(
          static_cast<unsigned int>(material.GetZ()));
      mInteraction.elementFrac.push_back(1.0f);
    } else {
      for (std::size_t i = 0; i < nElements; i++) {
        mInteraction.elementZ.push_back(
            static_cast<unsigned int>(elements->at(i)->GetZ()));
        mInteraction.elementFrac.push_back(static_cast<float>(fraction[i]));
      }
    }
  }

  assert(step.GetTrack() != nullptr);
  const G4Track& track = *step.GetTrack();
  const std::size_t trackId = track.GetTrackID();
  if (!eventStore().materialTracks.contains(trackId - 1)) {
    Acts::RecordedMaterialTrack rmTrack;
    const Acts::Vector3 vertex = convertPosition(track.GetVertexPosition());
    const Acts::Vector3 direction =
        convertDirection(track.GetMomentumDirection());
    rmTrack.first = {vertex, direction};
    rmTrack.second.materialInteractions.push_back(mInteraction);
    eventStore().materialTracks[trackId - 1] = rmTrack;
  } else {
    eventStore()
        .materialTracks[trackId - 1]
        .second.materialInteractions.push_back(mInteraction);
  }
}

}  // namespace ActsExamples::Geant4
