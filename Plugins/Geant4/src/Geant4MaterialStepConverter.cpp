// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Geant4/Geant4MaterialStepConverter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"

#include <algorithm>
#include <cstddef>

#include <CLHEP/Units/SystemOfUnits.h>
#include <G4Material.hh>
#include <G4Step.hh>
#include <G4Track.hh>

namespace {

constexpr double convertLengthToActs = Acts::UnitConstants::mm / CLHEP::mm;
constexpr double convertDensityToActs =
    (Acts::UnitConstants::g / Acts::UnitConstants::mm3) /
    (CLHEP::gram / CLHEP::mm3);

Acts::Vector3 convertPosition(const G4ThreeVector& g4vec) {
  return Acts::Vector3(g4vec[0], g4vec[1], g4vec[2]) * convertLengthToActs;
}

Acts::Vector3 convertDirection(const G4ThreeVector& g4vec) {
  return Acts::Vector3(g4vec[0], g4vec[1], g4vec[2]);
}

}  // namespace

namespace ActsPlugins {

std::optional<Acts::MaterialInteraction> Geant4MaterialStepConverter::convert(
    const G4Step& step) const {
  const G4Material* materialPtr = step.GetPreStepPoint()->GetMaterial();
  if (materialPtr == nullptr) {
    return std::nullopt;
  }
  const G4Material& material = *materialPtr;

  if (std::ranges::find(m_cfg.excludeMaterials, material.GetName()) !=
      m_cfg.excludeMaterials.end()) {
    return std::nullopt;
  }

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

  return mInteraction;
}

std::pair<Acts::Vector3, Acts::Vector3> Geant4MaterialStepConverter::trackStart(
    const G4Track& track) {
  return {convertPosition(track.GetVertexPosition()),
          convertDirection(track.GetMomentumDirection())};
}

}  // namespace ActsPlugins
