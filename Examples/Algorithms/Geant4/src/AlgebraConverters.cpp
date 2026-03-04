// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/AlgebraConverters.hpp"

#include "Acts/Definitions/Units.hpp"

#include "CLHEP/Units/SystemOfUnits.h"

namespace {
constexpr double convertLength = CLHEP::mm / Acts::UnitConstants::mm;
constexpr double convertTime = Acts::UnitConstants::ns / CLHEP::ns;
constexpr double convertEnergy = Acts::UnitConstants::GeV / CLHEP::GeV;
}  // namespace

namespace ActsExamples::Geant4 {
Acts::Vector3 convertPosition(const G4ThreeVector& g4vec) {
  return Acts::Vector3(g4vec[0] * convertLength, g4vec[1] * convertLength,
                       g4vec[2] * convertLength);
};

Acts::Vector4 convertPosition(const G4ThreeVector& g4vec, const double time) {
  return Acts::Vector4(g4vec[0] * convertLength, g4vec[1] * convertLength,
                       g4vec[2] * convertLength, time * convertTime);
}

Acts::Vector4 convertMomentum(const G4ThreeVector& g4vec, const double energy) {
  return Acts::Vector4{convertEnergy * g4vec[0], convertEnergy * g4vec[1],
                       convertEnergy * g4vec[2], convertEnergy * energy};
}

G4ThreeVector convertPosition(const Acts::Vector3& actsVec) {
  return G4ThreeVector(actsVec[0] / convertLength, actsVec[1] / convertLength,
                       actsVec[2] / convertLength);
}
Acts::Vector3 convertDirection(const G4ThreeVector& g4vec) {
  return Acts::Vector3{g4vec[0], g4vec[1], g4vec[2]};
}
G4ThreeVector convertDirection(const Acts::Vector3& actsVec) {
  return G4ThreeVector{actsVec[0], actsVec[1], actsVec[2]};
}

}  // namespace ActsExamples::Geant4
