// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4/AlgebraConverters.hpp"

#include "ActsExamples/Geant4/UnitConversion.hpp"

namespace ActsExamples {

Acts::Vector3 Geant4::convertPosition(const G4ThreeVector& g4vec) {
  return Acts::Vector3(g4vec[0] * convertLengthToActs,
                       g4vec[1] * convertLengthToActs,
                       g4vec[2] * convertLengthToActs);
}

Acts::Vector4 Geant4::convertPosition(const G4ThreeVector& g4vec,
                                      const double time) {
  return Acts::Vector4(
      g4vec[0] * convertLengthToActs, g4vec[1] * convertLengthToActs,
      g4vec[2] * convertLengthToActs, time * convertTimeToActs);
}

Acts::Vector4 Geant4::convertMomentum(const G4ThreeVector& g4vec,
                                      const double energy) {
  return Acts::Vector4{
      convertEnergyToActs * g4vec[0], convertEnergyToActs * g4vec[1],
      convertEnergyToActs * g4vec[2], convertEnergyToActs * energy};
}

G4ThreeVector Geant4::convertPosition(const Acts::Vector3& actsVec) {
  return G4ThreeVector(actsVec[0] / convertLengthToActs,
                       actsVec[1] / convertLengthToActs,
                       actsVec[2] / convertLengthToActs);
}

Acts::Vector3 Geant4::convertDirection(const G4ThreeVector& g4vec) {
  return Acts::Vector3{g4vec[0], g4vec[1], g4vec[2]};
}

G4ThreeVector Geant4::convertDirection(const Acts::Vector3& actsVec) {
  return G4ThreeVector{actsVec[0], actsVec[1], actsVec[2]};
}

}  // namespace ActsExamples
