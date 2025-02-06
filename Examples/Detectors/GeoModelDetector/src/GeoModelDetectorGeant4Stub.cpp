// SPDX-PackageName: "ACTS"
// SPDX-FileCopyrightText: 2016 CERN
// SPDX-License-Identifier: MPL-2.0

#include "ActsExamples/GeoModelDetector/GeoModelDetector.hpp"

namespace ActsExamples {

std::unique_ptr<G4VUserDetectorConstruction>
GeoModelDetector::buildGeant4DetectorConstruction(
    const Geant4ConstructionOptions& /*options*/) const {
  throw std::runtime_error("Geant4 is not enabled");
}

}  // namespace ActsExamples
