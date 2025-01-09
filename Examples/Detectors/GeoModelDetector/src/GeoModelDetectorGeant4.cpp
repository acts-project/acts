// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GeoModelDetector/GeoModelDetector.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelGeant4DetectorConstruction.hpp"

namespace ActsExamples {

std::unique_ptr<G4VUserDetectorConstruction>
GeoModelDetector::buildGeant4DetectorConstruction(
    const Geant4ConstructionOptions& options) const {
  return std::make_unique<GeoModelGeant4DetectorConstruction>(m_geoModel,
                                                              options);
}

}  // namespace ActsExamples
