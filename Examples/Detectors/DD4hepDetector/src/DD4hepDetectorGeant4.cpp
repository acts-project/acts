// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/DD4hepDetector/DDG4DetectorConstruction.hpp"
#include "ActsExamples/Geant4/RegionCreator.hpp"

namespace ActsExamples {

std::unique_ptr<G4VUserDetectorConstruction>
DD4hepDetector::buildGeant4DetectorConstruction(
    std::vector<std::shared_ptr<RegionCreator>> regionCreators) {
  if (m_detector == nullptr) {
    buildDD4hepGeometry();
  }

  return std::make_unique<DDG4DetectorConstruction>(shared_from_this(),
                                                    std::move(regionCreators));
}

}  // namespace ActsExamples
