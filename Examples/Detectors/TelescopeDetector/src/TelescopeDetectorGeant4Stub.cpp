// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"

namespace ActsExamples {

std::unique_ptr<G4VUserDetectorConstruction>
TelescopeDetector::buildGeant4DetectorConstruction(
    const Geant4ConstructionOptions& /*options*/) const {
  throw std::runtime_error("Geant4 is not enabled");
}

}  // namespace ActsExamples
