// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/DD4hepDetector/DDG4DetectorConstruction.hpp"

#include <G4VUserDetectorConstruction.hh>

namespace ActsExamples {

std::unique_ptr<G4VUserDetectorConstruction>
DD4hepDetectorBase::buildGeant4DetectorConstruction(
    const Geant4ConstructionOptions& options) const {
  return std::make_unique<DDG4DetectorConstruction>(m_detector, options);
}

}  // namespace ActsExamples
