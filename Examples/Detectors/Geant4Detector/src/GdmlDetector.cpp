// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4Detector/GdmlDetector.hpp"

#include "ActsExamples/Geant4Detector/GdmlDetectorConstruction.hpp"

#include <G4GDMLParser.hh>

namespace ActsExamples {

GdmlDetector::GdmlDetector(const Config& cfg)
    : Detector(Acts::getDefaultLogger("GdmlDetector", cfg.logLevel)),
      m_cfg(cfg) {}

std::unique_ptr<G4VUserDetectorConstruction>
GdmlDetector::buildGeant4DetectorConstruction(
    const Geant4ConstructionOptions& options) const {
  return std::make_unique<GdmlDetectorConstruction>(m_cfg.path, options);
}

}  // namespace ActsExamples
