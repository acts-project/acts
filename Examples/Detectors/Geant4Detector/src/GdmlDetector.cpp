// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Geant4Detector/GdmlDetector.hpp"

#include "ActsExamples/DetectorCommons/DetectorBase.hpp"
#include "ActsExamples/Geant4Detector/GdmlDetectorConstruction.hpp"

#include <G4GDMLParser.hh>

namespace ActsExamples {

GdmlDetectorFactory::GdmlDetectorFactory(const Config& cfg)
    : DetectorFactoryBase(
          Acts::getDefaultLogger("GdmlDetectorFactory", cfg.logLevel)),
      m_cfg(cfg) {}

std::shared_ptr<DetectorBase> GdmlDetectorFactory::buildDetector() const {
  return std::make_shared<GdmlDetector>(m_cfg);
}

GdmlDetector::GdmlDetector(const GdmlDetectorFactory::Config& config)
    : PreConstructedDetector(), m_cfg(config) {}

std::unique_ptr<G4VUserDetectorConstruction>
GdmlDetector::buildGeant4DetectorConstruction(
    const Geant4ConstructionOptions& options) const {
  return std::make_unique<GdmlDetectorConstruction>(m_cfg.path, options);
}

}  // namespace ActsExamples
