// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/GeoModelDetector/GeoModelDetector.hpp"

#include "Acts/Plugins/GeoModel/GeoModelReader.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "ActsExamples/DetectorCommons/DetectorBase.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelGeant4DetectorConstruction.hpp"

#include <G4GDMLParser.hh>

namespace ActsExamples {

GeoModelDetectorFactory::GeoModelDetectorFactory(const Config& cfg)
    : DetectorFactoryBase(
          Acts::getDefaultLogger("GeoModelDetectorFactory", cfg.logLevel)),
      m_cfg(cfg) {}

std::shared_ptr<DetectorBase> GeoModelDetectorFactory::buildDetector() const {
  return std::make_shared<GeoModelDetector>(
      Acts::GeoModelReader::readFromDb(m_cfg.path));
}

GeoModelDetector::GeoModelDetector(Acts::GeoModelTree geoModel)
    : PreConstructedDetector(), m_geoModel(std::move(geoModel)) {}

std::unique_ptr<G4VUserDetectorConstruction>
GeoModelDetector::buildGeant4DetectorConstruction(
    const Geant4ConstructionOptions& options) const {
  return std::make_unique<GeoModelGeant4DetectorConstruction>(m_geoModel,
                                                              options);
}

}  // namespace ActsExamples
