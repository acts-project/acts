// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/DetectorBase.hpp"

#include <memory>

namespace ActsExamples {

struct GeoModelDetectorFactory : public DetectorFactoryBase {
  struct Config {
    std::string path;

    /// Logging level of the child tools
    Acts::Logging::Level logLevel = Acts::Logging::INFO;
  };

  explicit GeoModelDetectorFactory(const Config& cfg);

  std::shared_ptr<DetectorBase> buildDetector() const override;

 private:
  Config m_cfg;
};

class GeoModelDetector : public PreConstructedDetector {
 public:
  explicit GeoModelDetector(Acts::GeoModelTree geoModel);

  std::unique_ptr<G4VUserDetectorConstruction> buildGeant4DetectorConstruction(
      const Geant4ConstructionOptions& options) const override;

 private:
  Acts::GeoModelTree m_geoModel;
};

}  // namespace ActsExamples
