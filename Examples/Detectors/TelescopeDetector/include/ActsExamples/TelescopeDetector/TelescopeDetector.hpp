// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/DetectorBase.hpp"

#include <array>
#include <memory>
#include <vector>

namespace ActsExamples {

class TelescopeDetectorFactory : public DetectorFactoryBase {
 public:
  struct Config {
    std::vector<double> positions{{0, 30, 60, 120, 150, 180}};
    std::vector<double> stereos{{0, 0, 0, 0, 0, 0}};
    std::array<double, 2> offsets{{0, 0}};
    std::array<double, 2> bounds{{25, 100}};
    double thickness{80 * Acts::UnitConstants::um};
    int surfaceType{0};
    int binValue{2};
    std::shared_ptr<const Acts::IMaterialDecorator> materialDecorator;
    Acts::Logging::Level logLevel{Acts::Logging::WARNING};
  };

  explicit TelescopeDetectorFactory(const Config& cfg);

  std::shared_ptr<DetectorBase> buildDetector() const override;

 private:
  Config m_cfg;
};

class TelescopeDetector : public PreConstructedDetector {
 public:
  TelescopeDetector(
      Acts::GeometryContext geometryContext,
      std::vector<std::shared_ptr<const Acts::DetectorElementBase>>
          detectorStore,
      std::shared_ptr<const Acts::TrackingGeometry> gen1Geometry,
      std::shared_ptr<Acts::Experimental::Detector> gen2Geometry,
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>
          contextDecorators);

  std::unique_ptr<G4VUserDetectorConstruction> buildGeant4DetectorConstruction(
      const Geant4ConstructionOptions& options) const override;

 private:
  TelescopeDetectorFactory::Config m_cfg;
};

}  // namespace ActsExamples
