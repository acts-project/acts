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
#include "ActsExamples/DetectorCommons/Detector.hpp"

#include <array>
#include <memory>
#include <vector>

namespace ActsExamples::Telescope {

class TelescopeDetectorElement;
class TelescopeG4DetectorConstruction;

class TelescopeDetector : public ActsExamples::DetectorCommons::Detector {
 public:
  using TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>;
  using ContextDecorators =
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;

  using DetectorElement = ActsExamples::Telescope::TelescopeDetectorElement;

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

  explicit TelescopeDetector(const Config& cfg);

 private:
  Config m_cfg;

  void buildTrackingGeometry() final;
};

}  // namespace ActsExamples::Telescope
