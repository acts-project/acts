// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"

#include <array>
#include <memory>
#include <utility>
#include <vector>

namespace Acts {
class TrackingGeometry;
class IMaterialDecorator;
}  // namespace Acts

namespace ActsExamples {

class IContextDecorator;

struct TelescopeDetector {
  using ContextDecorators = std::vector<std::shared_ptr<IContextDecorator>>;
  using TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>;

  struct Config {
    std::vector<double> positions{{0, 30, 60, 120, 150, 180}};
    std::vector<double> stereos{{0, 0, 0, 0, 0, 0}};
    std::array<double, 2> offsets{{0, 0}};
    std::array<double, 2> bounds{{25, 100}};
    double thickness{80 * Acts::UnitConstants::um};
    int surfaceType{0};
    int axisDir{2};
  };

  Config config;
  /// The store of the detector elements (lifetime: job)
  std::vector<std::shared_ptr<const Acts::DetectorElementBase>> detectorStore;

  std::pair<TrackingGeometryPtr, ContextDecorators> finalize(
      const Config& cfg,
      const std::shared_ptr<const Acts::IMaterialDecorator>& mdecorator);
};

}  // namespace ActsExamples
