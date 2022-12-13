// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Geant4Detector/Geant4DetectorService.hpp"

#include <memory>
#include <vector>

namespace Acts {
class TrackingGeometry;
namespace Experimental {
class Detector;
}
}  // namespace Acts

namespace ActsExamples {
class IContextDecorator;

namespace Geant4 {
struct Geant4Detector {
  using ContextDecorators =
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;

  using TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>;

  using DetectorPtr = std::shared_ptr<Acts::Experimental::Detector>;

  std::pair<DetectorPtr, ContextDecorators> constructDetector(
      const Geant4DetectorService::Config& cfg);

  std::pair<TrackingGeometryPtr, ContextDecorators> constructTrackingGeometry(
      const Geant4DetectorService::Config& cfg);
};
}  // namespace Geant4
}  // namespace ActsExamples
