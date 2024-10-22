// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"

#include <memory>
#include <vector>

namespace Acts {
class TrackingGeometry;
class DetectorElementBase;
class IMaterialDecorator;
class Logger;
namespace Experimental {
class Detector;
}  // namespace Experimental
}  // namespace Acts

namespace ActsExamples {
class IContextDecorator;
}  // namespace ActsExamples

namespace ActsExamples::DetectorCommons {

class Detector {
 public:
  using TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>;
  using DetectorPtr = std::shared_ptr<Acts::Experimental::Detector>;
  using ContextDecorators =
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;

  using DetectorElement = Acts::DetectorElementBase;
  using DetectorStore = std::vector<std::shared_ptr<const DetectorElement>>;

  explicit Detector(std::unique_ptr<const Acts::Logger> logger);
  virtual ~Detector() = default;

  virtual std::tuple<TrackingGeometryPtr, ContextDecorators, DetectorStore>
  trackingGeometry();

  virtual std::tuple<DetectorPtr, ContextDecorators, DetectorStore> detector();

  virtual void drop();

 protected:
  TrackingGeometryPtr m_trackingGeometry;
  DetectorPtr m_detector;
  ContextDecorators m_contextDecorators;
  DetectorStore m_detectorStore;

  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }

  virtual void buildTrackingGeometry() = 0;

  virtual void buildDetector();
};

}  // namespace ActsExamples::DetectorCommons
