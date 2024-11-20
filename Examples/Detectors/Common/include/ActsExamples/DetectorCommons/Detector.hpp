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

#include <memory>
#include <vector>

namespace Acts {
class GeometryContext;
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
  using ContextDecorators =
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;

  using DetectorElement = Acts::DetectorElementBase;
  using DetectorStore = std::vector<std::shared_ptr<const DetectorElement>>;

  explicit Detector(std::unique_ptr<const Acts::Logger> logger);
  Detector(const Detector&) = delete;
  Detector(Detector&&);
  virtual ~Detector();
  Detector& operator=(const Detector&) = delete;
  Detector& operator=(Detector&&);

  virtual std::tuple<std::shared_ptr<const Acts::TrackingGeometry>,
                     ContextDecorators, DetectorStore>
  trackingGeometry();

  virtual std::tuple<std::shared_ptr<Acts::Experimental::Detector>,
                     ContextDecorators, DetectorStore>
  detector();

  virtual void drop();

 protected:
  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
  std::shared_ptr<Acts::Experimental::Detector> m_detector;
  ContextDecorators m_contextDecorators;
  DetectorStore m_detectorStore;

  std::unique_ptr<const Acts::Logger> m_logger;

  const Acts::Logger& logger() const { return *m_logger; }

  virtual Acts::GeometryContext buildGeometryContext() const;

  virtual void buildTrackingGeometry(const Acts::GeometryContext& gctx) = 0;

  virtual void buildDetector();
};

}  // namespace ActsExamples::DetectorCommons
