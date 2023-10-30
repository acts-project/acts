// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/DD4hep/DD4hepFieldAdapter.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"

#include <memory>
#include <tuple>
#include <utility>
#include <vector>

namespace dd4hep {
class Detector;
}  // namespace dd4hep

namespace Acts {
class TrackingGeometry;
class IMaterialDecorator;
namespace Experimental {
class Detector;
}  // namespace Experimental
}  // namespace Acts

namespace ActsExamples {
class IContextDecorator;
}  // namespace ActsExamples

namespace ActsExamples {
namespace DD4hep {

struct DD4hepDetector {
  using ContextDecorators =
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;

  using TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>;

  using DetectorPtr = std::shared_ptr<const Acts::Experimental::Detector>;

  using DetectorStore = Acts::DD4hepDetectorElement::Store;

  DD4hepDetector();
  DD4hepDetector(std::shared_ptr<DD4hepGeometryService> geometryService);
  ~DD4hepDetector();

  std::shared_ptr<DD4hepGeometryService> geometryService;

  /// @brief Create a Acts::TrackingGeometry from a DD4hep::Detector
  ///
  /// @param config The configuration of the geometry service
  /// @param mdecorator The material decorator
  ///
  /// @return a pair of a TrackingGeometry and a vector of context decorators
  std::pair<TrackingGeometryPtr, ContextDecorators> finalize(
      DD4hepGeometryService::Config config,
      std::shared_ptr<const Acts::IMaterialDecorator> mdecorator);

  std::shared_ptr<Acts::DD4hepFieldAdapter> field() const;
};

}  // namespace DD4hep
}  // namespace ActsExamples
