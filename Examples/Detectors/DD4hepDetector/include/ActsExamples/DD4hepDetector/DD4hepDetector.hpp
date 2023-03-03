// This file is part of the Acts project.
//
// Copyright (C) 2018-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"

#include <memory>
#include <vector>

namespace Acts {
class TrackingGeometry;
class IMaterialDecorator;
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

  std::shared_ptr<DD4hepGeometryService> geometryService;
  dd4hep::Detector* lcdd = nullptr;

  std::pair<TrackingGeometryPtr, ContextDecorators> finalize(
      DD4hepGeometryService::Config config,
      std::shared_ptr<const Acts::IMaterialDecorator> mdecorator);
};

}  // namespace DD4hep
}  // namespace ActsExamples
