// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/DD4hep/GeometryModuleLoader.hpp"

#include "Acts/Geometry/GeometryModuleLoader.hpp"

#include <DD4hep/Detector.h>

namespace Acts {

std::shared_ptr<TrackingGeometry> loadDD4hepGeometryModule(
    const std::filesystem::path& modulePath, const dd4hep::Detector& detector,
    const Logger& logger) {
  return ::Acts::detail::loadGeometryModuleImpl(
      modulePath, "dd4hep::Detector", &detector, logger);
}

}  // namespace Acts
