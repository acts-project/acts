// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! [Write DD4hep Module]
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "ActsPlugins/DD4hep/GeometryModuleHelper.hpp"

#include <stdexcept>

#include <DD4hep/Detector.h>

namespace {
std::unique_ptr<Acts::TrackingGeometry> buildMyGeometry(
    const dd4hep::Detector& detector, const Acts::Logger& logger) {
  (void)detector;
  ACTS_INFO("Building DD4hep geometry");
  // ... use detector to build and return a TrackingGeometry ...
  return nullptr;
}
}  // namespace

ACTS_DEFINE_DD4HEP_GEOMETRY_MODULE(buildMyGeometry)
//! [Write DD4hep Module]
