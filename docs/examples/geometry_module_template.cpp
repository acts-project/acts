// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

//! [Write Plain Module]
#include "Acts/Geometry/GeometryModuleHelper.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include <stdexcept>

namespace {
std::unique_ptr<Acts::TrackingGeometry> buildMyGeometry(
    const Acts::Logger& logger) {
  ACTS_INFO("Building my geometry");
  // ... construct and return a TrackingGeometry ...
  return nullptr;
}
}  // namespace

ACTS_DEFINE_GEOMETRY_MODULE(buildMyGeometry)
//! [Write Plain Module]
