// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/DD4hep/GeometryModuleHelper.hpp"

#include <memory>
#include <stdexcept>

namespace {

std::unique_ptr<Acts::TrackingGeometry> buildDD4hepGeometryModule(
    const dd4hep::Detector& /*detector*/, const Acts::Logger& /*logger*/) {
  throw std::runtime_error("DD4hep geometry module stub - not implemented");
}

}  // namespace

ACTS_DEFINE_DD4HEP_GEOMETRY_MODULE(buildDD4hepGeometryModule)
