// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryModuleHelpers.hpp"

#include "GeometryModuleCommon.hpp"

namespace {

void* buildGeometryModule(const ActsGeometryModuleRequestV1* request) {
  (void)request;
  return ActsDownstream::buildTinyTrackingGeometryHandle().release();
}

void destroyGeometryModule(void* handle) {
  ActsDownstream::destroyTinyTrackingGeometryHandle(handle);
}

}  // namespace

ACTS_GEOMETRY_MODULE_DEFINE_V1(buildGeometryModule, destroyGeometryModule)

