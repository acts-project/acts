// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryModule.h"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include <exception>
#include <functional>
#include <memory>

#ifndef ACTS_GEOMETRY_MODULE_ABI_TAG
#error \
    "ACTS_GEOMETRY_MODULE_ABI_TAG must be provided via CMake (use acts_add_geometry_module)."
#endif

namespace Acts::detail {
using BuildFunction = std::unique_ptr<TrackingGeometry> (*)(const Logger&);
const ActsGeometryModuleV1* getGeometryModule(const char* module_abi_tag,
                                              BuildFunction buildFunc);
}  // namespace Acts::detail

#define ACTS_DEFINE_GEOMETRY_MODULE(build_function)                      \
  extern "C" const ActsGeometryModuleV1* acts_geometry_module_v1(void) { \
    return Acts::detail::getGeometryModule(ACTS_GEOMETRY_MODULE_ABI_TAG, \
                                           build_function);              \
  }