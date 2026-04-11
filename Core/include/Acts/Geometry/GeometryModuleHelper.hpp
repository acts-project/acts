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

#include <memory>

namespace Acts::detail {
using BuildFunction = std::unique_ptr<TrackingGeometry> (*)(const Logger&);
const ActsGeometryModuleV1* getGeometryModule(const char* module_abi_tag,
                                              const char* user_data_type,
                                              BuildFunction buildFunc);
// Low-level shared helper: accepts a raw build function pointer matching the
// C ABI struct's .build field. Handles static struct init and destroy.
const ActsGeometryModuleV1* getGeometryModuleFromRaw(
    const char* module_abi_tag, const char* user_data_type,
    void* (*buildFunc)(const void*, const void*));
}  // namespace Acts::detail

// Internal — do not use directly.
#define ACTS_IMPL_GEOMETRY_MODULE_ENTRY(get_module_expr)                 \
  extern "C" const ActsGeometryModuleV1* acts_geometry_module_v1(void) { \
    return (get_module_expr);                                            \
  }

// Emit a clear error only when the macro is actually expanded without the tag,
// rather than at include time — Acts internal code includes this header too.
#ifdef ACTS_GEOMETRY_MODULE_ABI_TAG
#define ACTS_DEFINE_GEOMETRY_MODULE(build_function)                \
  ACTS_IMPL_GEOMETRY_MODULE_ENTRY(Acts::detail::getGeometryModule( \
      ACTS_GEOMETRY_MODULE_ABI_TAG, nullptr, (build_function)))
#else
#define ACTS_DEFINE_GEOMETRY_MODULE(build_function)                  \
  static_assert(false,                                               \
                "ACTS_GEOMETRY_MODULE_ABI_TAG must be provided via " \
                "CMake (use acts_add_geometry_module).")
#endif
