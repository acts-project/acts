// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryModuleHelper.hpp"

#include "Acts/Geometry/GeometryModule.h"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <exception>
#include <memory>

namespace Acts::detail {

const ActsGeometryModuleV1* getGeometryModuleFromRaw(
    const char* module_abi_tag, void* (*buildFunc)(const void*, const void*)) {
  static const auto s_module = [module_abi_tag,
                                buildFunc]() -> ActsGeometryModuleV1 {
    return {
        .module_abi_tag = module_abi_tag,
        .build = buildFunc,
        .destroy =
            [](void* handle) noexcept {
              if (handle == nullptr) {
                return;
              }
              try {
                delete static_cast<TrackingGeometry*>(handle);
              } catch (...) {
              }
            },
    };
  }();

  return &s_module;
}

const ActsGeometryModuleV1* getGeometryModule(const char* module_abi_tag,
                                              BuildFunction buildFunc) {
  static BuildFunction s_buildFunc = buildFunc;

  return getGeometryModuleFromRaw(
      module_abi_tag,
      [](const void* /*userData*/, const void* loggerPtr) noexcept -> void* {
        if (loggerPtr == nullptr) {
          return nullptr;
        }
        const auto& logger = *static_cast<const Logger*>(loggerPtr);
        try {
          return s_buildFunc(logger).release();
        } catch (const std::exception& e) {
          ACTS_ERROR("Failed to build geometry module: " << e.what());
          return nullptr;
        } catch (...) {
          ACTS_ERROR("Failed to build geometry module");
          return nullptr;
        }
      });
}

}  // namespace Acts::detail
