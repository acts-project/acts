// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryModuleHelpers.hpp"

#include "Acts/Geometry/GeometryModule.h"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <exception>
#include <functional>
#include <memory>

namespace Acts::detail {

namespace {}  // namespace

const ActsGeometryModuleV1* getGeometryModule(const char* module_abi_tag,
                                              BuildFunction buildFunc) {
  static const auto s_module = [module_abi_tag,
                                buildFunc]() -> ActsGeometryModuleV1 {
    static BuildFunction s_buildFunc = buildFunc;

    return {
        .module_abi_tag = module_abi_tag,
        .build = [](const void* userData,
                    const void* loggerPtr) noexcept -> void* {
          if (loggerPtr == nullptr) {
            // Logger can't be null!
            return nullptr;
          }

          const auto& logger = *static_cast<const Logger*>(loggerPtr);

          try {
            static_cast<void>(userData);
            return s_buildFunc(logger).release();
          } catch (const std::exception&) {
            ACTS_ERROR("Failed to build geometry module");
            return nullptr;
          } catch (...) {
            ACTS_ERROR("Failed to build geometry module");
            return nullptr;
          }
        },
        .destroy =
            [](void* handle) noexcept {
              if (handle == nullptr) {
                return;
              }
              try {
                std::cout << "Destroying geometry module" << std::endl;
                delete static_cast<TrackingGeometry*>(handle);
              } catch (...) {
              }
            },
    };
  }();

  return &s_module;
}

}  // namespace Acts::detail
