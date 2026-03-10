// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/DD4hep/GeometryModuleHelper.hpp"

#include "Acts/Geometry/GeometryModuleHelper.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

namespace ActsPlugins::DD4hep::detail {

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

          const auto& logger = *static_cast<const Acts::Logger*>(loggerPtr);

          try {
            if (userData == nullptr) {
              throw std::invalid_argument("DD4hep detector is null");
            }

            const auto& detector =
                *static_cast<const dd4hep::Detector*>(userData);
            return s_buildFunc(detector, logger).release();
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
                delete static_cast<Acts::TrackingGeometry*>(handle);
              } catch (...) {
              }
            },
    };
  }();

  return &s_module;
}

}  // namespace ActsPlugins::DD4hep::detail