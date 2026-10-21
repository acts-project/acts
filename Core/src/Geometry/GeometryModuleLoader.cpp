// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryModuleLoader.hpp"

#include "Acts/Geometry/GeometryModule.h"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include <cstring>
#include <format>
#include <stdexcept>

#include <dlfcn.h>

#ifndef ACTS_GEOMETRY_MODULE_ABI_TAG
#error \
    "ACTS_GEOMETRY_MODULE_ABI_TAG must be provided by CMake when building ActsCore."
#endif

namespace {

using GeometryModuleEntryPointV1 = const ActsGeometryModuleV1* (*)(void);

std::shared_ptr<void> openSharedLibrary(const std::filesystem::path& path) {
  if (!std::filesystem::exists(path)) {
    throw std::runtime_error(
        std::format("Geometry module file does not exist: {}", path.string()));
  }

  void* rawHandle = ::dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (rawHandle == nullptr) {
    const char* error = ::dlerror();
    throw std::runtime_error(std::format(
        "Failed to load geometry module '{}': {}", path.string(), error));
  }
  return std::shared_ptr<void>(rawHandle, [](void* handle) {
    if (handle != nullptr) {
      ::dlclose(handle);
    }
  });
}

GeometryModuleEntryPointV1 resolveEntrypointV1(
    const std::filesystem::path& path, const std::shared_ptr<void>& library) {
  ::dlerror();
  void* symbol = ::dlsym(library.get(), "acts_geometry_module_v1");
  if (const char* error = ::dlerror(); error != nullptr) {
    throw std::runtime_error(
        std::format("Failed to resolve acts_geometry_module_v1 in '{}': {}",
                    path.string(), error));
  }
  if (symbol == nullptr) {
    throw std::runtime_error(std::format(
        "Entry point acts_geometry_module_v1 resolved to nullptr in '{}'",
        path.string()));
  }
  return reinterpret_cast<GeometryModuleEntryPointV1>(symbol);
}

const char* geometryModuleHostAbiTag() noexcept {
  return ACTS_GEOMETRY_MODULE_ABI_TAG;
}

}  // namespace

namespace Acts::detail {

std::shared_ptr<TrackingGeometry> loadGeometryModuleImpl(
    const std::filesystem::path& modulePath, const char* expectedAbiTag,
    const char* expectedUserDataType, const void* userData,
    const Logger& logger) {
  auto library = openSharedLibrary(modulePath);
  auto entryPoint = resolveEntrypointV1(modulePath, library);
  const ActsGeometryModuleV1* descriptor = entryPoint();
  if (descriptor == nullptr) {
    throw std::runtime_error("Geometry module descriptor is null");
  }
  if (descriptor->module_abi_tag == nullptr || descriptor->build == nullptr ||
      descriptor->destroy == nullptr) {
    throw std::runtime_error("Geometry module descriptor is incomplete");
  }
  if (expectedAbiTag == nullptr) {
    throw std::runtime_error("Expected geometry module ABI tag is null");
  }
  if (std::strcmp(descriptor->module_abi_tag, expectedAbiTag) != 0) {
    throw std::runtime_error(std::format(
        "Geometry module ABI mismatch: module='{}' host='{}' path='{}'",
        descriptor->module_abi_tag, expectedAbiTag, modulePath.string()));
  }

  // Validate user_data_type: both sides must agree on whether userData is
  // needed and what type it is.
  if (const char* actualType = descriptor->user_data_type;
      !((expectedUserDataType == nullptr && actualType == nullptr) ||
        (expectedUserDataType != nullptr && actualType != nullptr &&
         std::strcmp(actualType, expectedUserDataType) == 0))) {
    if (actualType == nullptr) {
      throw std::runtime_error(std::format(
          "Geometry module '{}' does not require user data; "
          "use loadGeometryModule(path, logger) without extra context",
          modulePath.string()));
    } else if (expectedUserDataType == nullptr) {
      throw std::runtime_error(std::format(
          "Geometry module '{}' requires user data of type '{}'; "
          "use the appropriate typed loader (e.g. loadDD4hepGeometryModule)",
          modulePath.string(), actualType));
    } else {
      throw std::runtime_error(
          std::format("Geometry module '{}' user_data_type mismatch: "
                      "expected '{}', module declares '{}'",
                      modulePath.string(), expectedUserDataType, actualType));
    }
  }

  void* rawHandle = descriptor->build(userData, &logger);
  if (rawHandle == nullptr) {
    throw std::runtime_error("Geometry module build returned null handle");
  }

  auto destroyFn = descriptor->destroy;
  return std::shared_ptr<TrackingGeometry>(
      static_cast<TrackingGeometry*>(rawHandle),
      [destroyFn, library = std::move(library)](TrackingGeometry* geometry) {
        destroyFn(static_cast<void*>(geometry));
      });
}

}  // namespace Acts::detail

namespace Acts {

std::shared_ptr<TrackingGeometry> loadGeometryModule(
    const std::filesystem::path& modulePath, const Logger& logger) {
  return detail::loadGeometryModuleImpl(modulePath, geometryModuleHostAbiTag(),
                                        nullptr, nullptr, logger);
}

}  // namespace Acts
