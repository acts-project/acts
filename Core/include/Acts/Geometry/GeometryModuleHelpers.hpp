// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryModule.h"

#include <exception>

#ifndef ACTS_GEOMETRY_MODULE_ABI_TAG
#error "ACTS_GEOMETRY_MODULE_ABI_TAG must be provided via CMake (link Acts::Core or use acts_add_geometry_module)."
#endif

#define ACTS_GEOMETRY_MODULE_DEFINE_V1(BuildFunction, DestroyFunction)         \
  namespace {                                                                   \
  void* acts_geometry_module_build_v1(void* user_data) noexcept {               \
    try {                                                                       \
      return BuildFunction(user_data);                                          \
    } catch (const std::exception&) {                                           \
      return nullptr;                                                           \
    } catch (...) {                                                             \
      return nullptr;                                                           \
    }                                                                           \
  }                                                                             \
                                                                                \
  void acts_geometry_module_destroy_v1(void* handle) noexcept {                 \
    if (handle == nullptr) {                                                    \
      return;                                                                   \
    }                                                                           \
    try {                                                                       \
      DestroyFunction(handle);                                                  \
    } catch (...) {                                                             \
    }                                                                           \
  }                                                                             \
                                                                                \
  const ActsGeometryModuleV1 g_actsGeometryModuleV1 = {                         \
      ACTS_GEOMETRY_MODULE_ABI_TAG,                                             \
      &acts_geometry_module_build_v1,                                           \
      &acts_geometry_module_destroy_v1};                                        \
  }                                                                             \
                                                                                \
  extern "C" const ActsGeometryModuleV1* acts_geometry_module_v1(void) {       \
    return &g_actsGeometryModuleV1;                                             \
  }
