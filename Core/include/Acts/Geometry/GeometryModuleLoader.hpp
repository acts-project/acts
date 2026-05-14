// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#if !defined(__unix__) && !defined(__APPLE__)
#error \
    "Runtime geometry modules are only supported on Unix-like systems (Linux, macOS)."
#endif

#include "Acts/Utilities/Logger.hpp"

#include <filesystem>
#include <memory>

namespace Acts {

class TrackingGeometry;

/// Load a module shared library, validate ABI compatibility, build and return
/// the tracking geometry. The returned deleter keeps the module loaded until
/// the geometry is destroyed. Throws if the module requires user data (e.g.
/// a DD4hep module) — use the appropriate typed loader instead.
/// @param modulePath Path to the geometry module shared library.
/// @param logger Logger instance used by the module loader.
/// @return Shared pointer to the loaded tracking geometry.
std::shared_ptr<TrackingGeometry> loadGeometryModule(
    const std::filesystem::path& modulePath,
    const Logger& logger = getDummyLogger());

namespace detail {
/// Low-level loader used by typed wrappers (e.g. loadDD4hepGeometryModule).
/// Validates that the module's ABI tag matches \a expectedAbiTag.
/// Validates that the module's user_data_type matches \a expectedUserDataType
/// (nullptr means the module must declare no user data requirement).
std::shared_ptr<TrackingGeometry> loadGeometryModuleImpl(
    const std::filesystem::path& modulePath, const char* expectedAbiTag,
    const char* expectedUserDataType, const void* userData,
    const Logger& logger);
}  // namespace detail

}  // namespace Acts
