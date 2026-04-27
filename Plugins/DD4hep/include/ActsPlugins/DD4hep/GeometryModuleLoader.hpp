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

namespace dd4hep {
class Detector;
}

namespace Acts {

class TrackingGeometry;

/// Load a DD4hep geometry module shared library with the given detector,
/// validate ABI compatibility, build and return the tracking geometry.
/// The returned deleter keeps the module loaded until the geometry is
/// destroyed.
/// @param modulePath Path to the geometry module shared library.
/// @param detector DD4hep detector instance passed to the module.
/// @param logger Logger instance used by the module loader.
/// @return Shared pointer to the loaded tracking geometry.
std::shared_ptr<TrackingGeometry> loadDD4hepGeometryModule(
    const std::filesystem::path& modulePath, const dd4hep::Detector& detector,
    const Logger& logger = getDummyLogger());

}  // namespace Acts
