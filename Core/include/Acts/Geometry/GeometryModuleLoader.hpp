// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryModule.h"

#include <filesystem>
#include <memory>

namespace Acts {

class TrackingGeometry;

/// Returns the host ABI tag expected by the runtime module loader.
const char* geometryModuleHostAbiTag() noexcept;

/// Load a module shared library, validate ABI compatibility, build and return
/// the tracking geometry. The returned deleter keeps the module loaded until
/// the geometry is destroyed.
std::shared_ptr<TrackingGeometry> loadGeometryModule(
    const std::filesystem::path& modulePath);

}  // namespace Acts
