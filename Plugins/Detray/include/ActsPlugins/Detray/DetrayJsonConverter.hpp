// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsPlugins/Json/TrackingGeometryJsonConverter.hpp"

namespace ActsPlugins {

/// Kind tag of @ref DetraySurfaceGridPolicy in the tracking geometry JSON
inline constexpr const char* kDetraySurfaceGridPolicyJsonKind =
    "DetraySurfaceGrid";

/// Teach a tracking geometry JSON converter the navigation policies of the
/// Detray plugin, so geometries carrying them can be written and read back
/// @param config The converter configuration to extend
void addDetrayNavigationPolicies(
    Acts::TrackingGeometryJsonConverter::Config& config);

}  // namespace ActsPlugins
