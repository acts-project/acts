// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

class TrackingVolume;

namespace MultiLayerNavigationJsonConverter {

/// Convert a MultiLayerNavigationPolicy to Json
///
/// Serializes the grid axis configuration, casts, and bin expansion.
/// Surfaces are not included — they belong to the tracking volume.
///
/// @param policy The policy to convert
/// @return The resulting json object
nlohmann::json toJson(
    const Experimental::MultiLayerNavigationPolicy& policy);

/// Reconstruct a MultiLayerNavigationPolicy from Json
///
/// Reconstructs the empty typed grid from the serialized axis configuration,
/// then delegates to the policy constructor which fills the grid from the
/// volume's surfaces.
///
/// @param j The json object produced by toJson
/// @param gctx The geometry context
/// @param volume The tracking volume whose surfaces populate the grid
/// @param logger A logging instance
/// @return The reconstructed policy
Experimental::MultiLayerNavigationPolicy fromJson(
    const nlohmann::json& j, const GeometryContext& gctx,
    const TrackingVolume& volume, const Logger& logger);

}  // namespace MultiLayerNavigationJsonConverter
}  // namespace Acts
