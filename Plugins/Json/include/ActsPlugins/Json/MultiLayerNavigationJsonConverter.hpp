// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/INavigationPolicy.hpp"
#include "Acts/Navigation/MultiLayerNavigationPolicy.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

#include <nlohmann/json.hpp>

namespace Acts {

class TrackingGeometryJsonConverter;
class TrackingVolume;

namespace MultiLayerNavigationJsonConverter {

/// Convert a MultiLayerNavigationPolicy to JSON, including its kind tag.
///
/// Serializes the grid axis configuration, casts, and bin expansion.
/// Surfaces are not included — they belong to the tracking volume.
///
/// @param policy The policy to convert
/// @param converter The tracking geometry converter (unused, required by dispatcher)
/// @return The resulting JSON object
nlohmann::json encodeMultiLayerNavigationPolicy(
    const Experimental::MultiLayerNavigationPolicy& policy,
    const TrackingGeometryJsonConverter& converter);

/// Reconstruct a MultiLayerNavigationPolicy from JSON.
///
/// Reconstructs the empty typed grid from the serialized axis configuration,
/// then delegates to the policy constructor which fills the grid from the
/// volume's surfaces.
///
/// @param encoded The JSON object produced by encodeMultiLayerNavigationPolicy
/// @param gctx The geometry context
/// @param converter The tracking geometry converter (unused, required by dispatcher)
/// @param volume The tracking volume whose surfaces populate the grid
/// @param logger A logging instance
/// @return Pointer to the reconstructed policy
std::unique_ptr<INavigationPolicy> decodeMultiLayerNavigationPolicy(
    const nlohmann::json& encoded, const GeometryContext& gctx,
    const TrackingGeometryJsonConverter& converter,
    const TrackingVolume& volume, const Logger& logger);

}  // namespace MultiLayerNavigationJsonConverter
}  // namespace Acts
