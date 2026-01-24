// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{
class IVolumeMaterial;
class TrackingVolume;

/// Convert tracking volume and material pair to JSON
/// @param j Destination JSON object
/// @param volume Source tracking volume and material pair to convert
void to_json(
    nlohmann::json& j,
    const std::pair<const Acts::TrackingVolume*,
                    std::shared_ptr<const Acts::IVolumeMaterial>>& volume);

/// Convert TrackingVolume to JSON
/// @param j Destination JSON object
/// @param volume Source TrackingVolume to convert
void to_json(nlohmann::json& j, const Acts::TrackingVolume& volume);

/// @}
}  // namespace Acts
