// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <nlohmann/json.hpp>

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.
namespace Acts {
class IVolumeMaterial;
class TrackingVolume;

using TrackingVolumePointer = const Acts::TrackingVolume*;
using TrackingVolumeAndMaterial =
    std::pair<const Acts::TrackingVolume*,
              std::shared_ptr<const Acts::IVolumeMaterial>>;

/// Conversion of a pair of tracking volume and material used for the material
/// mapping
void to_json(nlohmann::json& j, const TrackingVolumeAndMaterial& volume);

/// Conversion of a const pointer on a tracking volume used to write the
/// geometry
void to_json(nlohmann::json& j, const TrackingVolumePointer& volume);

/// Conversion of a tracking volume
void to_json(nlohmann::json& j, const Acts::TrackingVolume& volume);

}  // namespace Acts
