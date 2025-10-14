// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/ProtoDetector.hpp"
#include "ActsPlugins/Json/UtilitiesJsonConverter.hpp"

#include <nlohmann/json.hpp>

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.

namespace Acts {
struct ProtoDetector;
struct ProtoVolume;

/// Convert ProtoVolume to JSON
/// @param j Destination JSON object
/// @param pv Source ProtoVolume to convert
void to_json(nlohmann::json& j, const ProtoVolume& pv);

/// Convert JSON to ProtoVolume
/// @param j Source JSON object
/// @param pv Destination ProtoVolume to populate
void from_json(const nlohmann::json& j, ProtoVolume& pv);

/// Convert ProtoDetector to JSON
/// @param j Destination JSON object
/// @param pd Source ProtoDetector to convert
void to_json(nlohmann::json& j, const ProtoDetector& pd);

/// Convert JSON to ProtoDetector
/// @param j Source JSON object
/// @param pd Destination ProtoDetector to populate
void from_json(const nlohmann::json& j, ProtoDetector& pd);

}  // namespace Acts
