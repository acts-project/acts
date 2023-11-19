// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/ProtoDetector.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"

#include <nlohmann/json.hpp>

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.

namespace Acts {
struct ProtoDetector;
struct ProtoVolume;

void to_json(nlohmann::json& j, const ProtoVolume& pv);

void from_json(const nlohmann::json& j, ProtoVolume& pv);

void to_json(nlohmann::json& j, const ProtoDetector& pd);

void from_json(const nlohmann::json& j, ProtoDetector& pd);

}  // namespace Acts
