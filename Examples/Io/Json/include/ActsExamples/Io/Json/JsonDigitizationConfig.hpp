// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Json/GeometryHierarchyMapJsonConverter.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"

#include <nlohmann/json.hpp>

namespace ActsExamples {

void to_json(nlohmann::json& j, const ParameterSmearingConfig& psc);

void from_json(const nlohmann::json& j, ParameterSmearingConfig& psc);

void to_json(nlohmann::json& j, const GeometricDigitizationConfig& gdc);

void from_json(const nlohmann::json& j, GeometricDigitizationConfig& gdc);

void to_json(nlohmann::json& j, const SmearingConfig& sdc);

void from_json(const nlohmann::json& j, SmearingConfig& sdc);

void to_json(nlohmann::json& j, const DigitizationConfig& dc);

void from_json(const nlohmann::json& j, DigitizationConfig& dc);

using DigiConfigContainer = Acts::GeometryHierarchyMap<DigitizationConfig>;
using DigiConfigConverter =
    Acts::GeometryHierarchyMapJsonConverter<DigitizationConfig>;

}  // namespace ActsExamples
