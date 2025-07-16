// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

void to_json(nlohmann::json& j, const Direction& direction);

void from_json(const nlohmann::json& j, Direction& direction);

}  // namespace Acts
