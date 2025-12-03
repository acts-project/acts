// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/GridJsonConverter.hpp"
#include "Acts/Navigation/IndexGridNavigationPolicy.hpp"

#include <nlohmann/json.hpp>

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.
namespace Acts::IndexGridNavigationJsonConverter {

/// Convert an IndexGridNavigationPolicy for regular cylinder to Json
nlohmann::json toJson(const Experimental::RegularCylinderIndexGridNavigationPolicy& policy);

}  // namespace Acts::IndexGridNavigationJsonConverter