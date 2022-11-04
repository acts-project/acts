// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.
namespace Acts {
namespace Experimental {

class Detector;

/// @brief Convert to json format - specific for detray
///
/// @param detector the detector instance
/// @param gctx is the gometry context
/// @param detray is a special flag to indicate detray compatible writing
/// @param logLevel a screen oputput logging level
///
/// @return a json object
nlohmann::json toJson(const Acts::Experimental::Detector& detector,
            const Acts::GeometryContext& gctx = Acts::GeometryContext(),
            bool detray = false,
            Logging::Level logLevel = Logging::INFO);

/// Converstion to Surface from jsonn
///
/// @param j the read-in json object
/// @param logLevel a screen oputput logging level
///
/// @return a shared_ptr to a Detector
std::shared_ptr<Acts::Experimental::Detector> detectorFromJson(
    const nlohmann::json& j,
    Logging::Level logLevel = Logging::INFO);

}  // namespace Experimental
}  // namespace Acts
