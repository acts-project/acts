// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Experimental/NavigationDelegates.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Utilities/Logger.hpp"

// Custom Json encoder/decoders. Naming is mandated by nlohmann::json and thus
// can not match our naming guidelines.
namespace Acts {
namespace Experimental {

class Portal;
class DetectorVolume;

/// @brief Convert to json format
/// @param portal the detector portal instance
/// @param volumes is the volumes for the link association
/// @param gctx the geometry context
/// @param detray is a flag indicating if detray mode is switched on
/// @param logLevel the lscreen logging level
///
/// @note will use the default geometry context
nlohmann::json toJson(const Portal& portal,
                      const std::vector<const DetectorVolume*>& volumes = {},
                      const GeometryContext& gctx = Acts::GeometryContext(),
                      bool detray = false,
                      Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// Convert the volume link
///
/// @param delegate is the link implementation
/// @param volumes is the volumes for the link association
/// @param logLevel the lscreen logging level
///
nlohmann::json toJson(const IDelegateImpl& delegate,
                      const std::vector<const DetectorVolume*>& volumes,
                      Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// Converstion to Portal from jsonn
///
/// @param j the read-in json object
/// @param gctx the geometry context
/// @param logLevel the lscreen logging level
///
/// @return a shared_ptr to a Portal
std::shared_ptr<Portal> portalFromJson(
    const nlohmann::json& j,
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

/// Attach volume from json
///
/// @param j the read-in json object
/// @param portal the portal to which we attach the volume link
/// @param volumes already created volumes
/// @param logLevel the lscreen logging level
///
void attachVolumeLinkFromJson(
    const nlohmann::json& j, std::shared_ptr<Portal>& portal,
    const std::vector<std::shared_ptr<DetectorVolume>>& volumes,
    Acts::Logging::Level logLevel = Acts::Logging::INFO);

}  // namespace Experimental
}  // namespace Acts
