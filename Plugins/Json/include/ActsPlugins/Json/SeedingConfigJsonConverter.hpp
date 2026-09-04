// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding/CylindricalSpacePointGrid.hpp"
#include "Acts/Seeding/DoubletSeedFinder.hpp"
#include "Acts/Seeding/TripletSeedFinder.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

/// @addtogroup json_plugin
/// @{

/// Convert SeedConfirmationRangeConfig to JSON
/// @param j Destination JSON object
/// @param config Source SeedConfirmationRangeConfig to convert
void to_json(nlohmann::json& j, const SeedConfirmationRangeConfig& config);

/// Convert JSON to SeedConfirmationRangeConfig
/// @param j Source JSON object
/// @param config Destination SeedConfirmationRangeConfig to populate
void from_json(const nlohmann::json& j, SeedConfirmationRangeConfig& config);

/// @}
}  // namespace Acts

namespace Acts::Experimental {

/// @ingroup json_plugin
/// @{

/// Convert DoubletSeedFinder::Config to JSON
/// @param j Destination JSON object
/// @param config Source DoubletSeedFinder::Config to convert
void to_json(nlohmann::json& j, const DoubletSeedFinder::Config& config);
/// Convert DoubletSeedFinder::DerivedConfig to JSON
/// @param j Destination JSON object
/// @param config Source DoubletSeedFinder::DerivedConfig to convert
void to_json(nlohmann::json& j, const DoubletSeedFinder::DerivedConfig& config);
/// Convert TripletSeedFinder::Config to JSON
/// @param j Destination JSON object
/// @param config Source TripletSeedFinder::Config to convert
void to_json(nlohmann::json& j, const TripletSeedFinder::Config& config);
/// Convert TripletSeedFinder::DerivedConfig to JSON
/// @param j Destination JSON object
/// @param config Source TripletSeedFinder::DerivedConfig to convert
void to_json(nlohmann::json& j, const TripletSeedFinder::DerivedConfig& config);
/// Convert BroadTripletSeedFilter::Config to JSON
/// @param j Destination JSON object
/// @param config Source BroadTripletSeedFilter::Config to convert
void to_json(nlohmann::json& j, const BroadTripletSeedFilter::Config& config);
/// Convert CylindricalSpacePointGrid::Config to JSON
/// @param j Destination JSON object
/// @param config Source CylindricalSpacePointGrid::Config to convert
void to_json(nlohmann::json& j,
             const CylindricalSpacePointGrid::Config& config);

/// Convert JSON to DoubletSeedFinder::Config
/// @param j Source JSON object
/// @param config Destination DoubletSeedFinder::Config to populate
void from_json(const nlohmann::json& j, DoubletSeedFinder::Config& config);
/// Convert JSON to DoubletSeedFinder::DerivedConfig
/// @param j Source JSON object
/// @param config Destination DoubletSeedFinder::DerivedConfig to populate
void from_json(const nlohmann::json& j,
               DoubletSeedFinder::DerivedConfig& config);
/// Convert JSON to TripletSeedFinder::Config
/// @param j Source JSON object
/// @param config Destination TripletSeedFinder::Config to populate
void from_json(const nlohmann::json& j, TripletSeedFinder::Config& config);
/// Convert JSON to TripletSeedFinder::DerivedConfig
/// @param j Source JSON object
/// @param config Destination TripletSeedFinder::DerivedConfig to populate
void from_json(const nlohmann::json& j,
               TripletSeedFinder::DerivedConfig& config);
/// Convert JSON to BroadTripletSeedFilter::Config
/// @param j Source JSON object
/// @param config Destination BroadTripletSeedFilter::Config to populate
void from_json(const nlohmann::json& j, BroadTripletSeedFilter::Config& config);
/// Convert JSON to CylindricalSpacePointGrid::Config
/// @param j Source JSON object
/// @param config Destination CylindricalSpacePointGrid::Config to populate
void from_json(const nlohmann::json& j,
               CylindricalSpacePointGrid::Config& config);

/// @}

}  // namespace Acts::Experimental
