// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Seeding2/BroadTripletSeedFilter.hpp"
#include "Acts/Seeding2/BroadTripletSeedFinder.hpp"
#include "Acts/Seeding2/CylindricalSpacePointGrid2.hpp"
#include "Acts/Seeding2/DoubletSeedFinder.hpp"

#include <nlohmann/json.hpp>

namespace Acts {

void to_json(nlohmann::json& j, const SeedConfirmationRangeConfig& config);

void from_json(const nlohmann::json& j, SeedConfirmationRangeConfig& config);

}  // namespace Acts

namespace Acts::Experimental {

void to_json(nlohmann::json& j, const DoubletSeedFinder::Config& config);
void to_json(nlohmann::json& j, const DoubletSeedFinder::DerivedConfig& config);
void to_json(nlohmann::json& j, const BroadTripletSeedFinder::Options& options);
void to_json(nlohmann::json& j,
             const BroadTripletSeedFinder::TripletCuts& cuts);
void to_json(nlohmann::json& j,
             const BroadTripletSeedFinder::DerivedTripletCuts& cuts);
void to_json(nlohmann::json& j, const BroadTripletSeedFilter::Config& config);
void to_json(nlohmann::json& j,
             const CylindricalSpacePointGrid2::Config& config);

void from_json(const nlohmann::json& j, DoubletSeedFinder::Config& config);
void from_json(const nlohmann::json& j,
               DoubletSeedFinder::DerivedConfig& config);
void from_json(const nlohmann::json& j,
               BroadTripletSeedFinder::Options& options);
void from_json(const nlohmann::json& j,
               BroadTripletSeedFinder::TripletCuts& cuts);
void from_json(const nlohmann::json& j,
               BroadTripletSeedFinder::DerivedTripletCuts& cuts);
void from_json(const nlohmann::json& j, BroadTripletSeedFilter::Config& config);
void from_json(const nlohmann::json& j,
               CylindricalSpacePointGrid2::Config& config);

}  // namespace Acts::Experimental
