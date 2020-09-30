// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// System include(s)
#include <vector>

// SYCL plugin include(s)
#include "Acts/Plugins/Sycl/Seeding/DeviceExperimentCuts.hpp"
#include "Acts/Plugins/Sycl/Seeding/detail/Types.hpp"
#include "Acts/Plugins/Sycl/Utilities/QueueWrapper.hpp"

namespace Acts::Sycl {

/// @brief Seedfinding algorithm implemented in SYCL.
///
/// @param[in] wrappedQueue is a wrapper object of the SYCL queue
/// @param[in] seedfinderConfig includes the required configuration
/// parameters for the algorithm
/// @param[in] deviceCuts is an experiment specific object with customizable
/// seed weight altering and seed cutting member functions
/// @param[in] {bottom, middle, top} SPs are arrays of simplified internal space
/// point structures of {bottom, middle, top} space points
/// @param[out] seeds holds of the generated seed indices and weight
void createSeedsForGroupSycl(
    const QueueWrapper& wrappedQueue,
    const detail::DeviceSeedfinderConfig& seedfinderConfig,
    const DeviceExperimentCuts& deviceCuts,
    const std::vector<detail::DeviceSpacePoint>& bottomSPs,
    const std::vector<detail::DeviceSpacePoint>& middleSPs,
    const std::vector<detail::DeviceSpacePoint>& topSPs,
    std::vector<std::vector<detail::SeedData>>& seeds);
}  // namespace Acts::Sycl
