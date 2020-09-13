// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Sycl/Seeding/DeviceExperimentCuts.hpp"
#include "Acts/Plugins/Sycl/Seeding/detail/Types.hpp"
#include "Acts/Plugins/Sycl/Utilities/QueueWrapper.hpp"

#include <vector>

namespace Acts::Sycl {

/// @brief Seedfinding algorithm implemented in SYCL.
///
/// @param [in] wrappedQueue is a wrapper object of the SYCL queue
/// @param [in] seedfinderConfig includes the required configruation
/// parameters for the algorithm
/// @param [out] seeds consists of the generated seeds, which can also
/// be empty if none were found
void createSeedsForGroupSycl(
    const QueueWrapper& wrappedQueue,
    const detail::DeviceSeedfinderConfig& seedfinderConfig,
    const DeviceExperimentCuts& deviceCuts,
    const std::vector<detail::DeviceSpacePoint>& bottomSPs,
    const std::vector<detail::DeviceSpacePoint>& middleSPs,
    const std::vector<detail::DeviceSpacePoint>& topSPs,
    std::vector<std::vector<detail::SeedData>>& seeds);
}  // namespace Acts::Sycl
