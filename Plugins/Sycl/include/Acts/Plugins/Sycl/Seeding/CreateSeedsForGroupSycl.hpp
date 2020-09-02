// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Sycl/Seeding/DeviceExperimentCuts.hpp"
#include "Acts/Plugins/Sycl/Seeding/detail/Types.h"

#include <memory>
#include <vector>

/// Forward declaration of incomplete type cl::sycl::queue
inline namespace cl {
namespace sycl {
class queue;
}
};  // namespace cl

namespace Acts::Sycl {

/// @brief Seedfinding algorithm implemented with SYCL to offload computations
/// to GPUs.
void createSeedsForGroupSycl(
    const std::shared_ptr<cl::sycl::queue>& q,
    const detail::DeviceSeedfinderConfig& configData,
    const DeviceExperimentCuts& deviceCuts,
    const std::vector<detail::DeviceSpacePoint>& bottomSPs,
    const std::vector<detail::DeviceSpacePoint>& middleSPs,
    const std::vector<detail::DeviceSpacePoint>& topSPs,
    std::vector<std::vector<detail::SeedData>>& seeds);

}  // namespace Acts::Sycl