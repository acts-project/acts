// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s).
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Plugins/Sycl/Seeding/detail/Types.h"
#include "Acts/Plugins/Sycl/Seeding/DeviceExperimentCuts.hpp"

inline namespace cl{
  namespace sycl{
    class queue;
  }
};

namespace Acts::Sycl {

/// @brief Seedfinding algorithm implemented with SYCL to offload computations
/// to GPUs.

void offloadComputations( cl::sycl::queue* q,
                          const detail::deviceSeedfinderConfig& configData,
                          const DeviceExperimentCuts& deviceCuts,
                          const std::vector<detail::deviceSpacePoint>& bottomSPs,
                          const std::vector<detail::deviceSpacePoint>& middleSPs,
                          const std::vector<detail::deviceSpacePoint>& topSPs,
                          std::vector<std::vector<detail::SeedData>>& seeds);

/// @brief This function creates the SYCL queue object.
///
/// SYCL implementation details are hidden from this class, SYCL is only
/// linked to the one translation unit containing the actual implementation.
/// Because creating a queue is expensive, we only do it once.
///
/// @return A pointer to the queue.
cl::sycl::queue* createQueue(const std::string &);

template <typename external_spacepoint_t>
class Seedfinder {
  public:
  Seedfinder(Acts::SeedfinderConfig<external_spacepoint_t> config,
            Acts::Sycl::DeviceExperimentCuts cuts,
            const std::string &device_name_substring = "");

  ~Seedfinder() = default;
  Seedfinder() = delete;
  Seedfinder(const Seedfinder<external_spacepoint_t>&) = delete;
  Seedfinder<external_spacepoint_t>& operator=(
    const Seedfinder<external_spacepoint_t>&) = delete;

  /// Create all seeds from the space points in the three iterators.
  /// Can be used to parallelize the seed creation
  /// @param bottom group of space points to be used as innermost SP in a seed.
  /// @param middle group of space points to be used as middle SP in a seed.
  /// @param top group of space points to be used as outermost SP in a seed.
  /// Ranges must return pointers.
  /// Ranges must be separate objects for each parallel call.
  /// @return vector in which all found seeds for this group are stored.
  template <typename sp_range_t>
  std::vector<Seed<external_spacepoint_t> > createSeedsForGroup(
    sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const;

 private:

  Acts::SeedfinderConfig<external_spacepoint_t> m_config;
  Acts::Sycl::DeviceExperimentCuts m_deviceCuts;
  cl::sycl::queue* m_queue;
};

} // namespace Acts::Sycl

// Include the template implementation.
#include "Acts/Plugins/Sycl/Seeding/Seedfinder.ipp"