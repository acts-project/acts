// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Plugins/Cuda/Cuda.hpp"
#include "Acts/Plugins/Cuda/Seeding/Kernels.cuh"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"

#include <array>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Acts {

template <typename external_spacepoint_t>
class Seedfinder<external_spacepoint_t, Acts::Cuda> {
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////

 public:
  Seedfinder(Acts::SeedfinderConfig<external_spacepoint_t> config);

  ~Seedfinder() = default;
  /**    @name Disallow default instantiation, copy, assignment */
  //@{
  Seedfinder() = delete;
  Seedfinder(const Seedfinder<external_spacepoint_t, Acts::Cuda>&) = delete;
  Seedfinder<external_spacepoint_t, Acts::Cuda>& operator=(
      const Seedfinder<external_spacepoint_t, Acts::Cuda>&) = delete;
  //@}

  /// Create all seeds from the space points in the three iterators.
  /// Can be used to parallelize the seed creation
  /// @param bottomSPs group of space points to be used as innermost SP in a seed.
  /// @param middleSPs group of space points to be used as middle SP in a seed.
  /// @param topSPs group of space points to be used as outermost SP in a seed.
  /// Ranges must return pointers.
  /// Ranges must be separate objects for each parallel call.
  /// @return vector in which all found seeds for this group are stored.
  template <typename sp_range_t>
  std::vector<Seed<external_spacepoint_t> > createSeedsForGroup(
      sp_range_t bottomSPs, sp_range_t middleSPs, sp_range_t topSPs) const;

 private:
  Acts::SeedfinderConfig<external_spacepoint_t> m_config;
};

}  // namespace Acts

#include "Acts/Plugins/Cuda/Seeding/Seedfinder.ipp"
