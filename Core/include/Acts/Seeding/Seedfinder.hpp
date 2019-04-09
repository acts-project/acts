// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/SeedfinderConfig.hpp"
#include "Acts/Seeding/SeedfinderState.hpp"

#include <array>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace Acts {
struct LinCircle
{
  float Zo;
  float cotTheta;
  float iDeltaR;
  float Er;
  float U;
  float V;
};
template <typename external_spacepoint_t>
class Seedfinder
{
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////

public:
  /// The only constructor. Requires a config object.
  /// @param config the configuration for the Seedfinder
  Seedfinder(const Acts::SeedfinderConfig<external_spacepoint_t> config);
  ~Seedfinder() = default;
  /**    @name Disallow default instantiation, copy, assignment */
  //@{
  Seedfinder()                                 = delete;
  Seedfinder(const Seedfinder<external_spacepoint_t>&) = delete;
  Seedfinder<external_spacepoint_t>&
  operator=(const Seedfinder<external_spacepoint_t>&)
      = delete;
  //@}

  /// @brief Create a new space point grid, fill all space points from input parameter
  /// into the grid and save grid in the state.
  /// @param spBegin begin iterator to retrieve all input space points
  /// @param spEnd end iterator where to stop adding space points
  /// @param covTool covariance tool which is applied to each space point
  /// @param bottomBinFinder IBinFinder tool to store in the returned state
  /// @param topBinFinder IBinFinder tool to store in the returned state
  /// @return the state object containing all space points and container for
  /// found seeds
  template <typename spacepoint_iterator_t>
  Acts::SeedfinderState<external_spacepoint_t>
  initState(spacepoint_iterator_t spBegin,
            spacepoint_iterator_t spEnd,
            std::function<Acts::Vector2D(const external_spacepoint_t&,
                                         float,
                                         float,
                                         float)>          covTool,
            std::shared_ptr<Acts::IBinFinder<external_spacepoint_t>> bottomBinFinder,
            std::shared_ptr<Acts::IBinFinder<external_spacepoint_t>> topBinFinder) const;

  /// Create all seeds from the grid bin referenced by "it"
  /// This method can be called in parallel.
  /// @param it caches and updates the current space point and its
  /// neighbors. Iterator must be separate object for each parallel call.
  /// @param state the state object in which all found seeds are stored.
  /// state object can be shared between all parallel calls
  void
  createSeedsForRegion(SeedfinderStateIterator<external_spacepoint_t>  it,
                       Acts::SeedfinderState<external_spacepoint_t>& state) const;

private:
  void
  transformCoordinates(std::vector<const InternalSpacePoint<external_spacepoint_t>*>& vec,
                       const InternalSpacePoint<external_spacepoint_t>&               spM,
                       bool                    bottom,
                       std::vector<LinCircle>& linCircleVec) const;

  Acts::SeedfinderConfig<external_spacepoint_t> m_config;
};

}  // end of Acts namespace

#include "Acts/Seeding/Seedfinder.ipp"
