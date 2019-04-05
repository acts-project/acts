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
#include "Acts/Seeding/SeedmakerConfig.hpp"
#include "Acts/Seeding/SeedmakerState.hpp"

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
template <typename SpacePoint>
class New_Seedmaker
{
  ///////////////////////////////////////////////////////////////////
  // Public methods:
  ///////////////////////////////////////////////////////////////////

public:
  /// The only constructor. Requires a config object.
  /// @param config the configuration for the Seedmaker
  New_Seedmaker(const Acts::SeedmakerConfig<SpacePoint> config);
  ~New_Seedmaker() = default;
  /**    @name Disallow default instantiation, copy, assignment */
  //@{
  New_Seedmaker()                                 = delete;
  New_Seedmaker(const New_Seedmaker<SpacePoint>&) = delete;
  New_Seedmaker<SpacePoint>&
  operator=(const New_Seedmaker<SpacePoint>&)
      = delete;
  //@}

  /// converter function from internal seeds to external seeds. avoids
  /// templating all functions.
  /// @param intSeed the internal seed to be converted taken from a state
  /// @param inputSP the vector that has been used to call initState to generate
  /// intSeed
  /// @return the Seed return type with constructor Seed(SpacePoint, SpacePoint,
  /// Spacepoint, float z)
  template <typename Seed>
  std::unique_ptr<Seed>
  nextSeed(Acts::SeedmakerState<SpacePoint>* state);

  /// Create a new space point grid, fill all space points from input parameter
  /// into the grid and save grid in the state.
  /// @param spVec the unordered vector containing all input space points
  /// @param state the state of the object in which the space point grid will
  /// be stored
  template <typename SpacePointIterator>
  std::shared_ptr<Acts::SeedmakerState<SpacePoint>>
  initState(SpacePointIterator spBegin,
            SpacePointIterator spEnd,
            std::function<Acts::Vector2D(const SpacePoint* const,
                                         float,
                                         float,
                                         float)>          covTool,
            std::shared_ptr<Acts::IBinFinder<SpacePoint>> bottomBinFinder,
            std::shared_ptr<Acts::IBinFinder<SpacePoint>> topBinFinder) const;

  /// Create all seeds from the grid bin referenced by "it"
  void
  createSeedsForRegion(
      SeedingStateIterator<SpacePoint>                  it,
      std::shared_ptr<Acts::SeedmakerState<SpacePoint>> state) const;

private:
  void
  transformCoordinates(std::vector<const InternalSpacePoint<SpacePoint>*>& vec,
                       const InternalSpacePoint<SpacePoint>*               spM,
                       bool                    bottom,
                       std::vector<LinCircle>& linCircleVec) const;

  Acts::SeedmakerConfig<SpacePoint> m_config;
};

}  // end of Acts namespace

#include "Acts/Seeding/New_Seedmaker.ipp"
