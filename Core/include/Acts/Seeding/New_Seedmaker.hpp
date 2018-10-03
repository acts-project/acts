// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/SeedmakerConfig.hpp"
#include "Acts/Seeding/SeedmakerState.hpp"

#include <list>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <array>
#include <memory>

namespace Acts {
  struct LinCircle{
    float Zo;
    float cotTheta;
    float iDeltaR;
    float Er;
    float U;
    float V;
  };
  class New_Seedmaker 
    {
      ///////////////////////////////////////////////////////////////////
      // Public methods:
      ///////////////////////////////////////////////////////////////////
      
    public:
/// The only constructor. Requires a config object.
/// @param config the configuration for the Seedmaker
      New_Seedmaker(const Acts::SeedmakerConfig& config);
      ~New_Seedmaker() = default;

/// converter function from internal seeds to external seeds. avoids templating all functions.
/// @param intSeed the internal seed to be converted taken from a state
/// @param inputSP the vector that has been used to call initState to generate intSeed
/// @return the Seed return type with constructor Seed(SpacePoint, SpacePoint, Spacepoint, float z)
      template <typename Seed, typename SpacePoint>
      std::unique_ptr<Seed> nextSeed(Acts::SeedmakerState* state , std::vector<const SpacePoint*>& inputSP);

/// Create a new space point grid, fill all space points from input parameter
/// into the grid and save grid in the state.
/// @param spVec the unordered vector containing all input space points
/// @param state the state of the object in which the space point grid will
/// be stored
      template <typename SpacePoint>
      std::shared_ptr<Acts::SeedmakerState>
      initState
        (const std::vector<const SpacePoint*>& spVec,
         std::function<Acts::Vector2D(const SpacePoint*,float,float,float)> covTool,
         std::shared_ptr<Acts::IBinFinder> bottomBinFinder,
         std::shared_ptr<Acts::IBinFinder> topBinFinder) const;


/// Create all seed from the space point referenced by it
      void
      createSeedsForSP ( SeedingStateIterator it,
                         std::shared_ptr<Acts::SeedmakerState> state) const;

    private:
              /**    @name Disallow default instantiation, copy, assignment */
      //@{
      New_Seedmaker()                                = delete;
      New_Seedmaker(const New_Seedmaker&)            = delete;
      New_Seedmaker &operator=(const New_Seedmaker&) = delete;
      //@}

      // Private methods
//      void
//      createSeedsInRegion (std::vector<std::unique_ptr<const InternalSpacePoint > >& currentBin,
//                            std::set<size_t > bottomBins,
//                            std::set<size_t > topBins,
//                            std::shared_ptr<Acts::SeedmakerState> state) const ;

      void transformCoordinates (std::vector< const InternalSpacePoint* >& vec,
                                   const InternalSpacePoint* spM,
                                   bool bottom,
                                   std::vector<LinCircle>& linCircle) const ;

        Acts::SeedmakerConfig m_config;
     };


  } // end of Acts namespace

#include "Acts/Seeding/New_Seedmaker.ipp"
