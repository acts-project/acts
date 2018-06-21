// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/SPForSeed.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/SeedmakerConfig.hpp"
#include "Acts/Seeding/SeedmakerState.hpp"
#include "Acts/Seeding/SpacePointConcept.hpp"

#include <list>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <array>
#include <memory>

namespace Acts {
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

/// Create a new space point grid, fill all space points from input parameter
/// into the grid and save grid in the state.
/// @param spVec the unordered vector containing all input space points
/// @param state the state of the object in which the space point grid will
/// be stored
      std::shared_ptr<Acts::SeedmakerState>
      initState
        ( std::vector<const Acts::concept::AnySpacePoint<>*> spVec) const;


/// Create all seed from the space points passed in createSpacePointGrid
/// Specify number of seeds
      void
      createSeeds(std::shared_ptr<Acts::SeedmakerState> state) const ;

    private:
              /**    @name Disallow default instantiation, copy, assignment */
      //@{
      New_Seedmaker()                                = delete;
      New_Seedmaker(const New_Seedmaker&)            = delete;
      New_Seedmaker &operator=(const New_Seedmaker&) = delete;
      //@}

      // Private methods
      void
      createSeedsInRegion (std::vector<std::shared_ptr<SPForSeed > > currentBin,
                            std::set<size_t > bottomBins,
                            std::set<size_t > topBins,
                            std::shared_ptr<Acts::SeedmakerState> state) const ;

      void transformCoordinates (std::vector<std::shared_ptr<SPForSeed> >& vec,
                                   std::shared_ptr<SPForSeed> spM,
                                   bool bottom,
                                   std::vector<LinCircle>& linCircle) const ;

        Acts::SeedmakerConfig m_config;
     };

  } // end of Acts namespace

///////////////////////////////////////////////////////////////////
// Object-function for curvature seeds comparison
///////////////////////////////////////////////////////////////////

class comCurvature  {
public:
  
  bool operator ()
  (const std::pair<float,Acts::SPForSeed*>& i1, 
   const std::pair<float,Acts::SPForSeed*>& i2)
  {
    return i1.first < i2.first;
  }
};

#include "Acts/Seeding/New_Seedmaker.ipp"
