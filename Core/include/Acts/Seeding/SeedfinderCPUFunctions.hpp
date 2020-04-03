// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"

namespace Acts{

  template< typename external_spacepoint_t, typename sp_range_t >
  class SeedfinderCPUFunctions {

  public: 
    
    static std::vector<const InternalSpacePoint<external_spacepoint_t>*>
    searchDoublet(bool isBottom, sp_range_t& SPs,
		  const InternalSpacePoint<external_spacepoint_t>& spM,
		  const SeedfinderConfig<external_spacepoint_t>& config);

    static void transformCoordinates(std::vector<const InternalSpacePoint<external_spacepoint_t>*>& vec,
				     const InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
				     std::vector<LinCircle>& linCircleVec);


    static std::vector<std::pair< float, std::unique_ptr<const InternalSeed<external_spacepoint_t>>>>
    searchTriplet(const InternalSpacePoint<external_spacepoint_t>& spM,
		  const std::vector<const InternalSpacePoint<external_spacepoint_t>*>& compatBottomSP,
		  const std::vector<const InternalSpacePoint<external_spacepoint_t>*>& compatTopSP,
		  const std::vector<LinCircle>& linCircleBottom,
		  const std::vector<LinCircle>& linCircleTop,
		  const SeedfinderConfig<external_spacepoint_t>& config);

    
  private:
    
  };
  
} // namespace Acts

#include "Acts/Seeding/SeedfinderCPUFunctions.ipp"
