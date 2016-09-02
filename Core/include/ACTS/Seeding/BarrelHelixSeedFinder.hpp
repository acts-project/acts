// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SEEDING_BARRELHELIXSEEDFINDER_HPP
#define ACTS_SEEDING_BARRELHELIXSEEDFINDER_HPP

#include "ACTS/Seeding/SpacePoint.hpp"
#include "ACTS/Seeding/TrackSeed.hpp"
#include "ACTS/Seeding/detail/geometry.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {
namespace Seeding {

  struct HelixSeedConfig
  {
    double deltaPhi01 = 0.2;
    double deltaPhi12 = 0.2;
    double maxLambda  = M_PI_2;
  };

  template <typename Identifier>
  void
  findHelixSeeds(const HelixSeedConfig&               cfg,
                 const BarrelSpacePoints<Identifier>& barrel0,
                 const BarrelSpacePoints<Identifier>& barrel1,
                 const BarrelSpacePoints<Identifier>& barrel2,
                 TrackSeeds3<Identifier>&             seeds);

}  // namespace Seeding
}  // namespace Acts

template <typename Identifier>
inline void
Acts::Seeding::findHelixSeeds(const HelixSeedConfig&               cfg,
                              const BarrelSpacePoints<Identifier>& barrel0,
                              const BarrelSpacePoints<Identifier>& barrel1,
                              const BarrelSpacePoints<Identifier>& barrel2,
                              TrackSeeds3<Identifier>&             seeds)
{
  for (const auto& p0 : barrel0.points) {
    for (const auto& p1 : barrel1.rangePhiDelta(p0.phi(), cfg.deltaPhi01)) {
      auto phi01 = (p1.position() - p0.position()).phi();
      auto theta01 = (p1.position() - p0.position()).theta();
      if (cfg.maxLambda < std::abs(M_PI_2 - theta01)) continue;

      for (const auto& p2 : barrel2.rangePhiDelta(p1.phi(), cfg.deltaPhi12)) {
        auto theta12 = (p2.position() - p1.position()).theta();
        if (cfg.maxLambda < std::abs(M_PI_2 - theta12)) continue;

        double kappa = detail::calcCircleCurvature(p0, p1, p2);
        // use first doublet to estimate direction
        seeds.emplace_back(phi01, theta01, kappa, p0, p1, p2);
      }
    }
  }
}

#endif  // ACTS_SEEDING_BARRELHELIXSEEDFINDER_HPP
