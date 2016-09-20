// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SEEDING_BARRELSEEDFINDER_HPP
#define ACTS_SEEDING_BARRELSEEDFINDER_HPP

#include "ACTS/Seeding/SpacePoint.hpp"
#include "ACTS/Seeding/TrackSeed.hpp"
#include "ACTS/Seeding/detail/geometry.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {
namespace Seeding {

  struct HelixSeedConfig
  {
    double rangePhi1     = 0.2;  // search range in phi at layer 1
    double rangePhi2     = 0.2;  // search range in phi at layer 2
    double maxDeltaTheta = 0.1;  // cut on difference in theta between doublets
  };

  /// Find 3-point seeds assuming helix tracks with a combinatorial algorithm.
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
    for (const auto& p1 : barrel1.rangeDeltaPhi(p0.phi(), cfg.rangePhi1)) {
      Acts::Vector3D at2
          = detail::calcLineCircleIntersection(p0, p1, barrel2.radius);

      for (const auto& p2 : barrel2.rangeDeltaPhi(at2.phi(), cfg.rangePhi2)) {
        auto theta01 = (p1.position() - p0.position()).theta();
        auto theta12 = (p2.position() - p1.position()).theta();

        if (cfg.maxDeltaTheta < std::abs(theta12 - theta01)) continue;

        // use initial doublet to define direction
        // TODO add curvature corrections to angles
        auto phi01 = (p1.position() - p0.position()).phi();
        auto kappa = detail::calcCircleCurvature(p0, p1, p2);
        seeds.emplace_back(phi01, theta01, kappa, p0, p1, p2);
      }
    }
  }
}

#endif  // ACTS_SEEDING_BARRELSEEDFINDER_HPP
