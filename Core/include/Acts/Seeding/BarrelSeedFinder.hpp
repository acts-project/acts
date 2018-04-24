// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include "Acts/Seeding/SpacePoint.hpp"
#include "Acts/Seeding/TrackSeed.hpp"
#include "Acts/Seeding/detail/geometry.hpp"

namespace Acts {
namespace Seeding {

  struct HelixSeedConfig
  {
    double rangePhi1     = 0.2;  // search range in phi at layer 1
    double rangePhi2     = 0.2;  // search range in phi at layer 2
    double maxDeltaTheta = 0.1;  // cut on difference in theta between doublets
  };

  /// Find 3-point seeds with a combinatorial algorithm.
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
      Vector3D d01     = p1.position() - p0.position();
      double   theta01 = d01.theta();
      // Acts::Vector3D at2
      //     = detail::calcLineCircleIntersection(p0, d01, barrel2.radius);
      for (const auto& p2 : barrel2.rangeDeltaPhi(p1.phi(), cfg.rangePhi2)) {
        Vector3D d12     = p2.position() - p1.position();
        double   theta12 = d12.theta();

        if (cfg.maxDeltaTheta < std::abs(theta12 - theta01)) continue;

        double kappa = detail::calcCircleCurvature(d01, d12);
        // initial direction correction due to curvature, use
        //   chord = 2 * radius * sin(propagation angle / 2)
        // and assume sin(x) = x
        double phi01 = d01.phi() - d01.head<2>().norm() * kappa / 2;
        // track parameters defined at the first space point
        seeds.emplace_back(phi01, theta01, kappa, p0, p1, p2);
      }
    }
  }
}