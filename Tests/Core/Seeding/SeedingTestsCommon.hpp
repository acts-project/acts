// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>

#include "Acts/Seeding/SpacePoint.hpp"
#include "Acts/Seeding/TrackSeed.hpp"

/// Construct barrel layer w/ n points equidistant in phi.
///
/// Points go from 0.5 * delta to (n + 0.5) * delta
Acts::Seeding::BarrelSpacePoints<size_t>
makeBarrel(double radius, size_t nPoints)
{
  Acts::Seeding::BarrelSpacePoints<size_t> barrel(radius);

  for (size_t i = 0; i < nPoints; ++i) {
    // try to avoid floating point corner cases, e.g. 2*M_PI equivalent 0
    double phi = -M_PI + (2 * M_PI * (i + 0.5) / nPoints);
    barrel.points.emplace_back(
        radius * std::cos(phi), radius * std::sin(phi), 0, i);
  }
  barrel.sort();
  return barrel;
}

template <typename Identifier>
void
print(const std::vector<Acts::Seeding::SpacePoint<Identifier>>& points)
{
  for (const auto& point : points) {
    std::cout << point.identifier() << ' ' << point << '\n';
  }
}

template <typename Identifier>
void
print(const Acts::Seeding::BarrelSpacePoints<Identifier>& barrel)
{
  print(barrel.points);
}

template <typename Identifier, size_t N>
void
print(const std::vector<Acts::Seeding::TrackSeed<Identifier, N>>& seeds)
{
  for (const auto& seed : seeds) {
    std::cout << seed << '\n';
  }
}
