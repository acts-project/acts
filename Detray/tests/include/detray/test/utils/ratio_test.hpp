// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Detray tests include(s).
#include "detray/geometry/shapes/cuboid3D.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"

// Detray test include(s)
#include "detray/test/framework/types.hpp"

// System include(s)
#include <atomic>
#include <random>
#include <vector>

namespace detray::test {

template <typename shape_t>
inline std::vector<point3> generate_regular_points(const std::size_t,
                                                   const std::vector<scalar> &);

/// This method creates a number of regular spaced points in 3D cartesian
/// coordinates within a cube
template <>
inline std::vector<point3> generate_regular_points<cuboid3D>(
    const std::size_t n_points, const std::vector<scalar> &world) {
  const std::size_t n_steps{static_cast<std::size_t>(std::cbrt(n_points))};

  const scalar world_x{world.at(0)};
  const scalar world_y{(world.size() > 1 ? world[1] : world[0])};
  const scalar world_z{(world.size() > 2 ? world[2] : world[0])};
  const scalar dx{world_x / static_cast<scalar>(n_steps)};
  const scalar dy{world_y / static_cast<scalar>(n_steps)};
  const scalar dz{world_z / static_cast<scalar>(n_steps)};

  std::vector<point3> points{};
  points.reserve(n_points);

  for (std::size_t ix = 0u; ix < n_steps; ++ix) {
    scalar x{-0.5f * world_x + static_cast<scalar>(ix) * dx};

    for (std::size_t iy = 0u; iy < n_steps; ++iy) {
      scalar y{-0.5f * world_y + static_cast<scalar>(iy) * dy};

      for (std::size_t iz = 0u; iz < n_steps; ++iz) {
        scalar z{-0.5f * world_z + static_cast<scalar>(iz) * dz};

        points.push_back({x, y, z});
      }
    }
  }

  return points;
}

/// This method creates a number of regular spaced points in 3D cartesian
/// coordinates positioned on a cylinder surface of radius @param r
template <>
inline std::vector<point3> generate_regular_points<cylinder2D>(
    const std::size_t n_points, const std::vector<scalar> &world) {
  const std::size_t n_steps{static_cast<std::size_t>(std::sqrt(n_points))};

  const scalar world_r{world.at(0)};
  const scalar world_z{world.at(1)};
  const scalar dphi{2.f * constant<scalar>::pi / static_cast<scalar>(n_steps)};
  const scalar dz{world_z / static_cast<scalar>(n_steps)};

  std::vector<point3> points{};
  points.reserve(n_points);

  for (std::size_t iphi = 0u; iphi < n_steps; ++iphi) {
    scalar phi{-constant<scalar>::pi + static_cast<scalar>(iphi) * dphi};

    for (std::size_t iz = 0u; iz < n_steps; ++iz) {
      scalar z{-0.5f * world_z + static_cast<scalar>(iz) * dz};

      points.push_back({world_r * std::cos(phi), world_r * std::sin(phi), z});
    }
  }

  return points;
}

/// This method creates a number of random points in 3D cartesian coordinates
inline std::vector<point3> generate_random_points(
    const scalar world, const std::size_t n_points,
    const std::uint32_t seed = std::mt19937::default_seed) {
  std::vector<point3> points{};
  points.reserve(n_points);

  std::random_device rd;
  std::mt19937 mt(rd());
  mt.seed(seed);
  std::uniform_real_distribution<scalar> dist(-0.5f * world, 0.5f * world);

  for (std::size_t i = 0u; i < n_points; ++i) {
    points.push_back({dist(mt), dist(mt), dist(mt)});
  }

  return points;
}

/// Calculate the ratio of points that pass the given predicate
template <typename predicate_t, typename... Args>
inline scalar ratio_test(const std::vector<point3> &points, Args... args) {
  std::atomic_ulong yes = 0u;

  for (const auto &p : points) {
    if (predicate_t{}(p, std::forward<Args>(args)...)) {
      ++yes;
    }
  }

  return static_cast<scalar>(yes) / static_cast<scalar>(points.size());
}

}  // namespace detray::test
