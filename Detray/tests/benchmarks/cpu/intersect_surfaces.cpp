// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray core include(s).
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/tracks/ray.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/types.hpp"

// Detray test include(s)
#include "detray/test/common/track_generators.hpp"
#include "detray/test/utils/planes_along_direction.hpp"

// Google Benchmark include(s)
#include <benchmark/benchmark.h>

// Use the detray:: namespace implicitly.
using namespace detray;

using bench_algebra = benchmarks::algebra;

using ray_generator_t = uniform_track_generator<detail::ray<bench_algebra>>;

static const unsigned int theta_steps = 100u;
static const unsigned int phi_steps = 100u;

static const dvector<benchmarks::scalar> dists = {1.f, 2.f, 3.f, 4.f, 5.f,
                                                  6.f, 7.f, 8.f, 9.f, 10.f};

/// This benchmark runs intersection with the planar intersector
void BM_INTERSECT_PLANES(benchmark::State &state) {
  unsigned int sfhit = 0u;
  unsigned int sfmiss = 0u;

  auto [plane_descs, transforms] = test::planes_along_direction(
      dists, vector::normalize(benchmarks::vector3{1.f, 1.f, 1.f}));
  constexpr mask<rectangle2D, bench_algebra> rect{0u, 10.f, 20.f};

  // Iterate through uniformly distributed momentum directions
  auto ray_generator = ray_generator_t{};
  ray_generator.config().theta_steps(theta_steps).phi_steps(phi_steps);

  for (auto _ : state) {
    benchmark::DoNotOptimize(sfhit);
    benchmark::DoNotOptimize(sfmiss);

    // Iterate through uniformly distributed momentum directions
    for (const auto ray : ray_generator) {
      for (const auto &desc : plane_descs) {
        auto pi = ray_intersector<rectangle2D, bench_algebra>{};
        auto is = pi(ray, desc, rect, transforms[desc.transform()]);

        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);
        if (is.is_inside()) {
          ++sfhit;
        } else {
          ++sfmiss;
        }
        benchmark::ClobberMemory();
      }
    }
  }
}

BENCHMARK(BM_INTERSECT_PLANES)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

namespace {

enum class mask_id : unsigned int {
  e_rectangle2D = 0,
  e_cylinder2D = 1,
  e_conc_cylinder3 = 2,
};

enum class material_id : unsigned int {
  e_material_slab = 0,
};

}  // namespace

using mask_link_t = dtyped_index<mask_id, dindex>;
using material_link_t = dtyped_index<material_id, dindex>;

using surface_desc_t = surface_descriptor<mask_link_t, material_link_t>;
using intersection_t = intersection2D<surface_desc_t, bench_algebra>;

/// This benchmark runs intersection with the cylinder intersector
void BM_INTERSECT_CYLINDERS(benchmark::State &state) {
  using cylinder_mask = mask<cylinder2D, bench_algebra>;

  unsigned int sfhit = 0u;
  unsigned int sfmiss = 0u;
  dvector<cylinder_mask> cylinders;

  for (benchmarks::scalar r : dists) {
    cylinders.push_back(cylinder_mask{0u, r, -10.f, 10.f});
  }

  benchmarks::transform3 trf{};

  mask_link_t mask_link{mask_id::e_cylinder2D, 0};
  material_link_t material_link{material_id::e_material_slab, 0};
  surface_desc_t sf_desc(0u, mask_link, material_link, 0u,
                         surface_id::e_sensitive);

  // Iterate through uniformly distributed momentum directions
  auto ray_generator = ray_generator_t{};
  ray_generator.config().theta_steps(theta_steps).phi_steps(phi_steps);

  for (auto _ : state) {
    benchmark::DoNotOptimize(sfhit);
    benchmark::DoNotOptimize(sfmiss);

    // Iterate through uniformly distributed momentum directions
    for (const auto ray : ray_generator) {
      for (const auto &cylinder : cylinders) {
        auto ci = ray_intersector<cylinder2D, bench_algebra>{};
        auto inters = ci(ray, sf_desc, cylinder, trf);

        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);
        for (const auto &sfi : inters) {
          if (sfi.is_inside()) {
            ++sfhit;
          } else {
            ++sfmiss;
          }
        }
        benchmark::ClobberMemory();
      }
    }
  }
}

BENCHMARK(BM_INTERSECT_CYLINDERS)
#ifdef DETRAY_BENCHMARKS_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with a specialized portal cylinder
/// intersector
void BM_INTERSECT_PORTAL_CYLINDERS(benchmark::State &state) {
  using cylinder_mask = mask<concentric_cylinder2D, bench_algebra>;

  unsigned int sfhit = 0u;
  unsigned int sfmiss = 0u;
  dvector<cylinder_mask> cylinders;

  for (benchmarks::scalar r : dists) {
    cylinders.push_back(cylinder_mask{0u, r, -10.f, 10.f});
  }

  benchmarks::transform3 trf{};

  mask_link_t mask_link{mask_id::e_cylinder2D, 0u};
  material_link_t material_link{material_id::e_material_slab, 0u};
  surface_desc_t sf_desc(0u, mask_link, material_link, 0u,
                         surface_id::e_sensitive);

  // Iterate through uniformly distributed momentum directions
  auto ray_generator = ray_generator_t{};
  ray_generator.config().theta_steps(theta_steps).phi_steps(phi_steps);

  for (auto _ : state) {
    benchmark::DoNotOptimize(sfhit);
    benchmark::DoNotOptimize(sfmiss);

    // Iterate through uniformly distributed momentum directions
    for (const auto ray : ray_generator) {
      for (const auto &cylinder : cylinders) {
        auto cpi = ray_intersector<concentric_cylinder2D, bench_algebra>{};
        auto is = cpi(ray, sf_desc, cylinder, trf);

        benchmark::DoNotOptimize(sfhit);
        benchmark::DoNotOptimize(sfmiss);
        if (is.is_inside()) {
          ++sfhit;
        } else {
          ++sfmiss;
        }
        benchmark::ClobberMemory();
      }
    }
  }
}

BENCHMARK(BM_INTERSECT_PORTAL_CYLINDERS)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);
