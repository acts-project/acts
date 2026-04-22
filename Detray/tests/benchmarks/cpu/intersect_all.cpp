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

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/navigation/detail/intersection_kernel.hpp"
#include "detray/navigation/intersection/intersection_config.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/track_generators.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// Google Benchmark include(s)
#include <benchmark/benchmark.h>

// System include(s)
#include <iostream>
#include <map>
#include <string>

// Use the detray:: namespace implicitly.
using namespace detray;

using bench_algebra = benchmarks::algebra;

using trk_generator_t =
    uniform_track_generator<free_track_parameters<bench_algebra>>;

constexpr unsigned int theta_steps{100u};
constexpr unsigned int phi_steps{100u};

// This test runs intersection with all surfaces of the TrackML detector
void BM_INTERSECT_ALL(benchmark::State &state) {
  // Detector configuration
  vecmem::host_memory_resource host_mr;
  toy_det_config<benchmarks::scalar> toy_cfg{};
  toy_cfg.n_edc_layers(7u);
  auto [d, names] = build_toy_detector<bench_algebra>(host_mr, toy_cfg);

  using detector_t = decltype(d);
  using scalar_t = dscalar<bench_algebra>;
  using sf_desc_t = typename detector_t::surface_type;

  detector_t::geometry_context geo_context;

  const auto &transforms = d.transform_store(geo_context);

  std::size_t hits{0u};
  std::size_t missed{0u};
  std::size_t n_surfaces{0u};
  benchmarks::point3 origin{0.f, 0.f, 0.f};
  std::vector<intersection2D<sf_desc_t, bench_algebra>> intersections{};

  // Iterate through uniformly distributed momentum directions
  auto trk_generator = trk_generator_t{};
  trk_generator.config()
      .theta_steps(theta_steps)
      .phi_steps(phi_steps)
      .origin(origin);

  // Intersector configuration
  intersection::config intr_cfg{.min_mask_tolerance = 1.f * unit<float>::um,
                                .max_mask_tolerance = 1.f * unit<float>::mm};
  constexpr const scalar_t external_mask_tol{0.f};

  for (auto _ : state) {
    for (const auto track : trk_generator) {
      // Loop over all surfaces in detector
      for (const sf_desc_t &sf_desc : d.surfaces()) {
        const auto sf = geometry::surface{d, sf_desc};
        sf.template visit_mask<
            detail::intersection_initialize<ray_intersector>>(
            intersections, detail::ray(track), sf_desc, transforms, geo_context,
            intr_cfg, external_mask_tol);

        ++n_surfaces;
      }
      benchmark::DoNotOptimize(hits);
      benchmark::DoNotOptimize(missed);

      hits += intersections.size();
      missed += n_surfaces - intersections.size();

      n_surfaces = 0u;
      intersections.clear();
    }
  }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
  DETRAY_INFO_HOST("[detray] hits / missed / total = "
                   << hits << " / " << missed << " / " << hits + missed);
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_ALL)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);
