// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Algebra include(s).
// #include "detray/plugins/algebra/array_definitions.hpp"
// #include "detray/plugins/algebra/eigen_definitions.hpp"
#include "algebra/vc_aos.hpp"
#include "algebra/vc_soa.hpp"

// Detray core include(s).
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/navigation/intersection/ray_intersector.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/logging.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/types.hpp"

// Detray test include(s)
#include "detray/test/common/track_generators.hpp"
#include "detray/test/utils/planes_along_direction.hpp"

// Google Benchmark include(s)
#include <benchmark/benchmark.h>

// System include(s)
#include <algorithm>

using namespace detray;

static constexpr unsigned int theta_steps{2000u};
static constexpr unsigned int phi_steps{2000u};
static constexpr unsigned int n_surfaces{16u};

/// Linear algebra implementation using SoA memory layout
using algebra_v = detray::vc_soa<benchmarks::scalar>;

/// Linear algebra implementation using AoS memory layout
// using algebra_array = detray::array<benchmarks::scalar>;
using algebra_vc_aos = detray::vc_aos<benchmarks::scalar>;
// using algebra_eigen = detray::eigen<benchmarks::scalar>;

using algebra_s = algebra_vc_aos;

// Size of an SoA batch
constexpr std::size_t simd_size{dscalar<algebra_v>::size()};

using ray_t = detail::ray<algebra_s>;

namespace {

enum class mask_id : unsigned int {
  e_rectangle2D = 0,
  e_cylinder2D = 1,
  e_conc_cylinder3 = 2,
};

enum class material_id : unsigned int {
  e_material_slab = 0,
};

// Helper type definitions.
using mask_link_t = dtyped_index<mask_id, dindex>;
using material_link_t = dtyped_index<material_id, dindex>;

using surface_desc_t = surface_descriptor<mask_link_t, material_link_t>;

/// Generate a number of test rays
std::vector<ray_t> generate_rays() {
  using ray_generator_t = uniform_track_generator<ray_t>;

  // Iterate through uniformly distributed momentum directions
  auto ray_generator = ray_generator_t{};
  ray_generator.config().theta_steps(theta_steps).phi_steps(phi_steps);

  std::vector<ray_t> rays;
  std::ranges::copy(ray_generator, std::back_inserter(rays));

  return rays;
}

/// Generate the translation distances to place the surfaces
template <concepts::algebra algebra_t>
dvector<dscalar<algebra_t>> get_dists(std::size_t n) {
  using scalar_t = dscalar<algebra_t>;

  dvector<benchmarks::scalar> dists;

  for (std::size_t i = 1u; i <= n; ++i) {
    dists.push_back(static_cast<scalar_t>(i));
  }

  return dists;
}

/// Specialization for hthe SOA memory layout (need n/simd_size samples)
template <>
dvector<dscalar<algebra_v>> get_dists<algebra_v>(std::size_t n) {
  using scalar_t = dscalar<algebra_v>;
  using value_t = typename algebra_v::value_type;

  dvector<scalar_t> dists;
  dists.resize(static_cast<std::size_t>(std::ceil(n / simd_size)));
  for (std::size_t i = 0u; i < dists.size(); ++i) {
    dists[i] = scalar_t::IndexesFromZero() +
               scalar_t(static_cast<value_t>(i)) * simd_size + scalar_t(1.f);
  }

  return dists;
}

}  // namespace

/// This benchmark runs intersection with the planar intersector
void BM_INTERSECT_PLANES_AOS(benchmark::State& state) {
  using mask_t = mask<rectangle2D, algebra_s, std::uint_least16_t>;

  auto dists = get_dists<algebra_s>(n_surfaces);
  auto [plane_descs, transforms] = test::planes_along_direction<algebra_s>(
      dists, dvector3D<algebra_s>{1.f, 1.f, 1.f});

  constexpr mask_t rect{0u, 100.f, 200.f};
  std::vector<mask_t> masks(plane_descs.size(), rect);

  const auto rays = generate_rays();
  const auto pi = ray_intersector<rectangle2D, algebra_s>{};

#ifdef DETRAY_BENCHMARK_PRINTOUTS
  std::size_t hit{0u};
  std::size_t miss{0u};
#endif

  for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    hit = 0u;
    miss = 0u;
#endif

    // Iterate through uniformly distributed momentum directions
    for (const auto& ray : rays) {
      for (std::size_t i = 0u; i < plane_descs.size(); ++i) {
        auto is = pi(ray, plane_descs[i], masks[i], transforms[i]);

#ifdef DETRAY_BENCHMARK_PRINTOUTS
        if (is.is_inside()) {
          ++hit;
        } else {
          ++miss;
        }
#endif

        benchmark::DoNotOptimize(is);
      }
    }
  }
#ifdef DETRAY_BENCHMARK_PRINTOUTS
  DETRAY_INFO_HOST(mask_t::shape::name
                   << " AoS: hit/miss ... " << hit << " / " << miss
                   << " (total: " << rays.size() * masks.size() << ")");
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_PLANES_AOS)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the planar intersector
void BM_INTERSECT_PLANES_SOA(benchmark::State& state) {
  using mask_t = mask<rectangle2D, algebra_v, std::uint_least16_t>;
  using vector3_t = dvector3D<algebra_v>;

  auto dists = get_dists<algebra_v>(n_surfaces);
  auto [plane_descs, transforms] =
      test::planes_along_direction<algebra_v>(dists, vector3_t{1.f, 1.f, 1.f});

  std::vector<mask_t> masks{};
  for (std::size_t i = 0u; i < plane_descs.size(); ++i) {
    masks.emplace_back(0u, 100.f, 200.f);
  }

  const auto rays = generate_rays();
  const auto pi = ray_intersector<rectangle2D, algebra_v>{};

#ifdef DETRAY_BENCHMARK_PRINTOUTS
  std::size_t hit{0u};
  std::size_t miss{0u};
#endif

  for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    hit = 0u;
    miss = 0u;
#endif

    // Iterate through uniformly distributed momentum directions
    for (const auto& ray : rays) {
      for (std::size_t i = 0u; i < plane_descs.size(); ++i) {
        auto is = pi(ray, plane_descs[i], masks[i], transforms[i]);

        benchmark::DoNotOptimize(is);

#ifdef DETRAY_BENCHMARK_PRINTOUTS
        hit += is.is_inside().count();
        miss += simd_size - is.is_inside().count();
#endif
      }
    }
  }

#ifdef DETRAY_BENCHMARK_PRINTOUTS
  DETRAY_INFO_HOST(mask_t::shape::name
                   << " SoA: hit/miss ... " << hit << " / " << miss
                   << " (total: " << rays.size() * masks.size() * simd_size
                   << ")");
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_PLANES_SOA)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the cylinder intersector
void BM_INTERSECT_CYLINDERS_AOS(benchmark::State& state) {
  using transform3_t = dtransform3D<algebra_s>;
  using scalar_t = dscalar<algebra_s>;

  using mask_t = mask<cylinder2D, algebra_s, std::uint_least16_t>;

  std::vector<mask_t> masks;
  for (const scalar_t r : get_dists<algebra_s>(n_surfaces)) {
    masks.emplace_back(0u, r, -100.f, 100.f);
  }

  transform3_t trf{};

  mask_link_t mask_link{mask_id::e_conc_cylinder3, 0u};
  material_link_t material_link{material_id::e_material_slab, 0u};
  surface_desc_t cyl_desc(0u, mask_link, material_link, 0u,
                          surface_id::e_sensitive);

  // Iterate through uniformly distributed momentum directions
  const auto rays = generate_rays();
  const auto cci = ray_intersector<cylinder2D, algebra_s>{};
#ifdef DETRAY_BENCHMARK_PRINTOUTS
  std::size_t hit{0u};
  std::size_t miss{0u};
#endif
  for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    hit = 0u;
    miss = 0u;
#endif

    // Iterate through uniformly distributed momentum directions
    for (const auto& ray : rays) {
      for (const auto& cylinder : masks) {
        auto is = cci(ray, cyl_desc, cylinder, trf);

        static_assert(is.size() == 2u, "Wrong number of solutions");
#ifdef DETRAY_BENCHMARK_PRINTOUTS
        for (const auto& i : is) {
          if (i.is_inside()) {
            ++hit;
          } else {
            ++miss;
          }
        }
#endif
        benchmark::DoNotOptimize(is);
      }
    }
  }
#ifdef DETRAY_BENCHMARK_PRINTOUTS
  DETRAY_INFO_HOST(mask_t::shape::name
                   << " AoS: hit/miss ... " << hit << " / " << miss
                   << " (total: " << rays.size() * masks.size() << ")");
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_CYLINDERS_AOS)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the cylinder intersector
void BM_INTERSECT_CYLINDERS_SOA(benchmark::State& state) {
  using transform3_t = dtransform3D<algebra_v>;
  using scalar_t = dscalar<algebra_v>;

  using mask_t = mask<cylinder2D, algebra_v, std::uint_least16_t>;

  std::vector<mask_t> masks;
  for (const scalar_t r : get_dists<algebra_v>(n_surfaces)) {
    masks.emplace_back(0u, r, -100.f, 100.f);
  }

  transform3_t trf{};

  mask_link_t mask_link{mask_id::e_conc_cylinder3, 0u};
  material_link_t material_link{material_id::e_material_slab, 0u};
  surface_desc_t cyl_desc(0u, mask_link, material_link, 0u,
                          surface_id::e_sensitive);

  // Iterate through uniformly distributed momentum directions
  const auto rays = generate_rays();
  const auto cci = ray_intersector<cylinder2D, algebra_v>{};
#ifdef DETRAY_BENCHMARK_PRINTOUTS
  std::size_t hit{0u};
  std::size_t miss{0u};
#endif
  for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    hit = 0u;
    miss = 0u;
#endif

    // Iterate through uniformly distributed momentum directions
    for (const auto& ray : rays) {
      for (const auto& cylinder : masks) {
        auto is = cci(ray, cyl_desc, cylinder, trf);

        static_assert(is.size() == 2u, "Wrong number of solutions");
#ifdef DETRAY_BENCHMARK_PRINTOUTS
        for (const auto& i : is) {
          hit += i.is_inside().count();
          miss += simd_size - i.is_inside().count();
        }
#endif
        benchmark::DoNotOptimize(is);
      }
    }
  }
#ifdef DETRAY_BENCHMARK_PRINTOUTS
  DETRAY_INFO_HOST(mask_t::shape::name
                   << " SoA: hit/miss ... " << hit << " / " << miss
                   << " (total: " << rays.size() * masks.size() * simd_size
                   << ")");
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_CYLINDERS_SOA)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the concentric cylinder intersector
void BM_INTERSECT_CONCENTRIC_CYLINDERS_AOS(benchmark::State& state) {
  using transform3_t = dtransform3D<algebra_s>;
  using scalar_t = dscalar<algebra_s>;

  using mask_t = mask<concentric_cylinder2D, algebra_s, std::uint_least16_t>;

  std::vector<mask_t> masks;
  for (const scalar_t r : get_dists<algebra_s>(n_surfaces)) {
    masks.emplace_back(0u, r, -100.f, 100.f);
  }

  transform3_t trf{};

  mask_link_t mask_link{mask_id::e_conc_cylinder3, 0u};
  material_link_t material_link{material_id::e_material_slab, 0u};
  surface_desc_t cyl_desc(0u, mask_link, material_link, 0u,
                          surface_id::e_sensitive);

  // Iterate through uniformly distributed momentum directions
  const auto rays = generate_rays();
  const auto cci = ray_intersector<concentric_cylinder2D, algebra_s>{};
#ifdef DETRAY_BENCHMARK_PRINTOUTS
  std::size_t hit{0u};
  std::size_t miss{0u};
#endif

  for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    hit = 0u;
    miss = 0u;
#endif
    // Iterate through uniformly distributed momentum directions
    for (const auto& ray : rays) {
      for (const auto& cylinder : masks) {
        auto is = cci(ray, cyl_desc, cylinder, trf);
#ifdef DETRAY_BENCHMARK_PRINTOUTS
        if (is.is_inside()) {
          ++hit;
        } else {
          ++miss;
        }
#endif
        benchmark::DoNotOptimize(is);
      }
    }
  }
#ifdef DETRAY_BENCHMARK_PRINTOUTS
  DETRAY_INFO_HOST(mask_t::shape::name
                   << " AoS: hit/miss ... " << hit << " / " << miss
                   << " (total: " << rays.size() * masks.size() << ")");
#endif  // DETRAY_BENCHMARK_PRINTOUTS
}

BENCHMARK(BM_INTERSECT_CONCENTRIC_CYLINDERS_AOS)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);

/// This benchmark runs intersection with the concentric cylinder intersector
void BM_INTERSECT_CONCENTRIC_CYLINDERS_SOA(benchmark::State& state) {
  using transform3_t = dtransform3D<algebra_v>;
  using scalar_t = dscalar<algebra_v>;

  using mask_t = mask<concentric_cylinder2D, algebra_v, std::uint_least16_t>;

  std::vector<mask_t> masks;
  for (const scalar_t r : get_dists<algebra_v>(n_surfaces)) {
    masks.emplace_back(0u, r, -100.f, 100.f);
  }

  transform3_t trf{};

  mask_link_t mask_link{mask_id::e_conc_cylinder3, 0u};
  material_link_t material_link{material_id::e_material_slab, 0u};
  surface_desc_t cyl_desc(0u, mask_link, material_link, 0u,
                          surface_id::e_sensitive);

  // Iterate through uniformly distributed momentum directions
  const auto rays = generate_rays();
  const auto cci = ray_intersector<concentric_cylinder2D, algebra_v>{};
#ifdef DETRAY_BENCHMARK_PRINTOUTS
  std::size_t hit{0u};
  std::size_t miss{0u};
#endif
  for (auto _ : state) {
#ifdef DETRAY_BENCHMARK_PRINTOUTS
    hit = 0u;
    miss = 0u;
#endif
    // Iterate through uniformly distributed momentum directions
    for (const auto& ray : rays) {
      for (const auto& cylinder : masks) {
        auto is = cci(ray, cyl_desc, cylinder, trf);
#ifdef DETRAY_BENCHMARK_PRINTOUTS
        hit += is.is_inside().count();
        miss += simd_size - is.is_inside().count();
#endif
        benchmark::DoNotOptimize(is);
      }
    }
  }
#ifdef DETRAY_BENCHMARK_PRINTOUTS
  DETRAY_INFO_HOST(mask_t::shape::name
                   << " SoA: hit/miss ... " << hit << " / " << miss
                   << " (total: " << rays.size() * masks.size() * simd_size
                   << ")");
#endif
}

BENCHMARK(BM_INTERSECT_CONCENTRIC_CYLINDERS_SOA)
#ifdef DETRAY_BENCHMARK_MULTITHREAD
    ->ThreadRange(1, benchmark::CPUInfo::Get().num_cpus)
#endif
    ->Unit(benchmark::kMillisecond);
