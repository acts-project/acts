// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/tracks/tracks.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/benchmark_base.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace detray::benchmarks {

/// Register a micro benchmark
///
/// @tparam benchmark_t the benchmark functor
/// @param bench_cfg basic benchmark configuration
/// @param suffix to the name of the benchmark
template <typename benchmark_t>
  requires std::derived_from<benchmark_t, benchmark_base>
inline void register_benchmark(
    const typename benchmark_t::configuration &bench_cfg,
    const std::string &suffix) {
  benchmark_t bench{bench_cfg};

  ::benchmark::RegisterBenchmark((bench.name() + suffix).c_str(), bench)
      ->UseRealTime()
      ->MeasureProcessCPUTime()
      ->ThreadPerCpu();
}

/// @returns the default track generation configuration for detray benchmarks
template <typename track_generator_t>
inline typename track_generator_t::configuration get_default_trk_gen_config(
    const std::vector<int> &n_tracks) {
  using track_t = typename track_generator_t::track_type;
  using scalar_t = dscalar<typename track_t::algebra_type>;

  int n_trks{*std::ranges::max_element(n_tracks)};

  // Generate tracks
  typename track_generator_t::configuration trk_cfg{};
  trk_cfg.n_tracks(static_cast<std::size_t>(n_trks));
  trk_cfg.randomize_charge(true);
  trk_cfg.phi_range(-constant<scalar_t>::pi, constant<scalar_t>::pi);
  trk_cfg.eta_range(-3.f, 3.f);
  trk_cfg.mom_range(1.f * unit<scalar_t>::GeV, 100.f * unit<scalar_t>::GeV);
  trk_cfg.origin(0.f, 0.f, 0.f);
  trk_cfg.origin_stddev(0.f, 0.f, 0.f);

  return trk_cfg;
}

/// Precompute the tracks
///
/// @param mr memory resource to allocate the track vector
/// @param cfg the configuration of the track generator
/// @param do_sort sort the tracks by theta angle
template <typename track_generator_t>
inline auto generate_tracks(vecmem::memory_resource *mr, track_generator_t &gen,
                            bool do_sort = true) {
  using track_t = typename track_generator_t::track_type;
  using scalar_t = dscalar<typename track_t::algebra_type>;

  // Track collection
  dvector<track_t> tracks(mr);

  // Iterate through uniformly distributed momentum directions
  for (auto track : gen) {
    // Put it into vector of trajectories
    tracks.push_back(track);
  }

  if (do_sort) {
    // Sort by theta angle
    const auto traj_comp = [](const auto &lhs, const auto &rhs) {
      constexpr auto pi_2{constant<scalar_t>::pi_2};
      return math::fabs(pi_2 - vector::theta(lhs.dir())) <
             math::fabs(pi_2 - vector::theta(rhs.dir()));
    };

    std::ranges::sort(tracks, traj_comp);
  }

  return tracks;
}

/// Generate as many samples of track states as there are entries in the
/// @param n_tracks vector using and externally provided track generator
/// @param gen
template <typename track_generator_t>
inline auto generate_track_samples(vecmem::memory_resource *mr,
                                   const std::vector<int> &n_tracks,
                                   track_generator_t &gen,
                                   bool do_sort = true) {
  using track_t = typename track_generator_t::track_type;

  std::vector<dvector<track_t>> track_samples{};
  track_samples.reserve(n_tracks.size());

  for (const int n : n_tracks) {
    gen.config().n_tracks(static_cast<std::size_t>(n));
    track_samples.push_back(generate_tracks(mr, gen, do_sort));
  }

  return track_samples;
}

/// Generate as many samples of track states as there are entries in the
/// @param n_tracks vector
template <typename track_generator_t>
inline auto generate_track_samples(
    vecmem::memory_resource *mr, const std::vector<int> &n_tracks,
    typename track_generator_t::configuration &cfg = {}, bool do_sort = true) {
  track_generator_t gen{cfg};
  return generate_track_samples(mr, n_tracks, gen, do_sort);
}

}  // namespace detray::benchmarks
