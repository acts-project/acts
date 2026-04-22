// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/tracks/tracks.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/benchmark_utils.hpp"

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

/// Which propagate function to run
enum class propagation_opt {
  e_unsync = 0,
  e_sync = 1,
};

/// Register a propagation benchmark of type @tparam benchmark_t
///
/// @tparam benchmark_t the propagation benchmark functor
/// @tparam propagator_t full propagator type
/// @tparam detector_t host detector type
/// @tparam bfield_t covfie magnetic field type
///
/// @param name name for the benchmark
/// @param bench_cfg basic benchmark configuration
/// @param prop_cfg propagation configuration
/// @param det the detector
/// @param bfield the covfie field
/// @param actor_states tuple that contains all actor states (same order as in
///                     actor_chain_t)
/// @param tracks the pre-computed test tracks
/// @param n_samples the number of track to run
template <template <typename, typename,
                    detray::benchmarks::propagation_opt> class benchmark_t,
          typename propagator_t, typename detector_t, typename bfield_bknd_t,
          detray::benchmarks::propagation_opt kOPT =
              detray::benchmarks::propagation_opt::e_unsync>
inline void register_benchmark(
    const std::string &name, benchmark_base::configuration &bench_cfg,
    propagation::config &prop_cfg, const detector_t &det, bfield_bknd_t &bfield,
    typename propagator_t::actor_chain_type::state_tuple *actor_states,
    std::vector<
        dvector<free_track_parameters<typename detector_t::algebra_type>>>
        &track_samples,
    const std::vector<int> &n_samples = {10000},
    vecmem::memory_resource *dev_mr = nullptr,
    const std::vector<int> &n_host_threads = {static_cast<int>(
        std::thread::hardware_concurrency())},
    int max_chunk_size = 1, int openmp_sched = 2) {
  using algebra_t = typename detector_t::algebra_type;
  using propagation_benchmark_t =
      benchmark_t<propagator_t, bfield_bknd_t, kOPT>;

  assert(track_samples.size() == n_samples.size());

  const std::size_t bench_range{
      math::max(n_samples.size(), n_host_threads.size())};
  for (std::size_t i = 0u; i < bench_range; ++i) {
    auto &tracks =
        track_samples.size() == 1u ? track_samples[0] : track_samples[i];
    int host_threads{n_host_threads.size() == 1u ? n_host_threads[0]
                                                 : n_host_threads[i]};
    const int n{n_samples.size() == 1u ? n_samples[0] : n_samples[i]};
    assert(static_cast<std::size_t>(n) <= tracks.size());

    bench_cfg.n_samples(n);

    typename propagation_benchmark_t::configuration prop_bm_cfg{bench_cfg};
    prop_bm_cfg.propagation() = prop_cfg;

    // Configure the benchmark
    propagation_benchmark_t prop_benchmark{prop_bm_cfg};

    std::string bench_name = prop_benchmark.config().name() + "_" + name + "_" +
                             std::to_string(n) + "_TRACKS";

    std::clog << bench_name << "\n" << bench_cfg;

    if constexpr (std::is_invocable_v<
                      decltype(prop_benchmark), ::benchmark::State &,
                      dvector<free_track_parameters<algebra_t>> *,
                      const detector_t *, const bfield_bknd_t *,
                      typename propagator_t::actor_chain_type::state_tuple *,
                      int, int, int>) {
      // Cpu benchmark
      ::benchmark::RegisterBenchmark(bench_name.c_str(), prop_benchmark,
                                     &tracks, &det, &bfield, actor_states,
                                     host_threads, max_chunk_size, openmp_sched)
          ->UseRealTime();
    } else {
      // Device benchmark
      ::benchmark::RegisterBenchmark(bench_name.c_str(), prop_benchmark, dev_mr,
                                     &tracks, &det, &bfield, actor_states)
          ->UseRealTime();
    }
  }
}

/// Register a propagation benchmark of type @tparam benchmark_t
///
/// @tparam benchmark_t the propagation benchmark functor
/// @tparam stepper_t the stepper to use for track parameter transport
/// @tparam actor_chain_t types of actors
template <template <typename, typename,
                    detray::benchmarks::propagation_opt> class benchmark_t,
          typename stepper_t, typename actor_chain_t, typename detector_t,
          typename bfield_bknd_t,
          detray::benchmarks::propagation_opt kOPT =
              detray::benchmarks::propagation_opt::e_unsync>
inline void register_benchmark(
    const std::string &name, benchmark_base::configuration &bench_cfg,
    propagation::config &prop_cfg, const detector_t &det, bfield_bknd_t &bfield,
    typename actor_chain_t::state_tuple *actor_states,
    std::vector<
        dvector<free_track_parameters<typename detector_t::algebra_type>>>
        &tracks,
    const std::vector<int> &n_samples = {10000},
    const std::vector<int> &n_host_threads = {static_cast<int>(
        std::thread::hardware_concurrency())},
    int max_chunk_size = 1, int openmp_sched = 2) {
  using propagator_t =
      propagator<stepper_t, caching_navigator<detector_t>, actor_chain_t>;
  register_benchmark<benchmark_t, propagator_t, detector_t, bfield_bknd_t,
                     kOPT>(name, bench_cfg, prop_cfg, det, bfield, actor_states,
                           tracks, n_samples, nullptr, n_host_threads,
                           max_chunk_size, openmp_sched);
}

}  // namespace detray::benchmarks
