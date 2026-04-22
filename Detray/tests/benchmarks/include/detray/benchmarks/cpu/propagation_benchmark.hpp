// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/benchmark_base.hpp"
#include "detray/benchmarks/propagation_benchmark.hpp"
#include "detray/benchmarks/propagation_benchmark_utils.hpp"

// Benchmark include
#include <benchmark/benchmark.h>

#ifdef _OPENMP
// openMP include
#include <omp.h>
#endif

// System include(s)
#include <algorithm>
#include <cassert>
#include <ranges>
#include <string>

namespace detray::benchmarks {

/// Host propagation benchmarks
template <typename propagator_t, typename bfield_t,
          detray::benchmarks::propagation_opt kOPT =
              detray::benchmarks::propagation_opt::e_unsync>
struct host_propagation_bm
    : public propagation_benchmark<
          typename propagator_t::detector_type::algebra_type> {
  /// Detector dependent types
  using algebra_t = typename propagator_t::detector_type::algebra_type;
  using base_type = propagation_benchmark<algebra_t>;

  using base_type::base_type;

  /// Prepare data and run benchmark loop
  inline void operator()(
      ::benchmark::State &state,
      const dvector<free_track_parameters<algebra_t>> *tracks,
      const typename propagator_t::detector_type *det, const bfield_t *bfield,
      const typename propagator_t::actor_chain_type::state_tuple
          *input_actor_states,
      [[maybe_unused]] const int n_threads,
      [[maybe_unused]] const int max_chunk_size,
      [[maybe_unused]] const int thread_schedule) const {
    using actor_chain_t = typename propagator_t::actor_chain_type;
    using actor_states_t = typename actor_chain_t::state_tuple;

    assert(tracks != nullptr);
    assert(det != nullptr);
    assert(bfield != nullptr);
    assert(input_actor_states != nullptr);

    const int n_samples{this->config().n_samples()};
    const int n_warmup{this->config().n_warmup()};

    assert(static_cast<std::size_t>(n_samples) <= tracks->size());

#ifdef _OPENMP
    // Set the number of threads for the openMP parallel regions
    omp_set_num_threads(n_threads);
    // Clamp chunk size to [1, max_chunk_size]
    int chunk_size{
        math::min(static_cast<int>(n_samples / n_threads), max_chunk_size)};
    chunk_size = math::max(chunk_size, 1);
    omp_set_schedule(static_cast<omp_sched_t>(thread_schedule), chunk_size);

    DETRAY_VERBOSE_HOST("No. tracks " << n_samples);
    DETRAY_VERBOSE_HOST("No. threads " << n_threads);
    DETRAY_VERBOSE_HOST("Schedule type " << thread_schedule);
    DETRAY_VERBOSE_HOST("Chunk size " << chunk_size);
#endif

    // Create propagator
    propagator_t p{this->propagation()};

    // Call the host propagation
    auto run_propagation = [&p, det, bfield, input_actor_states](
                               const free_track_parameters<algebra_t> &track) {
      // Fresh copy of actor states
      actor_states_t actor_states(*input_actor_states);
      // Tuple of references to pass to the propagator
      typename actor_chain_t::state_ref_tuple actor_state_refs =
          actor_chain_t::setup_actor_states(actor_states);

      typename propagator_t::stepper_type::magnetic_field_type bfield_view(
          *bfield);
      typename propagator_t::state p_state(track, bfield_view, *det);
      // Particle hypothesis
      auto &ptc = p_state.stepping().particle_hypothesis();
      p_state.set_particle(update_particle_hypothesis(ptc, track));

      // Run propagation
      if constexpr (kOPT == detray::benchmarks::propagation_opt::e_unsync) {
        ::benchmark::DoNotOptimize(p.propagate(p_state, actor_state_refs));
      } else if constexpr (kOPT ==
                           detray::benchmarks::propagation_opt::e_sync) {
        /* Do nothing for now */
      }
      assert(p.finished(p_state));
    };

    // Warm-up
    if (this->config().do_warmup()) {
      assert(n_warmup > 0);
      int stride{n_samples / n_warmup};
      stride = (stride == 0) ? 10 : stride;
      assert(stride > 0);

#pragma omp parallel for
      for (int i = 0; i < n_samples; i += stride) {
        // The track gets copied into the stepper state, so that the
        // original track sample vector remains unchanged
        run_propagation((*tracks)[static_cast<std::size_t>(i)]);
      }
    } else {
      DETRAY_WARN_HOST("Running host benchmarks without warmup");
    }

    // Run the benchmark

    // Calculate the propagation rate
    // @see
    // https://github.com/google/benchmark/blob/main/docs/user_guide.md#custom-counters
    std::size_t total_tracks = 0u;
    for (auto _ : state) {
#pragma omp parallel for
      for (int i = 0; i < n_samples; ++i) {
        run_propagation((*tracks)[static_cast<std::size_t>(i)]);
      }
      total_tracks += static_cast<std::size_t>(n_samples);
    }

    // Report throughput
    state.counters["TracksPropagated"] = benchmark::Counter(
        static_cast<double>(total_tracks), benchmark::Counter::kIsRate);
  }
};

}  // namespace detray::benchmarks
