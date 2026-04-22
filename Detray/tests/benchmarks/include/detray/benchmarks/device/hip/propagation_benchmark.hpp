// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/logging.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/benchmark_base.hpp"
#include "detray/benchmarks/propagation_benchmark.hpp"
#include "detray/benchmarks/propagation_benchmark_utils.hpp"
#include "detray/benchmarks/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/hip/device_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/hip/copy.hpp>

// Benchmark include
#include <benchmark/benchmark.h>

// System include(s)
#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <string>

namespace detray::benchmarks {

// Define propagator type
template <concepts::algebra algebra_t>
using empty_chain = actor_chain<>;

template <concepts::algebra algebra_t>
using default_chain = actor_chain<actor::parameter_updater<
    algebra_t, actor::pointwise_material_interactor<algebra_t>>>;

using const_field_t = bfield::const_bknd_t<benchmarks::scalar>;

template <typename metadata_t, typename bfield_t,
          template <typename> class actor_chain_t>
using hip_propagator_type =
    propagator<rk_stepper<covfie::field_view<bfield_t>,
                          typename detector<metadata_t>::algebra_type>,
               caching_navigator<detector<metadata_t>>,
               actor_chain_t<typename detector<metadata_t>::algebra_type>>;

/// Launch the propagation kernelfor benchmarking
///
/// @param cfg the propagation configuration
/// @param det_view the detector vecmem view
/// @param field_data the magentic field view (maybe an empty field)
/// @param tracks_data the track collection view
/// @param navigation_cache_view the navigation cache vecemem view
/// @param opt which propagation to run (sync vs. unsync)
template <typename propagator_t,
          detray::benchmarks::propagation_opt kOPT =
              detray::benchmarks::propagation_opt::e_unsync>
void run_propagation_kernel(
    const propagation::config &,
    typename propagator_t::detector_type::view_type,
    typename propagator_t::stepper_type::magnetic_field_type,
    typename propagator_t::actor_chain_type::state_tuple *,
    vecmem::data::vector_view<
        free_track_parameters<typename propagator_t::algebra_type>>,
    const int);

/// Allocate actor state blueprint on device
/// @note This only works if each actor state in the tuple is essentially POD
template <typename propagator_t>
typename propagator_t::actor_chain_type::state_tuple *setup_actor_states(
    typename propagator_t::actor_chain_type::state_tuple *);

/// Release actor state blueprint
template <typename propagator_t>
void release_actor_states(
    typename propagator_t::actor_chain_type::state_tuple *);

/// Device Propagation becnhmark
template <typename propagator_t, typename bfield_bknd_t,
          detray::benchmarks::propagation_opt kOPT =
              detray::benchmarks::propagation_opt::e_unsync>
struct hip_propagation_bm
    : public propagation_benchmark<
          typename propagator_t::detector_type::algebra_type> {
  /// Detector dependent types
  using algebra_t = typename propagator_t::detector_type::algebra_type;
  using base_type = propagation_benchmark<algebra_t>;

  using base_type::base_type;

  /// Prepare data and run benchmark loop
  inline void operator()(::benchmark::State &state,
                         vecmem::memory_resource *dev_mr,
                         dvector<free_track_parameters<algebra_t>> *tracks,
                         const typename propagator_t::detector_type *det,
                         const bfield_bknd_t *bfield,
                         typename propagator_t::actor_chain_type::state_tuple
                             *input_actor_states) const {
    assert(dev_mr != nullptr);
    assert(tracks != nullptr);
    assert(det != nullptr);
    assert(bfield != nullptr);
    assert(input_actor_states != nullptr);

    // Helper object for performing memory copies (to HIP devices)
    vecmem::hip::copy hip_cpy;

    const int n_samples{this->config().n_samples()};
    const int n_warmup{this->config().n_warmup()};

    assert(static_cast<std::size_t>(n_samples) <= tracks->size());

    // Copy the track collection to device
    auto track_buffer =
        detray::get_buffer(vecmem::get_data(*tracks), *dev_mr, hip_cpy);

    // Copy the detector to device and get its view
    auto det_buffer = detray::get_buffer(*det, *dev_mr, hip_cpy);
    auto det_view = detray::get_data(det_buffer);

    // Copy blueprint actor states to device
    auto *device_actor_state_ptr =
        setup_actor_states<propagator_t>(input_actor_states);

    // Do a small warm up run
    if (this->config().do_warmup()) {
      auto warmup_track_buffer =
          detray::get_buffer(vecmem::get_data(*tracks), *dev_mr, hip_cpy);

      run_propagation_kernel<propagator_t, kOPT>(
          this->propagation(), det_view, *bfield, device_actor_state_ptr,
          warmup_track_buffer, math::min(n_warmup, n_samples));
    } else {
      DETRAY_WARN_HOST(
          "Running HIP benchmarks without warmup is "
          "not recommended");
    }

    // Calculate the propagation rate
    // @see
    // https://github.com/google/benchmark/blob/main/docs/user_guide.md#custom-counters
    std::size_t total_tracks = 0u;
    for (auto _ : state) {
      // Launch the propagator test for GPU device
      run_propagation_kernel<propagator_t, kOPT>(
          this->propagation(), det_view, *bfield, device_actor_state_ptr,
          track_buffer, n_samples);

      total_tracks += static_cast<std::size_t>(n_samples);
    }

    // Report throughput
    state.counters["TracksPropagated"] = benchmark::Counter(
        static_cast<double>(total_tracks), benchmark::Counter::kIsRate);

    release_actor_states<propagator_t>(device_actor_state_ptr);
  }
};

}  // namespace detray::benchmarks
