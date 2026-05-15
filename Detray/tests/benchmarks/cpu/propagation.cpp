// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/actors/parameter_updater.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/cpu/propagation_benchmark.hpp"
#include "detray/benchmarks/types.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/build_wire_chamber.hpp"
#include "detray/test/common/track_generators.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <iostream>
#include <string>

using namespace detray;

int main(int argc, char** argv) {
  using toy_detector_t = detector<benchmarks::toy_metadata>;
  using bench_algebra = typename toy_detector_t::algebra_type;
  using scalar = dscalar<bench_algebra>;
  using vector3 = dvector3D<bench_algebra>;

  using free_track_parameters_t = free_track_parameters<bench_algebra>;
  using uniform_gen_t =
      detail::random_numbers<scalar, std::uniform_real_distribution<scalar>>;
  using track_generator_t =
      random_track_generator<free_track_parameters_t, uniform_gen_t>;

  using field_t = bfield::const_field_t<scalar>;
  using stepper_t = rk_stepper<typename field_t::view_t, bench_algebra>;
  using empty_chain_t = actor_chain<>;
  using default_chain = actor_chain<actor::parameter_updater<
      bench_algebra, actor::pointwise_material_interactor<bench_algebra>>>;

  vecmem::host_memory_resource host_mr;

  //
  // Configuration
  //

  // Constant magnetic field
  vector3 B{static_cast<scalar>(0.f), static_cast<scalar>(0.f),
            2.f * unit<scalar>::T};

  // Configure toy detector
  toy_det_config<scalar> toy_cfg{};
  toy_cfg.use_material_maps(false).n_brl_layers(4u).n_edc_layers(7u);

  std::clog << toy_cfg << std::endl;

  // Configure wire chamber
  wire_chamber_config<scalar> wire_chamber_cfg{};
  wire_chamber_cfg.half_z(500.f * unit<scalar>::mm);

  std::clog << wire_chamber_cfg << std::endl;

  // Configure propagation
  propagation::config prop_cfg{};
  prop_cfg.navigation.search_window = {3u, 3u};
  // No scattering in test tracks
  prop_cfg.navigation.estimate_scattering_noise = false;

  std::clog << prop_cfg << std::endl;

  // Benchmark config
  detray::benchmarks::benchmark_base::configuration bench_cfg{};

  std::vector<int> n_tracks{8 * 8,     16 * 16,   32 * 32,  64 * 64,
                            128 * 128, 256 * 256, 512 * 512};

  auto trk_cfg =
      detray::benchmarks::get_default_trk_gen_config<track_generator_t>(
          n_tracks);

  // Specific configuration for the random track generation
  trk_cfg.seed(detail::random_numbers<scalar>::default_seed());

  // Add additional tracks for warmup
  bench_cfg.n_warmup(static_cast<int>(
      std::ceil(0.1f * static_cast<float>(trk_cfg.n_tracks()))));
  bench_cfg.do_warmup(true);

  //
  // Prepare data
  //
  auto track_samples =
      detray::benchmarks::generate_track_samples<track_generator_t>(
          &host_mr, n_tracks, trk_cfg);

  const auto [toy_det, names] =
      build_toy_detector<bench_algebra>(host_mr, toy_cfg);
  const auto [wire_chamber, _] =
      build_wire_chamber<bench_algebra>(host_mr, wire_chamber_cfg);

  auto bfield = create_const_field<scalar>(B);

  dtuple<> empty_state{};

  actor::parameter_updater_state<bench_algebra> updater_state{prop_cfg};
  actor::pointwise_material_interactor<bench_algebra>::state interactor_state{};

  auto actor_states =
      detail::make_tuple<dtuple>(updater_state, interactor_state);

  //
  // Register benchmarks
  //
  std::clog << "Propagation Benchmarks\n"
            << "----------------------\n\n";

  prop_cfg.stepping.do_covariance_transport = true;
  detray::benchmarks::register_benchmark<
      detray::benchmarks::host_propagation_bm, stepper_t, default_chain>(
      "TOY_DETECTOR_W_COV_TRANSPORT", bench_cfg, prop_cfg, toy_det, bfield,
      &actor_states, track_samples, n_tracks);

  prop_cfg.stepping.do_covariance_transport = false;
  detray::benchmarks::register_benchmark<
      detray::benchmarks::host_propagation_bm, stepper_t, empty_chain_t>(
      "TOY_DETECTOR", bench_cfg, prop_cfg, toy_det, bfield, &empty_state,
      track_samples, n_tracks);

  prop_cfg.stepping.do_covariance_transport = true;
  detray::benchmarks::register_benchmark<
      detray::benchmarks::host_propagation_bm, stepper_t, default_chain>(
      "WIRE_CHAMBER_W_COV_TRANSPORT", bench_cfg, prop_cfg, wire_chamber, bfield,
      &actor_states, track_samples, n_tracks);

  prop_cfg.stepping.do_covariance_transport = false;
  detray::benchmarks::register_benchmark<
      detray::benchmarks::host_propagation_bm, stepper_t, empty_chain_t>(
      "WIRE_CHAMBER", bench_cfg, prop_cfg, wire_chamber, bfield, &empty_state,
      track_samples, n_tracks);

  // Run benchmarks
  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
