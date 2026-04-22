// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s)
#include "detray/navigation/caching_navigator.hpp"
#include "detray/propagator/actors.hpp"
#include "detray/propagator/rk_stepper.hpp"
#include "detray/tracks/tracks.hpp"
#include "detray/utils/type_list.hpp"

// Detray IO include(s)
#include "detray/io/frontend/detector_reader.hpp"

// Detray benchmark include(s)
#include "detray/benchmarks/device/hip/propagation_benchmark.hpp"
#include "detray/benchmarks/types.hpp"

// Detray test include(s)
#include "detray/test/common/bfield.hpp"
#include "detray/test/common/track_generators.hpp"

// Detray tools include(s)
#include "detray/options/detector_io_options.hpp"
#include "detray/options/parse_options.hpp"
#include "detray/options/propagation_options.hpp"
#include "detray/options/track_generator_options.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// System include(s)
#include <algorithm>
#include <string>

namespace po = boost::program_options;

using namespace detray;

int main(int argc, char** argv) {
  // Use the most general type to be able to read in all detector files
  using detector_t = detray::detector<benchmarks::default_metadata>;
  using bench_algebra = typename detector_t::algebra_type;
  using scalar = dscalar<bench_algebra>;
  using vector3 = dvector3D<bench_algebra>;

  using free_track_parameters_t = free_track_parameters<bench_algebra>;
  using uniform_gen_t =
      detail::random_numbers<scalar, std::uniform_real_distribution<scalar>>;
  using track_generator_t =
      random_track_generator<free_track_parameters_t, uniform_gen_t>;

  using field_bknd_t = bfield::const_bknd_t<scalar>;

  // Host and device memory resources
  vecmem::host_memory_resource host_mr;
  vecmem::hip::device_memory_resource dev_mr;

  // Constant magnetic field
  vector3 B{0.f, 0.f, 2.f * unit<scalar>::T};

  // Number of tracks in the different benchmark cases
  std::vector<int> n_tracks{10,     100,    500,     1000,   5000,
                            10'000, 50'000, 100'000, 250'000};

  //
  // Configuration
  //

  // Google benchmark specific options
  ::benchmark::Initialize(&argc, argv);

  // Specific options for this test
  po::options_description desc("\ndetray propagation benchmark options");

  desc.add_options()("context", po::value<dindex>(),
                     "Index of the geometry context")(
      "bknd_name", po::value<std::string>(), "Name of the Processor")(
      "sort_tracks", "Sort track samples by theta angle");

  // Configs to be filled
  detray::io::detector_reader_config reader_cfg{};
  track_generator_t::configuration trk_cfg{};
  propagation::config prop_cfg{};
  detray::benchmarks::benchmark_base::configuration bench_cfg{};

  // Read options from commandline
  po::variables_map vm = detray::options::parse_options(
      desc, argc, argv, reader_cfg, trk_cfg, prop_cfg);

  // Custom options
  bool do_sort{(vm.count("sort_tracks") != 0)};

  // The geometry context to be used
  detector_t::geometry_context gctx;
  if (vm.count("context")) {
    gctx = detector_t::geometry_context{vm["context"].as<dindex>()};
  }
  std::string proc_name{"unknown"};
  if (vm.count("bknd_name")) {
    proc_name = vm["bknd_name"].as<std::string>();
  }

  // String that describes the detector setup
  std::string setup_str{};
  auto add_delim = [](std::string& str) { str += ", "; };
  if (!vm.count("grid_file")) {
    setup_str += "no grids";
  }
  if (!vm.count("material_file")) {
    if (!setup_str.empty()) {
      add_delim(setup_str);
    }
    setup_str += "no mat.";
  }

  //
  // Prepare data
  //

  // Read the detector geometry
  reader_cfg.do_check(true);

  const auto [det, names] =
      detray::io::read_detector<detector_t>(host_mr, reader_cfg);
  const std::string& det_name = det.name(names);

  // Generate the track samples
  auto track_samples =
      detray::benchmarks::generate_track_samples<track_generator_t>(
          &host_mr, n_tracks, trk_cfg, do_sort);

  // Create a constant b-field
  auto bfield = create_const_field<scalar>(B);

  // Build actor states
  dtuple<> empty_state{};

  actor::pointwise_material_interactor<bench_algebra>::state interactor_state{};
  actor::parameter_updater_state<bench_algebra> updater_state{prop_cfg};

  auto actor_states =
      detail::make_tuple<dtuple>(interactor_state, updater_state);

  //
  // Register benchmarks
  //

  // Number of warmup tracks
  const int n_max_tracks{*std::ranges::max_element(n_tracks)};
  bench_cfg.n_warmup(
      static_cast<int>(std::ceil(0.1f * static_cast<float>(n_max_tracks))));

  if (prop_cfg.stepping.do_covariance_transport) {
    detray::benchmarks::register_benchmark<
        detray::benchmarks::hip_propagation_bm,
        detray::benchmarks::hip_propagator_type<
            detray::benchmarks::default_metadata, field_bknd_t,
            detray::benchmarks::default_chain>>(
        det_name + "_W_COV_TRANSPORT", bench_cfg, prop_cfg, det, bfield,
        &actor_states, track_samples, n_tracks, &dev_mr);
  } else {
    detray::benchmarks::register_benchmark<
        detray::benchmarks::hip_propagation_bm,
        detray::benchmarks::hip_propagator_type<
            detray::benchmarks::default_metadata, field_bknd_t,
            detray::benchmarks::empty_chain>>(det_name, bench_cfg, prop_cfg,
                                              det, bfield, &empty_state,
                                              track_samples, n_tracks, &dev_mr);

    if (!setup_str.empty()) {
      add_delim(setup_str);
    }
    setup_str += "no cov.";
  }

  // These fields are needed by the plotting scripts, even if undefined
  ::benchmark::AddCustomContext("Backend", HIP_PLATFORM_NAME);
  ::benchmark::AddCustomContext("Backend Name", proc_name);
  ::benchmark::AddCustomContext("Algebra-plugin",
                                detray::types::get_name<bench_algebra>());
  ::benchmark::AddCustomContext("Detector Setup", setup_str);

  // Run benchmarks
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();
}
