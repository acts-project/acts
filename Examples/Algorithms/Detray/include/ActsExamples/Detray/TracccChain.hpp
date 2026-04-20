// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <string>

#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/copy.hpp>

#include "traccc/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/cuda/ambiguity_resolution/greedy_ambiguity_resolution_algorithm.hpp"
#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/cuda/clusterization/measurement_sorting_algorithm.hpp"
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/cuda/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/cuda/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/cuda/seeding/triplet_seeding_algorithm.hpp"
#include "traccc/cuda/utils/make_magnetic_field.hpp"
#include "traccc/cuda/utils/stream.hpp"
#include "traccc/device/container_d2h_copy_alg.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/io/read_cells.hpp"
#include "traccc/io/read_detector.hpp"
#include "traccc/io/read_detector_description.hpp"
#include "traccc/io/read_magnetic_field.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/track_params_estimation_config.hpp"
#include "traccc/utils/propagation.hpp"

namespace ActsExamples {

/// Holds all GPU/host objects that must outlive the per-event calls.
/// Constructed once at algorithm initialisation, reused for every event.
struct TracccChain {
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::host_memory_resource cuda_host_mr;
  vecmem::cuda::device_memory_resource device_mr;
  traccc::memory_resource mr;
  traccc::cuda::stream stream;
  vecmem::cuda::async_copy copy;
  vecmem::copy host_copy;

  traccc::detector_design_description::host host_det_descr;
  traccc::detector_conditions_description::host host_det_cond;

  traccc::detector_design_description::buffer device_det_descr;
  traccc::detector_conditions_description::buffer device_det_cond;

  traccc::host_detector host_detector;
  traccc::detector_buffer device_detector;

  traccc::magnetic_field host_field;
  traccc::magnetic_field device_field;

  traccc::seedfinder_config seedfinder_cfg;
  traccc::seedfilter_config seedfilter_cfg;
  traccc::spacepoint_grid_config spacepoint_grid_cfg;
  traccc::track_params_estimation_config track_params_cfg;
  traccc::finding_config finding_cfg;
  traccc::fitting_config fitting_cfg;
  traccc::host::greedy_ambiguity_resolution_algorithm::config_type
      resolution_cfg;

  traccc::cuda::clusterization_algorithm ca_cuda;
  traccc::cuda::measurement_sorting_algorithm ms_cuda;
  traccc::cuda::silicon_pixel_spacepoint_formation_algorithm sf_cuda;
  traccc::cuda::triplet_seeding_algorithm sa_cuda;
  traccc::cuda::seed_parameter_estimation_algorithm tp_cuda;
  traccc::cuda::combinatorial_kalman_filter_algorithm finding_cuda;
  traccc::cuda::greedy_ambiguity_resolution_algorithm resolution_cuda;
  traccc::cuda::kalman_fitting_algorithm fitting_cuda;

  TracccChain(const std::string& detector_file,
              const std::string& digitization_file,
              const std::string& conditions_file,
              const std::string& material_file,
              const std::string& grid_file,
              const std::string& bfield_file);
};

/// Per-event output counters.
struct EventResult {
  std::size_t n_cells = 0;
  std::size_t n_measurements = 0;
  std::size_t n_spacepoints = 0;
  std::size_t n_seeds = 0;
  std::size_t n_found_tracks = 0;
  std::size_t n_resolved_tracks = 0;
  std::size_t n_fitted_tracks = 0;
};

EventResult processEvent(std::shared_ptr<TracccChain> chain,
                         const std::string& data_directory,
                         std::size_t event_id);

}  // namespace ActsExamples