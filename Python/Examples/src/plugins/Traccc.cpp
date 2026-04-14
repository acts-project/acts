// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Traccc/TracccSeqAlg.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ActsExamples;
using namespace ActsPython;

// ---------------------------------------------------------------------------
// Holds all objects that must outlive the algorithm object.
// Mirrors what seq_run() in traccc_seq_example sets up before the event loop.
// ---------------------------------------------------------------------------
struct TracccSeqAlg {
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::host_memory_resource cuda_host_mr;
  vecmem::cuda::device_memory_resource device_mr;
  traccc::memory_resource mr;
  traccc::cuda::stream stream;
  vecmem::cuda::async_copy copy;
  vecmem::copy host_copy;

  traccc::detector_design_description::host host_det_descr;
  traccc::detector_conditions_description::host host_det_cond;

  // These are default-constructible, so we can leave them uninitialised until
  // we know the sizes
  traccc::detector_design_description::buffer device_det_descr;
  traccc::detector_conditions_description::buffer device_det_cond;

  traccc::host_detector host_detector;
  traccc::detector_buffer device_detector;

  traccc::magnetic_field host_field;
  traccc::magnetic_field device_field;

  // These are default-constructed, but can be modified by the user before
  // processing events if desired
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

  TracccSeqAlg(const std::string& detector_file,
               const std::string& digitization_file,
               const std::string& conditions_file,
               const std::string& material_file, const std::string& grid_file,
               const std::string& bfield_file)
      : host_mr{},
        cuda_host_mr{},
        device_mr{},
        mr{device_mr, &cuda_host_mr},
        stream{},
        copy{stream.cudaStream()},
        host_copy{},
        host_det_descr{host_mr},
        host_det_cond{host_mr}
        // leave device buffers, detector, field uninitialised for now
        ,
        seedfinder_cfg{},
        seedfilter_cfg{},
        spacepoint_grid_cfg{seedfinder_cfg},
        track_params_cfg{},
        finding_cfg{},
        fitting_cfg{},
        resolution_cfg{},
        ca_cuda{mr, copy, stream, traccc::clustering_config{},
                traccc::getDefaultLogger("CudaClusteringAlg",
                                         traccc::Logging::INFO)},
        ms_cuda{mr, copy, stream,
                traccc::getDefaultLogger("CudaMeasSortingAlg",
                                         traccc::Logging::INFO)},
        sf_cuda{mr, copy, stream,
                traccc::getDefaultLogger("CudaSpFormationAlg",
                                         traccc::Logging::INFO)},
        sa_cuda{
            seedfinder_cfg,
            spacepoint_grid_cfg,
            seedfilter_cfg,
            mr,
            copy,
            stream,
            traccc::getDefaultLogger("CudaSeedingAlg", traccc::Logging::INFO)},
        tp_cuda{track_params_cfg, mr, copy, stream,
                traccc::getDefaultLogger("CudaTrackParEstAlg",
                                         traccc::Logging::INFO)},
        finding_cuda{
            finding_cfg, mr, copy, stream,
            traccc::getDefaultLogger("CudaFindingAlg", traccc::Logging::INFO)},
        resolution_cuda{resolution_cfg, mr, copy, stream,
                        traccc::getDefaultLogger("CudaAmbiguityResolutionAlg",
                                                 traccc::Logging::INFO)},
        fitting_cuda{
            fitting_cfg, mr, copy, stream,
            traccc::getDefaultLogger("CudaFittingAlg", traccc::Logging::INFO)} {
    // Step 1: read detector description (host objects)
    traccc::io::read_detector_description(
        host_det_descr, host_det_cond, detector_file, digitization_file,
        conditions_file, traccc::data_format::json);

    traccc::io::read_detector(host_detector, host_mr, detector_file,
                              material_file, grid_file);

    device_detector =
        traccc::buffer_from_host_detector(host_detector, device_mr, copy);
    stream.synchronize();

    traccc::io::read_magnetic_field(host_field, bfield_file);
    device_field = traccc::cuda::make_magnetic_field(host_field);

    // Step 2: compute sizes the same way the standalone does for device objects
    // This is necessary for modules with varying design descriptions (e.g.
    // ATLAS ITk)
    device_det_descr = traccc::detector_design_description::buffer{
        [&]() {
          std::vector<unsigned int> sizes(host_det_descr.size());
          for (std::size_t i = 0; i < host_det_descr.size(); ++i) {
            auto this_design = host_det_descr.at(i);
            sizes[i] = std::max(
                static_cast<unsigned int>(this_design.bin_edges_x().size()),
                static_cast<unsigned int>(this_design.bin_edges_y().size()));
          }
          return sizes;
        }(),
        device_mr, &host_mr, vecmem::data::buffer_type::resizable};

    device_det_cond = traccc::detector_conditions_description::buffer{
        static_cast<traccc::detector_conditions_description::buffer::size_type>(
            host_det_cond.size()),
        device_mr};

    // Step 3: upload to device
    copy.setup(device_det_descr)->wait();
    copy(vecmem::get_data(host_det_descr), device_det_descr)->wait();
    copy.setup(device_det_cond)->wait();
    copy(vecmem::get_data(host_det_cond), device_det_cond)->wait();
  }
};

// ---------------------------------------------------------------------------
// Result of processing one event — counters only
// ---------------------------------------------------------------------------
struct EventResult {
  std::size_t n_cells = 0;
  std::size_t n_measurements = 0;
  std::size_t n_spacepoints = 0;
  std::size_t n_seeds = 0;
  std::size_t n_found_tracks = 0;
  std::size_t n_resolved_tracks = 0;
  std::size_t n_fitted_tracks = 0;
};

PYBIND11_MODULE(ActsExamplesPythonBindingsTraccc, traccc) {

  ACTS_PYTHON_DECLARE_ALGORITHM(TracccSeqAlgorithm, traccc,
                                "TracccSeqAlgorithm", detectorFile,
                                digitizationFile, conditionsFile, materialFile,
                                gridFile, bfieldFile, dataDirectory);
}