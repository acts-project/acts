// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsExamples/Traccc/DetrayPropagator.hpp"
#include "ActsExamples/Traccc/DetrayStore.hpp"
#include "ActsPlugins/Covfie/FieldConversion.hpp"
#include "ActsPlugins/Detray/DetrayConversionUtils.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

// detray includes
#include <detray/core/detector.hpp>
#include <detray/io/frontend/detector_reader.hpp>
#include <detray/navigation/volume_graph.hpp>
#include <detray/propagator/line_stepper.hpp>
#include <detray/propagator/propagator.hpp>
#include <detray/propagator/rk_stepper.hpp>
#include <pybind11/pybind11.h>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

// traccc includes
#include <memory>
#include <string>

#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>

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

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPlugins;
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
  /// Define host detray store
  {
    py::class_<DetrayHostStore, std::shared_ptr<DetrayHostStore>>(
        traccc, "DetrayHostStore");

    traccc.def("readDetectorHost", [](const std::string& geometry,
                                      const std::string& materials,
                                      const std::string& grids) {
      auto mr = std::make_shared<vecmem::host_memory_resource>();
      auto reader_cfg = detray::io::detector_reader_config{};
      reader_cfg.add_file(geometry);
      if (!materials.empty())
        reader_cfg.add_file(materials);
      if (!grids.empty())
        reader_cfg.add_file(grids);
      auto [det, names] =
          detray::io::read_detector<DetrayHostDetector>(*mr, reader_cfg);
      return DetrayHostStore{std::move(mr), std::move(det)};
    });
  }

  /// Propagators
  {
    traccc.def(
        "createSlPropagatorHost",
        [](const std::shared_ptr<const DetrayHostStore>& detrayStore,
           bool sterile) {
          using DetrayLineStepper =
              detray::line_stepper<typename DetrayHostDetector::algebra_type>;
          using DP = DetrayPropagator<DetrayLineStepper, DetrayHostStore>;
          DP::Config cfg{detrayStore, sterile};
          return std::shared_ptr<PropagatorInterface>(
              std::make_shared<DP>(cfg));
        },
        "store"_a, "sterile"_a = false);

    traccc.def(
        "createConstBFieldPropagatorHost",
        [](const std::shared_ptr<const DetrayHostStore>& detrayStore,
           Covfie::ConstantField cfield, bool sterile) {
          using Stepper =
              detray::rk_stepper<Covfie::ConstantField::view_t,
                                 typename DetrayHostDetector::algebra_type>;
          using DP = DetrayPropagator<Stepper, DetrayHostStore,
                                      Covfie::ConstantField::view_t>;
          DP::Config cfg{detrayStore, sterile, cfield};
          return std::shared_ptr<PropagatorInterface>(
              std::make_shared<DP>(cfg));
        },
        "store"_a, "field"_a, "sterile"_a = false);
  }

  /// The full GPU sequence example (mirrors traccc_seq_example setup)
  {
    py::class_<TracccSeqAlg, std::shared_ptr<TracccSeqAlg>>(traccc,
                                                            "TracccSeqAlg");

    traccc.def(
        "createSeqAlg",
        [](const std::string& detector_file,
           const std::string& digitization_file,
           const std::string& conditions_file, const std::string& material_file,
           const std::string& grid_file, const std::string& bfield_file) {
          return std::make_shared<TracccSeqAlg>(
              detector_file, digitization_file, conditions_file, material_file,
              grid_file, bfield_file);
        },
        "detector_file"_a, "digitization_file"_a, "conditions_file"_a = "",
        "material_file"_a = "", "grid_file"_a = "", "bfield_file"_a);
  }

  /// EventResult — mirrors per-event counters from traccc_seq_example
  {
    py::class_<EventResult>(traccc, "EventResult")
        .def_readonly("n_cells", &EventResult::n_cells)
        .def_readonly("n_measurements", &EventResult::n_measurements)
        .def_readonly("n_spacepoints", &EventResult::n_spacepoints)
        .def_readonly("n_seeds", &EventResult::n_seeds)
        .def_readonly("n_found_tracks", &EventResult::n_found_tracks)
        .def_readonly("n_resolved_tracks", &EventResult::n_resolved_tracks)
        .def_readonly("n_fitted_tracks", &EventResult::n_fitted_tracks);
  }

  /// Process one event
  {
    traccc.def(
        "processEvent",
        [](std::shared_ptr<TracccSeqAlg> s, const std::string& data_directory,
           std::size_t event_id) -> EventResult {
          EventResult result;

          // --- Read cells ---
          traccc::edm::silicon_cell_collection::host cells{s->host_mr};
          static constexpr bool DEDUPLICATE = true;
          traccc::io::read_cells(
              cells, event_id, data_directory,
              traccc::getDefaultLogger("ReadCells", traccc::Logging::INFO),
              &s->host_det_cond, traccc::data_format::csv, DEDUPLICATE, false);
          result.n_cells = cells.size();

          // --- Upload cells to device ---
          traccc::edm::silicon_cell_collection::buffer cells_buf(
              static_cast<unsigned int>(cells.size()), s->mr.main);
          s->copy.setup(cells_buf)->wait();
          s->copy(vecmem::get_data(cells), cells_buf)->wait();

          // --- Clusterization ---
          auto unsorted_meas =
              s->ca_cuda(cells_buf, s->device_det_descr, s->device_det_cond);
          traccc::edm::measurement_collection::buffer meas_buf =
              s->ms_cuda(unsorted_meas);
          s->stream.synchronize();

          traccc::edm::measurement_collection::host meas_host{s->host_mr};
          s->copy(meas_buf, meas_host)->wait();
          result.n_measurements = meas_host.size();

          // --- Spacepoint formation ---
          traccc::edm::spacepoint_collection::buffer sp_buf =
              s->sf_cuda(s->device_detector, meas_buf);
          s->stream.synchronize();

          traccc::edm::spacepoint_collection::host sp_host{s->host_mr};
          s->copy(sp_buf, sp_host)->wait();
          result.n_spacepoints = sp_host.size();

          // --- Seeding ---
          traccc::edm::seed_collection::buffer seeds_buf = s->sa_cuda(sp_buf);
          s->stream.synchronize();

          traccc::edm::seed_collection::host seeds_host{s->host_mr};
          s->copy(seeds_buf, seeds_host)->wait();
          result.n_seeds = seeds_host.size();

          // --- Track parameter estimation ---
          traccc::bound_track_parameters_collection_types::buffer params_buf =
              s->tp_cuda(s->device_field, meas_buf, sp_buf, seeds_buf);
          s->stream.synchronize();

          // --- Track finding ---
          traccc::edm::track_container<traccc::default_algebra>::buffer
              track_candidates_buf = s->finding_cuda(
                  s->device_detector, s->device_field, meas_buf, params_buf);

          traccc::edm::track_collection<traccc::default_algebra>::host
              track_candidates_host{s->host_mr};
          s->copy(track_candidates_buf.tracks, track_candidates_host,
                  vecmem::copy::type::device_to_host)
              ->wait();
          result.n_found_tracks = track_candidates_host.size();

          // --- Ambiguity resolution ---
          traccc::edm::track_container<traccc::default_algebra>::buffer
              resolved_buf = s->resolution_cuda(track_candidates_buf);

          traccc::edm::track_collection<traccc::default_algebra>::host
              resolved_host{s->host_mr};
          s->copy(resolved_buf.tracks, resolved_host,
                  vecmem::copy::type::device_to_host)
              ->wait();
          result.n_resolved_tracks = resolved_host.size();

          // --- Track fitting ---
          traccc::edm::track_container<traccc::default_algebra>::buffer
              fitted_buf = s->fitting_cuda(s->device_detector, s->device_field,
                                           resolved_buf);
          s->stream.synchronize();

          traccc::edm::track_collection<traccc::default_algebra>::host
              fitted_host{s->host_mr};
          s->copy(fitted_buf.tracks, fitted_host,
                  vecmem::copy::type::device_to_host)
              ->wait();
          result.n_fitted_tracks = fitted_host.size();

          return result;
        },
        "state"_a, "data_directory"_a, "event_id"_a);
  }
}
