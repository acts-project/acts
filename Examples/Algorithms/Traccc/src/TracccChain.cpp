// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/TracccChain.hpp"

namespace ActsExamples {

TracccChain::TracccChain(const std::string& detector_file,
                         const std::string& digitization_file,
                         const std::string& conditions_file,
                         const std::string& material_file,
                         const std::string& grid_file,
                         const std::string& bfield_file)
    : host_mr{},
      cuda_host_mr{},
      device_mr{},
      mr{device_mr, &cuda_host_mr},
      stream{},
      copy{stream.cudaStream()},
      host_copy{},
      host_det_descr{host_mr},
      host_det_cond{host_mr},
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
      sa_cuda{seedfinder_cfg,
              spacepoint_grid_cfg,
              seedfilter_cfg,
              mr,
              copy,
              stream,
              traccc::getDefaultLogger("CudaSeedingAlg",
                                       traccc::Logging::INFO)},
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

  copy.setup(device_det_descr)->wait();
  copy(vecmem::get_data(host_det_descr), device_det_descr)->wait();
  copy.setup(device_det_cond)->wait();
  copy(vecmem::get_data(host_det_cond), device_det_cond)->wait();
}

EventResult processEvent(std::shared_ptr<TracccChain> chain,
                         const std::string& data_directory,
                         std::size_t event_id) {
  EventResult result;

  traccc::edm::silicon_cell_collection::host cells{chain->host_mr};
  static constexpr bool DEDUPLICATE = true;
  traccc::io::read_cells(
      cells, event_id, data_directory,
      traccc::getDefaultLogger("ReadCells", traccc::Logging::INFO),
      &chain->host_det_cond, traccc::data_format::csv, DEDUPLICATE, false);
  result.n_cells = cells.size();

  traccc::edm::silicon_cell_collection::buffer cells_buf(
      static_cast<unsigned int>(cells.size()), chain->mr.main);
  chain->copy.setup(cells_buf)->wait();
  chain->copy(vecmem::get_data(cells), cells_buf)->wait();

  auto unsorted_meas =
      chain->ca_cuda(cells_buf, chain->device_det_descr, chain->device_det_cond);
  traccc::edm::measurement_collection::buffer meas_buf =
      chain->ms_cuda(unsorted_meas);
  chain->stream.synchronize();

  traccc::edm::measurement_collection::host meas_host{chain->host_mr};
  chain->copy(meas_buf, meas_host)->wait();
  result.n_measurements = meas_host.size();

  traccc::edm::spacepoint_collection::buffer sp_buf =
      chain->sf_cuda(chain->device_detector, meas_buf);
  chain->stream.synchronize();

  traccc::edm::spacepoint_collection::host sp_host{chain->host_mr};
  chain->copy(sp_buf, sp_host)->wait();
  result.n_spacepoints = sp_host.size();

  traccc::edm::seed_collection::buffer seeds_buf = chain->sa_cuda(sp_buf);
  chain->stream.synchronize();

  traccc::edm::seed_collection::host seeds_host{chain->host_mr};
  chain->copy(seeds_buf, seeds_host)->wait();
  result.n_seeds = seeds_host.size();

  traccc::bound_track_parameters_collection_types::buffer params_buf =
      chain->tp_cuda(chain->device_field, meas_buf, sp_buf, seeds_buf);
  chain->stream.synchronize();

  traccc::edm::track_container<traccc::default_algebra>::buffer
      track_candidates_buf = chain->finding_cuda(
          chain->device_detector, chain->device_field, meas_buf, params_buf);

  traccc::edm::track_collection<traccc::default_algebra>::host
      track_candidates_host{chain->host_mr};
  chain->copy(track_candidates_buf.tracks, track_candidates_host,
              vecmem::copy::type::device_to_host)->wait();
  result.n_found_tracks = track_candidates_host.size();

  traccc::edm::track_container<traccc::default_algebra>::buffer resolved_buf =
      chain->resolution_cuda(track_candidates_buf);

  traccc::edm::track_collection<traccc::default_algebra>::host resolved_host{
      chain->host_mr};
  chain->copy(resolved_buf.tracks, resolved_host,
              vecmem::copy::type::device_to_host)->wait();
  result.n_resolved_tracks = resolved_host.size();

  traccc::edm::track_container<traccc::default_algebra>::buffer fitted_buf =
      chain->fitting_cuda(chain->device_detector, chain->device_field,
                          resolved_buf);
  chain->stream.synchronize();

  traccc::edm::track_collection<traccc::default_algebra>::host fitted_host{
      chain->host_mr};
  chain->copy(fitted_buf.tracks, fitted_host,
              vecmem::copy::type::device_to_host)->wait();
  result.n_fitted_tracks = fitted_host.size();

  return result;
}

}  // namespace ActsExamples