/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/examples/cuda/full_chain_algorithm.hpp"

#include "cuda_error_check.hpp"
#include "traccc/examples/cuda/tbb_await.hpp"

// Project include(s).
#include "traccc/cuda/utils/make_magnetic_field.hpp"
#include "traccc/seeding/detail/track_params_estimation_config.hpp"

// VecMem include(s).
#include <vecmem/utils/cuda/copy.hpp>

// CUDA include(s).
#include <cuda_runtime_api.h>

// System include(s).
#include <iostream>

namespace traccc::cuda {

namespace {
await_function_type get_await_function(await_strategy await_mode) {
  switch (await_mode) {
    case await_strategy::sync_event:
      return await_sync_event;
    case await_strategy::callback:
      return tbb_await_callback;
    default:
      throw std::invalid_argument("Unknown await strategy");
  }
}
}  // namespace

full_chain_algorithm::shared_data::shared_data(
    vecmem::memory_resource& host_mr,
    const detector_design_description::host& det_descr,
    const detector_conditions_description::host& det_cond,
    const magnetic_field& field, const host_detector* detector)
    : m_device_mr(),
      m_field(make_magnetic_field(field)),
      m_device_det_descr(
          [&]() {
            // number of elements in the detector design description
            std::vector<unsigned int> sizes(det_descr.size());
            for (std::size_t i = 0; i < det_descr.size(); ++i) {
              auto this_design = det_descr.at(i);
              // now for each element, set the size to the largest size of
              // that element across all modules
              sizes[i] = std::max(static_cast<unsigned int>(
                                      ((this_design.bin_edges_x()).size())),
                                  static_cast<unsigned int>(
                                      ((this_design.bin_edges_y()).size())));
            }
            return sizes;
          }(),
          m_device_mr, &host_mr, vecmem::data::buffer_type::resizable),
      m_device_det_cond(
          static_cast<detector_conditions_description::buffer::size_type>(
              det_cond.size()),
          m_device_mr),
      m_detector(detector) {
  // Copy the detector (description) to the device.
  vecmem::cuda::copy copy;
  copy.setup(m_device_det_descr)->wait();
  copy(vecmem::get_data(det_descr), m_device_det_descr)->wait();
  copy(vecmem::get_data(det_cond), m_device_det_cond)->wait();
  if (m_detector != nullptr) {
    m_device_detector =
        traccc::buffer_from_host_detector(*m_detector, m_device_mr, copy);
  }
}

full_chain_algorithm::full_chain_algorithm(
    vecmem::memory_resource& host_mr,
    const clustering_config& clustering_config,
    const seedfinder_config& finder_config,
    const spacepoint_grid_config& grid_config,
    const seedfilter_config& filter_config,
    const gbts_seedfinder_config& gbts_config,
    const track_params_estimation_config& track_params_estimation_config,
    const finding_algorithm::config_type& finding_config,
    const fitting_algorithm::config_type& fitting_config,
    const shared_data& data, std::unique_ptr<const traccc::Logger> logger,
    bool useGBTS, await_strategy await_mode)
    : messaging(logger->clone()),
      m_host_mr(host_mr),
      m_pinned_host_mr(),
      m_cached_pinned_host_mr(m_pinned_host_mr),
      m_vecmem_stream{},
      m_stream{m_vecmem_stream.stream()},
      m_device_mr(),
      m_cached_device_mr(m_device_mr),
      m_copy(m_stream.cudaStream()),
      m_await_function(get_await_function(await_mode)),
      m_field(data.m_field),
      m_device_det_descr(data.m_device_det_descr),
      m_device_det_cond(data.m_device_det_cond),
      m_detector(data.m_detector),
      m_device_detector(data.m_device_detector),
      m_clusterization({m_cached_device_mr, &m_cached_pinned_host_mr}, m_copy,
                       m_stream, clustering_config,
                       logger->cloneWithSuffix("ClusteringAlg"),
                       m_await_function),
      m_measurement_sorting({m_cached_device_mr, &m_cached_pinned_host_mr},
                            m_copy, m_stream,
                            logger->cloneWithSuffix("MeasSortingAlg")),
      m_spacepoint_formation(
          {m_cached_device_mr, &m_cached_pinned_host_mr}, m_copy, m_stream,
          logger->cloneWithSuffix("SpFormationAlg"), m_await_function),
      m_seeding(finder_config, grid_config, filter_config,
                {m_cached_device_mr, &m_cached_pinned_host_mr}, m_copy,
                m_stream, logger->cloneWithSuffix("SeedingAlg"),
                m_await_function),
      m_gbts_seeding(gbts_config,
                     {m_cached_device_mr, &m_cached_pinned_host_mr}, m_copy,
                     m_stream, logger->cloneWithSuffix("GbtsAlg")),
      m_track_parameter_estimation(
          track_params_estimation_config,
          {m_cached_device_mr, &m_cached_pinned_host_mr}, m_copy, m_stream,
          logger->cloneWithSuffix("TrackParEstAlg"), m_await_function),
      m_finding(finding_config, {m_cached_device_mr, &m_cached_pinned_host_mr},
                m_copy, m_stream, logger->cloneWithSuffix("TrackFindingAlg")),
      m_fitting(fitting_config, {m_cached_device_mr, &m_cached_pinned_host_mr},
                m_copy, m_stream, logger->cloneWithSuffix("TrackFittingAlg")),
      usingGBTS(useGBTS) {
  // Tell the user what device is being used.
  int device = 0;
  CUDA_ERROR_CHECK(cudaGetDevice(&device));
  cudaDeviceProp props;
  CUDA_ERROR_CHECK(cudaGetDeviceProperties(&props, device));
  std::cout << "Using CUDA device: " << props.name << " [id: " << device
            << ", bus: " << props.pciBusID << ", device: " << props.pciDeviceID
            << "]" << std::endl;
}

full_chain_algorithm::~full_chain_algorithm() = default;

full_chain_algorithm::output_type full_chain_algorithm::operator()(
    const edm::silicon_cell_collection::host& cells) const {
  // Create device copy of input collections
  edm::silicon_cell_collection::buffer cells_buffer(
      static_cast<unsigned int>(cells.size()), m_cached_device_mr);
  m_copy(vecmem::get_data(cells), cells_buffer)->ignore();

  // Run the clusterization (asynchronously).
  const auto unsorted_measurements =
      m_clusterization(cells_buffer, m_device_det_descr, m_device_det_cond);
  const measurement_sorting_algorithm::output_type measurements =
      m_measurement_sorting(unsorted_measurements);

  // If we have a Detray detector, run the seeding, track finding and fitting.
  if (m_detector != nullptr) {
    // Run the seed-finding (asynchronously).
    const spacepoint_formation_algorithm::output_type spacepoints =
        m_spacepoint_formation(m_device_detector, measurements);

    triplet_seeding_algorithm::output_type seeds;
    if (usingGBTS) {
      seeds = m_gbts_seeding(spacepoints, measurements);
    } else {
      seeds = m_seeding(spacepoints);
    }
    const seed_parameter_estimation_algorithm::output_type track_params =
        m_track_parameter_estimation(m_field, measurements, spacepoints, seeds);

    // Run the track finding (asynchronously).
    const finding_algorithm::output_type track_candidates =
        m_finding(m_device_detector, m_field, measurements, track_params);

    // Copy a limited amount of result data back to the host.
    const auto host_tracks =
        m_copy.to(track_candidates.tracks, m_cached_pinned_host_mr, nullptr,
                  vecmem::copy::type::device_to_host);
    output_type result{m_host_mr};
    vecmem::copy host_copy;
    host_copy(host_tracks, result)->wait();
    return result;

  }
  // If not, copy the measurements back to the host, and return a dummy
  // object.
  else {
    // Copy the measurements back to the host.
    edm::measurement_collection::host measurements_host(m_host_mr);
    m_copy(measurements, measurements_host)->wait();

    // Return an empty object.
    return output_type{m_host_mr};
  }
}

bound_track_parameters_collection_types::host full_chain_algorithm::seeding(
    const edm::silicon_cell_collection::host& cells) const {
  // Create device copy of input collections
  edm::silicon_cell_collection::buffer cells_buffer(
      static_cast<unsigned int>(cells.size()), m_cached_device_mr);
  m_copy(vecmem::get_data(cells), cells_buffer)->ignore();

  // Run the clusterization (asynchronously).
  const auto unsorted_measurements =
      m_clusterization(cells_buffer, m_device_det_descr, m_device_det_cond);
  const measurement_sorting_algorithm::output_type measurements =
      m_measurement_sorting(unsorted_measurements);

  // If we have a Detray detector, run the seeding, track finding and fitting.
  if (m_detector != nullptr) {
    // Run the seed-finding (asynchronously).
    const spacepoint_formation_algorithm::output_type spacepoints =
        m_spacepoint_formation(m_device_detector, measurements);

    triplet_seeding_algorithm::output_type seeds;
    if (usingGBTS) {
      seeds = m_gbts_seeding(spacepoints, measurements);
    } else {
      seeds = m_seeding(spacepoints);
    }
    const seed_parameter_estimation_algorithm::output_type track_params =
        m_track_parameter_estimation(m_field, measurements, spacepoints, seeds);

    // Copy a limited amount of result data back to the host.
    const auto host_seeds = m_copy.to(track_params, m_cached_pinned_host_mr,
                                      vecmem::copy::type::device_to_host);
    bound_track_parameters_collection_types::host result{&m_host_mr};
    vecmem::copy host_copy;
    host_copy(host_seeds, result)->wait();
    return result;

  }
  // If not, copy the measurements back to the host, and return a dummy
  // object.
  else {
    // Copy the measurements back to the host.
    edm::measurement_collection::host measurements_host(m_host_mr);
    m_copy(measurements, measurements_host)->wait();

    // Return an empty object.
    return {};
  }
}

}  // namespace traccc::cuda
