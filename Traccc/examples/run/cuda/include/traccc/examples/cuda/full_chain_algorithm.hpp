/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/clusterization/clustering_config.hpp"
#include "traccc/cuda/clusterization/clusterization_algorithm.hpp"
#include "traccc/cuda/clusterization/measurement_sorting_algorithm.hpp"
#include "traccc/cuda/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/cuda/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/cuda/gbts_seeding/gbts_seeding_algorithm.hpp"
#include "traccc/cuda/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/cuda/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/cuda/seeding/triplet_seeding_algorithm.hpp"
#include "traccc/cuda/utils/stream_wrapper.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/track_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/messaging.hpp"
#include "traccc/utils/propagation.hpp"

// Local includes(s).
#include "traccc/examples/await_strategy.hpp"

// VecMem include(s).
#include <vecmem/containers/vector.hpp>
#include <vecmem/memory/binary_page_memory_resource.hpp>
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/stream_wrapper.hpp>

// System include(s).
#include <memory>

namespace traccc::cuda {

/// Algorithm performing the full chain of track reconstruction
///
/// At least as much as is implemented in the project at any given moment.
///
class full_chain_algorithm
    : public algorithm<edm::track_collection<default_algebra>::host(
          const edm::silicon_cell_collection::host&)>,
      public messaging {
 public:
  /// @name Type declaration(s)
  /// @{
  /// Spacepoint formation algorithm type
  using spacepoint_formation_algorithm =
      traccc::cuda::silicon_pixel_spacepoint_formation_algorithm;
  /// Clustering algorithm type
  using clustering_algorithm = traccc::cuda::clusterization_algorithm;
  /// Track finding algorithm type
  using finding_algorithm = traccc::cuda::combinatorial_kalman_filter_algorithm;
  /// Track fitting algorithm type
  using fitting_algorithm = traccc::cuda::kalman_fitting_algorithm;

  /// @}

  /// Device data shared by all instances of the algorithm
  ///
  /// It must outlive every algorithm that was constructed with it.
  ///
  struct shared_data {
    /// Constructor
    ///
    /// @param host_mr The host memory resource for the detector description
    ///                buffer
    /// @param det_descr The detector design description
    /// @param det_cond The detector conditions description
    /// @param field The magnetic field
    /// @param detector The host detector, or @c nullptr
    ///
    shared_data(vecmem::memory_resource& host_mr,
                const detector_design_description::host& det_descr,
                const detector_conditions_description::host& det_cond,
                const magnetic_field& field, const host_detector* detector);

    /// Device memory resource
    vecmem::cuda::device_memory_resource m_device_mr;
    /// B field for the track finding and fitting
    magnetic_field m_field;
    /// Detector description buffer
    detector_design_description::buffer m_device_det_descr;
    /// Detector conditions description buffer
    detector_conditions_description::buffer m_device_det_cond;
    /// Host detector
    const host_detector* m_detector;
    detector_buffer m_device_detector;
  };

  /// Algorithm constructor
  ///
  /// @param mr The memory resource to use for the intermediate and result
  ///           objects
  ///
  full_chain_algorithm(
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
      bool useGBTS = false, await_strategy = await_strategy::sync_event);

  /// Algorithm destructor
  ~full_chain_algorithm();

  /// Reconstruct track parameters in the entire detector
  ///
  /// @param cells The cells for every detector module in the event
  /// @return The track parameters reconstructed
  ///
  output_type operator()(
      const edm::silicon_cell_collection::host& cells) const override;

  /// Reconstruct track seeds in the entire detector
  ///
  /// @param cells The cells for every detector module in the event
  /// @return The track seeds reconstructed
  ///
  bound_track_parameters_collection_types::host seeding(
      const edm::silicon_cell_collection::host& cells) const;

 private:
  /// Host memory resource
  vecmem::memory_resource& m_host_mr;
  /// Pinned host memory resource
  vecmem::cuda::host_memory_resource m_pinned_host_mr;
  /// Cached pinned host memory resource
  mutable vecmem::binary_page_memory_resource m_cached_pinned_host_mr;
  /// (vecmem) CUDA stream to use
  vecmem::cuda::stream_wrapper m_vecmem_stream;
  /// (traccc) CUDA stream to use
  stream_wrapper m_stream;
  /// Device memory resource
  vecmem::cuda::device_memory_resource m_device_mr;
  /// Device caching memory resource
  mutable vecmem::binary_page_memory_resource m_cached_device_mr;
  /// (Asynchronous) Memory copy object
  mutable vecmem::cuda::async_copy m_copy;
  /// The function for awaiting asynchronous operations
  await_function_type m_await_function;

  /// B field for the track finding and fitting
  const magnetic_field& m_field;

  /// Detector description buffer
  const detector_design_description::buffer& m_device_det_descr;
  /// Detector conditions description buffer
  const detector_conditions_description::buffer& m_device_det_cond;
  /// Host detector
  const host_detector* m_detector;
  const detector_buffer& m_device_detector;

  /// @name Sub-algorithms used by this full-chain algorithm
  /// @{

  /// Clusterization algorithm
  clusterization_algorithm m_clusterization;
  /// Measurement sorting algorithm
  measurement_sorting_algorithm m_measurement_sorting;
  /// Spacepoint formation algorithm
  spacepoint_formation_algorithm m_spacepoint_formation;
  /// Seeding algorithm
  triplet_seeding_algorithm m_seeding;
  /// Seeding with GBTS algorithm
  gbts_seeding_algorithm m_gbts_seeding;
  /// Track parameter estimation algorithm
  seed_parameter_estimation_algorithm m_track_parameter_estimation;

  /// Track finding algorithm
  finding_algorithm m_finding;
  /// Track fitting algorithm
  fitting_algorithm m_fitting;

  /// @}

  /// @name Algorithm configurations
  /// @{

  bool usingGBTS;

  /// @}
};  // class full_chain_algorithm

}  // namespace traccc::cuda
