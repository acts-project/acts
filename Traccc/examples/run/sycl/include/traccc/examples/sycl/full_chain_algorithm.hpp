/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/edm/silicon_cell_collection.hpp"
#include "traccc/edm/track_parameters.hpp"
#include "traccc/gbts_seeding/gbts_seeding_config.hpp"
#include "traccc/geometry/detector.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/geometry/detector_conditions_description.hpp"
#include "traccc/geometry/detector_design_description.hpp"
#include "traccc/geometry/host_detector.hpp"
#include "traccc/sycl/clusterization/clusterization_algorithm.hpp"
#include "traccc/sycl/clusterization/measurement_sorting_algorithm.hpp"
#include "traccc/sycl/finding/combinatorial_kalman_filter_algorithm.hpp"
#include "traccc/sycl/fitting/kalman_fitting_algorithm.hpp"
#include "traccc/sycl/gbts_seeding/gbts_seeding_algorithm.hpp"
#include "traccc/sycl/seeding/seed_parameter_estimation_algorithm.hpp"
#include "traccc/sycl/seeding/silicon_pixel_spacepoint_formation_algorithm.hpp"
#include "traccc/sycl/seeding/triplet_seeding_algorithm.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/propagation.hpp"

// Local includes(s).
#include "traccc/examples/await_strategy.hpp"

// VecMem include(s).
#include <vecmem/memory/binary_page_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/memory/sycl/device_memory_resource.hpp>
#include <vecmem/memory/sycl/host_memory_resource.hpp>
#include <vecmem/utils/sycl/async_copy.hpp>

// System include(s).
#include <memory>

namespace traccc::sycl {
namespace details {
/// Internal data type used by @c traccc::sycl::full_chain_algorithm
struct full_chain_algorithm_data;
}  // namespace details

/// Algorithm performing the full chain of track reconstruction
///
/// At least as much as is implemented in the project at any given moment.
///
class full_chain_algorithm
    : public algorithm<edm::track_collection<default_algebra>::host(
          const edm::silicon_cell_collection::host&)>,
      public messaging {
 public:
  /// @name (For now dummy...) Type declaration(s)
  /// @{

  /// Spacepoint formation algorithm type
  using spacepoint_formation_algorithm =
      traccc::sycl::silicon_pixel_spacepoint_formation_algorithm;
  /// Clustering algorithm type
  using clustering_algorithm = clusterization_algorithm;
  /// Track finding algorithm type
  using finding_algorithm = traccc::sycl::combinatorial_kalman_filter_algorithm;
  /// Track fitting algorithm type
  using fitting_algorithm = traccc::sycl::kalman_fitting_algorithm;

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
    /// Destructor
    ~shared_data();

    /// Private data object
    std::unique_ptr<details::full_chain_algorithm_data> m_data;
    /// Device memory resource
    vecmem::sycl::device_memory_resource m_device_mr;
    /// B field for the track finding and fitting
    magnetic_field m_field;
    /// Detector description buffer
    detector_design_description::buffer m_device_det_descr;
    detector_conditions_description::buffer m_device_det_cond;
    /// Host detector
    const host_detector* m_detector;
    /// Buffer holding the detector's payload on the device
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
      bool useGBTS = false,
      await_strategy await_mode = await_strategy::sync_event);

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
  /// Private data object
  std::unique_ptr<details::full_chain_algorithm_data> m_data;
  /// Host memory resource
  std::reference_wrapper<vecmem::memory_resource> m_host_mr;
  /// Pinned host memory resource
  vecmem::sycl::host_memory_resource m_pinned_host_mr;
  /// Cached pinned host memory resource
  mutable vecmem::binary_page_memory_resource m_cached_pinned_host_mr;
  /// Device memory resource
  vecmem::sycl::device_memory_resource m_device_mr;
  /// Device caching memory resource
  mutable vecmem::binary_page_memory_resource m_cached_device_mr;
  /// Memory copy object
  mutable vecmem::sycl::async_copy m_copy;
  /// The function for awaiting asynchronous operations
  await_function_type m_await_function;

  /// B field for the track finding and fitting
  const magnetic_field& m_field;

  /// Detector description buffer
  const detector_design_description::buffer& m_device_det_descr;
  const detector_conditions_description::buffer& m_device_det_cond;

  /// Host detector
  const host_detector* m_detector;
  /// Buffer holding the detector's payload on the device
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

  /// @}

  /// @name Algorithm configurations
  /// @{

  bool usingGBTS;

  /// @}
};  // class full_chain_algorithm

}  // namespace traccc::sycl
