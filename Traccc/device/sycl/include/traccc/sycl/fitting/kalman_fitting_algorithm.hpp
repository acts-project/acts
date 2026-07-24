/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// SYCL library include(s).
#include "traccc/sycl/utils/algorithm_base.hpp"
#include "traccc/sycl/utils/queue_wrapper.hpp"

// Project include(s).
#include "traccc/fitting/device/kalman_fitting_algorithm.hpp"

namespace traccc::sycl {

/// Kalman filter based track fitting algorithm using SYCL
class kalman_fitting_algorithm : public device::kalman_fitting_algorithm,
                                 public sycl::algorithm_base {
 public:
  /// Constructor with the algorithm's configuration
  ///
  /// @param config The configuration object
  /// @param mr     The memory resource(s) used by the algorithm
  /// @param copy   The copy object used by the algorithm
  /// @param queue  The SYCL queue used by the algorithm
  /// @param logger The logger used by the algorithm
  ///
  kalman_fitting_algorithm(
      const config_type& config, const traccc::memory_resource& mr,
      const vecmem::copy& copy, queue_wrapper& queue,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone());

 private:
  /// @name Function(s) implemented from @c device::kalman_fitting_algorithm
  /// @{

  /// Prepare a buffer with the index order with which to fit the tracks
  ///
  /// @param[in] tracks The tracks to be fitted
  /// @param[out] track_sort_keys Buffer storing temporary sorting keys
  /// @param[out] track_indices The buffer to write the fitting order into
  ///
  void prepare_track_fit_order(
      const edm::track_collection<default_algebra>::const_view& tracks,
      vecmem::data::vector_view<device::sort_key>& track_sort_keys,
      vecmem::data::vector_view<unsigned int>& track_indices) const override;

  /// Kernel to prepare the fitting payloads
  ///
  /// @param payload The payload for the kernel(s)
  ///
  void fit_prelude_kernel(
      const device::fit_prelude_payload& payload) const override;

  /// Function preparing the fitting payload
  ///
  /// @param det             The detector buffer to prepare the payload for
  /// @param field           The magnetic field to prepare the payload for
  /// @param n_surfaces      The number of surfaces for each track to be
  ///                        fitted
  /// @param payload         The (non-templated) payload for the kernel(s)
  ///
  /// @return The prepared payload for the fitting kernel(s)
  ///
  fit_payload prepare_fit_payload(
      const detector_buffer& det, const magnetic_field& field,
      const std::vector<unsigned int>& n_surfaces,
      const device::fit_payload& payload) const override;

  /// Function launching the "forward fitting" kernel(s)
  ///
  /// @param config The fitting configuration
  /// @param payload The payload for the fitting kernel(s)
  ///
  void fit_forward_kernel(const fitting_config& config,
                          const fit_payload& payload) const override;

  /// Function launching the "backward fitting" kernel(s)
  ///
  /// @param config The fitting configuration
  /// @param payload The payload for the fitting kernel(s)
  ///
  void fit_backward_kernel(const fitting_config& config,
                           const fit_payload& payload) const override;

  /// @}

};  // class kalman_fitting_algorithm

}  // namespace traccc::sycl
