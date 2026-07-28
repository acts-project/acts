/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/sycl/utils/algorithm_base.hpp"
#include "traccc/sycl/utils/await.hpp"

// Project include(s).
#include "traccc/seeding/device/seed_parameter_estimation_algorithm.hpp"

namespace traccc::sycl {

/// Algorithm for estimating track parameters for seeds using oneAPI/SYCL
struct seed_parameter_estimation_algorithm
    : public device::seed_parameter_estimation_algorithm,
      public sycl::algorithm_base {
 public:
  /// Constructor for track_params_estimation
  ///
  /// @param config The track parameter estimation configuration
  /// @param mr is a struct of memory resources (shared or host & device)
  /// @param copy The copy object to use for copying data between device
  ///             and host memory blocks
  /// @param queue The SYCL queue to work with
  /// @param logger The logger instance to use
  /// @param await_func The function used to synchronize events
  ///
  seed_parameter_estimation_algorithm(
      const track_params_estimation_config& config,
      const traccc::memory_resource& mr, const vecmem::copy& copy,
      queue_wrapper& queue,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone(),
      await_function_type await_func = await_sync_event);

 private:
  /// The function used to synchronize events.
  await_function_type m_await_func;

  /// @name Function(s) inherited from
  /// @c traccc::device::seed_parameter_estimation_algorithm
  /// @{

  /// Seed parameter estimation kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  void estimate_seed_params_kernel(
      const struct estimate_seed_params_kernel_payload& payload) const override;

  /// Synchronize an event related to asynchronous operations
  /// @param event The event to synchronize
  ///
  void await(vecmem::abstract_event& event) const override;

  /// @}

};  // struct seed_parameter_estimation_algorithm

}  // namespace traccc::sycl
