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
#include "traccc/seeding/device/triplet_seeding_algorithm.hpp"

namespace traccc::sycl {

/// Main algorithm for performing the track seeding using oneAPI/SYCL
class triplet_seeding_algorithm : public device::triplet_seeding_algorithm,
                                  public sycl::algorithm_base {
 public:
  /// Constructor for the seed finding algorithm
  ///
  /// @param finder_config The seed finding configuration
  /// @param grid_config The spacepoint grid configuration
  /// @param filter_config The seed filtering configuration
  /// @param mr is a struct of memory resources (shared or host & device)
  /// @param copy The copy object to use for copying data between device
  ///             and host memory blocks
  /// @param queue The SYCL queue to work with
  /// @param logger The logger instance to use
  /// @param await_func The function used to synchronize events
  ///
  triplet_seeding_algorithm(
      const seedfinder_config& finder_config,
      const spacepoint_grid_config& grid_config,
      const seedfilter_config& filter_config, const traccc::memory_resource& mr,
      const vecmem::copy& copy, queue_wrapper& queue,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone(),
      await_function_type await_func = await_sync_event);

 private:
  /// The function used to synchronize events.
  await_function_type m_await_func;

  /// @name Function(s) inherited from @c traccc::device::seeding_algorithm
  /// @{

  /// Spacepoint grid capacity counting kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  void count_grid_capacities_kernel(
      const count_grid_capacities_kernel_payload& payload) const override;

  /// Spacepoint grid population kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  void populate_grid_kernel(
      const populate_grid_kernel_payload& payload) const override;

  /// Doublet counting kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  void count_doublets_kernel(
      const count_doublets_kernel_payload& payload) const override;

  /// Doublet finding kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  void find_doublets_kernel(
      const find_doublets_kernel_payload& payload) const override;

  /// Triplet counting kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  void count_triplets_kernel(
      const count_triplets_kernel_payload& payload) const override;

  /// Triplet count reduction kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  void triplet_counts_reduction_kernel(
      const triplet_counts_reduction_kernel_payload& payload) const override;

  /// Triplet finding kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  void find_triplets_kernel(
      const find_triplets_kernel_payload& payload) const override;

  /// Triplet weight updater/filler kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  void update_triplet_weights_kernel(
      const update_triplet_weights_kernel_payload& payload) const override;

  /// Seed selection/filling kernel launcher
  ///
  /// @param payload The payload for the kernel
  ///
  void select_seeds_kernel(
      const select_seeds_kernel_payload& payload) const override;

  /// Synchronize an event related to asynchronous operations
  /// @param event The event to synchronize
  ///
  void await(vecmem::abstract_event& event) const override;

};  // class triplet_seeding_algorithm

}  // namespace traccc::sycl
