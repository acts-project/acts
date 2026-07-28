/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/alpaka/utils/algorithm_base.hpp"
#include "traccc/alpaka/utils/await.hpp"

// Project include(s).
#include "traccc/seeding/device/silicon_pixel_spacepoint_formation_algorithm.hpp"

namespace traccc::alpaka {

/// Algorithm forming space points out of measurements
///
/// This algorithm performs the local-to-global transformation of the 2D
/// measurements made on every detector module, into 3D spacepoint coordinates.
///
class silicon_pixel_spacepoint_formation_algorithm
    : public device::silicon_pixel_spacepoint_formation_algorithm,
      public alpaka::algorithm_base {
 public:
  /// Constructor for the spacepoint formation algorithm
  ///
  /// @param mr The memory resource(s) to use in the algorithm
  /// @param copy The copy object to use for copying data between device
  ///             and host memory blocks
  /// @param q The Alpaka queue to use
  /// @param logger The logger instance to use
  /// @param await_func The function used to synchronize events
  ///
  silicon_pixel_spacepoint_formation_algorithm(
      const traccc::memory_resource& mr, const vecmem::copy& copy,
      alpaka::queue& q,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone(),
      await_function_type await_func = await_sync_event);

 private:
  /// The function used to synchronize events.
  await_function_type m_await_func;

  /// @name Function(s) inherited from
  /// @c traccc::device::silicon_pixel_spacepoint_formation_algorithm
  /// @{

  /// Launch the spacepoint formation kernel
  ///
  /// @param payload The payload for the kernel
  ///
  void form_spacepoints_kernel(
      const form_spacepoints_kernel_payload& payload) const override;

  /// Synchronize an event related to asynchronous operations
  /// @param event The event to synchronize
  ///
  void await(vecmem::abstract_event& event) const override;

  /// @}

};  // class silicon_pixel_spacepoint_formation_algorithm

}  // namespace traccc::alpaka
