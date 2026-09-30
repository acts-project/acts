/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/alpaka/utils/await.hpp"
#include "traccc/alpaka/utils/queue.hpp"

// Project include(s).
#include "traccc/device/abstract_awaitable.hpp"

// System include(s).
#include <cstddef>
#include <functional>

namespace traccc::alpaka {

/// Base class for all Alpaka algorithms
///
/// Holding on to data that all Alpaka algorithms make use of.
///
class algorithm_base : public virtual device::abstract_awaitable {
 public:
  /// Constructor
  ///
  /// @param q The Alpaka queue to perform the operations in
  /// @param await_func The function to use for synchronizing async operations
  ///
  explicit algorithm_base(alpaka::queue& q,
                          await_function_type await_func = await_sync_event);

  /// Get the Alpaka queue of the algorithm
  alpaka::queue& queue() const;

  /// Get the preferred warp size of the device being used
  unsigned int warp_size() const;

  /// Synchronize an event related to asynchronous operations
  /// @param event The event to synchronize
  ///
  void await(vecmem::abstract_event& event) const override;

 private:
  /// The Alpaka queue to use
  std::reference_wrapper<alpaka::queue> m_queue;
  /// Preferred warp size of the device being used
  unsigned int m_warp_size;
  /// The function to use for synchronizing async operations
  await_function_type m_await_func;
};  // class algorithm_base

}  // namespace traccc::alpaka
