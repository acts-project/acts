/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// SYCL library include(s).
#include "traccc/sycl/utils/queue_wrapper.hpp"

// Project include(s).
#include "traccc/device/abstract_awaitable.hpp"

// Local include(s).
#include "traccc/sycl/utils/await.hpp"

// System include(s).
#include <functional>

namespace traccc::sycl {

/// Base class for all SYCL algorithms
///
/// Holding on to data that all SYCL algorithms make use of.
///
class algorithm_base : public virtual device::abstract_awaitable {
 public:
  /// Constructor for the algorithm base
  ///
  /// @param queue Queue to be used by the algorithm
  /// @param await_func The function to use for synchronizing async operations
  ///
  algorithm_base(queue_wrapper& queue,
                 await_function_type await_func = await_sync_event);

  /// Access the queue of the algorithm
  queue_wrapper& queue() const;

  /// Get the preferred warp (sub-group) size of the device being used
  unsigned int warp_size() const;

  /// Synchronize an event related to asynchronous operations
  /// @param event The event to synchronize
  ///
  void await(vecmem::abstract_event& event) const override;

 private:
  /// The SYCL queue to use
  std::reference_wrapper<queue_wrapper> m_queue;
  /// Preferred warp (sub-group) size of the device being used
  unsigned int m_warp_size;
  /// The function to use for synchronizing async operations
  await_function_type m_await_func;

};  // class algorithm_base

}  // namespace traccc::sycl
