/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/cuda/utils/await.hpp"
#include "traccc/cuda/utils/stream_wrapper.hpp"

// Project include(s).
#include "traccc/device/abstract_awaitable.hpp"

namespace traccc::cuda {

/// Base class for all CUDA algorithms
///
/// Holding on to data that all CUDA algorithms make use of.
///
class algorithm_base : public virtual device::abstract_awaitable {
 public:
  /// Constructor for the algorithm base
  ///
  /// @param str The CUDA stream to perform all operations on
  /// @param await_func The function to use for synchronizing async operations
  ///
  explicit algorithm_base(const stream_wrapper& str,
                          await_function_type await_func = await_sync_event);

  /// Get the CUDA stream of the algorithm
  const stream_wrapper& stream() const;
  /// Get the warp size of the GPU being used
  unsigned int warp_size() const;

  /// Synchronize an event related to asynchronous operations
  /// @param event The event to synchronize
  ///
  void await(vecmem::abstract_event& event) const override;

 private:
  /// The CUDA stream to use
  stream_wrapper m_stream;
  /// Warp size of the GPU being used
  unsigned int m_warp_size;
  /// The function to use for synchronizing async operations
  await_function_type m_await_func;

};  // class algorithm_base

}  // namespace traccc::cuda
