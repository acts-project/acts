/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/cuda/utils/stream_wrapper.hpp"

// Vecmem include(s).
#include <vecmem/utils/abstract_event.hpp>

namespace traccc::cuda {

/// Suspend execution of a TBB task until work on the CUDA stream has completed.
/// A callback-based mechanism is used to handle the completion notification.
///
/// @param event The event to synchronize, unused.
/// @param stream The CUDA stream.
///
/// @note Should be called only from a TBB task.
///
void tbb_await_callback(vecmem::abstract_event& /*event*/,
                        const stream_wrapper& stream);
}  // namespace traccc::cuda
