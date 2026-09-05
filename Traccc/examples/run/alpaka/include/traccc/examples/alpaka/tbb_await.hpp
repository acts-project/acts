/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/alpaka/utils/queue.hpp"

// Vecmem include(s).
#include <vecmem/utils/abstract_event.hpp>

namespace traccc::alpaka {

/// Suspend execution of a TBB task until work on the Alpaka queue has
/// completed. A callback-based mechanism is used to handle the completion
/// notification.
///
/// @param event The event to synchronize, unused.
/// @param queue The Alpaka queue.
///
/// @note Should be called only from a TBB task.
///
void tbb_await_callback(vecmem::abstract_event& /*event*/, const queue& queue);
}  // namespace traccc::alpaka
