/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/alpaka/utils/queue.hpp"

// VecMem include(s).
#include <vecmem/utils/abstract_event.hpp>

// System include(s).
#include <functional>

namespace traccc::alpaka {

/// Type of the function used to synchronize async operations.
using await_function_type =
    std::function<void(vecmem::abstract_event&, const queue&)>;

/// Synchronize an event by waiting for it to complete.
///
/// @param event The event to synchronize.
/// @param queue The Alpaka queue, unused.
///
void await_sync_event(vecmem::abstract_event& event, const queue& queue);

}  // namespace traccc::alpaka
