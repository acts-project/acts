/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/alpaka/utils/await.hpp"

// System include(s).
#include <concepts>

namespace traccc::alpaka {

void await_sync_event(vecmem::abstract_event& event, const queue& /*queue*/) {
  event.wait();
}

static_assert(
    std::constructible_from<await_function_type, decltype(await_sync_event)>,
    "await_sync_event should be of type await_function_type");

}  // namespace traccc::alpaka
