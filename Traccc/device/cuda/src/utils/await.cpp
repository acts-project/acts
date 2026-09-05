/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/cuda/utils/await.hpp"

// System include(s).
#include <concepts>

namespace traccc::cuda {
void await_sync_event(vecmem::abstract_event& event,
                      const stream_wrapper& /*stream*/) {
  event.wait();
}

static_assert(
    std::constructible_from<await_function_type, decltype(await_sync_event)>,
    "await_sync_event should be of type await_function_type");

}  // namespace traccc::cuda
