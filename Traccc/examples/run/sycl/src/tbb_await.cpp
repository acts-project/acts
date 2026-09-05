/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/examples/sycl/tbb_await.hpp"

// TBB include(s).
#include <tbb/task.h>

namespace traccc::sycl {
void tbb_await_callback(vecmem::abstract_event& /*event*/,
                        const queue_wrapper& queue) {
  tbb::task::suspend([&queue](auto suspend_point) {
    queue.enqueue_callback(
        [suspend_point]() { tbb::task::resume(suspend_point); });
    // Exceptions thrown by a function passed to tbb::task::suspend occur
    // before the actual suspension, so the exception does not cause a deadlock.
  });
}
}  // namespace traccc::sycl
