/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/examples/alpaka/tbb_await.hpp"

// System include(s).
#include <exception>

// TBB include(s).
#include <tbb/task.h>

namespace traccc::alpaka {

void tbb_await_callback(vecmem::abstract_event& /*event*/, const queue& queue) {
  auto suspend_point =
      tbb::task::suspend_point{};  // suspension point address must remain valid
  // when resumption callback is called
  tbb::task::suspend([&queue, &suspend_point](auto tag) {
    suspend_point = tag;
    try {
      queue.enqueue_callback(
          [&suspend_point]() { tbb::task::resume(suspend_point); });
    } catch (const std::exception& e) {
      // resume immediately in case of an error to avoid deadlock
      tbb::task::resume(suspend_point);
      throw;
    }
  });
}
}  // namespace traccc::alpaka
