// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/examples/alpaka/tbb_await.hpp"

// TBB include(s).
#include <tbb/task.h>

namespace traccc::alpaka {

void tbb_await_callback(vecmem::abstract_event& event, const queue& queue) {
  event.ignore();  // ignore the event, as it is not needed for resumption
  tbb::task::suspend([&queue](auto suspend_point) {
    queue.enqueue_callback(
        [suspend_point]() { tbb::task::resume(suspend_point); });
    // Exceptions thrown by a function passed to tbb::task::suspend occur
    // before the actual suspension, so the exception does not cause a deadlock.
  });
}
}  // namespace traccc::alpaka
