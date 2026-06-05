/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "get_queue.hpp"

// System include(s).
#include <cassert>

namespace traccc::alpaka::details {

Queue& get_queue(queue& q) {

    assert(q.alpakaQueue() != nullptr);
    return *(reinterpret_cast<Queue*>(q.alpakaQueue()));
}

const Queue& get_queue(const queue& q) {

    assert(q.alpakaQueue() != nullptr);
    return *(reinterpret_cast<const Queue*>(q.alpakaQueue()));
}

}  // namespace traccc::alpaka::details
