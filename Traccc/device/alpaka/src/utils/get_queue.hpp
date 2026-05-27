/**
 * traccc library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/alpaka/utils/queue.hpp"
#include "utils.hpp"

namespace traccc::alpaka::details {

/// Helper function for getting a @c Queue out of @c queue (non-const)
Queue& get_queue(queue& q);

/// Helper function for getting a @c Queue out of @c queue (const)
const Queue& get_queue(const queue& q);

}  // namespace traccc::alpaka::details
