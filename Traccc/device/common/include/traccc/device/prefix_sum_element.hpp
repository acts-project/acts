/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/utils/pair.hpp"

namespace traccc::device {

/// Type for the individual elements in a prefix sum vector
typedef traccc::pair<unsigned int, unsigned int> prefix_sum_element_t;

}  // namespace traccc::device
