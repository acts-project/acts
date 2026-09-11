/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/sycl/utils/queue_wrapper.hpp"

// SYCL include(s).
#include <sycl/sycl.hpp>

namespace traccc::sycl::details {

/// Helper function for getting a @c sycl::queue out of
/// @c traccc::sycl::queue_wrapper (non-const)
::sycl::queue& get_queue(traccc::sycl::queue_wrapper& queue);

/// Helper function for getting a @c sycl::queue out of
/// @c traccc::sycl::queue_wrapper (const)
const ::sycl::queue& get_queue(const traccc::sycl::queue_wrapper& queue);

}  // namespace traccc::sycl::details
