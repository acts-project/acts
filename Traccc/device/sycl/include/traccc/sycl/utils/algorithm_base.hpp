/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// SYCL library include(s).
#include "traccc/sycl/utils/queue_wrapper.hpp"

// System include(s).
#include <functional>

namespace traccc::sycl {

/// Base class for all SYCL algorithms
///
/// Holding on to data that all SYCL algorithms make use of.
///
class algorithm_base {

    public:
    /// Constructor for the algorithm base
    ///
    /// @param queue Queue to be used by the algorithm
    ///
    algorithm_base(queue_wrapper& queue);

    /// Access the queue of the algorithm
    queue_wrapper& queue() const;

    /// Get the preferred warp (sub-group) size of the device being used
    unsigned int warp_size() const;

    private:
    /// The SYCL queue to use
    std::reference_wrapper<queue_wrapper> m_queue;
    /// Preferred warp (sub-group) size of the device being used
    unsigned int m_warp_size;

};  // class algorithm_base

}  // namespace traccc::sycl
