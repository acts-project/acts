/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Project include(s).
#include "traccc/sycl/utils/queue_wrapper.hpp"

namespace traccc::sycl {

queue_wrapper::queue_wrapper(void* queue) : m_queue(queue) {}

void* queue_wrapper::queue() {
    return m_queue;
}

const void* queue_wrapper::queue() const {
    return m_queue;
}

}  // namespace traccc::sycl
