/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/alpaka/utils/algorithm_base.hpp"

#include "get_queue.hpp"

// Alpaka include(s).
#include <alpaka/alpaka.hpp>

namespace traccc::alpaka {

algorithm_base::algorithm_base(alpaka::queue& q)
    : m_queue(q),
      m_warp_size(static_cast<unsigned int>(::alpaka::getPreferredWarpSize(
          ::alpaka::getDev(details::get_queue(q))))) {}

alpaka::queue& algorithm_base::queue() const {

    return m_queue.get();
}

unsigned int algorithm_base::warp_size() const {

    return m_warp_size;
}

}  // namespace traccc::alpaka
