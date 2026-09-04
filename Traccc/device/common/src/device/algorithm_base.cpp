/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/device/algorithm_base.hpp"

namespace traccc::device {

algorithm_base::algorithm_base(const memory_resource& mr,
                               const vecmem::copy& copy)
    : m_mr{mr}, m_copy{copy} {}

const memory_resource& algorithm_base::mr() const {
  return m_mr;
}

const vecmem::copy& algorithm_base::copy() const {
  return m_copy.get();
}

}  // namespace traccc::device
