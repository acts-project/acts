/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/cuda/utils/algorithm_base.hpp"

#include "../utils/utils.hpp"

namespace traccc::cuda {

algorithm_base::algorithm_base(const stream_wrapper& str,
                               await_function_type await_func)
    : m_stream(str),
      m_warp_size(details::get_warp_size(str.device())),
      m_await_func(std::move(await_func)) {}

const stream_wrapper& algorithm_base::stream() const {
  return m_stream;
}

unsigned int algorithm_base::warp_size() const {
  return m_warp_size;
}

void algorithm_base::await(vecmem::abstract_event& event) const {
  m_await_func(event, m_stream);
}

}  // namespace traccc::cuda
