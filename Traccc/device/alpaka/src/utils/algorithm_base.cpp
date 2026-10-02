// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/alpaka/utils/algorithm_base.hpp"

#include "get_queue.hpp"

// Alpaka include(s).
#include <alpaka/alpaka.hpp>

namespace traccc::alpaka {

algorithm_base::algorithm_base(alpaka::queue& q, await_function_type await_func)
    : m_queue(q),
      m_warp_size(static_cast<unsigned int>(::alpaka::getPreferredWarpSize(
          ::alpaka::getDev(details::get_queue(q))))),
      m_await_func(std::move(await_func)) {}

alpaka::queue& algorithm_base::queue() const {
  return m_queue.get();
}

unsigned int algorithm_base::warp_size() const {
  return m_warp_size;
}

void algorithm_base::await(vecmem::abstract_event& event) const {
  m_await_func(event, m_queue.get());
}

}  // namespace traccc::alpaka
