// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "utils.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/device/concepts/thread_id.hpp"

// Alpaka include(s).
#include <alpaka/alpaka.hpp>

namespace traccc::alpaka::details {

/// An Alpaka thread identifier type
template <typename TAcc>
struct thread_id1 {
  ALPAKA_FN_INLINE ALPAKA_FN_ACC explicit thread_id1(const TAcc& acc)
      : m_acc(acc) {}

  unsigned int inline ALPAKA_FN_ACC getLocalThreadId() const {
    return static_cast<unsigned int>(
        ::alpaka::getIdx<::alpaka::Block, ::alpaka::Threads>(m_acc)[0u]);
  }

  unsigned int inline ALPAKA_FN_ACC getLocalThreadIdX() const {
    return getLocalThreadId();
  }

  unsigned int inline ALPAKA_FN_ACC getGlobalThreadId() const {
    return getLocalThreadId() + getBlockIdX() * getBlockDimX();
  }

  unsigned int inline ALPAKA_FN_ACC getGlobalThreadIdX() const {
    return getLocalThreadId() + getBlockIdX() * getBlockDimX();
  }

  unsigned int inline ALPAKA_FN_ACC getBlockIdX() const {
    return static_cast<unsigned int>(
        ::alpaka::getIdx<::alpaka::Grid, ::alpaka::Blocks>(m_acc)[0u]);
  }

  unsigned int inline ALPAKA_FN_ACC getBlockDimX() const {
    return static_cast<unsigned int>(
        ::alpaka::getWorkDiv<::alpaka::Block, ::alpaka::Threads>(m_acc)[0u]);
  }

  unsigned int inline ALPAKA_FN_ACC getGridDimX() const {
    return static_cast<unsigned int>(
        ::alpaka::getWorkDiv<::alpaka::Grid, ::alpaka::Blocks>(m_acc)[0u]);
  }

 private:
  const TAcc& m_acc;
};

/// Verify that @c traccc::alpaka::details::thread_id1 fulfills the
/// @c traccc::device::concepts::thread_id1 concept.
static_assert(traccc::device::concepts::thread_id1<thread_id1<Acc>>);

}  // namespace traccc::alpaka::details
