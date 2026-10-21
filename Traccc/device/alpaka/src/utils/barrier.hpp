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

namespace traccc::alpaka {

template <typename TAcc>
struct barrier {
  ALPAKA_FN_INLINE ALPAKA_FN_ACC barrier(const TAcc* acc) : m_acc(acc) {};

  ALPAKA_FN_ACC
  void blockBarrier() const { ::alpaka::syncBlockThreads(*m_acc); }

  ALPAKA_FN_ACC
  bool blockOr(bool predicate) const {
    return ::alpaka::syncBlockThreadsPredicate<::alpaka::BlockOr>(*m_acc,
                                                                  predicate);
  }

  ALPAKA_FN_ACC
  bool blockAnd(bool predicate) const {
    return ::alpaka::syncBlockThreadsPredicate<::alpaka::BlockAnd>(*m_acc,
                                                                   predicate);
  }

  ALPAKA_FN_ACC
  bool blockCount(int threadCount) const {
    return ::alpaka::syncBlockThreadsPredicate<::alpaka::BlockCount>(
        *m_acc, threadCount);
  }

 private:
  const TAcc* m_acc;
};

}  // namespace traccc::alpaka
