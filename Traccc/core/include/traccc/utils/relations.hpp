// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/silicon_cell_collection.hpp"

namespace traccc {
struct [[maybe_unused]] channel0_major_cell_order_relation {
  template <typename T1, typename T2>
  [[maybe_unused]] TRACCC_HOST_DEVICE bool operator()(
      const edm::silicon_cell<T1>& a, const edm::silicon_cell<T2>& b) const {
    if (a.module_index() == b.module_index()) {
      if (a.channel1() == b.channel1()) {
        return a.channel0() <= b.channel0();
      } else {
        return a.channel1() <= b.channel1();
      }
    } else {
      return true;
    }
  }
};
}  // namespace traccc
