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

struct [[maybe_unused]] cell_module_projection {
  template <typename T>
  TRACCC_HOST_DEVICE auto operator()(const edm::silicon_cell<T>& c) const {
    return c.module_index();
  }
};

}  // namespace traccc
