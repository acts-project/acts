// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
