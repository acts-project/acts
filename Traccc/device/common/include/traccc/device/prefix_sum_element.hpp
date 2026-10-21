// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "traccc/utils/pair.hpp"

namespace traccc::device {

/// Type for the individual elements in a prefix sum vector
typedef traccc::pair<unsigned int, unsigned int> prefix_sum_element_t;

}  // namespace traccc::device
