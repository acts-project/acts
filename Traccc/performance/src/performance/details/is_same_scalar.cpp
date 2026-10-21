// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Library include(s).
#include "traccc/performance/details/is_same_scalar.hpp"

// System include(s).
#include <cmath>

namespace traccc::details {

bool is_same_scalar(scalar lhs, scalar rhs, scalar unc) {
  // The difference of the two values is meant to be smaller than
  // their average times the uncertainty.
  return (std::abs(lhs - rhs) <=
          (unc * ((std::abs(lhs) + std::abs(rhs)) / 2.f)));
}

}  // namespace traccc::details
