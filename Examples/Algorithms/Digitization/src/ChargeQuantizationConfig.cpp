// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/ChargeQuantizationConfig.hpp"

#include <cmath>
#include <stdexcept>

namespace ActsExamples {

void ChargeQuantizationConfig::validate() const {
  if (bits < 1 || bits > 16 || !std::isfinite(offset) || !std::isfinite(step) ||
      step <= 0.) {
    throw std::invalid_argument(
        "Charge quantization requires 1..16 bits, finite offset and positive "
        "finite step");
  }
  const unsigned int maxCode = (1u << bits) - 1u;
  const double upper = offset + maxCode * step;
  if (!std::isfinite(upper) || upper <= offset ||
      (maxCode > 1 && (offset + step <= offset || upper - step >= upper))) {
    throw std::invalid_argument(
        "Charge quantization bins are not representable");
  }
}

double ChargeQuantizationConfig::quantize(double activation) const {
  if (!enabled) {
    return activation;
  }
  if (!std::isfinite(activation)) {
    throw std::invalid_argument("Cannot quantize a non-finite activation");
  }
  // Search the implicit uniform boundaries instead of dividing and flooring:
  // division can round an activation just below a boundary onto that boundary.
  unsigned int first = 1;
  unsigned int last = (1u << bits) - 1u;
  while (first < last) {
    const unsigned int middle = first + (last - first) / 2;
    if (activation < offset + middle * step) {
      last = middle;
    } else {
      first = middle + 1;
    }
  }
  return first;
}

}  // namespace ActsExamples
