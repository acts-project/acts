// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <cmath>

namespace Acts {
namespace TrackParameterHelpers {

template <typename parameters_t>
inline auto charge(const parameters_t &parameters, double absCharge) {
  assert(absCharge >= 0 && parameters.qop() != 0);
  return (parameters.qop() >= 0 ? +1 : -1) * absCharge;
}

template <typename parameters_t>
inline auto absoluteMomentum(const parameters_t &parameters, double absCharge) {
  auto q = charge(parameters, absCharge);
  return (q == 0 ? 1 : q) / parameters.qop();
}

template <typename parameters_t>
inline auto transverseMomentum(const parameters_t &parameters,
                               double absCharge) {
  return std::sin(parameters.theta()) * absoluteMomentum(parameters, absCharge);
}

template <typename parameters_t>
inline Vector3 momentum(const parameters_t &parameters, double absCharge) {
  return absoluteMomentum(parameters, absCharge) * parameters.direction();
}

}  // namespace TrackParameterHelpers
}  // namespace Acts
