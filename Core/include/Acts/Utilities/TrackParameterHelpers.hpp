// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"

#include <cmath>

namespace Acts {
namespace TrackParameterHelpers {

template <typename parameters_t>
inline auto charge(const parameters_t &parameters, double absQ) {
  assert(absQ >= 0);
  return std::copysign(parameters.qop(), absQ);
}

template <typename parameters_t>
inline auto absoluteMomentum(const parameters_t &parameters, double absQ) {
  assert(absQ >= 0);
  return absQ / std::abs(parameters.qop());
}

template <typename parameters_t>
inline auto transverseMomentum(const parameters_t &parameters, double absQ) {
  return std::sin(parameters.theta()) * absoluteMomentum(parameters, absQ);
}

template <typename parameters_t>
inline auto momentum(const parameters_t &parameters, double absQ) {
  return parameters.direction() * absoluteMomentum(parameters, absQ);
}

}  // namespace TrackParameterHelpers
}  // namespace Acts
