// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/CurvilinearTrackParameters.hpp"
#include "Acts/EventData/FreeTrackParameters.hpp"

#include <cmath>

namespace Acts {

template <typename parameters_t>
inline auto charge(const parameters_t &parameters, double absQ) {
  return std::copysign(parameters.template get<eBoundQOverP>(), absQ);
}

template <typename parameters_t>
inline auto absoluteMomentum(const parameters_t &parameters, double absQ) {
  return absQ / std::abs(parameters.template get<eBoundQOverP>());
}

template <typename parameters_t>
inline auto transverseMomentum(const parameters_t &parameters, double absQ) {
  return std::sin(parameters.template get<eBoundTheta>()) *
         absoluteMomentum(parameters, absQ);
}

template <typename parameters_t>
inline auto momentum(const parameters_t &parameters, double absQ) {
  return parameters.direction() * absoluteMomentum(parameters, absQ);
}

}  // namespace Acts
