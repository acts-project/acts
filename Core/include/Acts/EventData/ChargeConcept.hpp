// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/HashedString.hpp"

#include <any>
#include <type_traits>

#if defined(__cpp_concepts)
#include <concepts>

namespace Acts {

template <typename C>
concept ChargeConcept = requires(C c, float f, double d) {
  {C{f}};

  { c == c } -> std::same_as<bool>;
  { c != c } -> std::same_as<bool>;

  { c.absQ() } -> std::same_as<float>;
  { c.extractCharge(d) } -> std::same_as<float>;
  { c.extractMomentum(d) } -> std::same_as<double>;
  { c.qOverP(d, d) } -> std::same_as<double>;
};

}  // namespace Acts

#endif
