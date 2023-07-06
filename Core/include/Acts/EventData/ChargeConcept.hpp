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

#if defined(ACTS_CONCEPTS_SUPPORTED)
#include <concepts>

namespace Acts {

template <typename C>
concept ChargeConcept = requires(C c, float f, double d) {
  {C{f}};
  { c.absQ() } -> std::same_as<float>;

  { c.extractCharge(f) } -> std::convertible_to<float>;
  { c.extractCharge(d) } -> std::convertible_to<float>;

  { c.extractMomentum(f) } -> std::convertible_to<float>;
  { c.extractMomentum(d) } -> std::convertible_to<float>;

  { c.qOverP(f, f) } -> std::same_as<float>;
  { c.qOverP(d, d) } -> std::same_as<double>;

  { c == c } -> std::same_as<bool>;
  { c != c } -> std::same_as<bool>;
};

}  // namespace Acts

#endif
