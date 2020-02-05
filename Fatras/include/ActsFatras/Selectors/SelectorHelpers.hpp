// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <climits>

#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {

/// Select all objects with an extracted value equal or larger than the cut.
template <typename cast_t>
struct Min {
  double valMin = 0.;

  template <typename detector_t, typename T>
  bool operator()(const detector_t &, const T &thing) const {
    return (valMin <= cast_t()(thing));
  }
};

/// Select all objects with an extracted value below the cut.
template <typename cast_t>
struct Max {
  double valMax = std::numeric_limits<double>::max();

  template <typename detector_t, typename T>
  bool operator()(const detector_t &, const T &thing) const {
    return (cast_t()(thing) < valMax);
  }
};

/// Select all objects with an extracted value within the range.
///
/// The range is defined as the left, half-open interval within the cuts.
template <typename cast_t>
struct Range {
  double valMin = std::numeric_limits<double>::lowest();
  double valMax = std::numeric_limits<double>::max();

  template <typename detector_t, typename T>
  bool operator()(const detector_t &, const T &thing) const {
    const auto val = cast_t()(thing);
    return ((valMin <= val) and (val < valMax));
  }
};

}  // namespace ActsFatras
