// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <climits>

namespace ActsFatras {

// static selectors
template <typename cast_t>
struct Min {
  cast_t cast;
  double valMin = 0.;

  /// Return true for all particles with transverse momentum
  /// bigger than the specified minimum value
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    double val = cast(particle);
    return (val >= valMin);
  }
};

template <typename cast_t>
struct Max {
  cast_t cast;
  double valMax = std::numeric_limits<double>::max();

  /// Return true for all particles with transverse momentum
  /// bigger than the specified minimum value
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    double val = cast(particle);
    return (val <= valMax);
  }
};

template <typename cast_t>
struct Range {
  cast_t cast;
  double valMin = 0.;
  double valMax = std::numeric_limits<double>::max();

  /// Return true for all particles with transverse momentum
  /// within the specified range
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    double val = cast(particle);
    return (val >= valMin && val <= valMax);
  }
};

}  // namespace ActsFatras
