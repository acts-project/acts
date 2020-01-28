// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace ActsFatras {

struct X0Limit {
  /// Return true if the limit in X0 is reached
  /// @todo modify, return what's left from the detector: needs non-const
  /// detector
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &detector,
                  const particle_t &particle) const {
    return particle.pathInX0() +
               detector.thickness() / detector.material().X0() >=
           particle.limitInX0();
  }
};

struct L0Limit {
  /// Return true if the limit in X0 is reached
  /// @todo modify, return what's left from the detector: needs non-const
  /// detector
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &detector,
                  const particle_t &particle) const {
    return particle.pathInL0() +
               detector.thickness() / detector.material().L0() >=
           particle.limitInL0();
  }
};

}  // namespace ActsFatras
