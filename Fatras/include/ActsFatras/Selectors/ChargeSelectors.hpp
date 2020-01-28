// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace ActsFatras {

struct ChargedSelector {
  /// Return true for all particles with charge != 0.
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return (particle.q() * particle.q() > 0.);
  }
};

struct NeutralSelector {
  /// Return true for all particles with charge == 0.
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return (particle.q() == 0.);
  }
};

struct PositiveSelector {
  /// Return true for all particles with charge > 0.
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return (particle.q() > 0.);
  }
};

struct NegativeSelector {
  /// Return true for all particles with charge<> 0.
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &particle) const {
    return (particle.q() * particle.q() > 0.);
  }
};

}  // namespace ActsFatras
