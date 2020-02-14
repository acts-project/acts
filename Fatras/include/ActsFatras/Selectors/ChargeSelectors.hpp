// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {

/// Select neutral particles.
struct NeutralSelector {
  template <typename detector_t>
  bool operator()(const detector_t &, const Particle &particle) const {
    return (particle.charge() == Particle::Scalar(0));
  }
};

/// Select all charged particles.
struct ChargedSelector {
  template <typename detector_t>
  bool operator()(const detector_t &, const Particle &particle) const {
    return (particle.charge() != Particle::Scalar(0));
  }
};

/// Select positively charged particles.
struct PositiveSelector {
  template <typename detector_t>
  bool operator()(const detector_t &, const Particle &particle) const {
    return (Particle::Scalar(0) < particle.charge());
  }
};

/// Select negatively charged particles.
struct NegativeSelector {
  template <typename detector_t>
  bool operator()(const detector_t &, const Particle &particle) const {
    return (particle.charge() < Particle::Scalar(0));
  }
};

}  // namespace ActsFatras
