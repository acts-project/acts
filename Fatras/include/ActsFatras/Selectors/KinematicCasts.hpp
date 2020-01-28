// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

#include "Acts/Utilities/Helpers.hpp"

namespace ActsFatras {
namespace Casts {

/// The Eta cast operator
struct eta {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return Acts::VectorHelpers::eta(particle.momentum());
  }
};

/// The Eta cast operator
struct absEta {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return std::abs(Acts::VectorHelpers::eta(particle.momentum()));
  }
};

/// The Pt cast operator
struct pT {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return Acts::VectorHelpers::perp(particle.momentum());
  }
};

/// The P cast operator
struct p {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.momentum().norm();
  }
};

/// The E cast operator
struct E {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.E();
  }
};

/// The E cast operator
struct vR {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return Acts::VectorHelpers::perp(particle.position());
  }
};

/// The E cast operator
struct vZ {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return particle.position().z();
  }
};

/// The E cast operator
struct AbsVz {
  template <typename particle_t>
  double operator()(const particle_t &particle) const {
    return std::abs(particle.position().z());
  }
};

}  // namespace Casts
}  // namespace ActsFatras
