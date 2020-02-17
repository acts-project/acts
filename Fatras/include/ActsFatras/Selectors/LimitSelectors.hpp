// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialProperties.hpp"
#include "ActsFatras/EventData/Particle.hpp"

namespace ActsFatras {

/// Select particles whose X0 limit would be reached after material passage.
struct PathLimitX0 {
  bool operator()(const Particle &particle,
                  const Acts::MaterialProperties &slab) const {
    return particle.pathLimitX0() <
           (particle.pathInX0() + slab.thicknessInX0());
  }
};

/// Select particles whose L0 limit would be reached after material passage.
struct PathLimitL0 {
  bool operator()(const Particle &particle,
                  const Acts::MaterialProperties &slab) const {
    return particle.pathLimitL0() <
           (particle.pathInL0() + slab.thicknessInL0());
  }
};

}  // namespace ActsFatras
