// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

namespace ActsFatras {
namespace detail {

/// Cast function definition for bound parameters
using BoundCastFunction =
    std::function<double(const Acts::Vector2D&, const Acts::Vector3D&, double)>;

/// Cast operator to provide local 0
struct BoundLoc0Cast {
  double operator()(const Acts::Vector2D& lhit,
                    const Acts::Vector3D& /*ignored*/, double /*ignored*/) {
    return lhit[Acts::eBoundLoc0];
  }
};

/// Cast operator to provide local 1
struct BoundLoc1Cast {
  double operator()(const Acts::Vector2D& lhit,
                    const Acts::Vector3D& /*ignored*/, double /*ignored*/) {
    return lhit[Acts::eBoundLoc1];
  }
};

/// Cast operator to provide phi
struct BoundPhiCast {
  double operator()(const Acts::Vector2D& /*ignored*/,
                    const Acts::Vector3D& dir, double /*ignored*/) {
    return Acts::VectorHelpers::phi(dir);
  }
};

/// Cast operator to provide theta
struct BoundThetaCast {
  double operator()(const Acts::Vector2D& /*ignored*/,
                    const Acts::Vector3D& dir, double /*ignored*/) {
    return Acts::VectorHelpers::theta(dir);
  }
};

/// Cast operator to provide 0 for QoP (not smeared)
struct BoundQoPCast {
  double operator()(const Acts::Vector2D& /*ignored*/,
                    const Acts::Vector3D& /*ignored*/, double /*ignored*/) {
    return 0.;
  }
};

/// Cast operator to provide theta
struct BoundTimeCast {
  double operator()(const Acts::Vector2D& /*ignored*/,
                    const Acts::Vector3D& /*ignored*/, double ctime) {
    return ctime;
  }
};

static std::array<BoundCastFunction, Acts::eBoundSize> BoundCasts = {
    BoundLoc0Cast{},  BoundLoc1Cast{}, BoundPhiCast{},
    BoundThetaCast{}, BoundQoPCast{},  BoundTimeCast{}};

}  // namespace detail
}  // namespace ActsFatras