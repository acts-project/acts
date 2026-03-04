// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts {

class NavigationStream;
class AppendOnlyNavigationStream;
class NavigationPolicyState;
class Logger;

/// Struct that serves as the argument to the navigation delegate.
/// It is not supposed to be used as an lvalue.
struct NavigationArguments {
  /// Current position in 3D space for navigation
  Vector3 position{};
  /// Direction vector for navigation propagation
  Vector3 direction{};

  /// Boundary tolerance for surface intersection calculations
  BoundaryTolerance tolerance = BoundaryTolerance::None();
};

/// Central alias for the navigation delegate. This type is owning to support
/// (type-erased) navigation delegate chains (i.e. multiple policies).
using NavigationDelegate = Delegate<void(
    const GeometryContext&, const NavigationArguments&, NavigationPolicyState&,
    AppendOnlyNavigationStream&, const Logger&)>;

}  // namespace Acts
