// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Navigation/NavigationStream.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts {

class NavigationStream;
class Logger;

/// Struct that serves as the argument to the navigation delegate.
/// It is not supposed to be used as an lvalue.
struct NavigationArguments {
  Vector3 position;
  Vector3 direction;

  BoundaryTolerance tolerance = BoundaryTolerance::None();

  bool wantsPortals = true;
  bool wantsSurfaces = true;
};

/// Central alias for the navigation delegate. This type is owning to support
/// (type-erased) navigation delegate chains (i.e. multiple policies).
// @TODO: Add geometry context to navigation delegate signature
using NavigationDelegate = OwningDelegate<void(
    const NavigationArguments&, AppendOnlyNavigationStream&, const Logger&)>;

}  // namespace Acts
