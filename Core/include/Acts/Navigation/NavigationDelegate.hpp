// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts {

class NavigationStream;

struct NavigationArguments {
  NavigationStream& main;
  Vector3 position;
  Vector3 direction;

  BoundaryTolerance tolerance = BoundaryTolerance::None();
};

using NavigationDelegate = OwningDelegate<void(const NavigationArguments&)>;

}  // namespace Acts
