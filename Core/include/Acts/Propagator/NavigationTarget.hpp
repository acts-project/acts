// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/BoundaryTolerance.hpp"

#include <utility>

namespace Acts {

class Surface;

struct NavigationTarget {
  const Surface* surface = nullptr;
  std::uint8_t surfaceIntersectionIndex = 0;
  BoundaryTolerance boundaryTolerance = BoundaryTolerance::None();

  static NavigationTarget invalid() { return NavigationTarget(); }

  NavigationTarget(const Surface& surface_,
                   std::uint8_t surfaceIntersectionIndex_,
                   BoundaryTolerance boundaryTolerance_)
      : surface(&surface_),
        surfaceIntersectionIndex(surfaceIntersectionIndex_),
        boundaryTolerance(std::move(boundaryTolerance_)) {}

  bool isValid() const { return surface != nullptr; }

 private:
  NavigationTarget() = default;
};

}  // namespace Acts
