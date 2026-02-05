// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>

namespace Acts {

class GeometryContext;
class Surface;

/// Plain navigator options carrying geometry context and surfaces.
struct NavigatorPlainOptions {
  /// NavigatorPlainOptions with context
  /// @param gctx The geometry context
  explicit NavigatorPlainOptions(const GeometryContext &gctx)
      : geoContext(gctx) {}

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;

  /// Start surface for navigation
  const Surface *startSurface{};
  /// Target surface for navigation
  const Surface *targetSurface{};
};

}  // namespace Acts
