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
#include "Acts/Utilities/Delegate.hpp"

#include <memory>
#include <vector>

namespace Acts {

  class Surface;
  class DetectorElementBase;

  /// A possible implementation of an alignment context via the surface
  using AlignmentDelegate = Delegate<const Transform3*(const Surface&)>;

  /// A possible implementation of using a simple delegate driven
  /// alignmnet infrastructure
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param detElement is the detector element
  ///
  /// Client code that wants to use this proposed alginment context can
  /// call this metod from the `transform(const GeometryContext&)` method
  ///
  /// It unpacks the GeometryContext to check if it contains an AlignmentContext
  /// and then calls the included delegate. If no transform is found, a nullptr
  /// is returned.
  ///
  /// Client code that implements its on alignmnet context can savely ignore this
  /// method and simply overwrite the `transform(const GeometryContext&)` method
  ///
  /// @return a pointer to the transform if found, otherwise nullptr
  const Transform3* contextualTransform(const GeometryContext& gctx, const DetectorElementBase& detElement);

} // namespace Acts
