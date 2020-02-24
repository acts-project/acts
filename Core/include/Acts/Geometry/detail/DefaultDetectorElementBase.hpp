// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DefaultDetectorElementBase.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <memory>
#include <vector>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class Surface;

/// @class DetectorElementBase
///
/// This is the default base class for all tracking detector elements
/// with read-out relevant information. It provides the minimal interface
/// for the Acts proxy mechanism for surfaces, i.e. surfaces in the
/// Tracking geometry representing actual detection devices
///
class DetectorElementBase {
 public:
  /// Constructor
  DetectorElementBase() = default;

  /// Virtual Destructor
  virtual ~DetectorElementBase() = default;

  /// Return the transform for the Element proxy mechanism
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  virtual const Transform3D& transform(const GeometryContext& gctx) const = 0;

  /// Return surface representation
  virtual const Surface& surface() const = 0;

  /// Returns the thickness of the module
  /// @return double that indicates the thickness of the module
  virtual double thickness() const = 0;
};

}  // end of namespace Acts
