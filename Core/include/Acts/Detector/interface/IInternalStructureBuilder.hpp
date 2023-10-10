// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/DetectorComponents.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

namespace Acts {
namespace Experimental {

/// @brief This is the interface definition of internal structure
/// builders for DetectorVolume construction.
///
/// It is assumed that each builder returns a consistent set of
/// DetectorVolume internals, which in turn can be directly provided
/// to a DetectorVolume constructor.
class IInternalStructureBuilder {
 public:
  virtual ~IInternalStructureBuilder() = default;
  /// The interface definition for internal structure creation
  ///
  /// @param gctx the geometry context at the creation of the internal structure
  ///
  /// @return a consistent set of detector volume internals
  virtual InternalStructure construct(const GeometryContext& gctx) const = 0;
};

}  // namespace Experimental
}  // namespace Acts
