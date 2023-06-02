// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"

#include <memory>

namespace Acts {

namespace Experimental {

class Detector;

/// @brief This is the interface defintition of a detector builder
class IDetectorBuilder {
 public:
  virtual ~IDetectorBuilder() = default;
  /// The interface definition for a detector builder
  ///
  /// @param gctx the geometry context at the creation of the detector
  ///
  /// @return a shared Detector object
  virtual std::shared_ptr<Detector> construct(
      const GeometryContext& gctx) const = 0;
};

}  // namespace Experimental
}  // namespace Acts
