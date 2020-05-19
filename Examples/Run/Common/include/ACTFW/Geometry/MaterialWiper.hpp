// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialWiper.hpp, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Surfaces/Surface.hpp"

// @note This file will go into the acts-core
namespace Acts {

/// @class MaterialWiper
///
/// This decorator sets the nulls-material
///
class MaterialWiper : public IMaterialDecorator {
 public:
  /// Decorate a surface
  ///
  /// @param surface the non-cost surface that is decorated
  void decorate(Surface& surface) const final {
    surface.assignSurfaceMaterial(nullptr);
  }

  /// Decorate a TrackingVolume
  ///
  /// @param volume the non-cost volume that is decorated
  virtual void decorate(TrackingVolume& volume) const final {
    volume.assignVolumeMaterial(nullptr);
  }
};

}  // namespace Acts