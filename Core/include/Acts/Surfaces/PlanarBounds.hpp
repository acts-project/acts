// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlanarBounds.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include <vector>

#include "Acts/Surfaces/SurfaceBounds.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

/// forward declare rectangle bounds as boundary box
class RectangleBounds;

///
/// @class PlanarBounds
///
/// common base class for all bounds that are in a local x/y cartesian frame
///  - simply introduced to avoid wrong bound assigments to surfaces
///
class PlanarBounds : public SurfaceBounds
{
public:
  /// Return the vertices - or, the points of the extremas
  virtual std::vector<Vector2D>
  vertices() const = 0;

  // Bounding box parameters
  virtual const RectangleBounds&
  boundingBox() const = 0;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const override = 0;
};

}  // namespace