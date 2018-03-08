// This file is part of the ACTS project.
//
// Copyright (C) 2016-2017 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlanarBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_PLANARBOUNDS_H
#define ACTS_SURFACES_PLANARBOUNDS_H 1

#include <vector>

#include "ACTS/Surfaces/SurfaceBounds.hpp"
#include "ACTS/Utilities/VariantDataFwd.hpp"

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

  virtual variant_data
  toVariantData() const = 0;
};

}  // end of namespace

#endif  // ACTS_SURFACES_PLANARBOUNDS_H
