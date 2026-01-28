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
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceConcept.hpp"

#include <memory>
#include <string>

namespace Acts {

struct Polyhedron;
class LineBounds;

///  @class StrawSurface
///
///  Class for a StrawSurface in the TrackingGeometry
///  to describe dirft tube and straw like detectors.
///
/// @image html LineSurface.png
///
class StrawSurface : public LineSurface {
  friend class Surface;

 protected:
  /// @brief Implement all constructors from the mother class
  using LineSurface::LineSurface;

 public:
  /// Return the surface type
  /// @return Surface type identifier for straw surfaces
  SurfaceType type() const final { return Surface::Straw; }

  /// Return properly formatted class name for screen output */
  /// @return String representation of the surface type name
  std::string name() const final { return "Acts::StrawSurface"; }

  /// Return a Polyhedron for the surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param quarterSegments is the number of segments used to describe curved
  /// segments in a quarter of the phi range. If it is 1, then only the extrema
  /// points in phi are inserted next to the segment corners.
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& gctx,
                                      unsigned int quarterSegments) const final;
};

static_assert(SurfaceConcept<StrawSurface>,
              "StrawSurface does not fulfill SurfaceConcept");

}  // namespace Acts
