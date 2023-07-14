// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/PointSurface.hpp"

#include <limits>

namespace Acts {
namespace Test {

class PointSurfaceStub : public PointSurface {
 public:
  PointSurfaceStub() = delete;
  //
  PointSurfaceStub(const Translation3& htrans, double radius)
      : GeometryObject(), PointSurface(htrans, radius) { /* nop */
  }
  //
  PointSurfaceStub(const Translation3& htrans,
                   std::shared_ptr<const PointBounds> pbounds = nullptr)
      : GeometryObject(), PointSurface(htrans, std::move(pbounds)) { /*nop */
  }
  //
  PointSurfaceStub(std::shared_ptr<const PointBounds> pbounds,
                   const DetectorElementBase& detelement)
      : GeometryObject(),
        PointSurface(std::move(pbounds), detelement) { /* nop */
  }

  //
  PointSurfaceStub(const PointSurfaceStub& ls)
      : GeometryObject(), PointSurface(ls) { /* nop */
  }

  PointSurfaceStub& operator=(const PointSurfaceStub& ls) = default;

  //
  PointSurfaceStub(const GeometryContext& gctx, const PointSurfaceStub& ls,
                   const Translation3& t)
      : GeometryObject(), PointSurface(gctx, ls, t) { /* nop */
  }

  /// Return method for the Surface type to avoid dynamic casts
  SurfaceType type() const final { return Surface::Point; }

  /// Simply return true to show object exists and is callable
  bool constructedOk() const { return true; }

  using Surface::normal;

  /// Return a Polyhedron for the surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lseg is ignored for a perigee @note ignored
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& /*gctx*/,
                                      size_t /*lseg*/) const final {
    return Polyhedron({}, {}, {});
  }
};
}  // namespace Test
}  // namespace Acts
