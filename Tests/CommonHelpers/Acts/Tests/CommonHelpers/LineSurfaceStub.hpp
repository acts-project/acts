// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/LineSurface.hpp"

#include <limits>

namespace Acts {
namespace Test {

class LineSurfaceStub : public LineSurface {
 public:
  LineSurfaceStub() = delete;
  //
  LineSurfaceStub(const Transform3& htrans, double radius, double halfz)
      : GeometryObject(), LineSurface(htrans, radius, halfz) { /* nop */
  }
  //
  LineSurfaceStub(const Transform3& htrans,
                  std::shared_ptr<const LineBounds> lbounds = nullptr)
      : GeometryObject(), LineSurface(htrans, std::move(lbounds)) { /*nop */
  }
  //
  LineSurfaceStub(std::shared_ptr<const LineBounds> lbounds,
                  const DetectorElementBase& detelement)
      : GeometryObject(),
        LineSurface(std::move(lbounds), detelement) { /* nop */
  }

  //
  LineSurfaceStub(const LineSurfaceStub& ls)
      : GeometryObject(), LineSurface(ls) { /* nop */
  }

  LineSurfaceStub& operator=(const LineSurfaceStub& ls) = default;

  //
  LineSurfaceStub(const GeometryContext& gctx, const LineSurfaceStub& ls,
                  const Transform3& t)
      : GeometryObject(), LineSurface(gctx, ls, t) { /* nop */
  }

  /// Return method for the Surface type to avoid dynamic casts
  SurfaceType type() const final { return Surface::Straw; }

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
