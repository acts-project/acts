// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once

#include <limits>

#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
namespace Test {

class LineSurfaceStub : public LineSurface {
 public:
  LineSurfaceStub() = delete;
  //
  LineSurfaceStub(std::shared_ptr<const Transform3D> htrans, double radius,
                  double halfz)
      : GeometryObject(), LineSurface(htrans, radius, halfz) { /* nop */
  }
  //
  LineSurfaceStub(std::shared_ptr<const Transform3D> htrans,
                  std::shared_ptr<const LineBounds> lbounds = nullptr)
      : GeometryObject(), LineSurface(htrans, lbounds) { /*nop */
  }
  //
  LineSurfaceStub(std::shared_ptr<const LineBounds> lbounds,
                  const DetectorElementBase& detelement)
      : GeometryObject(), LineSurface(lbounds, detelement) { /* nop */
  }

  //
  LineSurfaceStub(const LineSurfaceStub& ls)
      : GeometryObject(), LineSurface(ls) { /* nop */
  }
  //
  LineSurfaceStub(const GeometryContext& gctx, const LineSurfaceStub& ls,
                  const Transform3D& t)
      : GeometryObject(), LineSurface(gctx, ls, t) { /* nop */
  }
  /// pure virtual functions of baseclass implemented here
  std::shared_ptr<LineSurfaceStub> clone(const GeometryContext& /*gctx*/,
                                         const Transform3D& /*unused*/) const {
    return nullptr;
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

 private:
  Surface* clone_impl(const GeometryContext& /*gctx*/,
                      const Transform3D& /*unused*/) const {
    return nullptr;
  }
};
}  // namespace Test
}  // namespace Acts
