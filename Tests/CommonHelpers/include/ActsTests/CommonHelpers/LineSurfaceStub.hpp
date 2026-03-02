// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace ActsTests {

class LineSurfaceStub : public Acts::LineSurface {
 public:
  LineSurfaceStub() = delete;

  LineSurfaceStub(const Acts::Transform3& htrans, double radius, double halfz)
      : Acts::GeometryObject(),
        Acts::LineSurface(htrans, radius, halfz) { /* nop */ }

  explicit LineSurfaceStub(
      const Acts::Transform3& htrans,
      std::shared_ptr<const Acts::LineBounds> lbounds = nullptr)
      : Acts::GeometryObject(),
        Acts::LineSurface(htrans, std::move(lbounds)) { /*nop */ }

  LineSurfaceStub(std::shared_ptr<const Acts::LineBounds> lbounds,
                  const Acts::SurfacePlacementBase& placement)
      : Acts::GeometryObject(),
        Acts::LineSurface(std::move(lbounds), placement) { /* nop */ }

  LineSurfaceStub(const LineSurfaceStub& ls)
      : Acts::GeometryObject(), Acts::LineSurface(ls) { /* nop */ }

  LineSurfaceStub& operator=(const LineSurfaceStub& ls) = default;

  LineSurfaceStub(const Acts::GeometryContext& gctx, const LineSurfaceStub& ls,
                  const Acts::Transform3& t)
      : Acts::GeometryObject(), Acts::LineSurface(gctx, ls, t) { /* nop */ }

  /// Return method for the Surface type to avoid dynamic casts
  SurfaceType type() const final { return Acts::Surface::Straw; }

  /// Simply return true to show object exists and is callable
  bool constructedOk() const { return true; }

  using Acts::Surface::normal;

  /// Return a Polyhedron for the surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param ingoredSegmeent is ignored for the srub
  ///
  /// @return A list of vertices and a face/facett description of it
  Acts::Polyhedron polyhedronRepresentation(
      const Acts::GeometryContext& /*gctx*/,
      unsigned int /*lseg*/) const final {
    return Acts::Polyhedron({}, {}, {});
  }
};

}  // namespace ActsTests
