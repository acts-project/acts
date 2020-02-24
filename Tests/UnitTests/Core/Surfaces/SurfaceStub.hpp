// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"  //to get s_noBounds
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
/// Surface derived class stub
class SurfaceStub : public Surface {
 public:
  SurfaceStub(std::shared_ptr<const Transform3D> htrans = nullptr)
      : GeometryObject(), Surface(htrans) {}
  SurfaceStub(const GeometryContext& gctx, const SurfaceStub& sf,
              const Transform3D& transf)
      : GeometryObject(), Surface(gctx, sf, transf) {}
  SurfaceStub(const DetectorElementBase& detelement)
      : GeometryObject(), Surface(detelement) {}

  ~SurfaceStub() override { /*nop */
  }

  /// Implicit constructor
  Surface* clone(const GeometryContext& /*gctx*/,
                 const Transform3D& /*shift = nullptr*/) const {
    return nullptr;
  }

  /// Return method for the Surface type to avoid dynamic casts
  SurfaceType type() const final { return Surface::Other; }

  /// Return method for the normal vector of the surface
  const Vector3D normal(const GeometryContext& gctx,
                        const Vector2D& /*lpos*/) const final {
    return normal(gctx);
  }

  const Vector3D normal(const GeometryContext& gctx,
                        const Vector3D&) const final {
    return normal(gctx);
  }

  const Vector3D normal(const GeometryContext& /*gctx*/) const final {
    return Vector3D{0., 0., 0.};
  }

  /// Return method for SurfaceBounds
  const SurfaceBounds& bounds() const final {
    return s_noBounds;  // need to improve this for meaningful test
  }

  /// Local to global transformation
  void localToGlobal(const GeometryContext& /*gctx*/, const Vector2D& /*lpos*/,
                     const Vector3D& /*gmom*/, Vector3D& /*gpos*/) const final {
    // nop
  }

  /// Global to local transformation
  bool globalToLocal(const GeometryContext& /*cxt*/, const Vector3D& /*gpos*/,
                     const Vector3D& /*gmom*/, Vector2D& lpos) const final {
    lpos = Vector2D{20., 20.};
    return true;
  }

  /// Calculation of the path correction for incident
  double pathCorrection(const GeometryContext& /*cxt*/,
                        const Vector3D& /*gpos*/,
                        const Vector3D& /*gmom*/) const final {
    return 0.0;
  }

  /// Straight line intersection schema from parameters
  Intersection intersectionEstimate(
      const GeometryContext& /*cxt*/, const Vector3D& /*gpos*/,
      const Vector3D& /*gdir*/, const BoundaryCheck& /*bcheck*/) const final {
    const Intersection is{Vector3D{1, 1, 1}, 20.,
                          Intersection::Status::reachable};
    return is;
  }

  /// Inherited from GeometryObject base
  const Vector3D binningPosition(const GeometryContext& /*txt*/,
                                 BinningValue /*bValue*/) const final {
    const Vector3D v{0.0, 0.0, 0.0};
    return v;
  }

  /// Return properly formatted class name
  std::string name() const final { return std::string("SurfaceStub"); }

  /// Simply return true to check a method can be called on a constructed object
  bool constructedOk() const { return true; }

  /// Return a Polyhedron for the surfaces
  Polyhedron polyhedronRepresentation(const GeometryContext& /*gctx*/,
                                      size_t /*lseg */) const final {
    std::vector<Vector3D> vertices;
    std::vector<std::vector<size_t>> faces;
    std::vector<std::vector<size_t>> triangularMesh;

    return Polyhedron(vertices, faces, triangularMesh);
  }

 private:
  /// the bounds of this surface
  std::shared_ptr<const PlanarBounds> m_bounds;

  SurfaceStub* clone_impl(const GeometryContext& /*gctx*/,
                          const Transform3D& /* shift */) const override {
    return nullptr;
  }
};
}  // namespace Acts
