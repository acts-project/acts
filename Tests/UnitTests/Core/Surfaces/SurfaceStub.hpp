// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/InfiniteBounds.hpp"  //to get s_noBounds
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {
/// Surface derived class stub
class SurfaceStub : public Surface {
 public:
  SurfaceStub(const Transform3& htrans = Transform3::Identity())
      : GeometryObject(), Surface(htrans) {}
  SurfaceStub(const GeometryContext& gctx, const SurfaceStub& sf,
              const Transform3& transf)
      : GeometryObject(), Surface(gctx, sf, transf) {}
  SurfaceStub(const DetectorElementBase& detelement)
      : GeometryObject(), Surface(detelement) {}

  ~SurfaceStub() override = default;

  /// Return method for the Surface type to avoid dynamic casts
  SurfaceType type() const final { return Surface::Other; }

  /// Return method for the normal vector of the surface
  Vector3 normal(const GeometryContext& gctx,
                 const Vector2& /*lpos*/) const final {
    return normal(gctx);
  }

  Vector3 normal(const GeometryContext& gctx,
                 const Vector3& /*position*/) const final {
    return normal(gctx);
  }

  Vector3 normal(const GeometryContext& /*gctx*/) const final {
    return Vector3{0., 0., 0.};
  }

  /// Return method for SurfaceBounds
  const SurfaceBounds& bounds() const final {
    return s_noBounds;  // need to improve this for meaningful test
  }

  /// Local to global transformation
  Vector3 localToGlobal(const GeometryContext& /*gctx*/,
                        const Vector2& /*lpos*/,
                        const Vector3& /*gmom*/) const final {
    return Vector3(0., 0., 0.);
  }

  /// Global to local transformation
  Result<Vector2> globalToLocal(const GeometryContext& /*cxt*/,
                                const Vector3& /*gpos*/,
                                const Vector3& /*gmom*/,
                                double /*tolerance*/) const final {
    return Result<Vector2>::success(Vector2{20., 20.});
  }

  /// Calculation of the path correction for incident
  double pathCorrection(const GeometryContext& /*cxt*/, const Vector3& /*gpos*/,
                        const Vector3& /*gmom*/) const final {
    return 0.0;
  }

  /// Inherited from GeometryObject base
  Vector3 binningPosition(const GeometryContext& /*txt*/,
                          BinningValue /*bValue*/) const final {
    const Vector3 v{0.0, 0.0, 0.0};
    return v;
  }

  /// Surface intersction
  SurfaceIntersection intersect(const GeometryContext& /*gctx*/,
                                const Vector3& /*position*/,
                                const Vector3& /*direction*/,
                                const BoundaryCheck& /*bcheck*/) const final {
    Intersection3D stubIntersection(Vector3(20., 0., 0.), 20.,
                                    Intersection3D::Status::reachable);
    return SurfaceIntersection(stubIntersection, this);
  }

  /// Return properly formatted class name
  std::string name() const final { return std::string("SurfaceStub"); }

  /// Simply return true to check a method can be called on a constructed object
  bool constructedOk() const { return true; }

  /// Return a Polyhedron for the surfaces
  Polyhedron polyhedronRepresentation(const GeometryContext& /*gctx*/,
                                      size_t /*lseg */) const final {
    std::vector<Vector3> vertices;
    std::vector<std::vector<size_t>> faces;
    std::vector<std::vector<size_t>> triangularMesh;

    return Polyhedron(vertices, faces, triangularMesh);
  }

  // Cartesian 3D to local bound derivative
  ActsMatrix<2, 3> localCartesianToBoundLocalDerivative(
      const GeometryContext& /*gctx*/,
      const Vector3& /*position*/) const final {
    return ActsMatrix<2, 3>::Identity();
  };

 private:
  /// the bounds of this surface
  std::shared_ptr<const PlanarBounds> m_bounds;
};
}  // namespace Acts
