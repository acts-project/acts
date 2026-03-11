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
#include "Acts/Surfaces/InfiniteBounds.hpp"  //to get s_noBounds
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/RegularSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceConcept.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace ActsTests {
/// Surface derived class stub
class SurfaceStub : public Acts::RegularSurface {
 public:
  explicit SurfaceStub(
      const Acts::Transform3& htrans = Acts::Transform3::Identity())
      : Acts::GeometryObject(), Acts::RegularSurface(htrans) {}
  SurfaceStub(const Acts::GeometryContext& gctx, const SurfaceStub& sf,
              const Acts::Transform3& transf)
      : Acts::GeometryObject(), Acts::RegularSurface(gctx, sf, transf) {}
  explicit SurfaceStub(const Acts::SurfacePlacementBase& detelement)
      : Acts::GeometryObject(), Acts::RegularSurface(detelement) {}

  ~SurfaceStub() override = default;

  /// Return method for the Surface type to avoid dynamic casts
  Acts::Surface::SurfaceType type() const final { return Acts::Surface::Other; }

  /// Return method for the normal vector of the surface
  Acts::Vector3 normal(const Acts::GeometryContext& /*gctx*/,
                       const Acts::Vector3& /*position*/) const final {
    return Acts::Vector3::Zero();
  }

  Acts::Vector3 normal(const Acts::GeometryContext& /*gctx*/,
                       const Acts::Vector2& /*lposition*/) const final {
    return Acts::Vector3::Zero();
  }

  using Acts::RegularSurface::normal;

  /// Return method for SurfaceBounds
  const Acts::SurfaceBounds& bounds() const final {
    return Acts::s_noBounds;  // need to improve this for meaningful test
  }

  /// Local to global transformation
  Acts::Vector3 localToGlobal(const Acts::GeometryContext& /*gctx*/,
                              const Acts::Vector2& /*lpos*/
  ) const final {
    return Acts::Vector3::Zero();
  }

  using Acts::RegularSurface::localToGlobal;

  /// Global to local transformation
  Acts::Result<Acts::Vector2> globalToLocal(
      const Acts::GeometryContext& /*cxt*/, const Acts::Vector3& /*gpos*/,
      double /*tolerance*/) const final {
    return Acts::Result<Acts::Vector2>::success(Acts::Vector2{20., 20.});
  }

  using Acts::RegularSurface::globalToLocal;

  /// Calculation of the path correction for incident
  double pathCorrection(const Acts::GeometryContext& /*cxt*/,
                        const Acts::Vector3& /*gpos*/,
                        const Acts::Vector3& /*gmom*/) const final {
    return 0.0;
  }

  /// Inherited from GeometryObject base
  Acts::Vector3 referencePosition(const Acts::GeometryContext& /*txt*/,
                                  Acts::AxisDirection /*bValue*/) const final {
    return Acts::Vector3::Zero();
  }

  /// Surface intersction
  Acts::MultiIntersection3D intersect(
      const Acts::GeometryContext& /*gctx*/, const Acts::Vector3& /*position*/,
      const Acts::Vector3& /*direction*/,
      const Acts::BoundaryTolerance& /*boundaryTolerance*/,
      const double /*tolerance*/) const final {
    Acts::Intersection3D stubIntersection(Acts::Vector3(20., 0., 0.), 20.,
                                          Acts::IntersectionStatus::reachable);
    return Acts::MultiIntersection3D(stubIntersection,
                                     Acts::Intersection3D::Invalid());
  }

  /// Return properly formatted class name
  std::string name() const final { return std::string("SurfaceStub"); }

  /// Simply return true to check a method can be called on a constructed object
  bool constructedOk() const { return true; }

  /// Return a Polyhedron for the surfaces
  Acts::Polyhedron polyhedronRepresentation(
      const Acts::GeometryContext& /*gctx*/,
      unsigned int /* ignored */) const final {
    std::vector<Acts::Vector3> vertices;
    std::vector<std::vector<std::size_t>> faces;
    std::vector<std::vector<std::size_t>> triangularMesh;

    return Acts::Polyhedron(vertices, faces, triangularMesh);
  }

  // Cartesian 3D to local bound derivative
  Acts::Matrix<2, 3> localCartesianToBoundLocalDerivative(
      const Acts::GeometryContext& /*gctx*/,
      const Acts::Vector3& /*position*/) const final {
    return Acts::Matrix<2, 3>::Identity();
  };

 private:
  /// the bounds of this surface
  std::shared_ptr<const Acts::PlanarBounds> m_bounds;
};

static_assert(Acts::RegularSurfaceConcept<SurfaceStub>,
              "SurfaceStub does not fulfill RegularSurfaceConcept");

}  // namespace ActsTests
