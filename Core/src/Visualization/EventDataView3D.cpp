// This file is part of the Acts project.
//
// Copyright (C) 2020-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/EventDataView3D.hpp"

#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Utilities/Helpers.hpp"

#include <cmath>
#include <utility>

namespace Acts {
class IVisualization3D;
}  // namespace Acts

void Acts::EventDataView3D::drawCovarianceCartesian(
    IVisualization3D& helper, const Vector2& lposition,
    const SquareMatrix2& covariance, const Transform3& transform,
    double locErrorScale, const ViewConfig& viewConfig) {
  auto [lambda0, lambda1, theta] = decomposeCovariance(covariance);

  std::vector<Vector3> ellipse = createEllipse(
      lambda0 * locErrorScale, lambda1 * locErrorScale, theta,
      viewConfig.nSegments, viewConfig.offset, lposition, transform);

  ellipse.push_back(transform *
                    Vector3(lposition.x(), lposition.y(), viewConfig.offset));
  auto faces = detail::FacesHelper::convexFaceMesh(ellipse, true);
  Polyhedron ellipseHedron(ellipse, faces.first, faces.second);
  Acts::GeometryView3D::drawPolyhedron(helper, ellipseHedron, viewConfig);
}

void Acts::EventDataView3D::drawCovarianceAngular(
    IVisualization3D& helper, const Vector3& position, const Vector3& direction,
    const ActsSquareMatrix<2>& covariance, double directionScale,
    double angularErrorScale, const ViewConfig& viewConfig) {
  auto [lambda0, lambda1, theta] = decomposeCovariance(covariance);

  // Anker point
  Vector3 anker = position + directionScale * direction;

  double dphi = VectorHelpers::phi(direction);
  double dtheta = VectorHelpers::theta(direction);

  Transform3 eplane(Translation3(anker) *
                    AngleAxis3(dphi, Vector3(0., 0., 1.)) *
                    AngleAxis3(dtheta, Vector3(0., 1., 0.)));

  // Now generate the ellipse points
  std::vector<Vector3> ellipse =
      createEllipse(angularErrorScale * directionScale * lambda0 * sin(dtheta),
                    angularErrorScale * directionScale * lambda1, theta,
                    viewConfig.nSegments, 0., {0., 0.}, eplane);

  std::vector<Vector3> coneTop = ellipse;
  coneTop.push_back(anker);
  auto coneTopFaces = detail::FacesHelper::convexFaceMesh(coneTop, true);
  Polyhedron coneTopHedron(coneTop, coneTopFaces.first, coneTopFaces.second);
  GeometryView3D::drawPolyhedron(helper, coneTopHedron, viewConfig);

  std::vector<Vector3> cone = ellipse;
  cone.push_back(position);
  // Force triangular
  ViewConfig coneViewConfig = viewConfig;
  coneViewConfig.triangulate = true;
  auto coneFaces = detail::FacesHelper::convexFaceMesh(cone, true);
  Polyhedron coneHedron(cone, coneFaces.first, coneFaces.second);
  GeometryView3D::drawPolyhedron(helper, coneHedron, coneViewConfig);
}
