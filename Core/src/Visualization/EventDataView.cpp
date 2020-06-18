// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Visualization/EventDataView.hpp"

void Acts::EventDataView::drawCovarianceCartesian(
    IVisualization& helper, const Vector2D& lposition,
    const ActsSymMatrixD<2>& covariance, const Transform3D& transform,
    double locErrorScale, const ViewConfig& viewConfig) {
  auto [lambda0, lambda1, theta] = decomposeCovariance(covariance);

  std::vector<Vector3D> ellipse = createEllipse(
      lambda0 * locErrorScale, lambda1 * locErrorScale, theta,
      viewConfig.nSegments, viewConfig.offset, lposition, transform);

  ellipse.push_back(transform *
                    Vector3D(lposition.x(), lposition.y(), viewConfig.offset));
  auto faces = detail::FacesHelper::convexFaceMesh(ellipse, true);
  Polyhedron ellipseHedron(ellipse, faces.first, faces.second);
  Acts::GeometryView::drawPolyhedron(helper, ellipseHedron, viewConfig);
}

void Acts::EventDataView::drawCovarianceAngular(
    IVisualization& helper, const Vector3D& position, const Vector3D& direction,
    const ActsSymMatrixD<2>& covariance, double directionScale,
    double angularErrorScale, const ViewConfig& viewConfig) {
  auto [lambda0, lambda1, theta] = decomposeCovariance(covariance);

  // Anker point
  Vector3D anker = position + directionScale * direction;

  double dphi = VectorHelpers::phi(direction);
  double dtheta = VectorHelpers::theta(direction);

  Transform3D eplane(Translation3D(anker) *
                     AngleAxis3D(dtheta, Vector3D(1., 0., 0.)) *
                     AngleAxis3D(dphi, Vector3D(0., 0., 1.)));

  // Now generate the ellipse points
  std::vector<Vector3D> ellipse =
      createEllipse(angularErrorScale * directionScale * tan(lambda0),
                    angularErrorScale * directionScale * tan(lambda1), theta,
                    viewConfig.nSegments, 0., {0., 0.}, eplane);

  std::vector<Vector3D> coneTop = ellipse;
  coneTop.push_back(anker);
  auto coneTopFaces = detail::FacesHelper::convexFaceMesh(coneTop, true);
  Polyhedron coneTopHedron(coneTop, coneTopFaces.first, coneTopFaces.second);
  GeometryView::drawPolyhedron(helper, coneTopHedron, viewConfig);

  std::vector<Vector3D> cone = ellipse;
  cone.push_back(position);
  // Force triangular
  ViewConfig coneViewConfig = viewConfig;
  coneViewConfig.triangulate = true;
  auto coneFaces = detail::FacesHelper::convexFaceMesh(cone, true);
  Polyhedron coneHedron(cone, coneFaces.first, coneFaces.second);
  GeometryView::drawPolyhedron(helper, coneHedron, coneViewConfig);
}
