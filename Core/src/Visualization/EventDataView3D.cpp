// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

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
      viewConfig.quarterSegments, viewConfig.offset, lposition, transform);

  ellipse.push_back(transform *
                    Vector3(lposition.x(), lposition.y(), viewConfig.offset));
  auto [faces, triangularMesh] =
      detail::FacesHelper::convexFaceMesh(ellipse, true);
  Polyhedron ellipseHedron(ellipse, faces, triangularMesh);
  Acts::GeometryView3D::drawPolyhedron(helper, ellipseHedron, viewConfig);
}

void Acts::EventDataView3D::drawCovarianceAngular(
    IVisualization3D& helper, const Vector3& position, const Vector3& direction,
    const SquareMatrix<2>& covariance, double directionScale,
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
  std::vector<Vector3> ellipse = createEllipse(
      angularErrorScale * directionScale * lambda0 * std::sin(dtheta),
      angularErrorScale * directionScale * lambda1, theta,
      viewConfig.quarterSegments, 0., {0., 0.}, eplane);

  std::vector<Vector3> coneTop = ellipse;
  coneTop.push_back(anker);
  auto [faces, triangularMesh] =
      detail::FacesHelper::convexFaceMesh(coneTop, true);
  Polyhedron coneTopHedron(coneTop, faces, triangularMesh);
  GeometryView3D::drawPolyhedron(helper, coneTopHedron, viewConfig);

  std::vector<Vector3> cone = ellipse;
  cone.push_back(position);
  // Force triangular
  ViewConfig coneViewConfig = viewConfig;
  coneViewConfig.triangulate = true;
  auto [facesCone, triangularMeshCone] =
      detail::FacesHelper::convexFaceMesh(cone, true);
  Polyhedron coneHedron(cone, facesCone, triangularMeshCone);
  GeometryView3D::drawPolyhedron(helper, coneHedron, coneViewConfig);
}

void Acts::EventDataView3D::drawTrack(IVisualization3D& helper,
                                      const AnyConstTrackProxy& track,
                                      const GeometryContext& gctx) {
  auto tparams = track.parameters();
  auto tphi = tparams[eBoundPhi];
  auto ttheta = tparams[eBoundTheta];
  Vector2 tlocpos{tparams[eBoundLoc0], tparams[eBoundLoc1]};
  Vector3 tlocdir{std::sin(ttheta) * std::cos(tphi),
                  std::sin(ttheta) * std::sin(tphi), std::cos(ttheta)};

  auto& rs = track.referenceSurface();
  auto tglobpos = rs.localToGlobal(gctx, tlocpos, tlocdir);

  auto previouspos = tglobpos;

  for (auto ts : track.trackStatesReversed()) {
    auto params = ts.parameters();
    auto phi = params[eBoundPhi];
    auto theta = params[eBoundTheta];
    Vector2 locpos{params[eBoundLoc0], params[eBoundLoc1]};
    Vector3 locdir{std::sin(theta) * std::cos(phi),
                   std::sin(theta) * std::sin(phi), std::cos(theta)};
    auto& s = ts.referenceSurface();
    auto currentpos = s.localToGlobal(gctx, locpos, locdir);

    helper.line(previouspos, currentpos);
    previouspos = currentpos;
  }
}

Acts::StripSpacePointBuilder::StripEnds getStripEnds(
    const Acts::GeometryContext& gctx, const Acts::Surface& surface,
    const Acts::detail::ConstDynamicMeasurement& measurement, const int idx) {
  const auto* bounds =
      dynamic_cast<const Acts::PlanarBounds*>(&surface.bounds());
  if (bounds == nullptr) {
    throw std::invalid_argument(
        "SpacePointMaker: Encountered non-planar surface");
  }
  const Acts::RectangleBounds& boundingBox = bounds->boundingBox();

  double negEnd;
  double posEnd;

  Acts::Vector2 localTop;
  Acts::Vector2 localBottom;

  Acts::Vector3 globalTop;
  Acts::Vector3 globalBottom;

  if (idx == Acts::BoundIndices::eBoundLoc0) {
    const double loc0 = measurement[0];
    negEnd = boundingBox.get(Acts::RectangleBounds::eMinY);
    posEnd = boundingBox.get(Acts::RectangleBounds::eMaxY);
    localTop = {loc0, posEnd};
    localBottom = {loc0, negEnd};
  }

  else {
    const double loc1 = measurement[0];
    negEnd = boundingBox.get(Acts::RectangleBounds::eMinX);
    posEnd = boundingBox.get(Acts::RectangleBounds::eMaxX);
    localTop = {posEnd, loc1};
    localBottom = {negEnd, loc1};
  }

  globalTop = surface.localToGlobal(gctx, localTop, Acts::Vector3::Zero());
  globalBottom =
      surface.localToGlobal(gctx, localBottom, Acts::Vector3::Zero());

  return Acts::StripSpacePointBuilder::StripEnds{globalTop, globalBottom};
}

void Acts::EventDataView3D::drawCluster(IVisualization3D& helper,
                                        const AnyConstTrackProxy& track,
                                        const GeometryContext& gctx) {
  for (auto ts : track.trackStatesReversed()) {
    if (!ts.hasCalibrated()) {
      continue;
    }
    auto meas = ts.effectiveCalibrated();
    auto& s = ts.referenceSurface();
    auto indices = ts.projectorSubspaceIndices();

    const auto n = ts.calibratedSize();

    std::optional<double> loc0;
    std::optional<double> loc1;

    for (std::size_t i = 0; i < n; i++) {
      if (indices[i] == Acts::BoundIndices::eBoundLoc0) {
        loc0 = meas[i];
      }

      else if (indices[i] == Acts::BoundIndices::eBoundLoc1) {
        loc1 = meas[i];
      }

      else {
        continue;
      }
    }

    if (n == 1) {
      auto stripEnds = getStripEnds(gctx, s, meas, indices[0]);
      auto globalpos = 0.5 * (stripEnds.top + stripEnds.bottom);
      helper.vertex(globalpos);
      continue;
    }

    else if (loc0.has_value() && loc1.has_value()) {
      Vector2 locpos{*loc0, *loc1};
      auto globalpos = s.localToGlobal(gctx, locpos, Vector3::Zero());
      helper.vertex(globalpos);
    }

    else {
      throw std::invalid_argument(
          "Measurement has to contain either loc0 or loc1 or both");
    }
  }
}
