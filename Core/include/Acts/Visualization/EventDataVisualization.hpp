// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Visualization/GeometryVisualization.hpp"
#include "Acts/Visualization/IVisualization.hpp"

#include <optional>

namespace Acts {

namespace EventDataVisualization {

/// Helper to find the egen values and corr angle
///
/// @param covariance The covariance matrix
static inline std::array<double, 3> decomposeCovariance(
    const ActsSymMatrixD<2>& covariance) {
  double c00 = covariance(eLOC_0, eLOC_0);
  double c01 = covariance(eLOC_0, eLOC_1);
  double c11 = covariance(eLOC_1, eLOC_1);

  double cdsq = std::pow((c00 - c11), 2) / 4.;
  double cosq = c01 * c01;

  // Calculate the eigen values w.r.t reference frame
  double lambda0 = (c00 + c11) / 2. + std::sqrt(cdsq + cosq);
  double lambda1 = (c00 + c11) / 2. - std::sqrt(cdsq + cosq);
  double theta = atan2(lambda0 - c00, c01);

  return {lambda0, lambda1, theta};
}

/// Helper mehod to draw the ellipse points
///
/// @param lambda0 The Eigenvalue in 0
/// @param lambda1 The Eigenvalue in 1
/// @param theta The angle between the x/y frame and EV frame
/// @param lseg The number of segments
/// @param outOfPlane The out of plane offset for visibility
/// @param lposition The local anker point of the ellipse
/// @param transform The transform to global
static inline std::vector<Vector3D> createEllipse(
    double lambda0, double lambda1, double theta, size_t lseg,
    double outOfPlane, const Vector2D& lposition = Vector2D(0., 0.),
    const Transform3D& transform = Transform3D::Identity()) {
  double ctheta = std::cos(theta);
  double stheta = std::sin(theta);

  double l1sq = std::sqrt(lambda0);
  double l2sq = std::sqrt(lambda1);

  // Now generate the ellipse points
  std::vector<Vector3D> ellipse;
  ellipse.reserve(lseg);
  double thetaStep = 2 * M_PI / lseg;
  for (size_t it = 0; it < lseg; ++it) {
    double phi = -M_PI + it * thetaStep;
    double cphi = std::cos(phi);
    double sphi = std::sin(phi);
    double x = lposition.x() + (l1sq * ctheta * cphi - l2sq * stheta * sphi);
    double y = lposition.y() + (l1sq * stheta * cphi + l2sq * ctheta * sphi);
    ellipse.push_back(transform * Vector3D(x, y, outOfPlane));
  }
  return ellipse;
}

/// Helper method to draw error ellipse
///
/// @param helper [in, out] The visualization helper
/// @param lposition The local position
/// @param covariance The covariance matrix
/// @param transform The reference Frame transform
/// @param nsigma The sigmas to be drawn
/// @param locErrorScale The local Error scale
/// @param color The draw color
/// @param outOfPlane The out of plane drawning option
static inline void drawCovarianceCartesian(
    IVisualization& helper, const Vector2D& lposition,
    const ActsSymMatrixD<2>& covariance, const Transform3D& transform,
    std::vector<int> nsigma = {3}, double locErrorScale = 1, size_t lseg = 72,
    const IVisualization::ColorType& color = {20, 120, 20},
    double outOfPlane = 0.1) {
  auto [lambda0, lambda1, theta] = decomposeCovariance(covariance);

  // Now generate the ellipse points
  std::vector<Vector3D> ellipse =
      createEllipse(lambda0 * locErrorScale, lambda1 * locErrorScale, theta,
                    lseg, outOfPlane, lposition, transform);

  ellipse.push_back(transform *
                    Vector3D(lposition.x(), lposition.y(), outOfPlane));
  auto faces = detail::FacesHelper::convexFaceMesh(ellipse, true);
  Polyhedron ellipseHedron(ellipse, faces.first, faces.second);
  ellipseHedron.draw(helper, false, color);
}

/// Helper method to draw error cone of a direction
///
/// @param helper [in, out] The visualization helper
/// @param position Where the cone originates from
/// @param direction The direction parameters
/// @param covariance The 2x2 covariance matrix for phi/theta
/// @param nsigma The sigmas to be drawn
/// @param directionScale The direction arror length
/// @param angularErrorScale The local Error scale
/// @param lseg The number of segments
/// @param color The draw color
static inline void drawCovarianceAngular(
    IVisualization& helper, const Vector3D& position, const Vector3D& direction,
    const ActsSymMatrixD<2>& covariance, std::vector<int> nsigma = {3},
    double directionScale = 1, double angularErrorScale = 1, size_t lseg = 72,
    const IVisualization::ColorType& color = {20, 120, 20}) {
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
                    lseg, 0., {0., 0.}, eplane);

  std::vector<Vector3D> coneTop = ellipse;
  coneTop.push_back(anker);
  auto coneTopFaces = detail::FacesHelper::convexFaceMesh(coneTop, true);
  Polyhedron coneTopHedron(coneTop, coneTopFaces.first, coneTopFaces.second);
  coneTopHedron.draw(helper, false, color);

  std::vector<Vector3D> cone = ellipse;
  cone.push_back(position);
  auto coneFaces = detail::FacesHelper::convexFaceMesh(cone, true);
  Polyhedron coneHedron(cone, coneFaces.first, coneFaces.second);
  coneHedron.draw(helper, true, color);
}

/// Helper method to draw bound parameters object
///
/// @param helper [in, out] The visualization helper
/// @param parameters The bound parameters to be drawn
/// @param gctx The geometry context for which it is drawn
/// @param momentumScale The scale of the momentum
/// @param locErrorScale  The scale of the local error
/// @param angularErrorScale The sclae of the angular error
/// @param drawParameterSurface The indicator whether to draw the surface
/// @param lseg The number of segments for a full arch (if needed)
/// @param pcolor the (optional) color of the parameters to be written
/// @param scolor the (optional) color of the surface to be written
/// @param outOfPlane The out of plane drawning option
template <typename parameters_t>
static inline void drawBoundParameters(
    IVisualization& helper, const parameters_t& parameters,
    const GeometryContext& gctx = GeometryContext(), double momentumScale = 1.,
    double locErrorScale = 1., double angularErrorScale = 1.,
    bool drawParameterSurface = true, size_t lseg = 72,
    const IVisualization::ColorType& pcolor = {20, 120, 20},
    const IVisualization::ColorType& scolor = {235, 198, 52},
    double outOfPlane = 0.1) {
  // First, if necessary, draw the surface
  if (drawParameterSurface) {
    GeometryVisualization::drawSurface(helper, parameters.referenceSurface(),
                                       gctx, Transform3D::Identity(), lseg,
                                       false, scolor);
  }

  // Draw the parameter shaft and cone
  auto position = parameters.position();
  auto direction = parameters.momentum().normalized();
  double p = parameters.momentum().norm();

  Vector3D startmod = parameters.covariance().has_value()
                          ? 0.25 * p * momentumScale * direction
                          : Vector3D(0., 0., 0.);

  GeometryVisualization::drawArrowForward(
      helper, position, position + p * momentumScale * direction, 0.025, 0.05,
      2., 72, pcolor);

  if (parameters.covariance().has_value()) {
    auto paramVec = parameters.parameters();
    auto lposition = paramVec.template block<2, 1>(0, 0);

    // Draw the local covariance
    const auto& covariance = *parameters.covariance();
    drawCovarianceCartesian(helper, lposition,
                            covariance.template block<2, 2>(0, 0),
                            parameters.referenceSurface().transform(gctx), {3},
                            locErrorScale, 72, pcolor, outOfPlane);

    drawCovarianceAngular(
        helper, parameters.position(), parameters.momentum().normalized(),
        covariance.template block<2, 2>(2, 2), {3}, 0.9 * p * momentumScale,
        angularErrorScale, 72, pcolor);
  }
}

/// Helper method to draw one trajectory stored in a MultiTrajectory object
///
/// @tparam source_link_t The source link type

/// @param helper [in, out] The visualization helper
/// @param multiTraj The MultiTrajectory storing the trajectory to be drawn
/// @param entryIndex The trajectory entry index
/// @param gctx The geometry context for which it is drawn
/// @param momentumScale The scale of the momentum
/// @param locErrorScale  The scale of the local error
/// @param angularErrorScale The sclae of the angular error
/// @param drawParameterSurface The indicator whether to draw the surface
/// @param drawMeasurement The indicator whether to draw the (calibrated)
/// measurement
/// @param drawPredictedParameters The indicator whether to draw the predicted
/// track parameters
/// @param drawFilteredParameters The indicator whether to draw the filtered
/// track parameters
/// @param drawSmoothedParameters The indicator whether to draw the smoothed
/// track parameters
/// @param lseg The number of segments for a full arch (if needed)
/// @param scolor the (optional) color of the surface to be written
/// @param mcolor the (optional) color of the (calibrated) measurement to be
/// written
/// @param ppcolor the (optional) color of the predicted track parameters to be
/// written
/// @param fpcolor the (optional) color of the filtered track parameters to be
/// written
/// @param spcolor the (optional) color of the smoothed track parameters to be
/// written
/// @param outOfPlanes The out of plane drawning option for measurment,
/// predicted parameter, filtered parameter and smoothed parameters.
template <typename source_link_t>
static inline void drawMultiTrajectory(
    IVisualization& helper,
    const Acts::MultiTrajectory<source_link_t>& multiTraj,
    const size_t& entryIndex, const GeometryContext& gctx = GeometryContext(),
    double momentumScale = 1., double locErrorScale = 1.,
    double angularErrorScale = 1., bool drawParameterSurface = true,
    bool drawMeasurement = true, bool drawPredictedParameters = true,
    bool drawFilteredParameters = true, bool drawSmoothedParameters = true,
    size_t lseg = 72, const IVisualization::ColorType& scolor = {235, 198, 52},
    const IVisualization::ColorType& mcolor = {218, 165, 32},
    const IVisualization::ColorType& ppcolor = {20, 120, 20},
    const IVisualization::ColorType& fpcolor = {255, 102, 0},
    const IVisualization::ColorType& spcolor = {204, 153, 255},
    std::array<double, 4> outOfPlanes = {-0.01, -0.02, -0.03, -0.04}) {
  // Visit the track states on the trajectory
  multiTraj.visitBackwards(entryIndex, [&](const auto& state) {
    // Only draw the measurement states
    if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
      return true;
    }

    // Use unit scaling factors for the first state
    // @Todo: add parameter for the first state error scaling
    if (state.index() == 0) {
      locErrorScale = 1;
      angularErrorScale = 1;
    }

    // First, if necessary, draw the surface
    if (drawParameterSurface) {
      GeometryVisualization::drawSurface(helper, state.referenceSurface(), gctx,
                                         Transform3D::Identity(), lseg, false,
                                         scolor);
    }

    // Second, if necessary and present, draw the calibrated measurement (only
    // draw 2D measurement here)
    // @Todo: how to draw 1D measurement?
    if (drawMeasurement and state.hasCalibrated() and
        state.calibratedSize() == 2) {
      const Vector2D& lposition = state.calibrated().template head<2>();
      ActsSymMatrixD<2> covariance =
          state.calibratedCovariance().template topLeftCorner<2, 2>();
      drawCovarianceCartesian(helper, lposition, covariance,
                              state.referenceSurface().transform(gctx), {3},
                              locErrorScale, 72, mcolor, outOfPlanes.at(0));
    }

    // Last, if necessary and present, draw the track parameters
    // (a) predicted track parameters
    if (drawPredictedParameters and state.hasPredicted()) {
      drawBoundParameters(
          helper,
          BoundParameters(gctx, state.predictedCovariance(), state.predicted(),
                          state.referenceSurface().getSharedPtr()),
          gctx, momentumScale, locErrorScale, angularErrorScale, false, 72,
          ppcolor, scolor, outOfPlanes.at(1));
    }
    // (b) filtered track parameters
    if (drawFilteredParameters and state.hasFiltered()) {
      drawBoundParameters(
          helper,
          BoundParameters(gctx, state.filteredCovariance(), state.filtered(),
                          state.referenceSurface().getSharedPtr()),
          gctx, momentumScale, locErrorScale, angularErrorScale, false, 72,
          fpcolor, scolor, outOfPlanes.at(2));
    }
    // (c) smoothed track parameters
    if (drawSmoothedParameters and state.hasSmoothed()) {
      drawBoundParameters(
          helper,
          BoundParameters(gctx, state.smoothedCovariance(), state.smoothed(),
                          state.referenceSurface().getSharedPtr()),
          gctx, momentumScale, locErrorScale, angularErrorScale, false, 72,
          spcolor, scolor, outOfPlanes.at(3));
    }
    return true;
  });
}

}  // namespace EventDataVisualization
}  // namespace Acts
