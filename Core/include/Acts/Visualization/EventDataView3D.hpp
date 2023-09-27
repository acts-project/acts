// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/Polyhedron.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/detail/FacesHelper.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Visualization/GeometryView3D.hpp"
#include "Acts/Visualization/IVisualization3D.hpp"
#include "Acts/Visualization/ViewConfig.hpp"

#include <optional>

namespace Acts {

static ViewConfig s_viewParameter = ViewConfig({0, 0, 255});
static ViewConfig s_viewMeasurement = ViewConfig({255, 102, 0});
static ViewConfig s_viewPredicted = ViewConfig({51, 204, 51});
static ViewConfig s_viewFiltered = ViewConfig({255, 255, 0});
static ViewConfig s_viewSmoothed = ViewConfig({0, 102, 255});

struct EventDataView3D {
  /// Helper to find the egen values and corr angle
  ///
  /// @param covariance The covariance matrix
  static inline std::array<double, 3> decomposeCovariance(
      const ActsSymMatrix<2>& covariance) {
    double c00 = covariance(eBoundLoc0, eBoundLoc0);
    double c01 = covariance(eBoundLoc0, eBoundLoc1);
    double c11 = covariance(eBoundLoc1, eBoundLoc1);

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
  /// @param offset The out of plane offset for visibility
  /// @param lposition The local anker point of the ellipse
  /// @param transform The transform to global
  static inline std::vector<Vector3> createEllipse(
      double lambda0, double lambda1, double theta, size_t lseg, double offset,
      const Vector2& lposition = Vector2(0., 0.),
      const Transform3& transform = Transform3::Identity()) {
    double ctheta = std::cos(theta);
    double stheta = std::sin(theta);

    double l1sq = std::sqrt(lambda0);
    double l2sq = std::sqrt(lambda1);

    // Now generate the ellipse points
    std::vector<Vector3> ellipse;
    ellipse.reserve(lseg);
    double thetaStep = 2 * M_PI / lseg;
    for (size_t it = 0; it < lseg; ++it) {
      double phi = -M_PI + it * thetaStep;
      double cphi = std::cos(phi);
      double sphi = std::sin(phi);
      double x = lposition.x() + (l1sq * ctheta * cphi - l2sq * stheta * sphi);
      double y = lposition.y() + (l1sq * stheta * cphi + l2sq * ctheta * sphi);
      ellipse.push_back(transform * Vector3(x, y, offset));
    }
    return ellipse;
  }

  /// Helper method to draw error ellipse
  ///
  /// @param helper [in, out] The visualization helper
  /// @param lposition The local position
  /// @param covariance The covariance matrix
  /// @param transform The reference Frame transform
  /// @param locErrorScale The local Error scale
  /// @param viewConfig The visualization parameters
  static void drawCovarianceCartesian(
      IVisualization3D& helper, const Vector2& lposition,
      const SymMatrix2& covariance, const Transform3& transform,
      double locErrorScale = 1, const ViewConfig& viewConfig = s_viewParameter);

  /// Helper method to draw error cone of a direction
  ///
  /// @param helper [in, out] The visualization helper
  /// @param position Where the cone originates from
  /// @param direction The direction parameters
  /// @param covariance The 2x2 covariance matrix for phi/theta
  /// @param directionScale The direction arror length
  /// @param angularErrorScale The local Error scale
  /// @param viewConfig The visualization parameters
  static void drawCovarianceAngular(
      IVisualization3D& helper, const Vector3& position,
      const Vector3& direction, const ActsSymMatrix<2>& covariance,
      double directionScale = 1, double angularErrorScale = 1,
      const ViewConfig& viewConfig = s_viewParameter);

  /// Helper method to draw bound parameters object
  ///
  /// @param helper [in, out] The visualization helper
  /// @param parameters The bound parameters to be drawn
  /// @param gctx The geometry context for which it is drawn
  /// @param momentumScale The scale of the momentum
  /// @param locErrorScale  The scale of the local error
  /// @param angularErrorScale The sclae of the angular error
  /// @param parConfig The visualization options for the parameter
  /// @param covConfig The visualization option for the covariance
  /// @param surfConfig The visualization option for the surface
  template <typename parameters_t>
  static inline void drawBoundTrackParameters(
      IVisualization3D& helper, const parameters_t& parameters,
      const GeometryContext& gctx = GeometryContext(),
      double momentumScale = 1., double locErrorScale = 1.,
      double angularErrorScale = 1.,
      const ViewConfig& parConfig = s_viewParameter,
      const ViewConfig& covConfig = s_viewParameter,
      const ViewConfig& surfConfig = s_viewSensitive) {
    if (surfConfig.visible) {
      GeometryView3D::drawSurface(helper, parameters.referenceSurface(), gctx,
                                  Transform3::Identity(), surfConfig);
    }

    // Draw the parameter shaft and cone
    auto position = parameters.position(gctx);
    auto direction = parameters.unitDirection();
    double p = parameters.absoluteMomentum();

    ViewConfig lparConfig = parConfig;
    lparConfig.lineThickness = 0.05;
    Vector3 parLength = p * momentumScale * direction;

    GeometryView3D::drawArrowBackward(
        helper, position, position + 0.5 * parLength, 100., 1.0, lparConfig);

    GeometryView3D::drawArrowForward(helper, position + 0.5 * parLength,
                                     position + parLength, 4., 2.5, lparConfig);

    if (parameters.covariance().has_value()) {
      auto paramVec = parameters.parameters();
      auto lposition = paramVec.template block<2, 1>(0, 0);

      // Draw the local covariance
      const auto& covariance = *parameters.covariance();
      drawCovarianceCartesian(helper, lposition,
                              covariance.template block<2, 2>(0, 0),
                              parameters.referenceSurface().transform(gctx),
                              locErrorScale, covConfig);

      drawCovarianceAngular(
          helper, position, direction, covariance.template block<2, 2>(2, 2),
          0.9 * p * momentumScale, angularErrorScale, covConfig);
    }
  }

  /// Helper method to draw one trajectory stored in a MultiTrajectory object
  ///
  /// @param helper [in, out] The visualization helper
  /// @param multiTraj The MultiTrajectory storing the trajectory to be drawn
  /// @param entryIndex The trajectory entry index
  /// @param gctx The geometry context for which it is drawn
  /// @param momentumScale The scale of the momentum
  /// @param locErrorScale  The scale of the local error
  /// @param angularErrorScale The sclae of the angular error
  /// @param surfaceConfig The visualization options for the surface
  /// @param measurementConfig The visualization options for the measurement
  /// @param predictedConfig The visualization options for the predicted
  /// measurement
  /// @param filteredConfig The visualization options for the filtered
  /// parameters
  /// @param smoothedConfig The visualization options for the smoothed
  /// parameters
  template <typename D>
  static void drawMultiTrajectory(
      IVisualization3D& helper, const Acts::MultiTrajectory<D>& multiTraj,
      const size_t& entryIndex, const GeometryContext& gctx = GeometryContext(),
      double momentumScale = 1., double locErrorScale = 1.,
      double angularErrorScale = 1.,
      const ViewConfig& surfaceConfig = s_viewSensitive,
      const ViewConfig& measurementConfig = s_viewMeasurement,
      const ViewConfig& predictedConfig = s_viewPredicted,
      const ViewConfig& filteredConfig = s_viewFiltered,
      const ViewConfig& smoothedConfig = s_viewSmoothed) {
    // @TODO: Refactor based on Track class

    // Visit the track states on the trajectory
    multiTraj.visitBackwards(entryIndex, [&](const auto& state) {
      // Only draw the measurement states
      if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
        return true;
      }

      // Use smaller scaling factors for the first state
      // @Todo: add parameter for the first state error scaling
      if (state.index() == 0) {
        locErrorScale = locErrorScale * 0.1;
        angularErrorScale = angularErrorScale * 0.1;
      }

      // First, if necessary, draw the surface
      if (surfaceConfig.visible) {
        GeometryView3D::drawSurface(helper, state.referenceSurface(), gctx,
                                    Transform3::Identity(), surfaceConfig);
      }

      // Second, if necessary and present, draw the calibrated measurement (only
      // draw 2D measurement here)
      // @Todo: how to draw 1D measurement?
      if (measurementConfig.visible and state.hasCalibrated() and
          state.calibratedSize() == 2) {
        const Vector2& lposition = state.template calibrated<2>();
        const SymMatrix2 covariance = state.template calibratedCovariance<2>();
        drawCovarianceCartesian(helper, lposition, covariance,
                                state.referenceSurface().transform(gctx),
                                locErrorScale, measurementConfig);
      }

      // Last, if necessary and present, draw the track parameters
      // (a) predicted track parameters
      if (predictedConfig.visible and state.hasPredicted()) {
        drawBoundTrackParameters(
            helper,
            BoundTrackParameters(state.referenceSurface().getSharedPtr(),
                                 state.predicted(),
                                 state.predictedCovariance()),
            gctx, momentumScale, locErrorScale, angularErrorScale,
            predictedConfig, predictedConfig, ViewConfig(false));
      }
      // (b) filtered track parameters
      if (filteredConfig.visible and state.hasFiltered()) {
        drawBoundTrackParameters(
            helper,
            BoundTrackParameters(state.referenceSurface().getSharedPtr(),
                                 state.filtered(), state.filteredCovariance()),
            gctx, momentumScale, locErrorScale, angularErrorScale,
            filteredConfig, filteredConfig, ViewConfig(false));
      }
      // (c) smoothed track parameters
      if (smoothedConfig.visible and state.hasSmoothed()) {
        drawBoundTrackParameters(
            helper,
            BoundTrackParameters(state.referenceSurface().getSharedPtr(),
                                 state.smoothed(), state.smoothedCovariance()),
            gctx, momentumScale, locErrorScale, angularErrorScale,
            smoothedConfig, smoothedConfig, ViewConfig(false));
      }
      return true;
    });
  }
};

}  // namespace Acts
