// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<double> Acts::ImpactPoint3dEstimator<
    input_track_t, propagator_t,
    propagator_options_t>::calculateDistance(const GeometryContext& gctx,
                                             const BoundParameters& trkParams,
                                             const Vector3D& vtxPos) const {
  Vector3D deltaR;
  Vector3D momDir;

  auto res = getDistanceAndMomentum(gctx, trkParams, vtxPos, deltaR, momDir);

  if (!res.ok()) {
    return res.error();
  }

  // Return distance
  return deltaR.norm();
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<std::unique_ptr<const Acts::BoundParameters>>
Acts::ImpactPoint3dEstimator<input_track_t, propagator_t,
                             propagator_options_t>::
    getParamsAtClosestApproach(const GeometryContext& gctx,
                               const BoundParameters& trkParams,
                               const Vector3D& vtxPos) const {
  Vector3D deltaR;
  Vector3D momDir;

  auto res = getDistanceAndMomentum(gctx, trkParams, vtxPos, deltaR, momDir);

  if (!res.ok()) {
    return res.error();
  }

  // Normalize deltaR
  deltaR.normalize();

  // corrected deltaR for small deviations from orthogonality
  Vector3D corrDeltaR = deltaR - (deltaR.dot(momDir)) * momDir;

  // get perpendicular direction vector
  Vector3D perpDir = momDir.cross(corrDeltaR);

  Transform3D thePlane;
  // rotation matrix
  thePlane.matrix().block(0, 0, 3, 1) = corrDeltaR;
  thePlane.matrix().block(0, 1, 3, 1) = perpDir;
  thePlane.matrix().block(0, 2, 3, 1) = momDir;
  // translation
  thePlane.matrix().block(0, 3, 3, 1) = vtxPos;

  std::shared_ptr<PlaneSurface> planeSurface =
      Surface::makeShared<PlaneSurface>(
          std::make_shared<Transform3D>(thePlane));

  // Do the propagation to linPointPos
  auto result =
      m_cfg.propagator->propagate(trkParams, *planeSurface, m_cfg.pOptions);
  if (result.ok()) {
    return std::move((*result).endParameters);
  } else {
    return result.error();
  }
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<double> Acts::ImpactPoint3dEstimator<input_track_t, propagator_t,
                                                  propagator_options_t>::
    getVertexCompatibility(const GeometryContext& gctx,
                           const BoundParameters* trkParams,
                           const Vector3D& vertexPos) const {
  if (trkParams == nullptr) {
    return VertexingError::EmptyInput;
  }

  // surface rotation
  RotationMatrix3D myRotation =
      trkParams->referenceSurface().transform(gctx).rotation();
  // Surface translation
  Vector3D myTranslation =
      trkParams->referenceSurface().transform(gctx).translation();

  // x and y direction of plane
  Vector3D xDirPlane = myRotation.col(0);
  Vector3D yDirPlane = myRotation.col(1);

  // transform vertex position in local plane reference frame
  Vector3D vertexLocPlane = vertexPos - myTranslation;

  // local x/y vertex position
  Vector2D vertexLocXY{vertexLocPlane.dot(xDirPlane),
                       vertexLocPlane.dot(yDirPlane)};

  // track covariance
  auto cov = trkParams->covariance();
  ActsSymMatrixD<2> myWeightXY = cov->block<2, 2>(0, 0).inverse();

  // 2-dim residual
  Vector2D myXYpos =
      Vector2D(trkParams->parameters()[eX], trkParams->parameters()[eY]) -
      vertexLocXY;

  // return chi2
  return myXYpos.dot(myWeightXY * myXYpos);
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<double> Acts::ImpactPoint3dEstimator<
    input_track_t, propagator_t,
    propagator_options_t>::performNewtonApproximation(const Vector3D& trkPos,
                                                      const Vector3D& vtxPos,
                                                      double phi, double theta,
                                                      double r) const {
  double sinNewPhi = -std::sin(phi);
  double cosNewPhi = std::cos(phi);

  int nIter = 0;
  bool hasConverged = false;

  double cotTheta = 1. / std::tan(theta);

  // start iteration until convergence or max iteration reached
  while (!hasConverged && nIter < m_cfg.maxIterations) {
    double x0 = trkPos.x();
    double y0 = trkPos.y();
    double z0 = trkPos.z();

    double xc = vtxPos.x();
    double yc = vtxPos.y();
    double zc = vtxPos.z();

    double derivative = (x0 - xc) * (-r * cosNewPhi) +
                        (y0 - yc) * r * sinNewPhi +
                        (z0 - zc - r * phi * cotTheta) * (-r * cotTheta);
    double secDerivative = r * (-(x0 - xc) * sinNewPhi - (y0 - yc) * cosNewPhi +
                                r * cotTheta * cotTheta);

    if (secDerivative < 0.) {
      return VertexingError::NumericFailure;
    }

    double deltaPhi = -derivative / secDerivative;

    phi += deltaPhi;
    sinNewPhi = -std::sin(phi);
    cosNewPhi = std::cos(phi);

    nIter += 1;

    if (std::abs(deltaPhi) < m_cfg.precision) {
      hasConverged = true;
    }
  }  // end while loop

  if (!hasConverged) {
    // max iterations reached but did not converge
    return VertexingError::NotConverged;
  }
  return phi;
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<void> Acts::ImpactPoint3dEstimator<input_track_t, propagator_t,
                                                propagator_options_t>::
    getDistanceAndMomentum(const GeometryContext& gctx,
                           const BoundParameters& trkParams,
                           const Vector3D& vtxPos, Vector3D& deltaR,
                           Vector3D& momDir) const {
  Vector3D trkSurfaceCenter = trkParams.referenceSurface().center(gctx);

  double d0 = trkParams.parameters()[ParID_t::eLOC_D0];
  double z0 = trkParams.parameters()[ParID_t::eLOC_Z0];
  double phi = trkParams.parameters()[ParID_t::ePHI];
  double theta = trkParams.parameters()[ParID_t::eTHETA];
  double qOvP = trkParams.parameters()[ParID_t::eQOP];

  double cotTheta = 1. / std::tan(theta);

  // get B-field z-component at current position
  double bZ = m_cfg.bField.getField(trkSurfaceCenter)[eZ];

  // The radius
  double r;
  // Curvature is infinite w/o b field
  if (bZ == 0. || std::abs(qOvP) < m_cfg.minQoP) {
    r = m_cfg.maxRho;
  } else {
    // signed(!) r
    r = std::sin(theta) * (1. / qOvP) / bZ;
  }

  Vector3D vec0 = trkSurfaceCenter + Vector3D(-(d0 - r) * std::sin(phi),
                                              (d0 - r) * std::cos(phi),
                                              z0 + r * phi * cotTheta);

  // Perform newton approximation method
  // this will change the value of phi
  auto res = performNewtonApproximation(vec0, vtxPos, phi, theta, r);
  if (!res.ok()) {
    return res.error();
  }
  // Set new phi value
  phi = *res;

  // Set momentum direction
  momDir = Vector3D(std::sin(theta) * std::cos(phi),
                    std::sin(phi) * std::sin(theta), std::cos(theta));

  // point of closest approach in 3D
  Vector3D pointCA3d =
      vec0 + r * Vector3D(-std::sin(phi), std::cos(phi), -cotTheta * phi);

  // Set deltaR
  deltaR = pointCA3d - vtxPos;

  return {};
}
