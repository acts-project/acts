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
Acts::Result<double>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    calculate3dDistance(const GeometryContext& gctx,
                        const BoundTrackParameters& trkParams,
                        const Vector3& vtxPos, State& state) const {
  Vector3 deltaR;
  Vector3 momDir;

  auto res =
      getDistanceAndMomentum(gctx, trkParams, vtxPos, deltaR, momDir, state);

  if (!res.ok()) {
    return res.error();
  }

  // Return distance
  return deltaR.norm();
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<Acts::BoundTrackParameters>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    estimate3DImpactParameters(const GeometryContext& gctx,
                               const Acts::MagneticFieldContext& mctx,
                               const BoundTrackParameters& trkParams,
                               const Vector3& vtxPos, State& state) const {
  Vector3 deltaR;
  Vector3 momDir;

  auto res =
      getDistanceAndMomentum(gctx, trkParams, vtxPos, deltaR, momDir, state);

  if (!res.ok()) {
    return res.error();
  }

  // Normalize deltaR
  deltaR.normalize();

  // corrected deltaR for small deviations from orthogonality
  Vector3 corrDeltaR = deltaR - (deltaR.dot(momDir)) * momDir;

  // get perpendicular direction vector
  Vector3 perpDir = momDir.cross(corrDeltaR);

  Transform3 thePlane;
  // rotation matrix
  thePlane.matrix().block(0, 0, 3, 1) = corrDeltaR;
  thePlane.matrix().block(0, 1, 3, 1) = perpDir;
  thePlane.matrix().block(0, 2, 3, 1) = momDir;
  // translation
  thePlane.matrix().block(0, 3, 3, 1) = vtxPos;

  std::shared_ptr<PlaneSurface> planeSurface =
      Surface::makeShared<PlaneSurface>(thePlane);

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);
  pOptions.direction = NavigationDirection::Backward;

  // Do the propagation to linPointPos
  auto result = m_cfg.propagator->propagate(trkParams, *planeSurface, pOptions);
  if (result.ok()) {
    return *result->endParameters;
  } else {
    return result.error();
  }
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<double>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    get3dVertexCompatibility(const GeometryContext& gctx,
                             const BoundTrackParameters* trkParams,
                             const Vector3& vertexPos) const {
  if (trkParams == nullptr) {
    return VertexingError::EmptyInput;
  }

  // surface rotation
  RotationMatrix3 myRotation =
      trkParams->referenceSurface().transform(gctx).rotation();
  // Surface translation
  Vector3 myTranslation =
      trkParams->referenceSurface().transform(gctx).translation();

  // x and y direction of plane
  Vector3 xDirPlane = myRotation.col(0);
  Vector3 yDirPlane = myRotation.col(1);

  // transform vertex position in local plane reference frame
  Vector3 vertexLocPlane = vertexPos - myTranslation;

  // local x/y vertex position
  Vector2 vertexLocXY{vertexLocPlane.dot(xDirPlane),
                      vertexLocPlane.dot(yDirPlane)};

  // track covariance
  if (not trkParams->covariance().has_value()) {
    return VertexingError::NoCovariance;
  }
  auto cov = trkParams->covariance();
  SymMatrix2 myWeightXY = cov->block<2, 2>(0, 0).inverse();

  // 2-dim residual
  Vector2 myXYpos =
      Vector2(trkParams->parameters()[eX], trkParams->parameters()[eY]) -
      vertexLocXY;

  // return chi2
  return myXYpos.dot(myWeightXY * myXYpos);
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<double> Acts::ImpactPointEstimator<
    input_track_t, propagator_t,
    propagator_options_t>::performNewtonApproximation(const Vector3& trkPos,
                                                      const Vector3& vtxPos,
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
Acts::Result<void>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    getDistanceAndMomentum(const GeometryContext& gctx,
                           const BoundTrackParameters& trkParams,
                           const Vector3& vtxPos, Vector3& deltaR,
                           Vector3& momDir, State& state) const {
  Vector3 trkSurfaceCenter = trkParams.referenceSurface().center(gctx);

  double d0 = trkParams.parameters()[BoundIndices::eBoundLoc0];
  double z0 = trkParams.parameters()[BoundIndices::eBoundLoc1];
  double phi = trkParams.parameters()[BoundIndices::eBoundPhi];
  double theta = trkParams.parameters()[BoundIndices::eBoundTheta];
  double qOvP = trkParams.parameters()[BoundIndices::eBoundQOverP];

  double sinTheta = std::sin(theta);

  double cotTheta = 1. / std::tan(theta);

  // get B-field z-component at current position
  auto fieldRes = m_cfg.bField->getField(trkSurfaceCenter, state.fieldCache);
  if (!fieldRes.ok()) {
    return fieldRes.error();
  }
  double bZ = (*fieldRes)[eZ];

  // The radius
  double r = 0;
  // Curvature is infinite w/o b field
  if (bZ == 0. || std::abs(qOvP) < m_cfg.minQoP) {
    r = m_cfg.maxRho;
  } else {
    // signed(!) r
    r = sinTheta * (1. / qOvP) / bZ;
  }

  Vector3 vec0 = trkSurfaceCenter + Vector3(-(d0 - r) * std::sin(phi),
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

  double cosPhi = std::cos(phi);
  double sinPhi = std::sin(phi);

  // Set momentum direction
  momDir = Vector3(sinTheta * cosPhi, sinPhi * sinTheta, std::cos(theta));

  // point of closest approach in 3D
  Vector3 pointCA3d = vec0 + r * Vector3(-sinPhi, cosPhi, -cotTheta * phi);

  // Set deltaR
  deltaR = pointCA3d - vtxPos;

  return {};
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<Acts::ImpactParametersAndSigma>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    estimateImpactParameters(const BoundTrackParameters& track,
                             const Vertex<input_track_t>& vtx,
                             const GeometryContext& gctx,
                             const Acts::MagneticFieldContext& mctx) const {
  // estimating the d0 and its significance by propagating the trajectory state
  // towards
  // the vertex position. By this time the vertex should NOT contain this
  // trajectory anymore
  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtx.position());

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);
  pOptions.direction = NavigationDirection::Backward;

  // Do the propagation to linPoint
  auto result = m_cfg.propagator->propagate(track, *perigeeSurface, pOptions);

  if (!result.ok()) {
    return result.error();
  }

  const auto& propRes = *result;
  const auto& params = propRes.endParameters->parameters();
  const double d0 = params[BoundIndices::eBoundLoc0];
  const double z0 = params[BoundIndices::eBoundLoc1];
  const double phi = params[BoundIndices::eBoundPhi];
  const double theta = params[BoundIndices::eBoundTheta];

  const double sinPhi = std::sin(phi);
  const double sinTheta = std::sin(theta);
  const double cosPhi = std::cos(phi);
  const double cosTheta = std::cos(theta);

  SymMatrix2 vrtXYCov = vtx.covariance().template block<2, 2>(0, 0);

  // Covariance of perigee parameters after propagation to perigee surface
  if (not propRes.endParameters->covariance().has_value()) {
    return VertexingError::NoCovariance;
  }
  const auto& perigeeCov = *(propRes.endParameters->covariance());

  Vector2 d0JacXY(-sinPhi, cosPhi);

  ImpactParametersAndSigma newIPandSigma;

  newIPandSigma.IPd0 = d0;
  double d0_PVcontrib = d0JacXY.transpose() * (vrtXYCov * d0JacXY);
  if (d0_PVcontrib >= 0) {
    newIPandSigma.sigmad0 =
        std::sqrt(d0_PVcontrib + perigeeCov(BoundIndices::eBoundLoc0,
                                            BoundIndices::eBoundLoc0));
    newIPandSigma.PVsigmad0 = std::sqrt(d0_PVcontrib);
  } else {
    newIPandSigma.sigmad0 = std::sqrt(
        perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc0));
    newIPandSigma.PVsigmad0 = 0;
  }

  SymMatrix2 covPerigeeZ0Theta;
  covPerigeeZ0Theta(0, 0) =
      perigeeCov(BoundIndices::eBoundLoc1, BoundIndices::eBoundLoc1);
  covPerigeeZ0Theta(0, 1) =
      perigeeCov(BoundIndices::eBoundLoc1, BoundIndices::eBoundTheta);
  covPerigeeZ0Theta(1, 0) =
      perigeeCov(BoundIndices::eBoundTheta, BoundIndices::eBoundLoc1);
  covPerigeeZ0Theta(1, 1) =
      perigeeCov(BoundIndices::eBoundTheta, BoundIndices::eBoundTheta);

  double vtxZZCov = vtx.covariance()(eZ, eZ);

  Vector2 z0JacZ0Theta(sinTheta, z0 * cosTheta);

  if (vtxZZCov >= 0) {
    newIPandSigma.IPz0SinTheta = z0 * sinTheta;
    newIPandSigma.sigmaz0SinTheta = std::sqrt(
        z0JacZ0Theta.transpose() * (covPerigeeZ0Theta * z0JacZ0Theta) +
        sinTheta * vtxZZCov * sinTheta);

    newIPandSigma.PVsigmaz0SinTheta = std::sqrt(sinTheta * vtxZZCov * sinTheta);
    newIPandSigma.IPz0 = z0;
    newIPandSigma.sigmaz0 =
        std::sqrt(vtxZZCov + perigeeCov(BoundIndices::eBoundLoc1,
                                        BoundIndices::eBoundLoc1));
    newIPandSigma.PVsigmaz0 = std::sqrt(vtxZZCov);
  } else {
    // Remove contribution from PV
    newIPandSigma.IPz0SinTheta = z0 * sinTheta;
    double sigma2z0sinTheta =
        (z0JacZ0Theta.transpose() * (covPerigeeZ0Theta * z0JacZ0Theta));
    newIPandSigma.sigmaz0SinTheta = std::sqrt(sigma2z0sinTheta);
    newIPandSigma.PVsigmaz0SinTheta = 0;

    newIPandSigma.IPz0 = z0;
    newIPandSigma.sigmaz0 = std::sqrt(
        perigeeCov(BoundIndices::eBoundLoc1, BoundIndices::eBoundLoc1));
    newIPandSigma.PVsigmaz0 = 0;
  }

  return newIPandSigma;
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<std::pair<double, double>>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    getLifetimesSignOfTrack(const BoundTrackParameters& track,
                            const Vertex<input_track_t>& vtx,
                            const Acts::Vector3& direction,
                            const GeometryContext& gctx,
                            const MagneticFieldContext& mctx) const {
  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtx.position());

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);
  pOptions.direction = NavigationDirection::Backward;

  // Do the propagation to the perigeee
  auto result = m_cfg.propagator->propagate(track, *perigeeSurface, pOptions);

  if (!result.ok()) {
    return result.error();
  }

  const auto& propRes = *result;
  const auto& params = propRes.endParameters->parameters();
  const double d0 = params[BoundIndices::eBoundLoc0];
  const double z0 = params[BoundIndices::eBoundLoc1];
  const double phi = params[BoundIndices::eBoundPhi];
  const double theta = params[BoundIndices::eBoundTheta];

  double vs = std::sin(std::atan2(direction[1], direction[0]) - phi) * d0;
  double eta = -std::log(std::tan(theta / 2.));
  double dir_eta = VectorHelpers::eta(direction);

  double zs = (dir_eta - eta) * z0;

  std::pair<double, double> vszs;

  vszs.first = vs >= 0. ? 1. : -1.;
  vszs.second = zs >= 0. ? 1. : -1.;

  return vszs;
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<double>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    get3DLifetimeSignOfTrack(const BoundTrackParameters& track,
                             const Vertex<input_track_t>& vtx,
                             const Acts::Vector3& direction,
                             const GeometryContext& gctx,
                             const MagneticFieldContext& mctx) const {
  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtx.position());

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);
  pOptions.direction = NavigationDirection::Backward;

  // Do the propagation to the perigeee
  auto result = m_cfg.propagator->propagate(track, *perigeeSurface, pOptions);

  if (!result.ok()) {
    return result.error();
  }

  const auto& propRes = *result;
  const auto& params = propRes.endParameters;
  const Vector3 trkpos = params->position(gctx);
  const Vector3 trkmom = params->momentum();

  double sign =
      (direction.cross(trkmom)).dot(trkmom.cross(vtx.position() - trkpos));

  return sign >= 0. ? 1. : -1.;
}
