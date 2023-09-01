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
#include "Acts/Vertexing/VertexingError.hpp"

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<double>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    calculate3DDistance(const GeometryContext& gctx,
                        const BoundTrackParameters& trkParams,
                        const Vector3& vtxPos, State& state) const {
  auto res = getDistanceAndMomentum<3>(gctx, trkParams, vtxPos, state);

  if (!res.ok()) {
    return res.error();
  }

  // Return distance
  return res.value().first.norm();
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<Acts::BoundTrackParameters>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    estimate3DImpactParameters(const GeometryContext& gctx,
                               const Acts::MagneticFieldContext& mctx,
                               const BoundTrackParameters& trkParams,
                               const Vector3& vtxPos, State& state) const {
  auto res = getDistanceAndMomentum<3>(gctx, trkParams, vtxPos, state);

  if (!res.ok()) {
    return res.error();
  }

  // Vector pointing from vertex to 3D PCA
  Vector3 deltaR = res.value().first;

  // Get corresponding unit vector
  deltaR.normalize();

  // Momentum direction at vtxPos
  Vector3 momDir = res.value().second;

  // To understand why deltaR and momDir are not orthogonal, let us look at the
  // x-y-plane. Since we computed the 3D PCA, the 2D distance between the vertex
  // and the PCA is not necessarily minimal (see Fig. 4.2 in the reference). As
  // a consequence, the momentum and the vector connecting the vertex and the
  // PCA are not orthogonal to each other.
  Vector3 orthogonalDeltaR = deltaR - (deltaR.dot(momDir)) * momDir;

  // Vector perpendicular to momDir and orthogonalDeltaR
  Vector3 perpDir = momDir.cross(orthogonalDeltaR);

  // Cartesian coordinate system with:
  // -) origin at the vertex position
  // -) z-axis in momentum direction
  // -) x-axis approximately in direction of the 3D PCA (slight deviations
  // because it was modified to make if orthogonal to momDir)
  // -) y-axis is calculated to be orthogonal to x- and z-axis
  // The transformation is represented by a 4x4 matrix with 0 0 0 1 in the last
  // column.
  Transform3 coordinateSystem;
  // First three columns correspond to coordinate system axes
  coordinateSystem.matrix().block<3, 1>(0, 0) = orthogonalDeltaR;
  coordinateSystem.matrix().block<3, 1>(0, 1) = perpDir;
  coordinateSystem.matrix().block<3, 1>(0, 2) = momDir;
  // Fourth column corresponds to origin of the coordinate system
  coordinateSystem.matrix().block<3, 1>(0, 3) = vtxPos;

  // Surface with normal vector in direction of the z axis of coordinateSystem
  std::shared_ptr<PlaneSurface> planeSurface =
      Surface::makeShared<PlaneSurface>(coordinateSystem);

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);
  auto intersection = planeSurface->intersect(gctx, trkParams.position(gctx),
                                              trkParams.direction(), false);
  pOptions.direction =
      Direction::fromScalarZeroAsPositive(intersection.intersection.pathLength);

  // Propagate to the surface; intersection corresponds to an estimate of the 3D
  // PCA. If deltaR and momDir were orthogonal the calculation would be exact.
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
    get3DVertexCompatibility(const GeometryContext& gctx,
                             const BoundTrackParameters* trkParams,
                             const Vector3& vertexPos) const {
  if (trkParams == nullptr) {
    return VertexingError::EmptyInput;
  }

  // surface rotation
  RotationMatrix3 surfaceAxes =
      trkParams->referenceSurface().transform(gctx).rotation();
  // Surface translation
  Vector3 surfaceOrigin =
      trkParams->referenceSurface().transform(gctx).translation();

  // x and y direction of plane
  Vector3 xDir = surfaceAxes.col(0);
  Vector3 yDir = surfaceAxes.col(1);

  // transform vertex position in local plane reference frame
  Vector3 vertexLocPlane = vertexPos - surfaceOrigin;

  // local x/y vertex position
  Vector2 vertexLocXY{vertexLocPlane.dot(xDir), vertexLocPlane.dot(yDir)};

  // track covariance
  if (not trkParams->covariance().has_value()) {
    return VertexingError::NoCovariance;
  }
  auto cov = trkParams->covariance();
  SquareMatrix2 weightXY = cov->block<2, 2>(0, 0).inverse();

  // 2-dim residual
  Vector2 xyRes =
      Vector2(trkParams->parameters()[eX], trkParams->parameters()[eY]) -
      vertexLocXY;

  // return chi2
  return xyRes.dot(weightXY * xyRes);
}

template <typename input_track_t, typename propagator_t,
          typename propagator_options_t>
Acts::Result<double> Acts::ImpactPointEstimator<
    input_track_t, propagator_t,
    propagator_options_t>::performNewtonOptimization(const Vector3& helixCenter,
                                                     const Vector3& vtxPos,
                                                     double phi, double theta,
                                                     double rho) const {
  double sinPhi = std::sin(phi);
  double cosPhi = std::cos(phi);

  int nIter = 0;
  bool hasConverged = false;

  double cotTheta = 1. / std::tan(theta);

  double xO = helixCenter.x();
  double yO = helixCenter.y();
  double zO = helixCenter.z();

  double xVtx = vtxPos.x();
  double yVtx = vtxPos.y();
  double zVtx = vtxPos.z();

  // Iterate until convergence is reached or the maximum amount of iterations is
  // exceeded
  while (!hasConverged && nIter < m_cfg.maxIterations) {
    double derivative = rho * ((xVtx - xO) * cosPhi + (yVtx - yO) * sinPhi +
                               (zVtx - zO + rho * phi * cotTheta) * cotTheta);
    double secDerivative = rho * (-(xVtx - xO) * sinPhi + (yVtx - yO) * cosPhi +
                                  rho * cotTheta * cotTheta);

    if (secDerivative < 0.) {
      return VertexingError::NumericFailure;
    }

    double deltaPhi = -derivative / secDerivative;

    phi += deltaPhi;
    sinPhi = std::sin(phi);
    cosPhi = std::cos(phi);

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
template <unsigned int nDim>
Acts::Result<std::pair<Acts::ActsVector<nDim>, Acts::Vector3>>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    getDistanceAndMomentum(const GeometryContext& gctx,
                           const BoundTrackParameters& trkParams,
                           const ActsVector<nDim>& vtxPos, State& state,
                           const ActsScalar& massHypothesis,
                           const ActsScalar& chargeHypothesis) const {
  if (nDim != 3 and nDim != 4) {
    throw std::invalid_argument(
        "The number of dimensions N must be either 3 or 4 but was set to " +
        std::to_string(nDim) + ".");
  }

  // Reference point R
  Vector3 refPoint = trkParams.referenceSurface().center(gctx);

  // Perigee parameters (parameters of 2D PCA)
  double d0 = trkParams.parameters()[BoundIndices::eBoundLoc0];
  double z0 = trkParams.parameters()[BoundIndices::eBoundLoc1];
  double phiP = trkParams.parameters()[BoundIndices::eBoundPhi];
  double theta = trkParams.parameters()[BoundIndices::eBoundTheta];
  double qOvP = trkParams.parameters()[BoundIndices::eBoundQOverP];
  double tP = trkParams.parameters()[BoundIndices::eBoundTime];
  // Functions of Perigee parameters for later use
  double sinTheta = std::sin(theta);
  double cotTheta = 1. / std::tan(theta);

  // Set optimization variable phi to the angle at the 2D PCA as a first guess.
  // Note that phi corresponds to phiV in the reference.
  double phi = phiP;

  // B-field z-component at the reference position.
  // Note that we assume a constant B field here!
  auto fieldRes = m_cfg.bField->getField(refPoint, state.fieldCache);
  if (!fieldRes.ok()) {
    return fieldRes.error();
  }
  double bZ = (*fieldRes)[eZ];

  // Signed radius of the helix on which the particle moves
  double rho = 0;
  // Curvature is infinite w/o b field
  if (bZ == 0. || std::abs(qOvP) < m_cfg.minQoP) {
    rho = m_cfg.maxRho;
  } else {
    rho = sinTheta * (1. / qOvP) / bZ;
  }

  // Position of the helix center.
  // We can set the z-position to a convenient value since it is not fixed by
  // the Perigee parameters. Note that phi = phiP because we did not start the
  // optimization yet.
  Vector3 helixCenter =
      refPoint + Vector3(-(d0 - rho) * std::sin(phi),
                         (d0 - rho) * std::cos(phi), z0 + rho * phi * cotTheta);

  // Use Newton optimization method to iteratively change phi until we arrive at
  // the 3D PCA
  auto res = performNewtonOptimization(helixCenter, vtxPos.template head<3>(),
                                       phi, theta, rho);
  if (!res.ok()) {
    return res.error();
  }
  // Set new phi value
  phi = *res;

  double cosPhi = std::cos(phi);
  double sinPhi = std::sin(phi);

  // Momentum direction at the 3D PCA.
  // Note that we have thetaV = thetaP = theta since the polar angle does not
  // change in a constant B field.
  Vector3 momDir =
      Vector3(cosPhi * sinTheta, sinPhi * sinTheta, std::cos(theta));

  // 3D PCA (point P' in the reference)
  ActsVector<nDim> pca;
  pca.template head<3>() =
      helixCenter + rho * Vector3(-sinPhi, cosPhi, -cotTheta * phi);

  if (nDim == 4) {
    // TODO: get charge and mass hypothesis from track parameters once this is
    // possible
    ActsScalar p = std::abs(chargeHypothesis / qOvP);
    // Speed in units of c
    ActsScalar beta = p / std::hypot(p, massHypothesis);

    pca[3] = tP - rho / (beta * sinTheta) * (phi - phiP);
  }
  // Vector pointing from the vertex position to the 3D PCA
  ActsVector<nDim> deltaR = pca - vtxPos;

  return std::make_pair(deltaR, momDir);
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
  pOptions.direction = Direction::Backward;

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

  SquareMatrix2 vrtXYCov = vtx.covariance().template block<2, 2>(0, 0);

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

  SquareMatrix2 covPerigeeZ0Theta;
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
    getLifetimeSignOfTrack(const BoundTrackParameters& track,
                           const Vertex<input_track_t>& vtx,
                           const Acts::Vector3& direction,
                           const GeometryContext& gctx,
                           const MagneticFieldContext& mctx) const {
  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtx.position());

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);
  pOptions.direction = Direction::Backward;

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
  pOptions.direction = Direction::Backward;

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