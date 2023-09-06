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
template <unsigned int nDim>
Acts::Result<double>
Acts::ImpactPointEstimator<input_track_t, propagator_t, propagator_options_t>::
    getVertexCompatibility(const GeometryContext& gctx,
                           const BoundTrackParameters* trkParams,
                           const ActsVector<nDim>& vertexPos) const {
  if (nDim != 3 and nDim != 4) {
    throw std::invalid_argument(
        "The number of dimensions N must be either 3 or 4 but was set to " +
        std::to_string(nDim) + ".");
  }

  if (trkParams == nullptr) {
    return VertexingError::EmptyInput;
  }

  // Retrieve weight matrix of the track's local x- and y-coordinate. For this,
  // the covariance needs to be set.
  if (not trkParams->covariance().has_value()) {
    return VertexingError::NoCovariance;
  }
  auto covMat = trkParams->covariance();
  ActsSquareMatrix<nDim - 1> subCovMat;
  subCovMat.template block<2, 2>(0, 0) = covMat->block<2, 2>(0, 0);
  if (nDim == 4) {
    subCovMat.template block<2, 1>(0, 2) = covMat->block<2, 1>(0, eBoundTime);
    subCovMat.template block<1, 2>(2, 0) = covMat->block<1, 2>(eBoundTime, 0);
    subCovMat(2, 2) = covMat.value()(eBoundTime, eBoundTime);
  }
  ActsSquareMatrix<nDim - 1> weight = subCovMat.inverse();

  std::cout << "\n covMat:\n" << covMat.value();
  std::cout << "\n subCovMat:\n" << subCovMat;

  // Orientation of the surface (i.e., axes of the corresponding coordinate
  // system)
  RotationMatrix3 surfaceAxes =
      trkParams->referenceSurface().transform(gctx).rotation();
  // Origin of the surface coordinate system
  Vector3 surfaceOrigin =
      trkParams->referenceSurface().transform(gctx).translation();

  // x- and y-axis of the surface coordinate system
  Vector3 xAxis = surfaceAxes.col(0);
  Vector3 yAxis = surfaceAxes.col(1);

  // Vector pointing from the surface origin to the vertex position
  // TODO: The vertex should always be at the surfaceOrigin since the
  // track parameters should be obtained by estimate3DImpactParameters.
  // Therefore, originToVertex should always be 0, which is currently not the
  // case.
  Vector3 originToVertex = vertexPos.template head<3>() - surfaceOrigin;

  // x-, y-, and possibly time-coordinate of the vertex and the track in the
  // surface coordinate system
  ActsVector<nDim - 1> localVertexCoords;
  localVertexCoords.template head<2>() =
      Vector2(originToVertex.dot(xAxis), originToVertex.dot(yAxis));

  ActsVector<nDim - 1> localTrackCoords;
  localTrackCoords.template head<2>() =
      Vector2(trkParams->parameters()[eX], trkParams->parameters()[eY]);

  // Fill time coordinates if we check the 4D vertex compatibility
  if (nDim == 4) {
    localVertexCoords(2) = vertexPos(3);
    localTrackCoords(2) = trkParams->parameters()[eTime];
  }

  // residual
  ActsVector<nDim - 1> residual = localTrackCoords - localVertexCoords;

  // return chi2
  return residual.dot(weight * residual);
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
    getImpactParameters(const BoundTrackParameters& track,
                        const Vertex<input_track_t>& vtx,
                        const GeometryContext& gctx,
                        const Acts::MagneticFieldContext& mctx,
                        bool calculateTimeIP) const {
  const std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtx.position());

  // Create propagator options
  propagator_options_t pOptions(gctx, mctx);
  auto intersection = perigeeSurface->intersect(gctx, track.position(gctx),
                                                track.direction(), false);
  pOptions.direction =
      Direction::fromScalarZeroAsPositive(intersection.intersection.pathLength);

  // Do the propagation to linPoint
  auto result = m_cfg.propagator->propagate(track, *perigeeSurface, pOptions);

  if (!result.ok()) {
    return result.error();
  }

  const auto& propRes = *result;

  // Check if the covariance matrix of the Perigee parameters exists
  if (not propRes.endParameters->covariance().has_value()) {
    return VertexingError::NoCovariance;
  }

  // Extract Perigee parameters and corresponding covariance matrix
  const auto& params = propRes.endParameters->parameters();
  const double d0 = params[BoundIndices::eBoundLoc0];
  const double z0 = params[BoundIndices::eBoundLoc1];
  const auto& perigeeCov = *(propRes.endParameters->covariance());

  // Vertex variances
  // TODO: By looking at sigmaD0 and sigmaZ0 we neglect the offdiagonal terms
  // (i.e., we approximate the vertex as a sphere rather than an ellipsoid).
  // Using the full covariance matrix might furnish better results.
  double vtxVarX = vtx.covariance()(eX, eX);
  double vtxVarY = vtx.covariance()(eY, eY);
  double vtxVarZ = vtx.covariance()(eZ, eZ);

  ImpactParametersAndSigma ipAndSigma;

  ipAndSigma.d0 = d0;
  // Variance of the vertex extent in the x-y-plane
  double vtxVar2DExtent = vtxVarX + vtxVarY;
  // TODO: vtxVar2DExtent, vtxVarZ, and vtxVarT should always be >= 0. We need
  // to throw an error here once
  // https://github.com/acts-project/acts/issues/2231 is resolved.
  if (vtxVar2DExtent > 0) {
    ipAndSigma.sigmaD0 =
        std::sqrt(vtxVar2DExtent + perigeeCov(BoundIndices::eBoundLoc0,
                                              BoundIndices::eBoundLoc0));
  } else {
    ipAndSigma.sigmaD0 = std::sqrt(
        perigeeCov(BoundIndices::eBoundLoc0, BoundIndices::eBoundLoc0));
  }

  ipAndSigma.z0 = z0;
  if (vtxVarZ > 0) {
    ipAndSigma.sigmaZ0 =
        std::sqrt(vtxVarZ + perigeeCov(BoundIndices::eBoundLoc1,
                                       BoundIndices::eBoundLoc1));
  } else {
    ipAndSigma.sigmaZ0 = std::sqrt(
        perigeeCov(BoundIndices::eBoundLoc1, BoundIndices::eBoundLoc1));
  }

  if (calculateTimeIP) {
    ipAndSigma.deltaT = std::abs(vtx.time() - params[BoundIndices::eBoundTime]);
    double vtxVarT = vtx.fullCovariance()(eTime, eTime);
    if (vtxVarT > 0) {
      ipAndSigma.sigmaDeltaT =
          std::sqrt(vtxVarT + perigeeCov(BoundIndices::eBoundTime,
                                         BoundIndices::eBoundTime));
    } else {
      ipAndSigma.sigmaDeltaT = std::sqrt(
          perigeeCov(BoundIndices::eBoundTime, BoundIndices::eBoundTime));
    }
  }

  return ipAndSigma;
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