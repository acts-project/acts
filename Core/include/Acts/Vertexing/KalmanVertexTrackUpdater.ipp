// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

template <typename input_track_t>
void Acts::KalmanVertexTrackUpdater::update(TrackAtVertex<input_track_t>& track,
                                            const Vertex<input_track_t>& vtx) {
  const Vector4 vtxPos = vtx.fullPosition();

  // Get the linearized track
  const LinearizedTrack& linTrack = track.linearizedState;

  // Check if linearized state exists
  if (linTrack.covarianceAtPCA.determinant() == 0.) {
    // Track has no linearized state, returning w/o update
    return;
  }

  // Retrieve linTrack information
  const ActsMatrix<eBoundSize, 4> posJac = linTrack.positionJacobian;
  const ActsMatrix<eBoundSize, 3> momJac = linTrack.momentumJacobian;
  const BoundVector trkParams = linTrack.parametersAtPCA;
  const BoundSquareMatrix trkParamWeight = linTrack.weightAtPCA;

  // Calculate S matrix
  ActsSquareMatrix<3> sMat =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  const BoundVector residual = linTrack.constantTerm;

  // Refit track momentum
  Vector3 newTrkMomentum = sMat * momJac.transpose() * trkParamWeight *
                           (trkParams - residual - posJac * vtxPos);

  // Refit track parameters
  BoundVector newTrkParams(BoundVector::Zero());

  // Get phi and theta and correct for possible periodicity changes
  const auto correctedPhiTheta =
      Acts::detail::normalizePhiTheta(newTrkMomentum(0), newTrkMomentum(1));
  newTrkParams(BoundIndices::eBoundPhi) = correctedPhiTheta.first;     // phi
  newTrkParams(BoundIndices::eBoundTheta) = correctedPhiTheta.second;  // theta
  newTrkParams(BoundIndices::eBoundQOverP) = newTrkMomentum(2);        // qOverP

  // Vertex covariance and weight matrices
  const SquareMatrix4 vtxCov = vtx.fullCovariance();
  const SquareMatrix4 vtxWeight = vtxCov.inverse();

  // Cross covariance matrix between the vertex position and the refitted track
  // momentum
  const ActsMatrix<4, 3> crossCovVP =
      -vtxCov * posJac.transpose() * trkParamWeight * momJac * sMat;

  KalmanVertexUpdater::MatrixCache matrixCache;

  // Now determine the smoothed chi2 of the track in the following
  KalmanVertexUpdater::updatePosition<input_track_t>(
      vtx, linTrack, track.trackWeight, -1, matrixCache);

  // Corresponding weight matrix
  const SquareMatrix4& reducedVtxWeight = matrixCache.newVertexWeight;

  // Difference in positions
  Vector4 posDiff = vtx.fullPosition() - matrixCache.newVertexPos;

  // Get smoothed params
  BoundVector smoothedParams =
      trkParams -
      (residual + posJac * vtx.fullPosition() + momJac * newTrkMomentum);

  // New chi2 to be set later
  double chi2 = posDiff.dot(reducedVtxWeight * posDiff) +
                smoothedParams.dot(trkParamWeight * smoothedParams);

  Acts::BoundMatrix trkCov = detail::calculateTrackCovariance(
      sMat, crossCovVP, vtxWeight, vtxCov, newTrkParams);

  // Create new refitted parameters
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtx.position());

  BoundTrackParameters refittedPerigee =
      BoundTrackParameters(perigeeSurface, newTrkParams, std::move(trkCov),
                           track.fittedParams.particleHypothesis());

  // Set new properties
  track.fittedParams = refittedPerigee;
  track.chi2Track = chi2;
  track.ndf = 2 * track.trackWeight;

  return;
}

inline Acts::BoundMatrix
Acts::KalmanVertexTrackUpdater::detail::calculateTrackCovariance(
    const SquareMatrix3& sMat, const ActsMatrix<4, 3>& crossCovVP,
    const SquareMatrix4& vtxWeight, const SquareMatrix4& vtxCov,
    const BoundVector& newTrkParams) {
  // Now new momentum covariance
  SquareMatrix3 momCov = sMat + crossCovVP.transpose() * vtxWeight * crossCovVP;

  // Covariance matrix of the "free" track parameters, i.e., x, y, z, phi,
  // theta, q/p, and t. Note that the parameters are not actually free: Free
  // parametrization would correspond to x, y, z, t, p_x, p_y, p_z, and q/p.
  ActsSquareMatrix<7> freeTrkCov(ActsSquareMatrix<7>::Zero());

  freeTrkCov.block<3, 3>(0, 0) = vtxCov.template block<3, 3>(0, 0);
  freeTrkCov.block<3, 3>(0, 3) = crossCovVP.template block<3, 3>(0, 0);
  freeTrkCov.block<3, 3>(3, 0) =
      (crossCovVP.template block<3, 3>(0, 0)).transpose();
  freeTrkCov.block<3, 3>(3, 3) = momCov;
  freeTrkCov.block<3, 1>(0, 6) = vtxCov.template block<3, 1>(0, 3);
  freeTrkCov.block<1, 3>(6, 0) = vtxCov.template block<1, 3>(3, 0);
  freeTrkCov.block<3, 1>(3, 6) =
      (crossCovVP.template block<1, 3>(3, 0)).transpose();
  freeTrkCov.block<1, 3>(6, 3) = crossCovVP.template block<1, 3>(3, 0);
  freeTrkCov(6, 6) = vtxCov(3, 3);

  // TODO the Jacobian is not correct, this will be fixed once we make the
  // transition fittedParams -> fittedMomentum
  ActsMatrix<5, 6> trkJac(ActsMatrix<5, 6>::Zero());

  // First row
  trkJac(0, 0) = -std::sin(newTrkParams[2]);
  trkJac(0, 1) = std::cos(newTrkParams[2]);

  double tanTheta = std::tan(newTrkParams[3]);

  // Second row
  trkJac(1, 0) = -trkJac(0, 1) / tanTheta;
  trkJac(1, 1) = trkJac(0, 0) / tanTheta;

  trkJac.block<4, 4>(1, 2) = ActsMatrix<4, 4>::Identity();

  // Covariance matrix of the free track parameters, i.e., x, y, z, phi, theta,
  // q/p, and t
  BoundMatrix boundTrkCov(BoundMatrix::Identity());
  boundTrkCov.block<5, 5>(0, 0) =
      (trkJac * (freeTrkCov.block<6, 6>(0, 0) * trkJac.transpose()));

  boundTrkCov.block<5, 1>(0, 5) = freeTrkCov.block<5, 1>(0, 6);
  boundTrkCov.block<1, 5>(5, 0) = freeTrkCov.block<1, 5>(6, 0);
  boundTrkCov(5, 5) = freeTrkCov(6, 6);

  return boundTrkCov;
}
