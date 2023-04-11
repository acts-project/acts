// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

template <typename input_track_t>
void Acts::KalmanVertexTrackUpdater::update(TrackAtVertex<input_track_t>& track,
                                            const Vertex<input_track_t>& vtx) {
  const Vector3 vtxPos = vtx.fullPosition().template head<3>();

  // Get the linearized track
  const LinearizedTrack& linTrack = track.linearizedState;

  // Check if linearized state exists
  if (linTrack.covarianceAtPCA.determinant() == 0.) {
    // Track has no linearized state, returning w/o update
    return;
  }

  // Retrieve linTrack information
  const ActsMatrix<5, 3> posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  const ActsMatrix<5, 3> momJac = linTrack.momentumJacobian.block<5, 3>(0, 0);
  const ActsVector<5> trkParams = linTrack.parametersAtPCA.head<5>();
  const ActsSymMatrix<5> trkParamWeight =
      linTrack.weightAtPCA.block<5, 5>(0, 0);

  // Calculate S matrix
  ActsSymMatrix<3> sMat =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  const ActsVector<5> residual = linTrack.constantTerm.head<5>();

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
  const SymMatrix3 vtxCov = vtx.fullCovariance().template block<3, 3>(0, 0);
  const SymMatrix3 vtxWeight = vtxCov.inverse();

  // New track covariance matrix
  const SymMatrix3 newTrkCov =
      -vtxCov * posJac.transpose() * trkParamWeight * momJac * sMat;

  KalmanVertexUpdater::MatrixCache matrixCache;

  // Now determine the smoothed chi2 of the track in the following
  KalmanVertexUpdater::updatePosition<input_track_t>(
      vtx, linTrack, track.trackWeight, -1, matrixCache);

  // Corresponding weight matrix
  const SymMatrix3& reducedVtxWeight = matrixCache.newVertexWeight;

  // Difference in positions
  Vector3 posDiff = vtx.position() - matrixCache.newVertexPos;

  // Get smoothed params
  ActsVector<5> smParams =
      trkParams - (residual + posJac * vtx.fullPosition().template head<3>() +
                   momJac * newTrkMomentum);

  // New chi2 to be set later
  double chi2 = posDiff.dot(reducedVtxWeight * posDiff) +
                smParams.dot(trkParamWeight * smParams);

  // Not yet 4d ready. This can be removed together will all head<> statements,
  // once time is consistently introduced to vertexing
  ActsMatrix<4, 3> newFullTrkCov(ActsMatrix<4, 3>::Zero());
  newFullTrkCov.block<3, 3>(0, 0) = newTrkCov;

  SymMatrix4 vtxFullWeight(SymMatrix4::Zero());
  vtxFullWeight.block<3, 3>(0, 0) = vtxWeight;

  SymMatrix4 vtxFullCov(SymMatrix4::Zero());
  vtxFullCov.block<3, 3>(0, 0) = vtxCov;

  Acts::BoundMatrix fullPerTrackCov = detail::createFullTrackCovariance(
      sMat, newFullTrkCov, vtxFullWeight, vtxFullCov, newTrkParams);

  // Create new refitted parameters
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtx.position());

  BoundTrackParameters refittedPerigee = BoundTrackParameters(
      perigeeSurface, newTrkParams, std::move(fullPerTrackCov));

  // Set new properties
  track.fittedParams = refittedPerigee;
  track.chi2Track = chi2;
  track.ndf = 2 * track.trackWeight;

  return;
}

inline Acts::BoundMatrix
Acts::KalmanVertexTrackUpdater::detail::createFullTrackCovariance(
    const SymMatrix3& sMat, const ActsMatrix<4, 3>& newTrkCov,
    const SymMatrix4& vtxWeight, const SymMatrix4& vtxCov,
    const BoundVector& newTrkParams) {
  // Now new momentum covariance
  ActsSymMatrix<3> momCov =
      sMat + (newTrkCov.block<3, 3>(0, 0)).transpose() *
                 (vtxWeight.block<3, 3>(0, 0) * newTrkCov.block<3, 3>(0, 0));

  // Full (x,y,z,phi, theta, q/p) covariance matrix
  // To be made 7d again after switching to (x,y,z,phi, theta, q/p, t)
  ActsSymMatrix<6> fullTrkCov(ActsSymMatrix<6>::Zero());

  fullTrkCov.block<3, 3>(0, 0) = vtxCov.block<3, 3>(0, 0);
  fullTrkCov.block<3, 3>(0, 3) = newTrkCov.block<3, 3>(0, 0);
  fullTrkCov.block<3, 3>(3, 0) = (newTrkCov.block<3, 3>(0, 0)).transpose();
  fullTrkCov.block<3, 3>(3, 3) = momCov;

  // Combined track jacobian
  ActsMatrix<5, 6> trkJac(ActsMatrix<5, 6>::Zero());

  // First row
  trkJac(0, 0) = -std::sin(newTrkParams[2]);
  trkJac(0, 1) = std::cos(newTrkParams[2]);

  double tanTheta = std::tan(newTrkParams[3]);

  // Second row
  trkJac(1, 0) = -trkJac(0, 1) / tanTheta;
  trkJac(1, 1) = trkJac(0, 0) / tanTheta;

  trkJac.block<4, 4>(1, 2) = ActsMatrix<4, 4>::Identity();

  // Full perigee track covariance
  BoundMatrix fullPerTrackCov(BoundMatrix::Identity());
  fullPerTrackCov.block<5, 5>(0, 0) =
      (trkJac * (fullTrkCov * trkJac.transpose()));

  return fullPerTrackCov;
}
