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
  // TODO we could use `linTrack.weightAtPCA` but only if we would use time
  const ActsSquareMatrix<5> trkParamWeight =
      linTrack.covarianceAtPCA.block<5, 5>(0, 0).inverse();

  // Calculate S matrix
  ActsSquareMatrix<3> sMat =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  const ActsVector<5> residual = linTrack.constantTerm.head<5>();

  // Refit track momentum
  Vector3 newTrkMomentum = sMat * momJac.transpose() * trkParamWeight *
                           (trkParams - residual - posJac * vtxPos);

  // Get phi and theta and correct for possible periodicity changes
  const auto correctedPhiTheta =
      Acts::detail::normalizePhiTheta(newTrkMomentum(0), newTrkMomentum(1));
  newTrkMomentum(0) = correctedPhiTheta.first;
  newTrkMomentum(1) = correctedPhiTheta.second;

  // Vertex covariance and weight matrices
  const SquareMatrix3 vtxCov = vtx.fullCovariance().template block<3, 3>(0, 0);
  const SquareMatrix3 vtxWeight = vtxCov.inverse();

  // New track covariance matrix
  const SquareMatrix3 newTrkCov =
      -vtxCov * posJac.transpose() * trkParamWeight * momJac * sMat;

  KalmanVertexUpdater::MatrixCache matrixCache;

  // Now determine the smoothed chi2 of the track in the following
  KalmanVertexUpdater::updatePosition<input_track_t>(
      vtx, linTrack, track.weight, -1, matrixCache);

  // Corresponding weight matrix
  const SquareMatrix3& reducedVtxWeight = matrixCache.newVertexWeight;

  // Difference in positions
  Vector3 posDiff = vtx.position() - matrixCache.newVertexPos;

  // Get smoothed params
  ActsVector<5> smParams =
      trkParams - (residual + posJac * vtx.fullPosition().template head<3>() +
                   momJac * newTrkMomentum);

  // New chi2 to be set later
  double chi2 = posDiff.dot(reducedVtxWeight * posDiff) +
                smParams.dot(trkParamWeight * smParams);

  // Fitted momentum and its covariance matrix
  ActsSquareMatrix<3> momCov =
      sMat +
      (newTrkCov).transpose() * (vtxWeight.block<3, 3>(0, 0) * newTrkCov);
  FittedMomentum fittedMom(newTrkMomentum, momCov);

  // Set new properties
  track.fittedMomentum = fittedMom;
  track.chi2 = chi2;
  track.ndf = 2 * track.weight;

  return;
}
