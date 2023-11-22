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
  // Check if linearized state exists
  if (!track.isLinearized) {
    throw std::invalid_argument("TrackAtVertex object must be linearized.");
  }

  // Extract vertex position and covariance
  // \tilde{x_n}
  const Vector3& vtxPos = vtx.position();
  // C_n
  const SquareMatrix3& vtxCov = vtx.covariance();

  // Get the linearized track
  const LinearizedTrack& linTrack = track.linearizedState;

  // Retrieve variables from the track linearization. The comments indicate the
  // corresponding symbol used in the Ref. (1).
  // A_k
  const ActsMatrix<5, 3> posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  // B_k
  const ActsMatrix<5, 3> momJac = linTrack.momentumJacobian.block<5, 3>(0, 0);
  // p_k
  const ActsVector<5> trkParams = linTrack.parametersAtPCA.head<5>();
  // TODO we could use `linTrack.weightAtPCA` but only if we would use time
  // G_k
  const ActsSquareMatrix<5> trkParamWeight =
      linTrack.covarianceAtPCA.block<5, 5>(0, 0).inverse();
  // c_k
  const ActsVector<5> constTerm = linTrack.constantTerm.head<5>();

  // Cache object filled with zeros
  KalmanVertexUpdater::Cache cache;

  // Calculate the update of the vertex position when the track is removed. This
  // might be unintuitive, but it is needed to compute a symmetric chi2.
  KalmanVertexUpdater::calculateUpdate<input_track_t>(
      vtx, linTrack, track.trackWeight, -1, cache);

  // Refit track momentum with the final vertex position
  Vector3 newTrkMomentum = cache.wMat * momJac.transpose() * trkParamWeight *
                           (trkParams - constTerm - posJac * vtxPos);

  // Track parameters, impact parameters are set to 0 and momentum corresponds
  // to newTrkMomentum. TODO: Make transition fitterParams -> fittedMomentum.
  BoundVector newTrkParams(BoundVector::Zero());

  // Get phi and theta and correct for possible periodicity changes
  const auto correctedPhiTheta =
      Acts::detail::normalizePhiTheta(newTrkMomentum(0), newTrkMomentum(1));
  newTrkParams(BoundIndices::eBoundPhi) = correctedPhiTheta.first;     // phi
  newTrkParams(BoundIndices::eBoundTheta) = correctedPhiTheta.second;  // theta
  newTrkParams(BoundIndices::eBoundQOverP) = newTrkMomentum(2);        // qOverP

  // E_k^n
  const SquareMatrix3 crossCovVP =
      -vtxCov * posJac.transpose() * trkParamWeight * momJac * cache.wMat;

  // Difference in positions. cache.newVertexPos corresponds to \tilde{x_k^{n*}} in Ref. (1).
  Vector3 posDiff = vtxPos - cache.newVertexPos;

  // r_k^n
  ActsVector<5> paramDiff =
      trkParams - (constTerm + posJac * vtxPos + momJac * newTrkMomentum);

  // New chi2 to be set later
  double chi2 = posDiff.dot(cache.newVertexWeight * posDiff) +
                paramDiff.dot(trkParamWeight * paramDiff);

  Acts::BoundMatrix fullPerTrackCov = detail::createFullTrackCovariance(
      cache.wMat, crossCovVP, vtxCov, newTrkParams);

  // Create new refitted parameters
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtxPos);

  BoundTrackParameters refittedPerigee = BoundTrackParameters(
      perigeeSurface, newTrkParams, std::move(fullPerTrackCov),
      track.fittedParams.particleHypothesis());

  // Set new properties
  track.fittedParams = refittedPerigee;
  track.chi2Track = chi2;
  track.ndf = 2 * track.trackWeight;

  return;
}

inline Acts::BoundMatrix
Acts::KalmanVertexTrackUpdater::detail::createFullTrackCovariance(
    const SquareMatrix3& wMat, const SquareMatrix3& crossCovVP,
    const SquareMatrix3& vtxCov, const BoundVector& newTrkParams) {
  // D_k^n
  ActsSquareMatrix<3> momCov =
      wMat + crossCovVP.transpose() * vtxCov.inverse() * crossCovVP;

  // Full (x,y,z,phi, theta, q/p) covariance matrix
  // To be made 7d again after switching to (x,y,z,phi, theta, q/p, t)
  ActsSquareMatrix<6> fullTrkCov(ActsSquareMatrix<6>::Zero());

  fullTrkCov.block<3, 3>(0, 0) = vtxCov;
  fullTrkCov.block<3, 3>(0, 3) = crossCovVP;
  fullTrkCov.block<3, 3>(3, 0) = crossCovVP.transpose();
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
