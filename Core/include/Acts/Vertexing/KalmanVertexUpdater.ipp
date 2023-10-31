// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>

template <typename input_track_t>
void Acts::KalmanVertexUpdater::updateVertexWithTrack(
    Vertex<input_track_t>& vtx, TrackAtVertex<input_track_t>& trk) {
  detail::update<input_track_t>(vtx, trk, 1);
}

template <typename input_track_t>
void Acts::KalmanVertexUpdater::detail::update(
    Vertex<input_track_t>& vtx, TrackAtVertex<input_track_t>& trk, int sign) {
  const double& trackWeight = trk.trackWeight;

  // Set up cache where entire content is set to 0
  Cache cache;

  // Calculate update and save result in cache
  calculateUpdate(vtx, trk.linearizedState, trackWeight, sign, cache);

  // Get fit quality parameters wrt to old vertex
  std::pair fitQuality = vtx.fitQuality();
  double chi2 = fitQuality.first;
  double ndf = fitQuality.second;

  // Chi2 of the track parameters
  double trkChi2 =
      detail::trackParametersChi2<input_track_t>(trk.linearizedState, cache);

  // Update of the chi2 of the vertex position
  double vtxPosChi2Update =
      detail::vertexPositionChi2Update<input_track_t>(vtx, cache);

  // Calculate new chi2
  chi2 += sign * (vtxPosChi2Update + trackWeight * trkChi2);

  // Calculate ndf
  ndf += sign * trackWeight * 2.;

  // Updating the vertex
  vtx.setFullPosition(cache.newVertexPos);
  vtx.setFullCovariance(cache.newVertexCov);
  vtx.setFitQuality(chi2, ndf);

  if (sign == 1) {
    // Update track
    trk.chi2Track = trkChi2;
    trk.ndf = 2 * trackWeight;
  }
  // Remove trk from current vertex by setting its weight to 0
  else if (sign == -1) {
    trk.trackWeight = 0.;
  } else {
    throw std::invalid_argument(
        "Sign for adding/removing track must be +1 (add) or -1 (remove).");
  }
}

template <typename input_track_t>
void Acts::KalmanVertexUpdater::updateTrack(TrackAtVertex<input_track_t>& track,
                                            const Vertex<input_track_t>& vtx) {
  // Check if linearized state exists
  if (!track.isLinearized) {
    throw std::invalid_argument("TrackAtVertex object must be linearized.");
  }

  // Get the linearized track
  const LinearizedTrack& linTrack = track.linearizedState;

  // Retrieve variables from the track linearization. The comments indicate the
  // corresponding symbol used in the reference.
  // A_k
  const ActsMatrix<eBoundSize, 4>& posJac = linTrack.positionJacobian;
  // B_k
  const ActsMatrix<eBoundSize, 3>& momJac = linTrack.momentumJacobian;
  // p_k
  const BoundVector& trkParams = linTrack.parametersAtPCA;
  // G_k
  const BoundSquareMatrix& trkParamWeight = linTrack.weightAtPCA;

  // Set up cache where entire content is set to 0
  Cache cache;

  // Calculate update when removing track and save result in cache. Note that
  // the track is not really removed, this is just a way of computing a
  // symmetric chi2 (see reference).
  calculateUpdate<input_track_t>(vtx, linTrack, track.trackWeight, -1, cache);

  // 4D vertex position after the fit ...
  const Vector4& vtxPos = vtx.fullPosition();
  // ... and corresponding covariance
  const SquareMatrix4& vtxCov = vtx.fullCovariance();

  // Refit track parameters
  BoundVector newTrkParams(BoundVector::Zero());
  newTrkParams(BoundIndices::eBoundPhi) = cache.newTrackMomentum(0);
  newTrkParams(BoundIndices::eBoundTheta) = cache.newTrackMomentum(1);
  newTrkParams(BoundIndices::eBoundQOverP) = cache.newTrackMomentum(2);

  // Cross covariance matrix between the vertex position and the refitted track
  // momentum
  const ActsMatrix<4, 3> crossCovVP =
      -vtxCov * posJac.transpose() * trkParamWeight * momJac * cache.wMat;

  // Difference in vertex position before and after the track removal
  Vector4 posDiff = vtxPos - cache.newVertexPos;

  // r_k
  BoundVector paramDiff = trkParams - cache.linearizedTrackParameters;

  // Updated chi2 of the track (includes contribution from the vertex, see Ref.
  // (1))
  double chi2 = posDiff.dot(cache.newVertexWeight * posDiff) +
                paramDiff.dot(trkParamWeight * paramDiff);

  Acts::BoundMatrix trkCov = detail::calculateTrackCovariance(
      cache.wMat, crossCovVP, vtxCov, newTrkParams);

  // Create new refitted parameters
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtx.position());

  BoundTrackParameters refittedParams =
      BoundTrackParameters(perigeeSurface, newTrkParams, std::move(trkCov),
                           track.fittedParams.particleHypothesis());

  // Set new properties
  track.fittedParams = refittedParams;
  track.chi2Track = chi2;
  track.ndf = 2 * track.trackWeight;

  return;
}

template <typename input_track_t>
void Acts::KalmanVertexUpdater::calculateUpdate(
    const Acts::Vertex<input_track_t>& vtx,
    const Acts::LinearizedTrack& linTrack, const double& trackWeight, int sign,
    Cache& cache) {
  // Retrieve variables from the track linearization. The comments indicate the
  // corresponding symbol used in Ref. (1).
  // A_k
  const ActsMatrix<eBoundSize, 4>& posJac = linTrack.positionJacobian;
  // B_k
  const ActsMatrix<eBoundSize, 3>& momJac = linTrack.momentumJacobian;
  // p_k
  const BoundVector& trkParams = linTrack.parametersAtPCA;
  // c_k
  const BoundVector& constTerm = linTrack.constantTerm;
  // G_k
  // Note that, when removing a track, G_k -> - G_k, see Ref. (1).
  // Further note that, as we use the weighted formalism, the track weight
  // matrix (i.e., the inverse track covariance matrix) should be multiplied
  // with the track weight from the AMVF formalism. Here, we choose to
  // consider these two multiplicative factors directly in the updates of
  // newVertexWeight and newVertexPos.
  const BoundSquareMatrix& trkParamWeight = linTrack.weightAtPCA;

  // Retrieve current position of the vertex and its current weight matrix
  const Vector4& oldVtxPos = vtx.fullPosition();
  // C_{k-1}^-1
  cache.oldVertexWeight = (vtx.fullCovariance()).inverse();

  // W_k
  cache.wMat = (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // G_k^B = G_k - G_k*B_k*W_k*B_k^(T)*G_k
  BoundSquareMatrix gBMat =
      trkParamWeight - trkParamWeight * momJac * cache.wMat *
                           momJac.transpose() * trkParamWeight;

  // C_k^-1
  cache.newVertexWeight = cache.oldVertexWeight + sign * trackWeight *
                                                      posJac.transpose() *
                                                      gBMat * posJac;
  // C_k
  cache.newVertexCov = cache.newVertexWeight.inverse();

  // \tilde{x_k}
  cache.newVertexPos =
      cache.newVertexCov * (cache.oldVertexWeight * oldVtxPos +
                            sign * trackWeight * posJac.transpose() * gBMat *
                                (trkParams - constTerm));

  // The following computations are independent of whether we add/remove the
  // track. While it does not make much sense to compute an updated track
  // momentum in the case where the track is removed, the computation is useful
  // when we do the smoothing (i.e., when we update the track with the final
  // estimation of the vertex position). Note that, for the moment, this is the
  // only case where we remove tracks.
  // A_k * \tilde{x_k}
  const BoundVector posJacVtxPos = posJac * cache.newVertexPos;

  // \tilde{q_k}
  Vector3 newTrkMom = cache.wMat * momJac.transpose() * trkParamWeight *
                      (trkParams - constTerm - posJacVtxPos);

  // Correct phi and theta for possible periodicity changes
  const auto correctedPhiTheta =
      Acts::detail::normalizePhiTheta(newTrkMom(0), newTrkMom(1));
  newTrkMom(0) = correctedPhiTheta.first;   // phi
  newTrkMom(1) = correctedPhiTheta.second;  // theta

  cache.newTrackMomentum = newTrkMom;

  // Updated linearizedq track parameters \tilde{p_k}
  cache.linearizedTrackParameters =
      constTerm + posJacVtxPos + momJac * cache.newTrackMomentum;
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::vertexPositionChi2Update(
    const Vertex<input_track_t>& oldVtx, const Cache& cache) {
  Vector4 posDiff = cache.newVertexPos - oldVtx.fullPosition();

  // Calculate and return corresponding chi2
  return posDiff.transpose() * (cache.oldVertexWeight * posDiff);
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::trackParametersChi2(
    const LinearizedTrack& linTrack, const Cache& cache) {
  // Track properties
  const BoundVector& trkParams = linTrack.parametersAtPCA;
  const BoundSquareMatrix& trkParamWeight = linTrack.weightAtPCA;

  // Parameter difference
  BoundVector paramDiff = trkParams - cache.linearizedTrackParameters;

  // Return chi2
  return paramDiff.transpose() * (trkParamWeight * paramDiff);
}

inline Acts::BoundMatrix
Acts::KalmanVertexUpdater::detail::calculateTrackCovariance(
    const SquareMatrix3& wMat, const ActsMatrix<4, 3>& crossCovVP,
    const SquareMatrix4& vtxCov, const BoundVector& newTrkParams) {
  // New momentum covariance
  SquareMatrix3 momCov =
      wMat + crossCovVP.transpose() * vtxCov.inverse() * crossCovVP;

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
