// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
  double trackWeight = trk.trackWeight;

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
  vtx.setPosition(cache.newVertexPos);
  vtx.setCovariance(cache.newVertexCov);
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
void Acts::KalmanVertexUpdater::calculateUpdate(
    const Acts::Vertex<input_track_t>& vtx,
    const Acts::LinearizedTrack& linTrack, const double trackWeight,
    const int sign, Cache& cache) {
  // Retrieve variables from the track linearization. The comments indicate the
  // corresponding symbol used in Ref. (1).
  // A_k
  const ActsMatrix<5, 3> posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  // B_k
  const ActsMatrix<5, 3> momJac = linTrack.momentumJacobian.block<5, 3>(0, 0);
  // p_k
  const ActsVector<5> trkParams = linTrack.parametersAtPCA.head<5>();
  // c_k
  const ActsVector<5> constTerm = linTrack.constantTerm.head<5>();
  // TODO we could use `linTrack.weightAtPCA` but only if we would use time
  // G_k
  // Note that, when removing a track, G_k -> - G_k, see Ref. (1).
  // Further note that, as we use the weighted formalism, the track weight
  // matrix (i.e., the inverse track covariance matrix) should be multiplied
  // with the track weight from the AMVF formalism. Here, we choose to
  // consider these two multiplicative factors directly in the updates of
  // newVertexWeight and newVertexPos.
  const ActsSquareMatrix<5> trkParamWeight =
      linTrack.covarianceAtPCA.block<5, 5>(0, 0).inverse();

  // Retrieve current position of the vertex and its current weight matrix
  const Vector3& oldVtxPos = vtx.position();
  // C_{k-1}^-1
  cache.oldVertexWeight = (vtx.covariance()).inverse();

  // W_k
  cache.wMat = (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // G_k^B = G_k - G_k*B_k*W_k*B_k^(T)*G_k
  ActsSquareMatrix<5> gBMat =
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
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::vertexPositionChi2Update(
    const Vertex<input_track_t>& oldVtx, const Cache& cache) {
  Vector3 posDiff = cache.newVertexPos - oldVtx.position();

  // Calculate and return corresponding chi2
  return posDiff.transpose() * (cache.oldVertexWeight * posDiff);
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::trackParametersChi2(
    const LinearizedTrack& linTrack, const Cache& cache) {
  // A_k
  const ActsMatrix<5, 3> posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  // B_k
  const ActsMatrix<5, 3> momJac = linTrack.momentumJacobian.block<5, 3>(0, 0);
  // p_k
  const ActsVector<5> trkParams = linTrack.parametersAtPCA.head<5>();
  // c_k
  const ActsVector<5> constTerm = linTrack.constantTerm.head<5>();
  // TODO we could use `linTrack.weightAtPCA` but only if we would use time
  // G_k
  const ActsSquareMatrix<5> trkParamWeight =
      linTrack.covarianceAtPCA.block<5, 5>(0, 0).inverse();

  // A_k * \tilde{x_k}
  const ActsVector<5> posJacVtxPos = posJac * cache.newVertexPos;

  // \tilde{q_k}
  Vector3 newTrkMom = cache.wMat * momJac.transpose() * trkParamWeight *
                      (trkParams - constTerm - posJacVtxPos);

  // Correct phi and theta for possible periodicity changes
  // Commented out because of broken ATHENA tests.
  // TODO: uncomment
  /*
  const auto correctedPhiTheta =
      Acts::detail::normalizePhiTheta(newTrkMom(0), newTrkMom(1));
  newTrkMom(0) = correctedPhiTheta.first;   // phi
  newTrkMom(1) = correctedPhiTheta.second;  // theta
  */

  // \tilde{p_k}
  ActsVector<5> linearizedTrackParameters =
      constTerm + posJacVtxPos + momJac * newTrkMom;

  // r_k
  ActsVector<5> paramDiff = trkParams - linearizedTrackParameters;

  // Return chi2
  return paramDiff.transpose() * (trkParamWeight * paramDiff);
}
