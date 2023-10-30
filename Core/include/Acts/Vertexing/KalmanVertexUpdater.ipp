// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/TrackParametrization.hpp"
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

  MatrixCache matrixCache;

  updatePosition(vtx, trk.linearizedState, trackWeight, sign, matrixCache);

  // Get fit quality parameters wrt to old vertex
  std::pair fitQuality = vtx.fitQuality();
  double chi2 = fitQuality.first;
  double ndf = fitQuality.second;

  // Chi2 wrt to track parameters
  double trkChi2 = detail::trackParametersChi2<input_track_t>(
      trk.linearizedState, matrixCache);

  // Calculate new chi2
  chi2 += sign * (detail::vertexPositionChi2<input_track_t>(vtx, matrixCache) +
                  trackWeight * trkChi2);

  // Calculate ndf
  ndf += sign * trackWeight * 2.;

  // Updating the vertex
  vtx.setFullPosition(matrixCache.newVertexPos);
  vtx.setFullCovariance(matrixCache.newVertexCov);
  vtx.setFitQuality(chi2, ndf);

  if (sign == 1) {
    // Update track
    trk.chi2Track = trkChi2;
    trk.ndf = 2 * trackWeight;
  }
  // Remove trk from current vertex
  else if (sign == -1) {
    trk.trackWeight = 0;
  } else {
    throw std::invalid_argument(
        "Sign for adding/removing track must be +1 (add) or -1 (remove).");
  }
}

template <typename input_track_t>
void Acts::KalmanVertexUpdater::updatePosition(
    const Acts::Vertex<input_track_t>& vtx,
    const Acts::LinearizedTrack& linTrack, const double& trackWeight, int sign,
    MatrixCache& matrixCache) {
  // Retrieve variables from the track linearization. The comments indicate the
  // corresponding symbol used in the reference.
  // A_k
  const ActsMatrix<eBoundSize, 4>& posJac = linTrack.positionJacobian;
  // B_k
  const ActsMatrix<eBoundSize, 3>& momJac = linTrack.momentumJacobian;
  // p_k
  const BoundVector& trkParams = linTrack.parametersAtPCA;
  // c_k
  const BoundVector& constTerm = linTrack.constantTerm;
  // G_k
  const BoundSquareMatrix& trkParamWeight = linTrack.weightAtPCA;

  // Retrieve current position of the vertex and its current weight matrix
  const Vector4& oldVtxPos = vtx.fullPosition();
  matrixCache.oldVertexWeight = (vtx.fullCovariance()).inverse();

  // W_k
  matrixCache.wMat = (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // G_k^B = G_k - G_k*B_k*W_k*B_k^(T)*G_k
  BoundSquareMatrix gMat =
      trkParamWeight - trkParamWeight * momJac * matrixCache.wMat *
                           momJac.transpose() * trkParamWeight;

  // C_k^-1
  // TODO should the trackWeight not be multiplied into trkParamWeight?
  matrixCache.newVertexWeight =
      matrixCache.oldVertexWeight +
      sign * trackWeight * posJac.transpose() * gMat * posJac;
  // C_k
  matrixCache.newVertexCov = matrixCache.newVertexWeight.inverse();

  // New vertex position
  matrixCache.newVertexPos =
      matrixCache.newVertexCov * (matrixCache.oldVertexWeight * oldVtxPos +
                                  sign * trackWeight * posJac.transpose() *
                                      gMat * (trkParams - constTerm));
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::vertexPositionChi2(
    const Vertex<input_track_t>& oldVtx, const MatrixCache& matrixCache) {
  Vector4 posDiff = matrixCache.newVertexPos - oldVtx.fullPosition();

  // Calculate and return corresponding chi2
  return posDiff.transpose() * (matrixCache.oldVertexWeight * posDiff);
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::trackParametersChi2(
    const LinearizedTrack& linTrack, const MatrixCache& matrixCache) {
  // Track properties
  const ActsMatrix<eBoundSize, 4>& posJac = linTrack.positionJacobian;
  const ActsMatrix<eBoundSize, 3>& momJac = linTrack.momentumJacobian;
  const BoundVector& trkParams = linTrack.parametersAtPCA;
  const BoundVector& constTerm = linTrack.constantTerm;
  const BoundSquareMatrix& trkParamWeight = linTrack.weightAtPCA;

  // A_k * \tilde{x_k}
  const BoundVector posJacVtxPos = posJac * matrixCache.newVertexPos;

  // \tilde{q_k}
  Vector3 newTrackMomentum = matrixCache.wMat * momJac.transpose() *
                             trkParamWeight *
                             (trkParams - constTerm - posJacVtxPos);

  // Updated linearized track parameters \tilde{p_k}
  BoundVector linearizedTrkParams =
      constTerm + posJacVtxPos + momJac * newTrackMomentum;

  // Parameter difference
  BoundVector paramDiff = trkParams - linearizedTrkParams;

  // Return chi2
  return paramDiff.transpose() * (trkParamWeight * paramDiff);
}
