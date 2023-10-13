// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>

template <typename input_track_t, unsigned int nDimVertex>
void Acts::KalmanVertexUpdater::updateVertexWithTrack(
    Vertex<input_track_t>& vtx, TrackAtVertex<input_track_t>& trk) {
  static_assert(nDimVertex == 3 or nDimVertex == 4,
                "The dimension of the vertex must either be 3 (vertexing "
                "without time) or 4 (vertexing with time).");

  detail::update<input_track_t, nDimVertex>(vtx, trk, 1);
}

template <typename input_track_t, unsigned int nDimVertex>
void Acts::KalmanVertexUpdater::detail::update(
    Vertex<input_track_t>& vtx, TrackAtVertex<input_track_t>& trk, int sign) {
  static_assert(nDimVertex == 3 or nDimVertex == 4,
                "The dimension of the vertex must either be 3 (vertexing "
                "without time) or 4 (vertexing with time).");

  double trackWeight = trk.trackWeight;

  MatrixCache<nDimVertex> matrixCache;

  updatePosition<input_track_t, nDimVertex>(vtx, trk.linearizedState,
                                            trackWeight, sign, matrixCache);

  // Get fit quality parameters wrt to old vertex
  std::pair fitQuality = vtx.fitQuality();
  double chi2 = fitQuality.first;
  double ndf = fitQuality.second;

  // Chi2 wrt to track parameters
  double trkChi2 = detail::trackParametersChi2<input_track_t, nDimVertex>(
      trk.linearizedState, matrixCache);

  // Calculate new chi2
  chi2 += sign * (detail::vertexPositionChi2<input_track_t, nDimVertex>(
                      vtx, matrixCache) +
                  trackWeight * trkChi2);

  // Calculate ndf
  ndf += sign * trackWeight * 2.;

  // Updating the vertex
  vtx.setPosition(matrixCache.newVertexPos);
  vtx.setCovariance(matrixCache.newVertexCov);
  vtx.setFitQuality(chi2, ndf);

  // Updates track at vertex if already there
  // by removing it first and adding new one.
  // Otherwise just adds track to existing list of tracks at vertex
  if (sign > 0) {
    // Update track
    trk.chi2Track = trkChi2;
    trk.ndf = 2 * trackWeight;
  }
  // Remove trk from current vertex
  if (sign < 0) {
    trk.trackWeight = 0;
  }
}

template <typename input_track_t, unsigned int nDimVertex>
void Acts::KalmanVertexUpdater::updatePosition(
    const Acts::Vertex<input_track_t>& vtx,
    const Acts::LinearizedTrack& linTrack, double trackWeight, int sign,
    MatrixCache<nDimVertex>& matrixCache) {
  static_assert(nDimVertex == 3 or nDimVertex == 4,
                "The dimension of the vertex must either be 3 (vertexing "
                "without time) or 4 (vertexing with time).");

  // Number of track parameters. We have nDimVertex - 1 impact parameters (i.e.,
  // d0, z0, and, in the case of time vertexing, t0). Additionally, we have 3
  // parameters describing the momentum.
  constexpr unsigned int nParams = nDimVertex + 2;
  // Retrieve position and momentum Jacobian. posJack corresponds to A_k and
  // momJack to B_k in the reference and the comments below.
  const ActsMatrix<nParams, nDimVertex> posJac =
      linTrack.positionJacobian.block<nParams, nDimVertex>(0, 0);
  const ActsMatrix<nParams, 3> momJac =
      linTrack.momentumJacobian.block<nParams, 3>(0, 0);
  const ActsVector<nParams> trkParams =
      linTrack.parametersAtPCA.head<nParams>();
  const ActsVector<nParams> constTerm = linTrack.constantTerm.head<nParams>();
  const ActsSquareMatrix<nParams> trkParamWeight =
      linTrack.covarianceAtPCA.block<nParams, nParams>(0, 0).inverse();

  // Vertex to be updated
  const ActsVector<nDimVertex>& oldVtxPos =
      vtx.fullPosition().template head<nDimVertex>();
  matrixCache.oldVertexWeight = (vtx.covariance()).inverse();

  // W_k matrix
  matrixCache.momWeightInv =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // G_k^B = G_k - G_k*B_k*W_k*B_k^(T)*G_k^T
  ActsSquareMatrix<nParams> gMat =
      trkParamWeight -
      trkParamWeight *
          (momJac * (matrixCache.momWeightInv * momJac.transpose())) *
          trkParamWeight.transpose();

  // New vertex weight matrix
  matrixCache.newVertexWeight =
      matrixCache.oldVertexWeight +
      trackWeight * sign * posJac.transpose() * (gMat * posJac);
  matrixCache.newVertexCov = matrixCache.newVertexWeight.inverse();

  // New vertex position
  matrixCache.newVertexPos =
      matrixCache.newVertexCov * (matrixCache.oldVertexWeight * oldVtxPos +
                                  trackWeight * sign * posJac.transpose() *
                                      gMat * (trkParams - constTerm));
}

template <typename input_track_t, unsigned int nDimVertex>
double Acts::KalmanVertexUpdater::detail::vertexPositionChi2(
    const Vertex<input_track_t>& oldVtx,
    const MatrixCache<nDimVertex>& matrixCache) {
  static_assert(nDimVertex == 3 or nDimVertex == 4,
                "The dimension of the vertex must either be 3 (vertexing "
                "without time) or 4 (vertexing with time).");
  ActsVector<nDimVertex> posDiff =
      matrixCache.newVertexPos -
      oldVtx.fullPosition().template head<nDimVertex>();

  // Calculate and return corresponding chi2
  return posDiff.transpose() * (matrixCache.oldVertexWeight * posDiff);
}

template <typename input_track_t, unsigned int nDimVertex>
double Acts::KalmanVertexUpdater::detail::trackParametersChi2(
    const LinearizedTrack& linTrack,
    const MatrixCache<nDimVertex>& matrixCache) {
  static_assert(nDimVertex == 3 or nDimVertex == 4,
                "The dimension of the vertex must either be 3 (vertexing "
                "without time) or 4 (vertexing with time).");
  // Number of track parameters. We have nDimVertex - 1 impact parameters (i.e.,
  // d0, z0, and, in the case of time vertexing, t0). Additionally, we have 3
  // parameters describing the momentum.
  constexpr unsigned int nParams = nDimVertex + 2;

  const ActsMatrix<nParams, nDimVertex> posJac =
      linTrack.positionJacobian.block<nParams, nDimVertex>(0, 0);
  const ActsMatrix<nParams, 3> momJac =
      linTrack.momentumJacobian.block<nParams, 3>(0, 0);
  const ActsVector<nParams> trkParams =
      linTrack.parametersAtPCA.head<nParams>();
  const ActsVector<nParams> constTerm = linTrack.constantTerm.head<nParams>();
  const ActsSquareMatrix<nParams> trkParamWeight =
      linTrack.covarianceAtPCA.block<nParams, nParams>(0, 0).inverse();

  const ActsVector<nParams> jacVtx = posJac * matrixCache.newVertexPos;

  // Refitted track momentum
  Vector3 newTrackMomentum = matrixCache.momWeightInv * momJac.transpose() *
                             trkParamWeight * (trkParams - constTerm - jacVtx);

  // Refitted track parameters
  ActsVector<nParams> newTrkParams =
      constTerm + jacVtx + momJac * newTrackMomentum;

  // Parameter difference
  ActsVector<nParams> paramDiff = trkParams - newTrkParams;

  // Return chi2
  return paramDiff.transpose() * (trkParamWeight * paramDiff);
}
