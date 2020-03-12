// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <algorithm>
#include "Acts/Vertexing/VertexingError.hpp"

template <typename input_track_t>
void Acts::KalmanVertexUpdater::updateVertexWithTrack(
    Vertex<input_track_t>& vtx, TrackAtVertex<input_track_t>& trk) {
  detail::update<input_track_t>(vtx, trk, 1);
}

template <typename input_track_t>
void Acts::KalmanVertexUpdater::detail::update(
    Vertex<input_track_t>& vtx, TrackAtVertex<input_track_t>& trk, int sign) {
  double trackWeight = trk.trackWeight;

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

template <typename input_track_t>
void Acts::KalmanVertexUpdater::updatePosition(
    const Acts::Vertex<input_track_t>& vtx,
    const Acts::LinearizedTrack& linTrack, double trackWeight, int sign,
    MatrixCache& matrixCache) {
  // Retrieve linTrack information
  // TODO: To make 4-D compatible, remove block<> and head<> statements
  const auto posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  const auto momJac =
      linTrack.momentumJacobian.block<5, 3>(0, 0);  // B_k in comments below
  const auto trkParams = linTrack.parametersAtPCA.head<5>();
  const auto constTerm = linTrack.constantTerm.head<5>();
  const auto trkParamWeight = linTrack.weightAtPCA.block<5, 5>(0, 0);

  // Vertex to be updated
  const auto& oldVtxPos = vtx.position();
  matrixCache.oldVertexWeight = (vtx.covariance()).inverse();

  // W_k matrix
  matrixCache.momWeightInv =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // G_b = G_k - G_k*B_k*W_k*B_k^(T)*G_k^T
  auto gBmat = trkParamWeight -
               trkParamWeight *
                   (momJac * (matrixCache.momWeightInv * momJac.transpose())) *
                   trkParamWeight.transpose();

  // New vertex cov matrix
  matrixCache.newVertexWeight =
      matrixCache.oldVertexWeight +
      trackWeight * sign * posJac.transpose() * (gBmat * posJac);
  matrixCache.newVertexCov = matrixCache.newVertexWeight.inverse();

  // New vertex position
  matrixCache.newVertexPos =
      matrixCache.newVertexCov * (matrixCache.oldVertexWeight * oldVtxPos +
                                  trackWeight * sign * posJac.transpose() *
                                      gBmat * (trkParams - constTerm));
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::vertexPositionChi2(
    const Vertex<input_track_t>& oldVtx, const MatrixCache& matrixCache) {
  Vector3D posDiff = matrixCache.newVertexPos - oldVtx.position();

  // Calculate and return corresponding chi2
  return posDiff.transpose() * (matrixCache.oldVertexWeight * posDiff);
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::trackParametersChi2(
    const LinearizedTrack& linTrack, const MatrixCache& matrixCache) {
  // Track properties
  const auto posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  const auto momJac = linTrack.momentumJacobian.block<5, 3>(0, 0);
  const auto trkParams = linTrack.parametersAtPCA.head<5>();
  const auto constTerm = linTrack.constantTerm.head<5>();
  const auto trkParamWeight = linTrack.weightAtPCA.block<5, 5>(0, 0);

  const auto jacVtx = posJac * matrixCache.newVertexPos;

  // Refitted track momentum
  Vector3D newTrackMomentum = matrixCache.momWeightInv * momJac.transpose() *
                              trkParamWeight * (trkParams - constTerm - jacVtx);

  // Refitted track parameters
  auto newTrkParams = constTerm + jacVtx + momJac * newTrackMomentum;

  // Parameter difference
  auto paramDiff = trkParams - newTrkParams;

  // Return chi2
  return paramDiff.transpose() * (trkParamWeight * paramDiff);
}
