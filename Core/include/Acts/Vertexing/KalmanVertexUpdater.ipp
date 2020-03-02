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
    Vertex<input_track_t>* vtx, TrackAtVertex<input_track_t>& trk) {
  detail::update<input_track_t>(vtx, trk, 1);
}

template <typename input_track_t>
void Acts::KalmanVertexUpdater::updatePosition(
    const Acts::Vertex<input_track_t>* vtx,
    const Acts::LinearizedTrack& linTrack, double trackWeight, int sign,
    Vector3D& newVtxPos, ActsSymMatrixD<3>& newVtxCov,
    ActsSymMatrixD<3>& oldVtxWeight, ActsSymMatrixD<5>& trkParamWeight) {
  // Retrieve linTrack information
  // To make 4-D compatible, remove block<> and head<> statements
  const auto& posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  const auto& momJac =
      linTrack.momentumJacobian.block<5, 3>(0, 0);  // B_k in comments below
  const auto& trkParams = linTrack.parametersAtPCA.head<5>();
  const auto& constTerm = linTrack.constantTerm.head<5>();
  trkParamWeight = (linTrack.covarianceAtPCA.block<5, 5>(0, 0))
                       .inverse();  // G_k in comments below

  // Vertex to be updated
  const auto& oldVtxPos = vtx->position();
  oldVtxWeight = (vtx->covariance()).inverse();

  // W_k matrix
  ActsSymMatrixD<3> wMat =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // G_b = G_k - G_k*B_k*W_k*B_k^(T)*G_k^T
  auto gBmat = trkParamWeight - trkParamWeight *
                                    (momJac * (wMat * momJac.transpose())) *
                                    trkParamWeight.transpose();
  // New vertex cov matrix
  newVtxCov = (oldVtxWeight +
               trackWeight * sign * posJac.transpose() * (gBmat * posJac))
                  .inverse();

  // New vertex position
  newVtxPos = newVtxCov * (oldVtxWeight * oldVtxPos +
                           trackWeight * sign * posJac.transpose() * gBmat *
                               (trkParams - constTerm));
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::vertexPositionChi2(
    const Vertex<input_track_t>* oldVtx, const Vector3D& newVtxPos,
    const ActsSymMatrixD<3>& oldVertexWeight) {
  auto posDiff = newVtxPos - oldVtx->position();

  // Calculate and return corresponding chi2
  return posDiff.transpose() * (oldVertexWeight * posDiff);
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::trackParametersChi2(
    const Vector3D& vtxPos, const LinearizedTrack& linTrack,
    const ActsSymMatrixD<5>& trkParamWeight) {
  // Track properties
  const auto& posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  const auto& momJac = linTrack.momentumJacobian.block<5, 3>(0, 0);
  const auto& trkParams = linTrack.parametersAtPCA.head<5>();
  const auto& constTerm = linTrack.constantTerm.head<5>();

  // Calculate temp matrix S
  ActsSymMatrixD<3> matS =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  const auto jacVtx = posJac * vtxPos;

  // Refitted track momentum
  Vector3D newTrackMomentum = matS * momJac.transpose() * trkParamWeight *
                              (trkParams - constTerm - jacVtx);

  // Refitted track parameters
  auto newTrkParams = constTerm + jacVtx + momJac * newTrackMomentum;

  // Parameter difference
  auto paramDiff = trkParams - newTrkParams;

  // Return chi2
  return paramDiff.transpose() * (trkParamWeight * paramDiff);
}

template <typename input_track_t>
void Acts::KalmanVertexUpdater::detail::update(
    Vertex<input_track_t>* vtx, TrackAtVertex<input_track_t>& trk, int sign) {
  double trackWeight = trk.trackWeight;

  Vector3D newVertexPos = Vector3D::Zero();
  ActsSymMatrixD<3> newVertexCov = ActsSymMatrixD<3>::Zero();
  ActsSymMatrixD<3> oldVertexWeight = ActsSymMatrixD<3>::Zero();
  ActsSymMatrixD<5> trkParamWeight = ActsSymMatrixD<5>::Zero();

  updatePosition(vtx, trk.linearizedState, trackWeight, sign, newVertexPos,
                 newVertexCov, oldVertexWeight, trkParamWeight);

  // Get fit quality parameters wrt to old vertex
  std::pair fitQuality = vtx->fitQuality();
  double chi2 = fitQuality.first;
  double ndf = fitQuality.second;

  // Chi2 wrt to track parameters
  double trkChi2 = detail::trackParametersChi2<input_track_t>(
      newVertexPos, trk.linearizedState, trkParamWeight);

  // Calculate new chi2
  chi2 += sign * (detail::vertexPositionChi2<input_track_t>(vtx, newVertexPos,
                                                            oldVertexWeight) +
                  trackWeight * trkChi2);

  // Calculate ndf
  ndf += sign * trackWeight * 2.;

  // Updating the vertex
  vtx->setPosition(newVertexPos);
  vtx->setCovariance(newVertexCov);
  vtx->setFitQuality(chi2, ndf);

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
