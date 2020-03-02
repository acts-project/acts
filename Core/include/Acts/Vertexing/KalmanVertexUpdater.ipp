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
Acts::Result<void> Acts::KalmanVertexUpdater::updateVertexWithTrack(
    Vertex<input_track_t>* vtx, TrackAtVertex<input_track_t>& trk) {
  if (vtx == nullptr) {
    return VertexingError::EmptyInput;
  }

  auto res = detail::update<input_track_t>(vtx, trk, 1);

  if (!res.ok()) {
    return res.error();
  }

  return {};
}

template <typename input_track_t>
Acts::Result<Acts::Vertex<input_track_t>>
Acts::KalmanVertexUpdater::updatePosition(
    const Acts::Vertex<input_track_t>* vtx,
    const Acts::LinearizedTrack& linTrack, double trackWeight, int sign) {
  if (vtx == nullptr) {
    return VertexingError::EmptyInput;
  }

  // Retrieve linTrack information
  // To make 4-D compatible, remove block<> and head<> statements
  const auto& posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  const auto& momJac =
      linTrack.momentumJacobian.block<5, 3>(0, 0);  // B_k in comments below
  const auto& trkParams = linTrack.parametersAtPCA.head<5>();
  const auto& constTerm = linTrack.constantTerm.head<5>();
  const auto& trkParamWeight = (linTrack.covarianceAtPCA.block<5, 5>(0, 0))
                                   .inverse();  // G_k in comments below

  // Vertex to be updated
  const auto& oldVtxPos = vtx->position();
  const auto& oldVtxWeight = (vtx->covariance()).inverse().eval();

  // W_k matrix
  ActsSymMatrixD<3> wMat =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // G_b = G_k - G_k*B_k*W_k*B_k^(T)*G_k^T
  auto gBmat = trkParamWeight - trkParamWeight *
                                    (momJac * (wMat * momJac.transpose())) *
                                    trkParamWeight.transpose();
  // New vertex cov matrix
  auto newVtxCov = (oldVtxWeight +
                    trackWeight * sign * posJac.transpose() * (gBmat * posJac))
                       .inverse();
  // New vertex position
  auto newVtxPos = newVtxCov * (oldVtxWeight * oldVtxPos +
                                trackWeight * sign * posJac.transpose() *
                                    gBmat * (trkParams - constTerm));
  // Create return vertex with new position
  // and covariance, but w/o tracks
  Vertex<input_track_t> returnVertex;

  // Set position
  returnVertex.setPosition(newVtxPos);
  // Set cov
  returnVertex.setCovariance(newVtxCov);
  // Set fit quality
  returnVertex.setFitQuality(vtx->fitQuality().first, vtx->fitQuality().second);

  return returnVertex;
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::vertexPositionChi2(
    const Vertex<input_track_t>* oldVtx, const Vertex<input_track_t>* newVtx) {
  auto oldWeight =
      (oldVtx->fullCovariance().template block<3, 3>(0, 0)).inverse();
  auto posDiff =
      (newVtx->fullPosition() - oldVtx->fullPosition()).template head<3>();

  // Calculate and return corresponding chi2
  return posDiff.transpose() * (oldWeight * posDiff);
}

template <typename input_track_t>
double Acts::KalmanVertexUpdater::detail::trackParametersChi2(
    const Vertex<input_track_t>& vtx, const LinearizedTrack& linTrack) {
  const auto& vtxPos = vtx.fullPosition().template head<3>();

  // Track properties
  const auto& posJac = linTrack.positionJacobian.block<5, 3>(0, 0);
  const auto& momJac = linTrack.momentumJacobian.block<5, 3>(0, 0);
  const auto& trkParams = linTrack.parametersAtPCA.head<5>();
  const auto& constTerm = linTrack.constantTerm.head<5>();
  const auto& trkParamWeight =
      linTrack.covarianceAtPCA.block<5, 5>(0, 0).inverse();

  // Calculate temp matrix S
  ActsSymMatrixD<3> matS =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // Refitted track momentum
  Vector3D newTrackMomentum = matS * momJac.transpose() * trkParamWeight *
                              (trkParams - constTerm - posJac * vtxPos);

  // Refitted track parameters
  auto newTrkParams = constTerm + posJac * vtxPos + momJac * newTrackMomentum;

  // Parameter difference
  auto paramDiff = trkParams - newTrkParams;

  // Return chi2
  return paramDiff.transpose() * (trkParamWeight * paramDiff);
}

template <typename input_track_t>
Acts::Result<void> Acts::KalmanVertexUpdater::detail::update(
    Vertex<input_track_t>* vtx, TrackAtVertex<input_track_t>& trk, int sign) {
  double trackWeight = trk.trackWeight;

  auto res = updatePosition(vtx, trk.linearizedState, trackWeight, sign);

  if (!res.ok()) {
    return res.error();
  }

  Vertex<input_track_t> tempVtx = *res;

  // Get fit quality parameters wrt to old vertex
  std::pair fitQuality = vtx->fitQuality();
  double chi2 = fitQuality.first;
  double ndf = fitQuality.second;

  // Chi2 wrt to track parameters
  double trkChi2 =
      detail::trackParametersChi2<input_track_t>(tempVtx, trk.linearizedState);

  // Calculate new chi2
  chi2 += sign * (detail::vertexPositionChi2<input_track_t>(vtx, &tempVtx) +
                  trackWeight * trkChi2);

  // Calculate ndf
  ndf += sign * trackWeight * 2.;

  // Updating the vertex
  vtx->setFullPosition(tempVtx.fullPosition());
  vtx->setFullCovariance(tempVtx.fullCovariance());
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

  return {};
}
