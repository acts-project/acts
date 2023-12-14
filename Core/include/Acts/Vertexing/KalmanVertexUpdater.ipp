// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/VertexingError.hpp"

#include <algorithm>

template <typename input_track_t, unsigned int nDimVertex>
void Acts::KalmanVertexUpdater::updateVertexWithTrack(
    Vertex<input_track_t>& vtx, TrackAtVertex<input_track_t>& trk) {
  std::pair<double, double> fitQuality = vtx.fitQuality();
  detail::update<nDimVertex>(vtx.fullPosition(), vtx.fullCovariance(), fitQuality, trk, 1);
  vtx.setFitQuality(fitQuality);
}

template <unsigned int nDimVertex>
void Acts::KalmanVertexUpdater::detail::update(
    Vector4& vtxPos, SquareMatrix4& vtxCov,
    std::pair<double, double>& fitQuality,
    TrackAtVertexRef trk, int sign) {
  if constexpr (nDimVertex != 3 && nDimVertex != 4) {
    throw std::invalid_argument(
        "The vertex dimension must either be 3 (when fitting the spatial "
        "coordinates) or 4 (when fitting the spatial coordinates + time).");
  }

  double trackWeight = trk.trackWeight;

  // Set up cache where entire content is set to 0
  Cache<nDimVertex> cache;

  // Calculate update and save result in cache
  calculateUpdate(vtxPos, vtxCov, trk.linearizedState, trackWeight, sign, cache);

  // Get fit quality parameters wrt to old vertex
  double chi2;double ndf;
  std::tie(chi2, ndf) = fitQuality;

  // Chi2 of the track parameters
  double trkChi2 =
      detail::trackParametersChi2(trk.linearizedState, cache);

  // Update of the chi2 of the vertex position
  double vtxPosChi2Update =
      detail::vertexPositionChi2Update(vtxPos, cache);

  // Calculate new chi2
  chi2 += sign * (vtxPosChi2Update + trackWeight * trkChi2);

  // Calculate ndf
  ndf += sign * trackWeight * 2.;

  // Updating the vertex
  if constexpr (nDimVertex == 3) {
    vtxPos.head<3>() = cache.newVertexPos.template head<3>();
    vtxCov.template topLeftCorner<3, 3>() = cache.newVertexCov.template topLeftCorner<3, 3>();
  } else if constexpr (nDimVertex == 4) {
    vtxPos = cache.newVertexPos;
    vtxCov = cache.newVertexCov;
  }

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


template <unsigned int nDimVertex>
double Acts::KalmanVertexUpdater::detail::vertexPositionChi2Update(
    const Vector4& oldVtxPos, const Cache<nDimVertex>& cache) {
  ActsVector<nDimVertex> posDiff =
      cache.newVertexPos - oldVtxPos.template head<nDimVertex>();

  // Calculate and return corresponding chi2
  return posDiff.transpose() * (cache.oldVertexWeight * posDiff);
}

template <unsigned int nDimVertex>
double Acts::KalmanVertexUpdater::detail::trackParametersChi2(
    const LinearizedTrack& linTrack, const Cache<nDimVertex>& cache) {
  constexpr unsigned int nBoundParams = nDimVertex + 2;
  using ParameterVector = ActsVector<nBoundParams>;
  using ParameterMatrix = ActsSquareMatrix<nBoundParams>;
  // A_k
  const ActsMatrix<nBoundParams, nDimVertex> posJac =
      linTrack.positionJacobian.block<nBoundParams, nDimVertex>(0, 0);
  // B_k
  const ActsMatrix<nBoundParams, 3> momJac =
      linTrack.momentumJacobian.block<nBoundParams, 3>(0, 0);
  // p_k
  const ParameterVector trkParams =
      linTrack.parametersAtPCA.head<nBoundParams>();
  // c_k
  const ParameterVector constTerm = linTrack.constantTerm.head<nBoundParams>();
  // TODO we could use `linTrack.weightAtPCA` but only if we would always fit
  // time.
  // G_k
  const ParameterMatrix trkParamWeight =
      linTrack.covarianceAtPCA.block<nBoundParams, nBoundParams>(0, 0)
          .inverse();

  // A_k * \tilde{x_k}
  const ParameterVector posJacVtxPos = posJac * cache.newVertexPos;

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
  ParameterVector linearizedTrackParameters =
      constTerm + posJacVtxPos + momJac * newTrkMom;

  // r_k
  ParameterVector paramDiff = trkParams - linearizedTrackParameters;

  // Return chi2
  return paramDiff.transpose() * (trkParamWeight * paramDiff);
}
