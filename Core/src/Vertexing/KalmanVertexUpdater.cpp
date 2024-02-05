// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/KalmanVertexUpdater.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

namespace Acts {
namespace {
template <unsigned int nDimVertex>
void calculateUpdateImpl(const Vector4& vtxPos, const SquareMatrix4& vtxCov,
                         const Acts::LinearizedTrack& linTrack,
                         const double trackWeight, const int sign,
                         KalmanVertexUpdater::Cache<nDimVertex>& cache) {
  constexpr unsigned int nBoundParams = nDimVertex + 2;
  using ParameterVector = ActsVector<nBoundParams>;
  using ParameterMatrix = ActsSquareMatrix<nBoundParams>;
  // Retrieve variables from the track linearization. The comments indicate the
  // corresponding symbol used in Ref. (1).
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
  // Note that, when removing a track, G_k -> - G_k, see Ref. (1).
  // Further note that, as we use the weighted formalism, the track weight
  // matrix (i.e., the inverse track covariance matrix) should be multiplied
  // with the track weight from the AMVF formalism. Here, we choose to
  // consider these two multiplicative factors directly in the updates of
  // newVertexWeight and newVertexPos.
  const ParameterMatrix trkParamWeight =
      linTrack.covarianceAtPCA.block<nBoundParams, nBoundParams>(0, 0)
          .inverse();

  // Retrieve current position of the vertex and its current weight matrix
  const ActsVector<nDimVertex> oldVtxPos = vtxPos.template head<nDimVertex>();
  // C_{k-1}^-1
  cache.oldVertexWeight =
      (vtxCov.template block<nDimVertex, nDimVertex>(0, 0)).inverse();

  // W_k
  cache.wMat = (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // G_k^B = G_k - G_k*B_k*W_k*B_k^(T)*G_k
  ParameterMatrix gBMat = trkParamWeight - trkParamWeight * momJac *
                                               cache.wMat * momJac.transpose() *
                                               trkParamWeight;

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

template <unsigned int nDimVertex>
double vertexPositionChi2Update(
    const Vector4& oldVtxPos,
    const KalmanVertexUpdater::Cache<nDimVertex>& cache) {
  ActsVector<nDimVertex> posDiff =
      cache.newVertexPos - oldVtxPos.template head<nDimVertex>();

  // Calculate and return corresponding chi2
  return posDiff.transpose() * (cache.oldVertexWeight * posDiff);
}

template <unsigned int nDimVertex>
double trackParametersChi2(
    const LinearizedTrack& linTrack,
    const KalmanVertexUpdater::Cache<nDimVertex>& cache) {
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

  // \tilde{p_k}
  ParameterVector linearizedTrackParameters =
      constTerm + posJacVtxPos + momJac * newTrkMom;

  // r_k
  ParameterVector paramDiff = trkParams - linearizedTrackParameters;

  // Return chi2
  return paramDiff.transpose() * (trkParamWeight * paramDiff);
}

template <unsigned int nDimVertex>
void update(Vector4& vtxPos, SquareMatrix4& vtxCov,
            std::pair<double, double>& fitQuality, TrackAtVertexRef trk,
            int sign) {
  if constexpr (nDimVertex != 3 && nDimVertex != 4) {
    throw std::invalid_argument(
        "The vertex dimension must either be 3 (when fitting the spatial "
        "coordinates) or 4 (when fitting the spatial coordinates + time).");
  }

  double trackWeight = trk.trackWeight;

  // Set up cache where entire content is set to 0
  KalmanVertexUpdater::Cache<nDimVertex> cache;

  // Calculate update and save result in cache
  calculateUpdate(vtxPos, vtxCov, trk.linearizedState, trackWeight, sign,
                  cache);

  // Get fit quality parameters wrt to old vertex
  double chi2 = 0.;
  double ndf = 0.;
  std::tie(chi2, ndf) = fitQuality;

  // Chi2 of the track parameters
  double trkChi2 = trackParametersChi2(trk.linearizedState, cache);

  // Update of the chi2 of the vertex position
  double vtxPosChi2Update = vertexPositionChi2Update(vtxPos, cache);

  // Calculate new chi2
  chi2 += sign * (vtxPosChi2Update + trackWeight * trkChi2);

  // Calculate ndf
  ndf += sign * trackWeight * 2.;

  // Updating the vertex
  if constexpr (nDimVertex == 3) {
    vtxPos.setZero();
    vtxPos.head<3>() = cache.newVertexPos.template head<3>();
    vtxCov.setZero();
    vtxCov.template topLeftCorner<3, 3>() =
        cache.newVertexCov.template topLeftCorner<3, 3>();
  } else if constexpr (nDimVertex == 4) {
    vtxPos = cache.newVertexPos;
    vtxCov = cache.newVertexCov;
  }
  fitQuality = {chi2, ndf};

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

}  // namespace

void KalmanVertexUpdater::detail::calculateUpdate3(
    const Vector4& vtxPos, const SquareMatrix4& vtxCov,
    const Acts::LinearizedTrack& linTrack, const double trackWeight,
    const int sign, Cache<3>& cache) {
  calculateUpdateImpl<3>(vtxPos, vtxCov, linTrack, trackWeight, sign, cache);
}

void KalmanVertexUpdater::detail::calculateUpdate4(
    const Vector4& vtxPos, const SquareMatrix4& vtxCov,
    const Acts::LinearizedTrack& linTrack, const double trackWeight,
    const int sign, Cache<4>& cache) {
  calculateUpdateImpl<4>(vtxPos, vtxCov, linTrack, trackWeight, sign, cache);
}

void KalmanVertexUpdater::detail::updateVertexWithTrack(
    Vector4& vtxPos, SquareMatrix4& vtxCov,
    std::pair<double, double>& fitQuality, TrackAtVertexRef trk, int sign,
    unsigned int nDimVertex) {
  if (nDimVertex == 3) {
    update<3>(vtxPos, vtxCov, fitQuality, trk, sign);
  } else if (nDimVertex == 4) {
    update<4>(vtxPos, vtxCov, fitQuality, trk, sign);
  } else {
    throw std::invalid_argument(
        "The vertex dimension must either be 3 (when fitting the spatial "
        "coordinates) or 4 (when fitting the spatial coordinates + time).");
  }
}

}  // namespace Acts
