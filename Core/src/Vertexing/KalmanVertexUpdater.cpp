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

}  // namespace Acts
