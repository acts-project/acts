// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/KalmanVertexTrackUpdater.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"

namespace Acts {

namespace {

template <unsigned int nDimVertex>
Acts::BoundMatrix calculateTrackCovariance(
    const SquareMatrix3& wMat, const ActsMatrix<nDimVertex, 3>& crossCovVP,
    const ActsSquareMatrix<nDimVertex>& vtxCov,
    const BoundVector& newTrkParams) {
  // D_k^n
  ActsSquareMatrix<3> momCov =
      wMat + crossCovVP.transpose() * vtxCov.inverse() * crossCovVP;

  // Full x, y, z, phi, theta, q/p, and, optionally, t covariance matrix. Note
  // that we call this set of parameters "free" in the following even though
  // that is not quite correct (free parameters correspond to
  // x, y, z, t, px, py, pz)
  constexpr unsigned int nFreeParams = nDimVertex + 3;
  ActsSquareMatrix<nFreeParams> freeTrkCov(
      ActsSquareMatrix<nFreeParams>::Zero());

  freeTrkCov.template block<3, 3>(0, 0) = vtxCov.template block<3, 3>(0, 0);
  freeTrkCov.template block<3, 3>(0, 3) = crossCovVP.template block<3, 3>(0, 0);
  freeTrkCov.template block<3, 3>(3, 0) =
      crossCovVP.template block<3, 3>(0, 0).transpose();
  freeTrkCov.template block<3, 3>(3, 3) = momCov;

  // Fill time (cross-)covariances
  if constexpr (nFreeParams == 7) {
    freeTrkCov.template block<3, 1>(0, 6) = vtxCov.template block<3, 1>(0, 3);
    freeTrkCov.template block<3, 1>(3, 6) =
        crossCovVP.template block<1, 3>(3, 0).transpose();
    freeTrkCov.template block<1, 3>(6, 0) = vtxCov.template block<1, 3>(3, 0);
    freeTrkCov.template block<1, 3>(6, 3) =
        crossCovVP.template block<1, 3>(3, 0);
    freeTrkCov(6, 6) = vtxCov(3, 3);
  }

  // Jacobian relating "free" and bound track parameters
  constexpr unsigned int nBoundParams = nDimVertex + 2;
  ActsMatrix<nBoundParams, nFreeParams> freeToBoundJac(
      ActsMatrix<nBoundParams, nFreeParams>::Zero());

  // TODO: Jacobian is not quite correct
  // First row
  freeToBoundJac(0, 0) = -std::sin(newTrkParams[2]);
  freeToBoundJac(0, 1) = std::cos(newTrkParams[2]);

  double tanTheta = std::tan(newTrkParams[3]);

  // Second row
  freeToBoundJac(1, 0) = -freeToBoundJac(0, 1) / tanTheta;
  freeToBoundJac(1, 1) = freeToBoundJac(0, 0) / tanTheta;

  // Dimension of the part of the Jacobian that is an identity matrix
  constexpr unsigned int nDimIdentity = nFreeParams - 2;
  freeToBoundJac.template block<nDimIdentity, nDimIdentity>(1, 2) =
      ActsMatrix<nDimIdentity, nDimIdentity>::Identity();

  // Full perigee track covariance
  BoundMatrix boundTrackCov(BoundMatrix::Identity());
  boundTrackCov.block<nBoundParams, nBoundParams>(0, 0) =
      (freeToBoundJac * (freeTrkCov * freeToBoundJac.transpose()));

  return boundTrackCov;
}

template <unsigned int nDimVertex>
void updateImpl(TrackAtVertexRef track, const Vector4& vtxPosFull,
                const SquareMatrix4& vtxCovFull) {
  if constexpr (nDimVertex != 3 && nDimVertex != 4) {
    throw std::invalid_argument(
        "The vertex dimension must either be 3 (when fitting the spatial "
        "coordinates) or 4 (when fitting the spatial coordinates + time).");
  }

  using VertexVector = ActsVector<nDimVertex>;
  using VertexMatrix = ActsSquareMatrix<nDimVertex>;
  constexpr unsigned int nBoundParams = nDimVertex + 2;
  using ParameterVector = ActsVector<nBoundParams>;
  using ParameterMatrix = ActsSquareMatrix<nBoundParams>;
  // Check if linearized state exists
  if (!track.isLinearized) {
    throw std::invalid_argument("TrackAtVertex object must be linearized.");
  }

  // Extract vertex position and covariance
  // \tilde{x_n}
  const VertexVector vtxPos = vtxPosFull.template head<nDimVertex>();
  // C_n
  const VertexMatrix vtxCov =
      vtxCovFull.template block<nDimVertex, nDimVertex>(0, 0);

  // Get the linearized track
  const LinearizedTrack& linTrack = track.linearizedState;
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
  const ParameterMatrix trkParamWeight =
      linTrack.covarianceAtPCA.block<nBoundParams, nBoundParams>(0, 0)
          .inverse();

  // Cache object filled with zeros
  KalmanVertexUpdater::Cache<nDimVertex> cache;

  // Calculate the update of the vertex position when the track is removed. This
  // might be unintuitive, but it is needed to compute a symmetric chi2.
  KalmanVertexUpdater::calculateUpdate(vtxPosFull, vtxCovFull, linTrack,
                                       track.trackWeight, -1, cache);

  // Refit track momentum with the final vertex position
  Vector3 newTrkMomentum = cache.wMat * momJac.transpose() * trkParamWeight *
                           (trkParams - constTerm - posJac * vtxPos);

  // Track parameters, impact parameters are set to 0 and momentum corresponds
  // to newTrkMomentum. TODO: Make transition fitterParams -> fittedMomentum.
  BoundVector newTrkParams(BoundVector::Zero());

  // Get phi and theta and correct for possible periodicity changes
  const auto correctedPhiTheta =
      Acts::detail::normalizePhiTheta(newTrkMomentum(0), newTrkMomentum(1));
  newTrkParams(BoundIndices::eBoundPhi) = correctedPhiTheta.first;     // phi
  newTrkParams(BoundIndices::eBoundTheta) = correctedPhiTheta.second;  // theta
  newTrkParams(BoundIndices::eBoundQOverP) = newTrkMomentum(2);        // qOverP

  // E_k^n
  const ActsMatrix<nDimVertex, 3> crossCovVP =
      -vtxCov * posJac.transpose() * trkParamWeight * momJac * cache.wMat;

  // Difference in positions. cache.newVertexPos corresponds to \tilde{x_k^{n*}} in Ref. (1).
  VertexVector posDiff =
      vtxPos - cache.newVertexPos.template head<nDimVertex>();

  // r_k^n
  ParameterVector paramDiff =
      trkParams - (constTerm + posJac * vtxPos + momJac * newTrkMomentum);

  // New chi2 to be set later
  double chi2 =
      posDiff.dot(
          cache.newVertexWeight.template block<nDimVertex, nDimVertex>(0, 0) *
          posDiff) +
      paramDiff.dot(trkParamWeight * paramDiff);

  Acts::BoundMatrix newTrackCov = calculateTrackCovariance<nDimVertex>(
      cache.wMat, crossCovVP, vtxCov, newTrkParams);

  // Create new refitted parameters
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtxPos.template head<3>());

  BoundTrackParameters refittedPerigee =
      BoundTrackParameters(perigeeSurface, newTrkParams, std::move(newTrackCov),
                           track.fittedParams.particleHypothesis());

  // Set new properties
  track.fittedParams = refittedPerigee;
  track.chi2Track = chi2;
  track.ndf = 2 * track.trackWeight;

  return;
}

}  // namespace

void Acts::KalmanVertexTrackUpdater::update(TrackAtVertexRef track,
                                            const Vector4& vtxPosFull,
                                            const SquareMatrix4& vtxCovFull,
                                            unsigned int nDimVertex) {
  if (nDimVertex == 3) {
    updateImpl<3>(track, vtxPosFull, vtxCovFull);
  } else if (nDimVertex == 4) {
    updateImpl<4>(track, vtxPosFull, vtxCovFull);
  } else {
    throw std::invalid_argument(
        "The vertex dimension must either be 3 (when fitting the spatial "
        "coordinates) or 4 (when fitting the spatial coordinates + time).");
  }
}
}  // namespace Acts
