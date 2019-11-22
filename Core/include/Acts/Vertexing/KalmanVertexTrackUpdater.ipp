// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/detail/periodic.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

template <typename input_track_t>
Acts::Result<void> Acts::KalmanVertexTrackUpdater::update(
    const GeometryContext& gctx, TrackAtVertex<input_track_t>& track,
    const Vertex<input_track_t>* vtx) {
  if (vtx == nullptr) {
    return VertexingError::EmptyInput;
  }

  const SpacePointVector& vtxPos = vtx->fullPosition();

  // Get the linearized track
  const LinearizedTrack& linTrack = track.linearizedState;

  // Check if linearized state exists
  if (linTrack.covarianceAtPCA.determinant() == 0.) {
    // Track has no linearized state, returning w/o update
    return {};
  }

  // Retrieve linTrack information
  const SpacePointToBoundMatrix& posJac = linTrack.positionJacobian;
  const ActsMatrixD<BoundParsDim, 3>& momJac = linTrack.momentumJacobian;
  const BoundVector& trkParams = linTrack.parametersAtPCA;
  const BoundSymMatrix& trkParamWeight = linTrack.covarianceAtPCA.inverse();

  // Calculate S matrix
  ActsSymMatrixD<3> sMat =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  const BoundVector& residual = linTrack.constantTerm;

  // Refit track momentum
  Vector3D newTrkMomentum = sMat * momJac.transpose() * trkParamWeight *
                            (trkParams - residual - posJac * vtxPos);

  // Refit track parameters
  BoundVector newTrkParams(BoundVector::Zero());

  // Get phi and theta and correct for possible periodicity changes
  auto correctedPhiTheta =
      Acts::detail::ensureThetaBounds(newTrkMomentum(0), newTrkMomentum(1));

  newTrkParams(ParID_t::ePHI) = correctedPhiTheta.first;     // phi
  newTrkParams(ParID_t::eTHETA) = correctedPhiTheta.second;  // theta
  newTrkParams(ParID_t::eQOP) = newTrkMomentum(2);           // qOverP

  // Vertex covariance and weight matrices
  const SpacePointSymMatrix& vtxCov = vtx->fullCovariance();
  const SpacePointSymMatrix vtxWeight = vtxCov.inverse();

  // New track covariance matrix
  ActsMatrixD<SpacePointDim, 3> newTrkCov =
      -vtxCov * posJac.transpose() * trkParamWeight * momJac * sMat;

  // Now determine the smoothed chi2 of the track in the following
  // get updated position, this removes track from vtx
  auto res = KalmanVertexUpdater::updatePosition<input_track_t>(
      vtx, linTrack, track.trackWeight, -1);

  if (!res.ok()) {
    return res.error();
  }

  Vertex<input_track_t> reducedVtx = *res;

  // Corresponding weight matrix
  const SpacePointSymMatrix reducedVtxWeight =
      reducedVtx.fullCovariance().inverse();

  // Difference in positions
  SpacePointVector posDiff = vtx->fullPosition() - reducedVtx.fullPosition();

  // Get smoothed params
  BoundVector smParams = trkParams - (residual + posJac * vtx->fullPosition() +
                                      momJac * newTrkMomentum);

  // New chi2 to be set later
  double chi2 = posDiff.dot(reducedVtxWeight * posDiff) +
                smParams.dot(trkParamWeight * smParams);

  const BoundMatrix& fullPerTrackCov = detail::createFullTrackCovariance(
      sMat, newTrkCov, vtxWeight, vtxCov, newTrkParams);

  // Create new refitted parameters
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(VectorHelpers::position(vtxPos));

  BoundParameters refittedPerigee = BoundParameters(
      gctx, std::move(fullPerTrackCov), newTrkParams, perigeeSurface);

  // Set new properties
  track.fittedParams = refittedPerigee;
  track.chi2Track = chi2;
  track.ndf = 2 * track.trackWeight;

  return {};
}

Acts::BoundMatrix
Acts::KalmanVertexTrackUpdater::detail::createFullTrackCovariance(
    const ActsSymMatrixD<3>& sMat,
    const ActsMatrixD<SpacePointDim, 3>& newTrkCov,
    const SpacePointSymMatrix& vtxWeight, const SpacePointSymMatrix& vtxCov,
    const BoundVector& newTrkParams) {
  // Now new momentum covariance
  ActsSymMatrixD<3> momCov =
      sMat + newTrkCov.transpose() * (vtxWeight * newTrkCov);

  // Full (x,y,z,phi, theta, q/p, t) covariance matrix
  ActsSymMatrixD<7> fullTrkCov(ActsSymMatrixD<7>::Zero());

  fullTrkCov.block<4, 4>(0, 0) = vtxCov;
  fullTrkCov.block<4, 3>(0, 4) = newTrkCov;
  fullTrkCov.block<3, 4>(4, 0) = newTrkCov.transpose();
  fullTrkCov.block<3, 3>(4, 4) = momCov;

  // Combined track jacobian
  ActsMatrixD<BoundParsDim, 7> trkJac(ActsMatrixD<BoundParsDim, 7>::Zero());

  // First row
  trkJac(0, 0) = -std::sin(newTrkParams[2]);
  trkJac(0, 1) = std::cos(newTrkParams[2]);

  double tanTheta = std::tan(newTrkParams[3]);

  // Second row
  trkJac(1, 0) = -trkJac(0, 1) / tanTheta;
  trkJac(1, 1) = trkJac(0, 0) / tanTheta;

  trkJac.block<5, 5>(1, 2) = ActsSymMatrixD<5>::Identity();

  // Full perigee track covariance
  BoundMatrix fullPerTrackCov(trkJac * (fullTrkCov * trkJac.transpose()));

  return fullPerTrackCov;
}
