// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"

template <typename input_track_t>
void Acts::KalmanVertexTrackUpdator<input_track_t>::update(
    const GeometryContext& gctx, TrackAtVertex<input_track_t>& track,
    const Vertex<input_track_t>* vtx) const {
  const SpacePointVector& vtxPos = vtx->fullPosition();

  // Get the linearized track
  const LinearizedTrack& linTrack = track.linearizedState;

  // Check if linearized state exists
  if (linTrack.covarianceAtPCA.trace() == 0.) {
    // Track has no linearized state, returning w/o update
    return;
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
      correctPhiThetaPeriodicity(newTrkMomentum(0), newTrkMomentum(1));

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
  Vertex<input_track_t> reducedVtx =
      m_cfg.vtx_updator.updatePosition(vtx, linTrack, track.trackWeight, -1);

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

  // Second row
  trkJac(1, 0) = -std::cos(newTrkParams[2]) / std::tan(newTrkParams[3]);
  trkJac(1, 1) = -std::sin(newTrkParams[2]) / std::tan(newTrkParams[3]);

  trkJac.block<5, 5>(1, 2) = ActsSymMatrixD<5>::Identity();

  // Full perigee track covariance
  std::unique_ptr<BoundMatrix> fullPerTrackCov =
      std::make_unique<BoundMatrix>(trkJac * (fullTrkCov * trkJac.transpose()));

  // Create new refitted parameters
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(VectorHelpers::position(vtxPos));

  BoundParameters refittedPerigee = BoundParameters(
      gctx, std::move(fullPerTrackCov), newTrkParams, perigeeSurface);

  // Set new properties
  track.fittedParams = refittedPerigee;
  track.chi2Track = chi2;
  track.trackWeight =
      2 * track.trackWeight;  // TODO: where is this factor coming from?
}

template <typename input_track_t>
std::pair<double, double>
Acts::KalmanVertexTrackUpdator<input_track_t>::correctPhiThetaPeriodicity(
    double phiIn, double thetaIn) const {
  double tmpPhi = std::fmod(phiIn, 2 * M_PI);  // temp phi
  if (tmpPhi > M_PI) {
    tmpPhi -= 2 * M_PI;
  }
  if (tmpPhi < -M_PI && tmpPhi > -2 * M_PI) {
    tmpPhi += 2 * M_PI;
  }

  double tmpTht = std::fmod(thetaIn, 2 * M_PI);  // temp theta
  if (tmpTht < -M_PI) {
    tmpTht = std::abs(tmpTht + 2 * M_PI);
  } else if (tmpTht < 0) {
    tmpTht *= -1;
    tmpPhi += M_PI;
    tmpPhi = tmpPhi > M_PI ? tmpPhi - 2 * M_PI : tmpPhi;
  }
  if (tmpTht > M_PI) {
    tmpTht = 2 * M_PI - tmpTht;
    tmpPhi += M_PI;
    tmpPhi = tmpPhi > M_PI ? (tmpPhi - 2 * M_PI) : tmpPhi;
  }

  return std::pair<double, double>(tmpPhi, tmpTht);
}