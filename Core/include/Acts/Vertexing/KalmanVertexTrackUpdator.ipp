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
    const Vertex<input_track_t>& vtx) const {
  const Vector3D& vtxPos = vtx.position();
  // TODO: const ref or ref?
  LinearizedTrack& linTrack = track.linearizedState;

  // check if lin state exists // TODO: different check?
  if (linTrack.covarianceAtPCA.trace() == 0.) {
    if (track.trackWeight > m_cfg.maxWeight) {
      std::cout << "Current track does not have a LinearizedTrack. \n"
                << "Track returned not refitted." << std::endl;
    } else {
      std::cout << "Current track does not have a LinearizedTrack. \n"
                << "Track weight however smaller than " << m_cfg.maxWeight
                << "Track returned not refitted." << std::endl;
    }

    return;
  }

  // Retrieve linTrack information
  const auto& posJac = linTrack.positionJacobian;
  const auto& momJac = linTrack.momentumJacobian;
  const auto& trkParams = linTrack.parametersAtPCA;
  const auto& trkParamWeight = linTrack.covarianceAtPCA.inverse();

  // Calculate S matrix // TODO: check paper
  ActsSymMatrixD<3> sMat =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  const ActsVectorD<5>& residual = linTrack.constantTerm;

  // refit track momentum
  Vector3D newTrkMomentum = sMat * momJac.transpose() * trkParamWeight *
                            (trkParams - residual - posJac * vtxPos);

  // refit track parameters
  ActsVectorD<5> newTrkParams(ActsVectorD<5>::Zero());

  // Get phi and theta and correct for possible periodicity changes
  auto correctedPhiTheta =
      correctPhiThetaPeriodicity(newTrkMomentum(0), newTrkMomentum(1));

  newTrkParams(ParID_t::ePHI) = correctedPhiTheta.first;     // phi
  newTrkParams(ParID_t::eTHETA) = correctedPhiTheta.second;  // theta
  newTrkParams(ParID_t::eQOP) = newTrkMomentum(2);           // qOverP

  // vertex covariance and weight matrices
  const ActsSymMatrixD<3>& vtxCov = vtx.covariance();
  const ActsSymMatrixD<3> vtxWeight = vtxCov.inverse();

  // new track covariance matrix
  ActsSymMatrixD<3> newTrkCov =
      -vtxCov * posJac.transpose() * trkParamWeight * momJac * sMat;

  // now determine the smoothed chi2 of the track in the following
  // get updated position, this removes track from vtx
  Vertex<input_track_t> reducedVtx =
      m_cfg.vtx_updator.updatePosition(vtx, linTrack, track.trackWeight, -1);

  // corresponding weight matrix
  const ActsSymMatrixD<3> reducedVtxWeight = reducedVtx.covariance().inverse();

  // difference in positions
  Vector3D posDiff = vtx.position() - reducedVtx.position();

  // get smoothed params
  ActsVectorD<5> smParams = trkParams - (residual + posJac * vtx.position() +
                                         momJac * newTrkMomentum);

  // new chi2 to be set later
  double chi2 = posDiff.dot(reducedVtxWeight * posDiff) +
                smParams.dot(trkParamWeight * smParams);

  // now new momentum covariance
  ActsSymMatrixD<3> momCov =
      sMat + newTrkCov.transpose() * (vtxWeight * newTrkCov);

  // full (x,y,z,phi, theta, q/p) covariance matrix
  ActsSymMatrixD<6> fullTrkCov(ActsSymMatrixD<6>::Zero());

  fullTrkCov.block<3, 3>(0, 0) = vtxCov;
  fullTrkCov.block<3, 3>(0, 3) = newTrkCov;
  fullTrkCov.block<3, 3>(3, 0) = newTrkCov.transpose();
  fullTrkCov.block<3, 3>(3, 3) = momCov;

  // following is 'temporary hack', as implemented and used in athena
  ActsMatrixD<5, 6> trkJac(ActsMatrixD<5, 6>::Zero());

  trkJac(4, 5) = 1.;
  trkJac(3, 4) = 1.;
  trkJac(2, 3) = 1.;

  // first row
  trkJac(0, 0) = -std::sin(newTrkParams[2]);
  trkJac(0, 1) = std::cos(newTrkParams[2]);

  // second row
  trkJac(1, 0) = -std::cos(newTrkParams[2]) / tan(newTrkParams[3]);
  trkJac(1, 1) = -std::sin(newTrkParams[2]) / tan(newTrkParams[3]);
  trkJac(1, 2) = 1.;

  // full perigee track covariance
  std::unique_ptr<ActsSymMatrixD<5>> fullPerTrackCov =
      std::make_unique<ActsSymMatrixD<5>>(trkJac *
                                          (fullTrkCov * trkJac.transpose()));

  // TODO: why at vtxPos and not at updated position?
  std::shared_ptr<PerigeeSurface> perigeeSurface =
      Surface::makeShared<PerigeeSurface>(vtxPos);

  BoundParameters refittedPerigee = BoundParameters(
      gctx, std::move(fullPerTrackCov), newTrkParams, perigeeSurface);

  // set new properties
  track.fittedParams = refittedPerigee;
  track.chi2Track = chi2;
  track.trackWeight = 2 * track.trackWeight;  // TODO: is that factor correct?
                                              // If yes, why, if not, just do
                                              // not touch trackWeight
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