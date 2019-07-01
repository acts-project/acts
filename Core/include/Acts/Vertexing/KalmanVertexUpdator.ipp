// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename input_track_t>
void Acts::KalmanVertexUpdator<input_track_t>::addAndUpdate(
    Vertex<input_track_t>* vtx, TrackAtVertex<input_track_t> trk) const {
  update(vtx, trk, 1);
}

template <typename input_track_t>
Acts::Vertex<input_track_t>
Acts::KalmanVertexUpdator<input_track_t>::updatePosition(
    const Acts::Vertex<input_track_t>* vtx,
    const Acts::LinearizedTrack& linTrack, double trackWeight, int sign) const {
  // std::cout << "KalmanVertexUpdator::updatePosition: " << std::endl;

  // Retrieve linTrack information
  const SpacePointToBoundMatrix& posJac = linTrack.positionJacobian;
  const ActsMatrixD<BoundParsDim, 3>& momJac =
      linTrack.momentumJacobian;  // B_k in comments below
  const BoundVector& trkParams = linTrack.parametersAtPCA;
  const BoundVector& constTerm = linTrack.constantTerm;
  const BoundSymMatrix& trkParamWeight =
      linTrack.covarianceAtPCA.inverse();  // G_k in comments below

  // Vertex to be updated
  const SpacePointVector& oldVtxPos = vtx->fullPosition();
  const SpacePointSymMatrix& oldVtxWeight = vtx->fullCovariance().inverse();

  // W_k matrix
  ActsSymMatrixD<3> wMat =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // G_b = G_k - G_k*B_k*W_k*B_k^(T)*G_k^T
  BoundSymMatrix gBmat =
      trkParamWeight - trkParamWeight * (momJac * (wMat * momJac.transpose())) *
                           trkParamWeight.transpose();

  // new vertex cov matrix
  SpacePointSymMatrix newVtxCov =
      (oldVtxWeight +
       trackWeight * sign * posJac.transpose() * (gBmat * posJac))
          .inverse();

  // new vertex position
  SpacePointVector newVtxPos =
      newVtxCov *
      (oldVtxWeight * oldVtxPos + trackWeight * sign * posJac.transpose() *
                                      gBmat * (trkParams - constTerm));

  // create return vertex with new position
  // and covariance, but w/o tracks
  Vertex<input_track_t> returnVertex;

  // set position
  returnVertex.setFullPosition(newVtxPos);
  // set cov
  returnVertex.setFullCovariance(newVtxCov);
  // set fit quality
  returnVertex.setFitQuality(vtx->fitQuality().first, vtx->fitQuality().second);

  return returnVertex;
}

template <typename input_track_t>
double Acts::KalmanVertexUpdator<input_track_t>::vertexPositionChi2(
    const Vertex<input_track_t>* oldVtx,
    const Vertex<input_track_t>* newVtx) const {
  SpacePointSymMatrix oldWeight = oldVtx->fullCovariance().inverse();
  SpacePointVector posDiff = newVtx->fullPosition() - oldVtx->fullPosition();

  // calculate and return corresponding chi2
  return posDiff.transpose() * (oldWeight * posDiff);
}

template <typename input_track_t>
double Acts::KalmanVertexUpdator<input_track_t>::trackParametersChi2(
    const Vertex<input_track_t>& vtx, const LinearizedTrack& linTrack) const {
  const SpacePointVector& vtxPos = vtx.fullPosition();

  // track properties
  const SpacePointToBoundMatrix& posJac = linTrack.positionJacobian;
  const ActsMatrixD<BoundParsDim, 3>& momJac = linTrack.momentumJacobian;
  const BoundVector& trkParams = linTrack.parametersAtPCA;
  const BoundVector& constTerm = linTrack.constantTerm;
  const BoundSymMatrix& trkParamWeight = linTrack.covarianceAtPCA.inverse();

  // calculate temp matrix S
  ActsSymMatrixD<3> matS =
      (momJac.transpose() * (trkParamWeight * momJac)).inverse();

  // refitted track momentum
  Vector3D newTrackMomentum = matS * momJac.transpose() * trkParamWeight *
                              (trkParams - constTerm - posJac * vtxPos);

  // refitted track parameters
  auto newTrkParams = constTerm + posJac * vtxPos + momJac * newTrackMomentum;

  // parameter difference
  auto paramDiff = trkParams - newTrkParams;

  // return chi2
  return paramDiff.transpose() * (trkParamWeight * paramDiff);
}

template <typename input_track_t>
void Acts::KalmanVertexUpdator<input_track_t>::update(
    Vertex<input_track_t>* vtx, TrackAtVertex<input_track_t> trk,
    int sign) const {
  if (sign < 0) {
    // std::cout << "Error: Removal not supported yet. Returning." << std::endl;
    return;
  }

  double trackWeight = trk.trackWeight;

  Vertex<input_track_t> tempVtx =
      updatePosition(vtx, trk.linearizedState, trackWeight, sign);

  // get fit quality parameters wrt to old vertex
  std::pair fitQuality = vtx->fitQuality();
  double chi2 = fitQuality.first;
  double ndf = fitQuality.second;

  // chi2 wrt to track parameters
  double trkChi2 = trackParametersChi2(tempVtx, trk.linearizedState);

  // calculate new chi2
  chi2 += sign * (vertexPositionChi2(vtx, &tempVtx) + trackWeight * trkChi2);

  ndf += sign * trackWeight * 2;

  // updating the vertex
  vtx->setFullPosition(tempVtx.fullPosition());
  vtx->setFullCovariance(tempVtx.fullCovariance());
  vtx->setFitQuality(chi2, ndf);

  // add track to existing list of tracks at vertex
  if (sign > 1) {
    auto tracksAtVertex = vtx->tracks();
    // update track and add to list
    trk.chi2Track = trkChi2;
    trk.ndf = 2 * trackWeight;
    tracksAtVertex.push_back(trk);
    vtx->setTracksAtVertex(tracksAtVertex);
  }
}
