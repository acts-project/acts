// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename input_track_t>
double Acts::Chi2TrackCompatibilityEstimator<input_track_t>::getCompatibility(
    const GeometryContext& gctx, const BoundParameters& track,
    const Vector3D& vertexPos) const {
  // surface rotation
  RotationMatrix3D myRotation =
      track.referenceSurface().transform(gctx).rotation();
  // surface translation
  Vector3D myTranslation =
      track.referenceSurface().transform(gctx).translation();

  // x and y direction of plane
  Vector3D xDirPlane = myRotation.col(0);
  Vector3D yDirPlane = myRotation.col(1);

  // transform vertex position in local plane reference frame
  Vector3D vertexLocPlane = vertexPos - myTranslation;

  // local x/y vertex position
  Vector2D vertexLocXY{vertexLocPlane.dot(xDirPlane),
                       vertexLocPlane.dot(yDirPlane)};

  // track covariance
  auto cov = track.covariance();
  ActsSymMatrixD<2> myWeightXY = (*cov).block<2, 2>(0, 0).inverse();

  // TODO: is parameters()[eX] etc correct?
  Vector2D myXYpos =
      Vector2D(track.parameters()[eX], track.parameters()[eY]) - vertexLocXY;

  // return chi2
  return myXYpos.dot(myWeightXY * myXYpos);
}

template <typename input_track_t>
void Acts::Chi2TrackCompatibilityEstimator<input_track_t>::
    setTrackCompatibility(const GeometryContext& gctx,
                          TrackAtVertex<input_track_t>& trackAtVertex,
                          const Vector3D& vertexPos,
                          const std::function<BoundParameters(input_track_t)>
                              extractParameters) const {
  trackAtVertex.vertexCompatibility = getCompatibility(
      gctx, extractParameters(trackAtVertex.originalTrack), vertexPos);
}
