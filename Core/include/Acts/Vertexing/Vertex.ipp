// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename input_track_t>
Acts::Vertex<input_track_t>::Vertex(const Vector3D& position)

{
  m_position.setZero();
  VectorHelpers::position(m_position) = position;
}

template <typename input_track_t>
Acts::Vertex<input_track_t>::Vertex(const SpacePointVector& position)

{
  m_position = position;
}

template <typename input_track_t>
Acts::Vertex<input_track_t>::Vertex(
    const Vector3D& position, const ActsSymMatrixD<3>& covariance,
    const std::vector<TrackAtVertex<input_track_t>>& tracks)
    : m_tracksAtVertex(tracks) {
  m_position.setZero();
  VectorHelpers::position(m_position) = position;
  m_covariance.setZero();
  m_covariance.block<3, 3>(0, 0) = covariance;
}

template <typename input_track_t>
Acts::Vertex<input_track_t>::Vertex(
    const SpacePointVector& position, const SpacePointSymMatrix& covariance,
    const std::vector<TrackAtVertex<input_track_t>>& tracks)
    : m_position(position),
      m_covariance(covariance),
      m_tracksAtVertex(tracks) {}

template <typename input_track_t>
Acts::Vector3D Acts::Vertex<input_track_t>::position() const {
  return VectorHelpers::position(m_position);
}

template <typename input_track_t>
Acts::ParValue_t Acts::Vertex<input_track_t>::time() const {
  return VectorHelpers::time(m_position);
}

template <typename input_track_t>
const Acts::SpacePointVector& Acts::Vertex<input_track_t>::fullPosition()
    const {
  return m_position;
}

template <typename input_track_t>
Acts::ActsSymMatrixD<3> Acts::Vertex<input_track_t>::covariance() const {
  return m_covariance.block<3, 3>(0, 0);
}

template <typename input_track_t>
const Acts::SpacePointSymMatrix& Acts::Vertex<input_track_t>::fullCovariance()
    const {
  return m_covariance;
}

template <typename input_track_t>
const std::vector<Acts::TrackAtVertex<input_track_t>>&
Acts::Vertex<input_track_t>::tracks() const {
  return m_tracksAtVertex;
}

template <typename input_track_t>
std::pair<double, double> Acts::Vertex<input_track_t>::fitQuality() const {
  return std::pair<double, double>(m_chiSquared, m_numberDoF);
}

template <typename input_track_t>
void Acts::Vertex<input_track_t>::setPosition(const Vector3D& position,
                                              ParValue_t time) {
  m_position.setZero();
  VectorHelpers::position(m_position) = position;
  VectorHelpers::time(m_position) = time;
}

template <typename input_track_t>
void Acts::Vertex<input_track_t>::setFullPosition(
    const SpacePointVector& fullPosition) {
  m_position = fullPosition;
}

template <typename input_track_t>
void Acts::Vertex<input_track_t>::setTime(ParValue_t time) {
  VectorHelpers::time(m_position) = time;
}

template <typename input_track_t>
void Acts::Vertex<input_track_t>::setCovariance(
    const ActsSymMatrixD<3>& covariance) {
  m_covariance.setZero();
  m_covariance.block<3, 3>(0, 0) = covariance;
}

template <typename input_track_t>
void Acts::Vertex<input_track_t>::setFullCovariance(
    const SpacePointSymMatrix& covariance) {
  m_covariance = covariance;
}

template <typename input_track_t>
void Acts::Vertex<input_track_t>::setTracksAtVertex(
    const std::vector<TrackAtVertex<input_track_t>>& tracks) {
  m_tracksAtVertex = tracks;
}

template <typename input_track_t>
void Acts::Vertex<input_track_t>::setFitQuality(double chiSquared,
                                                double numberDoF) {
  m_chiSquared = chiSquared;
  m_numberDoF = numberDoF;
}

template <typename input_track_t>
void Acts::Vertex<input_track_t>::setFitQuality(
    std::pair<double, double> fitQuality) {
  m_chiSquared = fitQuality.first;
  m_numberDoF = fitQuality.second;
}
