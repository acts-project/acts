// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

inline Acts::Vertex::Vertex(const Vector3& position) {
  m_position.head<3>() = position;
  m_seedPosition.head<3>() = position;
}

inline Acts::Vertex::Vertex(const Vector4& position)
    : m_position(position), m_seedPosition(position) {}

inline Acts::Vertex::Vertex(const Vector3& position,
                            const SquareMatrix3& covariance,
                            std::vector<TrackAtVertex> tracks)
    : m_tracksAtVertex(std::move(tracks)) {
  m_position.head<3>() = position;
  m_seedPosition.head<3>() = position;
  m_covariance.block<3, 3>(ePos0, ePos0) = covariance;
}

inline Acts::Vertex::Vertex(const Vector4& position,
                            const SquareMatrix4& covariance,
                            std::vector<TrackAtVertex> tracks)
    : m_position(position),
      m_seedPosition(position),
      m_covariance(covariance),
      m_tracksAtVertex(std::move(tracks)) {}

inline Acts::Vector3 Acts::Vertex::position() const {
  return VectorHelpers::position(m_position);
}

inline Acts::ActsScalar Acts::Vertex::time() const {
  return m_position[eTime];
}

inline const Acts::Vector4& Acts::Vertex::fullPosition() const {
  return m_position;
}

inline Acts::Vector4& Acts::Vertex::fullPosition() {
  return m_position;
}

inline Acts::Vector4& Acts::Vertex::fullSeedPosition() {
  return m_seedPosition;
}

inline Acts::SquareMatrix3 Acts::Vertex::covariance() const {
  return m_covariance.block<3, 3>(ePos0, ePos0);
}

inline const Acts::SquareMatrix4& Acts::Vertex::fullCovariance() const {
  return m_covariance;
}

inline Acts::SquareMatrix4& Acts::Vertex::fullCovariance() {
  return m_covariance;
}

inline const std::vector<Acts::TrackAtVertex>& Acts::Vertex::tracks() const {
  return m_tracksAtVertex;
}

inline std::pair<double, double> Acts::Vertex::fitQuality() const {
  return std::pair<double, double>(m_chiSquared, m_numberDoF);
}

inline void Acts::Vertex::setPosition(const Vector3& position,
                                      ActsScalar time) {
  m_position[ePos0] = position[ePos0];
  m_position[ePos1] = position[ePos1];
  m_position[ePos2] = position[ePos2];
  m_position[eTime] = time;
}

inline void Acts::Vertex::setFullPosition(const Vector4& fullPosition) {
  m_position = fullPosition;
}

inline void Acts::Vertex::setTime(ActsScalar time) {
  m_position[eTime] = time;
}

inline void Acts::Vertex::setCovariance(const SquareMatrix3& covariance) {
  m_covariance.setZero();
  m_covariance.block<3, 3>(ePos0, ePos0) = covariance;
}

inline void Acts::Vertex::setFullCovariance(const SquareMatrix4& covariance) {
  m_covariance = covariance;
}

inline void Acts::Vertex::setTracksAtVertex(std::vector<TrackAtVertex> tracks) {
  m_tracksAtVertex = std::move(tracks);
}

inline void Acts::Vertex::setFitQuality(double chiSquared, double numberDoF) {
  m_chiSquared = chiSquared;
  m_numberDoF = numberDoF;
}

inline void Acts::Vertex::setFitQuality(std::pair<double, double> fitQuality) {
  m_chiSquared = fitQuality.first;
  m_numberDoF = fitQuality.second;
}
