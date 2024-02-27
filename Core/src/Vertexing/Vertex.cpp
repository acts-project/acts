// This file is part of the Acts project.
//
// Copyright (C) 2019-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/Vertex.hpp"

Acts::Vertex::Vertex(const Vector3& position) {
  m_position[ePos0] = position[ePos0];
  m_position[ePos1] = position[ePos1];
  m_position[ePos2] = position[ePos2];
}

Acts::Vertex::Vertex(const Vector4& position) : m_position(position) {}

Acts::Vertex::Vertex(const Vector3& position, const SquareMatrix3& covariance,
                     std::vector<TrackAtVertex> tracks)
    : m_tracksAtVertex(std::move(tracks)) {
  m_position[ePos0] = position[ePos0];
  m_position[ePos1] = position[ePos1];
  m_position[ePos2] = position[ePos2];
  m_covariance.block<3, 3>(ePos0, ePos0) = covariance;
}

Acts::Vertex::Vertex(const Vector4& position, const SquareMatrix4& covariance,
                     std::vector<TrackAtVertex> tracks)
    : m_position(position),
      m_covariance(covariance),
      m_tracksAtVertex(std::move(tracks)) {}

Acts::Vector3 Acts::Vertex::position() const {
  return VectorHelpers::position(m_position);
}

Acts::ActsScalar Acts::Vertex::time() const {
  return m_position[eTime];
}

const Acts::Vector4& Acts::Vertex::fullPosition() const {
  return m_position;
}

Acts::Vector4& Acts::Vertex::fullPosition() {
  return m_position;
}

Acts::SquareMatrix3 Acts::Vertex::covariance() const {
  return m_covariance.block<3, 3>(ePos0, ePos0);
}

const Acts::SquareMatrix4& Acts::Vertex::fullCovariance() const {
  return m_covariance;
}

Acts::SquareMatrix4& Acts::Vertex::fullCovariance() {
  return m_covariance;
}

const std::vector<Acts::TrackAtVertex>& Acts::Vertex::tracks() const {
  return m_tracksAtVertex;
}

std::pair<double, double> Acts::Vertex::fitQuality() const {
  return std::pair<double, double>(m_chiSquared, m_numberDoF);
}

void Acts::Vertex::setPosition(const Vector3& position, ActsScalar time) {
  m_position[ePos0] = position[ePos0];
  m_position[ePos1] = position[ePos1];
  m_position[ePos2] = position[ePos2];
  m_position[eTime] = time;
}

void Acts::Vertex::setFullPosition(const Vector4& fullPosition) {
  m_position = fullPosition;
}

void Acts::Vertex::setTime(ActsScalar time) {
  m_position[eTime] = time;
}

void Acts::Vertex::setCovariance(const SquareMatrix3& covariance) {
  m_covariance.setZero();
  m_covariance.block<3, 3>(ePos0, ePos0) = covariance;
}

void Acts::Vertex::setFullCovariance(const SquareMatrix4& covariance) {
  m_covariance = covariance;
}

void Acts::Vertex::setTracksAtVertex(std::vector<TrackAtVertex> tracks) {
  m_tracksAtVertex = std::move(tracks);
}

void Acts::Vertex::setFitQuality(double chiSquared, double numberDoF) {
  m_chiSquared = chiSquared;
  m_numberDoF = numberDoF;
}

void Acts::Vertex::setFitQuality(std::pair<double, double> fitQuality) {
  m_chiSquared = fitQuality.first;
  m_numberDoF = fitQuality.second;
}
