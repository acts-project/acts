// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename input_track_t>
Acts::Vertex<input_track_t>::Vertex(const Vector3D& position)
  : m_position(position)
{
}

template <typename input_track_t>
Acts::Vertex<input_track_t>::Vertex(
    const Vector3D&                            position,
    const ActsSymMatrixD<3>&                   covariance,
    std::vector<TrackAtVertex<input_track_t>>& tracks)
  : m_position(position)
  , m_covariance(covariance)
  , m_tracksAtVertex(std::move(tracks))
{
}

template <typename input_track_t>
const Acts::Vector3D&
Acts::Vertex<input_track_t>::position() const
{
  return m_position;
}

template <typename input_track_t>
const Acts::ActsSymMatrixD<3>&
Acts::Vertex<input_track_t>::covariance() const
{
  return m_covariance;
}

template <typename input_track_t>
const std::vector<Acts::TrackAtVertex<input_track_t>>&
Acts::Vertex<input_track_t>::tracks() const
{
  return m_tracksAtVertex;
}

template <typename input_track_t>
std::pair<double, double>
Acts::Vertex<input_track_t>::fitQuality() const
{
  return std::pair<double, double>(m_chiSquared, m_numberDoF);
}

template <typename input_track_t>
void
Acts::Vertex<input_track_t>::setPosition(const Vector3D& position)
{
  m_position = position;
}

template <typename input_track_t>
void
Acts::Vertex<input_track_t>::setCovariance(const ActsSymMatrixD<3>& covariance)
{
  m_covariance = covariance;
}

template <typename input_track_t>
void
Acts::Vertex<input_track_t>::setTracksAtVertex(
    const std::vector<TrackAtVertex<input_track_t>>& tracks)
{
  m_tracksAtVertex = std::move(tracks);
}

template <typename input_track_t>
void
Acts::Vertex<input_track_t>::setFitQuality(double chiSquared, double numberDoF)
{
  m_chiSquared = chiSquared;
  m_numberDoF  = numberDoF;
}
