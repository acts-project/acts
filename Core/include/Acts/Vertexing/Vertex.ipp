// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename InputTrack>
Acts::Vertex<InputTrack>::Vertex(const Vector3D& position)
  : m_position(position)
{
}

template <typename InputTrack>
Acts::Vertex<InputTrack>::Vertex(const Vector3D&          position,
                                 const ActsSymMatrixD<3>& covariance,
                                 std::vector<TrackAtVertex<InputTrack>>& tracks)
  : m_position(position)
  , m_covariance(covariance)
  , m_tracksAtVertex(std::move(tracks))
{
}

template <typename InputTrack>
const Acts::Vector3D&
Acts::Vertex<InputTrack>::position() const
{
  return m_position;
}

template <typename InputTrack>
const Acts::ActsSymMatrixD<3>&
Acts::Vertex<InputTrack>::covariance() const
{
  return m_covariance;
}

template <typename InputTrack>
const std::vector<Acts::TrackAtVertex<InputTrack>>&
Acts::Vertex<InputTrack>::tracks() const
{
  return m_tracksAtVertex;
}

template <typename InputTrack>
std::pair<double, double>
Acts::Vertex<InputTrack>::fitQuality() const
{
  return std::pair<double, double>(m_chiSquared, m_numberDoF);
}

template <typename InputTrack>
void
Acts::Vertex<InputTrack>::setPosition(const Vector3D& position)
{
  m_position = position;
}

template <typename InputTrack>
void
Acts::Vertex<InputTrack>::setCovariance(const ActsSymMatrixD<3>& covariance)
{
  m_covariance = covariance;
}

template <typename InputTrack>
void
Acts::Vertex<InputTrack>::setTracksAtVertex(
    const std::vector<TrackAtVertex<InputTrack>>& tracks)
{
  m_tracksAtVertex = std::move(tracks);
}

template <typename InputTrack>
void
Acts::Vertex<InputTrack>::setFitQuality(double chiSquared, double numberDoF)
{
  m_chiSquared = chiSquared;
  m_numberDoF  = numberDoF;
}
