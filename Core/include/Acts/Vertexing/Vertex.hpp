// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"

namespace Acts {

/// @class Vertex
///
/// @brief Class for storing vertex objects
///
/// @tparam InputTrack Track object type
///
template <typename InputTrack>
class Vertex
{

public:
  /// @brief Default constructor
  Vertex() = default;

  /// @brief Construct for vertex at given position, sets covariance to zero
  ///
  /// @param position Vertex position
  Vertex(const Vector3D& position);

  /// @brief Vertex constructor
  ///
  /// @param position Vertex position
  /// @param covariance Position covariance matrix
  /// @param tracks Vector of tracks associated with the vertex
  Vertex(const Vector3D&                         position,
         const ActsSymMatrixD<3>&                covariance,
         std::vector<TrackAtVertex<InputTrack>>& tracks);

  /// @return Returns 3-position
  const Vector3D&
  position() const;
  /// @return Returns position covariance
  const ActsSymMatrixD<3>&
  covariance() const;

  /// @return Returns vector of tracks associated with the vertex
  const std::vector<TrackAtVertex<InputTrack>>&
  tracks() const;

  /// @return Returns pair of (chi2, numberDoF)
  std::pair<double, double>
  fitQuality() const;

  /// @param position Vertex position
  void
  setPosition(const Vector3D& position);

  /// @param covariance Position covariance matrix
  void
  setCovariance(const ActsSymMatrixD<3>& covariance);

  /// @param tracks Vector of tracks at vertex
  void
  setTracksAtVertex(const std::vector<TrackAtVertex<InputTrack>>& tracks);

  /// @param chiSquared Chi2 of fit
  /// @param numberDoF Number of degrees of freedom
  void
  setFitQuality(double chiSquared, double numberDoF);

private:
  Vector3D          m_position   = Vector3D(0., 0., 0.);
  ActsSymMatrixD<3> m_covariance = ActsSymMatrixD<3>::Zero();
  std::vector<TrackAtVertex<InputTrack>> m_tracksAtVertex;
  double m_chiSquared = std::numeric_limits<double>::max();  // chi2 of the fit
  double m_numberDoF  = 0;  // number of degrees of freedom
};

}  // namespace Acts

#include "Acts/Vertexing/Vertex.ipp"
