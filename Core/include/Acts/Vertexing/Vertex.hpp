// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
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
/// @tparam input_track_t Track object type
///
template <typename input_track_t>
class Vertex {
 public:
  /// @brief Default constructor
  Vertex() = default;

  /// @brief Construct for vertex at given 3d-position, sets covariance to zero
  ///
  /// @param position Vertex position
  Vertex(const Vector3D& position);

  /// @brief Construct for vertex at given 4d-position, sets covariance to zero
  ///
  /// @param position Vertex position
  Vertex(const Vector4D& position);

  /// @brief Vertex constructor
  ///
  /// @param position Vertex position
  /// @param covariance Position covariance matrix
  /// @param tracks Vector of tracks associated with the vertex
  Vertex(const Vector3D& position, const SymMatrix3D& covariance,
         const std::vector<TrackAtVertex<input_track_t>>& tracks);

  /// @brief Vertex constructor
  ///
  /// @param position Full vertex position
  /// @param covariance 4x4 covariance matrix
  /// @param tracks Vector of tracks associated with the vertex
  Vertex(const Vector4D& position, const SymMatrix4D& covariance,
         const std::vector<TrackAtVertex<input_track_t>>& tracks);

  /// @return Returns 3-position
  Vector3D position() const;

  /// @return Returns time
  BoundScalar time() const;

  /// @return Returns 4-position
  const Vector4D& fullPosition() const;

  /// @return Returns position covariance
  SymMatrix3D covariance() const;

  /// @return Returns 4x4 covariance
  const SymMatrix4D& fullCovariance() const;

  /// @return Returns vector of tracks associated with the vertex
  const std::vector<TrackAtVertex<input_track_t>>& tracks() const;

  /// @return Returns pair of (chi2, numberDoF)
  std::pair<double, double> fitQuality() const;

  /// @brief Set position and time
  ///
  /// @param position Vertex position
  /// @param time The time
  void setPosition(const Vector3D& position, BoundScalar time = 0);

  /// @brief Set position and time
  ///
  /// @param fullPosition Vertex position and time
  void setFullPosition(const Vector4D& fullPosition);

  /// @brief Sets time
  ///
  /// @param time The time
  void setTime(BoundScalar time);

  /// @brief Sets 3x3 covariance
  ///
  /// @param covariance Position covariance matrix
  void setCovariance(const ActsSymMatrixD<3>& covariance);

  /// @brief Sets 4x4 covariance
  ///
  /// @param covariance The 4x4 covariance matrix
  void setFullCovariance(const SymMatrix4D& covariance);

  /// @param tracks Vector of tracks at vertex
  void setTracksAtVertex(
      const std::vector<TrackAtVertex<input_track_t>>& tracks);

  /// @param chiSquared Chi2 of fit
  /// @param numberDoF Number of degrees of freedom
  void setFitQuality(double chiSquared, double numberDoF);

  /// @param fitQuality pair of (chi2, numberDoF)
  void setFitQuality(std::pair<double, double> fitQuality);

 private:
  Vector4D m_position = Vector4D::Zero();
  SymMatrix4D m_covariance = SymMatrix4D::Zero();
  std::vector<TrackAtVertex<input_track_t>> m_tracksAtVertex;
  double m_chiSquared = 0.;  // chi2 of the fit
  double m_numberDoF = 0.;   // number of degrees of freedom
};

}  // namespace Acts

#include "Vertex.ipp"
