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
  Vertex(const SpacePointVector& position);

  /// @brief Vertex constructor
  ///
  /// @param position Vertex position
  /// @param covariance Position covariance matrix
  /// @param tracks Vector of tracks associated with the vertex
  Vertex(const Vector3D& position, const ActsSymMatrixD<3>& covariance,
         const std::vector<TrackAtVertex<input_track_t>>& tracks);

  /// @brief Vertex constructor
  ///
  /// @param position Full vertex position
  /// @param covariance 4x4 covariance matrix
  /// @param tracks Vector of tracks associated with the vertex
  Vertex(const SpacePointVector& position,
         const SpacePointSymMatrix& covariance,
         const std::vector<TrackAtVertex<input_track_t>>& tracks);

  /// @return Returns 3-position
  Vector3D position() const;

  /// @return Returns time
  ParValue_t time() const;

  /// @return Returns 4-position
  const SpacePointVector& fullPosition() const;

  /// @return Returns position covariance
  ActsSymMatrixD<3> covariance() const;

  /// @return Returns 4x4 covariance
  const SpacePointSymMatrix& fullCovariance() const;

  /// @return Returns vector of tracks associated with the vertex
  const std::vector<TrackAtVertex<input_track_t>>& tracks() const;

  /// @return Returns pair of (chi2, numberDoF)
  std::pair<double, double> fitQuality() const;

  /// @brief Set position and time
  ///
  /// @param position Vertex position
  /// @param time The time
  void setPosition(const Vector3D& position, ParValue_t time = 0);

  /// @brief Set position and time
  ///
  /// @param fullPosition Vertex position and time
  void setFullPosition(const SpacePointVector& fullPosition);

  /// @brief Sets time
  ///
  /// @param time The time
  void setTime(ParValue_t time);

  /// @brief Sets 3x3 covariance
  ///
  /// @param covariance Position covariance matrix
  void setCovariance(const ActsSymMatrixD<3>& covariance);

  /// @brief Sets 4x4 covariance
  ///
  /// @param covariance The 4x4 covariance matrix
  void setFullCovariance(const SpacePointSymMatrix& covariance);

  /// @param tracks Vector of tracks at vertex
  void setTracksAtVertex(
      const std::vector<TrackAtVertex<input_track_t>>& tracks);

  /// @param chiSquared Chi2 of fit
  /// @param numberDoF Number of degrees of freedom
  void setFitQuality(double chiSquared, double numberDoF);

  /// @param fitQuality pair of (chi2, numberDoF)
  void setFitQuality(std::pair<double, double> fitQuality);

 private:
  SpacePointVector m_position = SpacePointVector::Zero();
  SpacePointSymMatrix m_covariance = SpacePointSymMatrix::Zero();
  std::vector<TrackAtVertex<input_track_t>> m_tracksAtVertex;
  double m_chiSquared = 1e9;  // chi2 of the fit
  double m_numberDoF = 0.;    // number of degrees of freedom
};

}  // namespace Acts

#include "Vertex.ipp"
