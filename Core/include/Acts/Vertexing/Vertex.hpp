// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"

namespace Acts {

class Vertex
{

public:
  /// Default constructor
  Vertex();

  Vertex(const Acts::Vector3D&                              position,
         const Acts::ActsSymMatrixD<3>&                     covariance,
         std::vector<std::unique_ptr<Acts::TrackAtVertex>>& tracks);

  /// Return 3-position
  const Acts::Vector3D&
  position() const;
  /// Return covariance
  const Acts::ActsSymMatrixD<3>&
  covariance() const;

  const std::vector<std::unique_ptr<Acts::TrackAtVertex>>&
  tracks() const;

  /// Set 3-position
  void
  setPosition(const Acts::Vector3D& position);
  /// Set covariance
  void
  setCovariance(const Acts::ActsSymMatrixD<3>& covariance);
  /// Set tracks at vertex
  void
  setTracksAtVertex(std::vector<std::unique_ptr<Acts::TrackAtVertex>>& tracks);

private:
  Acts::Vector3D                                    m_position;
  Acts::ActsSymMatrixD<3>                           m_covariance;
  std::vector<std::unique_ptr<Acts::TrackAtVertex>> m_tracksAtVertex;
};

}  // namespace Acts
