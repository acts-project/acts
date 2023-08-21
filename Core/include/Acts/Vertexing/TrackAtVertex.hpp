// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"

#include <functional>

namespace Acts {

/// @struct FittedMomentum
///
/// @brief Vertex fitters return a vertex position and updated track momenta at said position
/// (+corresponding covariances). The updated track momenta and their
/// covariances are saved in the following struct.
struct FittedMomentum {
  FittedMomentum(const Vector3& mom, const ActsSymMatrix<3>& cov)
      : momentum(mom), covariance(cov) {}

  Vector3 momentum;
  ActsSymMatrix<3> covariance;
};

/// @struct TrackAtVertex
///
/// @brief Defines a track at vertex object
///
/// @tparam input_track_t Track object type
template <typename input_track_t>
struct TrackAtVertex {
  /// Deleted default constructor
  TrackAtVertex() = delete;

  /// @brief Constructor used before the vertex fit (i.e., when we don't know the fitted momentum yet)
  ///
  /// @param chi2PerTrack Chi2 of the track
  /// @param originalTrack Original perigee parameter
  TrackAtVertex(const input_track_t* originalTrack, double chi2PerTrack)
      : originalParams(originalTrack), chi2(chi2PerTrack) {}

  /// @brief Constructor used before the vertex fit (i.e., when we don't know the fitted momentum yet) with default chi2
  ///
  /// @param originalTrack Original perigee parameter
  TrackAtVertex(const input_track_t* originalTrack)
      : originalParams(originalTrack) {}

  /// @brief Constructor used when we know the momentum after the fit
  ///
  /// @param chi2PerTrack Chi2 of the track
  /// @param fittedMom updated momentum after the vertex fit
  /// @param originalTrack Original perigee parameter
  TrackAtVertex(const input_track_t* originalTrack,
                std::optional<FittedMomentum> fittedMom, double chi2PerTrack)
      : originalParams(originalTrack),
        fittedMomentum(std::move(fittedMom)),
        chi2(chi2PerTrack) {}

  /// @brief Constructor used when we know the momentum after the fit with default chi2
  ///
  /// @param fittedMom updated momentum after the vertex fit
  /// @param originalTrack Original perigee parameter
  TrackAtVertex(const input_track_t* originalTrack,
                std::optional<FittedMomentum> fittedMom)
      : originalParams(originalTrack), fittedMomentum(std::move(fittedMom)) {}

  /// Original track parameters
  const input_track_t* originalParams;

  /// Momentum after the vertex fit
  /// The full track parameters after the vertex fit are:
  /// d0 = 0,
  /// z0 = 0,
  /// phi = fittedMomentum.momentum(0),
  /// theta = fittedMomentum.momentum(1),
  /// qOverP = fittedMomentum.momentum(2).
  std::optional<FittedMomentum> fittedMomentum = std::nullopt;

  /// Chi2 of the track
  double chi2 = 0;

  /// Number degrees of freedom
  /// Note: Can be different from integer value since annealing can result in
  /// effective non-integer values.
  double ndf = 0;

  /// Value of the compatibility of the track to the actual vertex, based
  /// on the estimation of the 3d distance between the track and the vertex
  double vertexCompatibility = 0;

  /// Weight of track in fit
  double weight = 1;

  /// The linearized state of the track at vertex
  LinearizedTrack linearizedState;

  /// Is already linearized
  bool isLinearized = false;
};

}  // namespace Acts
