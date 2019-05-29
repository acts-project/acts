// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"

namespace Acts {

/// @class ImpactPoint3dEstimator
///
/// @brief Estimates point of closest approach in 3D
template <typename input_track_t>
class ImpactPoint3dEstimator {
 public:

  struct Config {
    Config(const bfield_t& bIn, const propagator_t& propagatorIn)
        : bField(bIn), propagator(propagatorIn) {}
    /// Magnetic field
    bfield_t bField;

    /// Propagator
    propagator_t propagator;

    /// Max. number of iterations in Newton method
    // TODO: correct value
    int maxIterations = 10;

    /// Desired precision in Newton method
    // TODO: correct value
    double precision = 1.e-4;
  };

  /// @brief Constructor used if input_track_t type == BoundParameters
  ///
  /// @param cfg Configuration object
  template <typename T = input_track_t,
            std::enable_if_t<std::is_same<T, BoundParameters>::value, int> = 0>
  ImpactPoint3dEstimator(const Config& cfg)
      : m_cfg(cfg), m_extractParameters([](T params) { return params; }) {}

  /// @brief Constructor for user-defined input_track_t type =! BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundParameters from input_track_t object
  ImpactPoint3dEstimator(const Config& cfg,
                         std::function<BoundParameters(input_track_t)> func)
      : m_cfg(cfg), m_extractParameters(func) {}


  /// @brief Calculates 3D distance between a track and a 3D point
  ///
  /// @param params Track parameters
  /// @param refPos Position to calculate distance to
  ///
  /// @return Distance
  double calculateDistance(const BoundParameters& params,
                           const Vector3D& refPos) const;

  /// @brief Creates track parameters bound to plane
  /// at point of closest approach in 3d to given
  /// reference position. The parameters and errors
  /// are defined on the plane intersecting the track
  /// at point of closest approach, with track ortogonal
  /// to the plane and center of the plane defined as the
  /// given reference point (vertex).
  ///
  /// @param geoCtx The geometry context
  /// @param trk Track at vertex
  /// @param refPos Reference position (vertex)
  ///
  /// @return New track params

  Result<BoundParameters> getParamsAtIP3d(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      const TrackAtVertex<input_track_t>& trk, const Vector3D& vtxPos) const;

 private:
  /// Configuration object
  const Config m_cfg;

  /// @brief Function to extract track parameters,
  /// input_track_t objects are BoundParameters by default, function to be
  /// overwritten to return BoundParameters for other input_track_t objects.
  ///
  /// @param input_track_t object to extract track parameters from
  const std::function<BoundParameters(input_track_t)> m_extractParameters;

  // TODO: add docs

  Result<void> performNewtonApproximation(const Vector3D& trkPos,
                                          const Vector3D& vtxPos,
                                          double& newPhi, double theta,
                                          double r) const;

};

}  // namespace Acts

#include "Acts/Vertexing/ImpactPoint3dEstimator.ipp"
