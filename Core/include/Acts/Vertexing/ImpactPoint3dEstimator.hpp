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
/// together with corresponding track parameters
template <typename bfield_t, typename input_track_t, typename propagator_t>
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

  /// @brief Constructor
  ///
  /// @param cfg Configuration object
  ImpactPoint3dEstimator(const Config& cfg) : m_cfg(cfg) {}


  /// @brief Calculates 3D distance between a track and a 3D point
  ///
  /// @param trkParams Track parameters
  /// @param refPos Position to calculate distance to
  ///
  /// @return Distance
  double calculateDistance(const BoundParameters& trkParams,
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
  /// @param trkParams Track parameters
  /// @param vtxPos Reference position (vertex)
  ///
  /// @return New track params
  Result<std::unique_ptr<const BoundParameters>> getParamsAtIP3d(
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      const BoundParameters& trkParams, const Vector3D& vtxPos) const;

 private:
  /// Configuration object
  const Config m_cfg;

  /// @brief: Performs a Newton approxation to retrieve a point
  /// of closest approach in 3D to a reference position
  ///
  /// @param trkPos Initial position
  /// @param vtxPos Reference position
  /// @param newPhi Phi along the helix which will be changed by
  ///        the Newton method
  /// @param theta Track theta
  /// @param r     Helix radius
  Result<void> performNewtonApproximation(const Vector3D& trkPos,
                                          const Vector3D& vtxPos,
                                          double& newPhi, double theta,
                                          double r) const;

};

}  // namespace Acts

#include "Acts/Vertexing/ImpactPoint3dEstimator.ipp"
