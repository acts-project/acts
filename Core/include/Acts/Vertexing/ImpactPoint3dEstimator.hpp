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

  /// @struct Configuration struct
  struct Config {
    /// @brief Configuration constructor
    ///
    /// @param bIn The magnetic field
    /// @param propagatorIn The propagator
    Config(const bfield_t& bIn, const propagator_t& propagatorIn)
        : bField(bIn), propagator(propagatorIn) {}
    /// Magnetic field
    bfield_t bField;

    /// Propagator
    propagator_t propagator;

    /// Max. number of iterations in Newton method
    int maxIterations = 20;

    /// Desired precision in deltaPhi in Newton method
    double precision = 1.e-10;
  };

  /// @brief Constructor
  ///
  /// @param cfg Configuration object
  ImpactPoint3dEstimator(const Config& cfg) : m_cfg(cfg) {}


  /// @brief Calculates 3D distance between a track and a 3D point
  ///
  /// @param gctx The geometry context
  /// @param trkParams Track parameters
  /// @param vtxPos Position to calculate distance to
  ///
  /// @return Distance
  Result<double> calculateDistance(const GeometryContext& gctx,
                                   const BoundParameters& trkParams,
                                   const Vector3D& vtxPos) const;

  /// @brief Creates track parameters bound to plane
  /// at point of closest approach in 3d to given
  /// reference position. The parameters and errors
  /// are defined on the plane intersecting the track
  /// at point of closest approach, with track orthogonal
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

  /// @brief Estimates the compatibility of a
  /// track to a vertex position based on the 3d
  /// distance between the track and the vertex
  ///
  /// @param gctx The Geometry context
  /// @param track Track parameters at point of closest
  /// approach in 3d as retrieved by getParamsAtIP3d
  /// @param vertexPos The vertex position
  ///
  /// @return The compatibility value
  double getVtxCompatibility(const GeometryContext& gctx,
                             const BoundParameters* trkParams,
                             const Vector3D& vertexPos) const;

 private:
  /// Configuration object
  const Config m_cfg;

  /// @brief Performs a Newton approximation to retrieve a point
  /// of closest approach in 3D to a reference position
  ///
  /// @param trkPos Initial position
  /// @param vtxPos Reference position
  /// @param phi Phi along the helix which will be changed by
  ///        the Newton method
  /// @param theta Track theta
  /// @param r     Helix radius
  ///
  /// @return New phi value
  Result<double> performNewtonApproximation(const Vector3D& trkPos,
                                            const Vector3D& vtxPos, double phi,
                                            double theta, double r) const;

  /// @brief Helper function to calculate relative
  /// distance between track and vtxPos and the
  /// direction of the momentum
  ///
  /// @param gctx The geometry context
  /// @param trkParams Track parameters
  /// @param vtxPos The vertex position
  /// @param deltaR Relative position between
  ///   track and vtxPos, to be determined by method
  /// @param momDir Momentum direction, to be
  ///   determined by method
  Result<void> getDistanceAndMomentum(const GeometryContext& gctx,
                                      const BoundParameters& trkParams,
                                      const Vector3D& vtxPos, Vector3D& deltaR,
                                      Vector3D& momDir) const;
};

}  // namespace Acts

#include "Acts/Vertexing/ImpactPoint3dEstimator.ipp"
