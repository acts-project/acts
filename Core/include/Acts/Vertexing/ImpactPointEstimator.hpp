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
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"

namespace Acts {

struct ImpactParametersAndSigma {
  double IPd0 = 0.;
  double IPz0 = 0.;
  double IPz0SinTheta = 0.;
  double sigmad0 = 0.;
  double sigmaz0 = 0.;
  double sigmaz0SinTheta = 0.;
  double PVsigmad0 = 0.;
  double PVsigmaz0 = 0.;
  double PVsigmaz0SinTheta = 0.;
};

/// @class ImpactPointEstimator
///
/// @brief Estimator for impact point calculations
template <typename input_track_t, typename propagator_t,
          typename propagator_options_t = PropagatorOptions<>>
class ImpactPointEstimator {
 public:
  /// State struct
  struct State {
    /// @brief The state constructor
    ///
    /// @param fieldCacheIn The magnetic field cache
    State(MagneticFieldProvider::Cache fieldCacheIn)
        : fieldCache(std::move(fieldCacheIn)) {}
    /// Magnetic field cache
    MagneticFieldProvider::Cache fieldCache;
  };

  struct Config {
    /// @brief Config constructor if magnetic field is present
    ///
    /// @param bIn The magnetic field
    /// @param prop The propagator
    Config(std::shared_ptr<const MagneticFieldProvider> bIn,
           std::shared_ptr<const propagator_t> prop)
        : bField(std::move(bIn)), propagator(std::move(prop)) {}

    /// @brief Config constructor without B field -> uses NullBField
    /// provided)
    ///
    /// @param prop The propagator
    Config(std::shared_ptr<propagator_t> prop)
        : bField{std::make_shared<NullBField>()}, propagator(std::move(prop)) {}

    /// Magnetic field
    std::shared_ptr<const MagneticFieldProvider> bField;
    /// Propagator
    std::shared_ptr<const propagator_t> propagator;
    /// Max. number of iterations in Newton method
    int maxIterations = 20;
    /// Desired precision in deltaPhi in Newton method
    double precision = 1.e-10;
    /// Minimum q/p value
    double minQoP = 1e-15;
    /// Maximum curvature value
    double maxRho = 1e+15;
  };

  /// @brief Constructor
  ///
  /// @param cfg Configuration object
  ImpactPointEstimator(const Config& cfg) : m_cfg(cfg) {}

  /// @brief Calculates 3D distance between a track and a 3D point
  ///
  /// @param gctx The geometry context
  /// @param trkParams Track parameters
  /// @param vtxPos Position to calculate distance to
  /// @param state The state object
  ///
  /// @return Distance
  Result<double> calculate3dDistance(const GeometryContext& gctx,
                                     const BoundTrackParameters& trkParams,
                                     const Vector3& vtxPos, State& state) const;

  /// @brief Creates track parameters bound to plane
  /// at point of closest approach in 3d to given
  /// reference position. The parameters and errors
  /// are defined on the plane intersecting the track
  /// at point of closest approach, with track orthogonal
  /// to the plane and center of the plane defined as the
  /// given reference point (vertex).
  ///
  /// @param gctx The geometry context
  /// @param mctx The magnetic field context
  /// @param trkParams Track parameters
  /// @param vtxPos Reference position (vertex)
  /// @param state The state object
  ///
  /// @return New track params
  Result<BoundTrackParameters> estimate3DImpactParameters(
      const GeometryContext& gctx, const Acts::MagneticFieldContext& mctx,
      const BoundTrackParameters& trkParams, const Vector3& vtxPos,
      State& state) const;

  /// @brief Estimates the compatibility of a
  /// track to a vertex position based on the 3d
  /// distance between the track and the vertex
  ///
  /// @param gctx The Geometry context
  /// @param trkParams Track parameters at point of closest
  /// approach in 3d as retrieved by estimate3DImpactParameters
  /// @param vertexPos The vertex position
  ///
  /// @return The compatibility value
  Result<double> get3dVertexCompatibility(const GeometryContext& gctx,
                                          const BoundTrackParameters* trkParams,
                                          const Vector3& vertexPos) const;

  /// @brief Estimates the impact parameters and their errors of a given
  /// track w.r.t. a vertex by propagating the trajectory state
  /// towards the vertex position.
  ///
  /// @param track Track to estimate IP from
  /// @param vtx Vertex the track belongs to
  /// @param gctx The geometry context
  /// @param mctx The magnetic field context
  Result<ImpactParametersAndSigma> estimateImpactParameters(
      const BoundTrackParameters& track, const Vertex<input_track_t>& vtx,
      const GeometryContext& gctx, const MagneticFieldContext& mctx) const;

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
  Result<double> performNewtonApproximation(const Vector3& trkPos,
                                            const Vector3& vtxPos, double phi,
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
  /// @param state The state object
  Result<void> getDistanceAndMomentum(const GeometryContext& gctx,
                                      const BoundTrackParameters& trkParams,
                                      const Vector3& vtxPos, Vector3& deltaR,
                                      Vector3& momDir, State& state) const;
};

}  // namespace Acts

#include "Acts/Vertexing/ImpactPointEstimator.ipp"
