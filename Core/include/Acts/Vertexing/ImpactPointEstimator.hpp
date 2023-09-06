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
  // Impact parameters ...
  double d0 = 0.;
  double z0 = 0.;
  // ... and their standard deviations wrt a vertex, e.g.:
  // sigmaD0 = sqrt(Var(X) + Var(Y) + Var(d0)),
  // where X and Y are the x- and y-coordinate of the vertex
  double sigmaD0 = 0.;
  double sigmaZ0 = 0.;
  // Absolute difference in time between the vertex and the track at the 2D PCA
  // ...
  std::optional<double> deltaT = std::nullopt;
  // ... and standard deviation wrt a vertex
  std::optional<double> sigmaDeltaT = std::nullopt;
};

/// @class ImpactPointEstimator
///
/// @brief Estimator for impact point calculations
/// A description of the underlying mathematics can be found here:
/// https://github.com/acts-project/acts/pull/2414
/// TODO: Upload reference at a better place
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
  /// @param vtxPos 3D position to calculate the distance to
  /// @param state The state object
  ///
  /// @return Distance
  Result<double> calculateDistance(const GeometryContext& gctx,
                                   const BoundTrackParameters& trkParams,
                                   const Vector3& vtxPos, State& state) const;

  /// @brief Creates track parameters bound to plane
  /// at point of closest approach in 3D to given
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

  /// @brief Estimates the compatibility of a track to a vertex based on their
  /// 3D or 4D distance and the track covariance.
  ///
  /// @param gctx The Geometry context
  /// @param trkParams Track parameters at point of closest
  /// approach in 3D as retrieved by estimate3DImpactParameters
  /// @param vertexPos The vertex position
  ///
  /// @return The compatibility value
  template <unsigned int nDim = 3>
  Result<double> getVertexCompatibility(
      const GeometryContext& gctx, const BoundTrackParameters* trkParams,
      const ActsVector<nDim>& vertexPos) const;

  /// @brief Calculate the distance between a track and a vertex by finding the
  /// corresponding 3D PCA. Returns also the momentum direction at the 3D PCA.
  ///
  /// @param gctx Geometry context
  /// @param trkParams Track parameters
  /// @param vtxPos Vertex position
  /// @param state The state object
  /// @param massHypothesis Mass hypothesis
  /// @param chargeHypothesis Charge hypothesis
  template <unsigned int nDim = 3>
  Result<std::pair<Acts::ActsVector<nDim>, Acts::Vector3>>
  getDistanceAndMomentum(const GeometryContext& gctx,
                         const BoundTrackParameters& trkParams,
                         const ActsVector<nDim>& vtxPos, State& state,
                         const ActsScalar& massHypothesis = 0.13957018,
                         const ActsScalar& chargeHypothesis = 1.) const;

  /// @brief Calculates the impact parameters of a track w.r.t. a vertex. The
  /// corresponding errors are approximated by summing the variances of the
  /// track and the vertex.
  ///
  /// @param track Track whose impact parameters are calculated
  /// @param vtx Vertex corresponding to the track
  /// @param gctx The geometry context
  /// @param mctx The magnetic field context
  /// @param calculateTimeIP If true, the difference in time is computed
  Result<ImpactParametersAndSigma> getImpactParameters(
      const BoundTrackParameters& track, const Vertex<input_track_t>& vtx,
      const GeometryContext& gctx, const MagneticFieldContext& mctx,
      bool calculateTimeIP = false) const;

  /// @brief Estimates the sign of the 2D and Z lifetime of a given track
  /// w.r.t. a vertex and a direction (e.g. a jet direction)
  /// by propagating the trajectory state towards the vertex position
  /// and computing the scalar product with the direction vector
  ///
  /// @param track Track to estimate the IP from
  /// @param vtx   Vertex the track belongs to
  /// @param direction   The direction
  /// @param gctx  The geometry context
  /// @param mctx  The magnetic field context
  ///
  /// @return A pair holding the sign for the 2D and Z lifetimes
  Result<std::pair<double, double>> getLifetimeSignOfTrack(
      const BoundTrackParameters& track, const Vertex<input_track_t>& vtx,
      const Acts::Vector3& direction, const GeometryContext& gctx,
      const MagneticFieldContext& mctx) const;

  /// @brief Estimates the sign of the 3D lifetime of a given track
  /// w.r.t. a vertex and a direction (e.g. a jet direction)
  ///
  /// @param track Track to estimate the IP from
  /// @param vtx   Vertex the track belongs to
  /// @param direction   The direction
  /// @param gctx  The geometry context
  /// @param mctx  The magnetic field context
  ///
  /// @return The value of the 3D lifetime
  Result<double> get3DLifetimeSignOfTrack(
      const BoundTrackParameters& track, const Vertex<input_track_t>& vtx,
      const Acts::Vector3& direction, const GeometryContext& gctx,
      const MagneticFieldContext& mctx) const;

 private:
  /// Configuration object
  const Config m_cfg;

  /// @brief Performs a Newton approximation to retrieve a point
  /// of closest approach in 3D to a reference position
  ///
  /// @param helixCenter Position of the helix center
  /// @param vtxPos Vertex position
  /// @param phi Azimuthal momentum angle
  /// @note Modifying phi corresponds to moving along the track. This function
  /// optimizes phi until we reach a 3D PCA.
  /// @param theta Polar momentum angle (constant along the track)
  /// @param rho Signed helix radius
  ///
  /// @return Phi value at 3D PCA
  Result<double> performNewtonOptimization(const Vector3& helixCenter,
                                           const Vector3& vtxPos, double phi,
                                           double theta, double rho) const;
};

}  // namespace Acts

#include "Acts/Vertexing/ImpactPointEstimator.ipp"
