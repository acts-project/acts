// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Tolerance.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/TrackParameterHelpers.hpp"
#include "Acts/EventData/TransformationHelpers.hpp"
#include "Acts/EventData/detail/PrintParameters.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <cassert>
#include <cmath>
#include <memory>

namespace Acts {

/// Track parameters bound to a reference surface for a single track.
///
/// @tparam particle_hypothesis_t Helper type to interpret the particle charge/momentum
///
/// This is intended as a user-facing data class that adds additional accessors
/// and charge/momentum interpretation on top of the pure parameters vector. All
/// parameters and their corresponding covariance matrix are stored in bound
/// parametrization. The specific definition of the local spatial parameters is
/// defined by the associated surface.
///
/// @note This class holds shared ownership on its reference surface.
template <class particle_hypothesis_t>
class GenericBoundTrackParameters {
 public:
  /// Type alias for bound track parameters vector
  using ParametersVector = BoundVector;
  /// Type alias for bound track covariance matrix
  using CovarianceMatrix = BoundMatrix;
  /// Type alias for particle hypothesis used in reconstruction
  using ParticleHypothesis = particle_hypothesis_t;

  /// Factory to construct from four-position, direction, absolute momentum, and
  /// charge.
  ///
  /// @param geoCtx Geometry context for the local-to-global transformation
  /// @param surface Reference surface the parameters are defined on
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored
  /// @param qOverP Charge over momentum
  /// @param cov Bound parameters covariance matrix
  /// @param particleHypothesis Particle hypothesis
  /// @param tolerance Tolerance used for globalToLocal
  ///
  /// @return Result containing the constructed bound track parameters or error
  /// @note The returned result indicates whether the free parameters could
  /// successfully be converted to on-surface parameters.
  static Result<GenericBoundTrackParameters> create(
      const GeometryContext& geoCtx, std::shared_ptr<const Surface> surface,
      const Vector4& pos4, const Vector3& dir, double qOverP,
      std::optional<CovarianceMatrix> cov,
      ParticleHypothesis particleHypothesis,
      double tolerance = s_onSurfaceTolerance) {
    Result<BoundVector> bound =
        transformFreeToBoundParameters(pos4.segment<3>(ePos0), pos4[eTime], dir,
                                       qOverP, *surface, geoCtx, tolerance);

    if (!bound.ok()) {
      return bound.error();
    }

    return GenericBoundTrackParameters{std::move(surface), std::move(*bound),
                                       std::move(cov),
                                       std::move(particleHypothesis)};
  }

  /// Construct from four-position, direction, and qOverP.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored.
  /// @param qOverP Charge over momentum
  /// @param cov Curvilinear bound parameters covariance matrix
  /// @param particleHypothesis Particle hypothesis
  /// @return Curvilinear bound track parameters
  static GenericBoundTrackParameters createCurvilinear(
      const Vector4& pos4, const Vector3& dir, double qOverP,
      std::optional<CovarianceMatrix> cov,
      ParticleHypothesis particleHypothesis) {
    return GenericBoundTrackParameters(
        CurvilinearSurface(pos4.segment<3>(ePos0), dir).surface(),
        transformFreeToCurvilinearParameters(pos4[eTime], dir, qOverP),
        std::move(cov), std::move(particleHypothesis));
  }

  /// Construct from four-position, angles, and qOverP.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param phi Transverse track direction angle
  /// @param theta Longitudinal track direction angle
  /// @param qOverP Charge over momentum
  /// @param cov Curvilinear bound parameters covariance matrix
  /// @param particleHypothesis Particle hypothesis
  /// @return Curvilinear bound track parameters
  static GenericBoundTrackParameters createCurvilinear(
      const Vector4& pos4, double phi, double theta, double qOverP,
      std::optional<CovarianceMatrix> cov,
      ParticleHypothesis particleHypothesis) {
    return GenericBoundTrackParameters(
        CurvilinearSurface(pos4.segment<3>(ePos0),
                           makeDirectionFromPhiTheta(phi, theta))
            .surface(),
        transformFreeToCurvilinearParameters(pos4[eTime], phi, theta, qOverP),
        std::move(cov), std::move(particleHypothesis));
  }

  /// Construct from a parameters vector on the surface and particle charge.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param params Bound parameters vector
  /// @param cov Bound parameters covariance matrix
  /// @param particleHypothesis Particle hypothesis
  ///
  /// In principle, only the charge magnitude is needed her to allow unambiguous
  /// extraction of the absolute momentum. The particle charge is required as
  /// an input here to be consistent with the other constructors below that
  /// that also take the charge as an input. The charge sign is only used in
  /// debug builds to check for consistency with the q/p parameter.
  GenericBoundTrackParameters(std::shared_ptr<const Surface> surface,
                              const ParametersVector& params,
                              std::optional<CovarianceMatrix> cov,
                              ParticleHypothesis particleHypothesis)
      : m_params(params),
        m_cov(std::move(cov)),
        m_surface(std::move(surface)),
        m_particleHypothesis(std::move(particleHypothesis)) {
    // TODO set `validateAngleRange` to `true` after fixing caller code
    assert(isBoundVectorValid(m_params, false) &&
           "Invalid bound parameters vector");
    assert(m_surface != nullptr && "Reference surface must not be null");
    normalizePhiTheta();
  }

  /// Converts a bound track parameter with a different hypothesis.
  /// @param other The other bound track parameters to convert from
  template <typename other_particle_hypothesis_t>
  explicit GenericBoundTrackParameters(
      const GenericBoundTrackParameters<other_particle_hypothesis_t>& other)
      : GenericBoundTrackParameters(
            other.referenceSurface().getSharedPtr(), other.parameters(),
            other.covariance(),
            ParticleHypothesis{other.particleHypothesis()}) {}

  /// Convert this track parameter object to the general type-erased one
  /// @return Type-erased bound track parameters
  GenericBoundTrackParameters<Acts::ParticleHypothesis> toBound() const {
    return GenericBoundTrackParameters<Acts::ParticleHypothesis>{*this};
  }

  /// Parameters vector.
  /// @return Mutable reference to the parameters vector
  ParametersVector& parameters() { return m_params; }
  /// Parameters vector.
  /// @return Const reference to the parameters vector
  const ParametersVector& parameters() const { return m_params; }
  /// Vector of spatial impact parameters (i.e., d0 and z0)
  /// @return Two-dimensional vector of spatial impact parameters
  Vector2 spatialImpactParameters() const { return m_params.head<2>(); }
  /// Vector of spatial and temporal impact parameters (i.e., d0, z0, and t)
  /// @return Three-dimensional vector of impact parameters
  Vector3 impactParameters() const {
    Vector3 ip;
    ip.template head<2>() = m_params.template head<2>();
    ip(2) = m_params(eBoundTime);
    return ip;
  }

  /// Optional covariance matrix.
  /// @return Mutable reference to the optional covariance matrix
  std::optional<CovarianceMatrix>& covariance() { return m_cov; }
  /// Optional covariance matrix.
  /// @return Const reference to the optional covariance matrix
  const std::optional<CovarianceMatrix>& covariance() const { return m_cov; }
  /// Covariance matrix of the spatial impact parameters (i.e., of d0 and z0)
  /// @return Optional 2x2 covariance matrix of spatial impact parameters
  std::optional<SquareMatrix<2>> spatialImpactParameterCovariance() const {
    if (!m_cov.has_value()) {
      return std::nullopt;
    }

    return m_cov.value().template topLeftCorner<2, 2>();
  }

  /// Covariance matrix of the spatial and temporal impact parameters (i.e., of
  /// d0, z0, and t)
  /// @return Optional 3x3 covariance matrix of impact parameters
  std::optional<SquareMatrix<3>> impactParameterCovariance() const {
    if (!m_cov.has_value()) {
      return std::nullopt;
    }

    SquareMatrix<3> ipCov;
    ipCov.template topLeftCorner<2, 2>() =
        m_cov.value().template topLeftCorner<2, 2>();
    ipCov.template block<2, 1>(0, 2) =
        m_cov.value().template block<2, 1>(0, eBoundTime);
    ipCov.template block<1, 2>(2, 0) =
        m_cov.value().template block<1, 2>(eBoundTime, 0);
    ipCov(2, 2) = m_cov.value()(eBoundTime, eBoundTime);
    return ipCov;
  }

  /// Access a single parameter value identified by its index.
  ///
  /// @tparam kIndex Track parameter index
  /// @return The parameter value at the specified index
  template <BoundIndices kIndex>
  double get() const {
    return m_params[kIndex];
  }

  /// Local spatial position two-vector.
  /// @return Two-dimensional local position vector on the reference surface
  Vector2 localPosition() const { return m_params.segment<2>(eBoundLoc0); }
  /// Space-time position four-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
  ///
  /// @return Four-dimensional position vector (x, y, z, t) in global coordinates
  ///
  /// This uses the associated surface to transform the local position on
  /// the surface to globalcoordinates. This requires a geometry context to
  /// select the appropriate transformation and might be a computationally
  /// expensive operation.
  Vector4 fourPosition(const GeometryContext& geoCtx) const {
    Vector4 pos4;
    pos4.segment<3>(ePos0) =
        m_surface->localToGlobal(geoCtx, localPosition(), direction());
    pos4[eTime] = m_params[eBoundTime];
    return pos4;
  }
  /// Spatial position three-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
  ///
  /// @return Three-dimensional position vector in global coordinates
  ///
  /// This uses the associated surface to transform the local position on
  /// the surface to globalcoordinates. This requires a geometry context to
  /// select the appropriate transformation and might be a computationally
  /// expensive operation.
  Vector3 position(const GeometryContext& geoCtx) const {
    return m_surface->localToGlobal(geoCtx, localPosition(), direction());
  }
  /// Time coordinate.
  /// @return The time coordinate value
  double time() const { return m_params[eBoundTime]; }

  /// Phi direction.
  /// @return The azimuthal angle phi in radians
  double phi() const { return m_params[eBoundPhi]; }
  /// Theta direction.
  /// @return The polar angle theta in radians
  double theta() const { return m_params[eBoundTheta]; }
  /// Charge over momentum.
  /// @return The charge over momentum ratio
  double qOverP() const { return m_params[eBoundQOverP]; }

  /// Unit direction three-vector, i.e. the normalized momentum
  /// three-vector.
  /// @return Normalized direction vector
  Vector3 direction() const {
    return makeDirectionFromPhiTheta(m_params[eBoundPhi],
                                     m_params[eBoundTheta]);
  }
  /// Absolute momentum.
  /// @return The absolute momentum magnitude
  double absoluteMomentum() const {
    return m_particleHypothesis.extractMomentum(m_params[eBoundQOverP]);
  }
  /// Transverse momentum.
  /// @return The transverse momentum magnitude
  double transverseMomentum() const {
    return std::sin(m_params[eBoundTheta]) * absoluteMomentum();
  }
  /// Momentum three-vector.
  /// @return Three-dimensional momentum vector
  Vector3 momentum() const { return absoluteMomentum() * direction(); }

  /// Particle electric charge.
  /// @return The particle electric charge
  double charge() const {
    return m_particleHypothesis.extractCharge(get<eBoundQOverP>());
  }

  /// Particle hypothesis.
  /// @return Reference to the particle hypothesis
  const ParticleHypothesis& particleHypothesis() const {
    return m_particleHypothesis;
  }

  /// Reference surface onto which the parameters are bound.
  /// @return Reference to the bound reference surface
  const Surface& referenceSurface() const { return *m_surface; }
  /// Reference frame in which the local error is defined.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
  ///
  /// @return The reference frame rotation matrix
  ///
  /// For planar surfaces, this is the transformation local-to-global
  /// rotation matrix. For non-planar surfaces, it is the local-to-global
  /// rotation matrix of the tangential plane at the track position.
  RotationMatrix3 referenceFrame(const GeometryContext& geoCtx) const {
    return m_surface->referenceFrame(geoCtx, position(geoCtx), momentum());
  }

  /// Reflect the parameters in place.
  void reflectInPlace() { m_params = reflectBoundParameters(m_params); }

  /// Reflect the parameters.
  /// @return Reflected parameters.
  GenericBoundTrackParameters<ParticleHypothesis> reflect() const {
    GenericBoundTrackParameters<ParticleHypothesis> reflected = *this;
    reflected.reflectInPlace();
    return reflected;
  }

 private:
  BoundVector m_params;
  std::optional<BoundMatrix> m_cov;
  /// reference surface
  std::shared_ptr<const Surface> m_surface;
  // TODO use [[no_unique_address]] once we switch to C++20
  ParticleHypothesis m_particleHypothesis;

  /// Ensure phi and theta angles are within bounds.
  void normalizePhiTheta() {
    auto [phi, theta] =
        detail::normalizePhiTheta(m_params[eBoundPhi], m_params[eBoundTheta]);
    m_params[eBoundPhi] = phi;
    m_params[eBoundTheta] = theta;
  }

  /// Compare two bound track parameters for bitwise equality.
  ///
  /// @note Comparing track parameters for bitwise equality is not a good
  /// idea.
  ///   Depending on the context you might want to compare only the
  ///   parameter values, or compare them for compatibility instead of
  ///   equality; you might also have different (floating point) thresholds
  ///   of equality in different contexts. None of that can be handled by
  ///   this operator. Users should think really hard if this is what they
  ///   want and we might decided that we will remove this in the future.
  friend bool operator==(const GenericBoundTrackParameters& lhs,
                         const GenericBoundTrackParameters& rhs) {
    return (lhs.m_params == rhs.m_params) && (lhs.m_cov == rhs.m_cov) &&
           (lhs.m_surface == rhs.m_surface) &&
           (lhs.m_particleHypothesis == rhs.m_particleHypothesis);
  }

  /// Print information to the output stream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const GenericBoundTrackParameters& tp) {
    detail::printBoundParameters(
        os, tp.referenceSurface(), tp.particleHypothesis(), tp.parameters(),
        tp.covariance().has_value() ? &tp.covariance().value() : nullptr);
    return os;
  }
};

}  // namespace Acts
