// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/detail/PrintParameters.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <type_traits>

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
  using Scalar = ActsScalar;
  using ParametersVector = BoundVector;
  using CovarianceMatrix = BoundSquareMatrix;
  using ParticleHypothesis = particle_hypothesis_t;

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
    assert(m_surface);
    normalizePhiTheta();
  }

  /// Converts a bound track parameter with a different hypothesis.
  template <typename other_particle_hypothesis_t>
  GenericBoundTrackParameters(
      const GenericBoundTrackParameters<other_particle_hypothesis_t>& other)
      : GenericBoundTrackParameters(other.referenceSurface().getSharedPtr(),
                                    other.parameters(), other.covariance(),
                                    other.particleHypothesis()) {}

  /// Factory to construct from four-position, direction, absolute momentum, and
  /// charge.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param geoCtx Geometry context for the local-to-global transformation
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored
  /// @param qOverP Charge over momentum
  /// @param cov Bound parameters covariance matrix
  /// @param particleHypothesis Particle hypothesis
  ///
  /// @note The returned result indicates whether the free parameters could
  /// successfully be converted to on-surface parameters.
  static Result<GenericBoundTrackParameters> create(
      std::shared_ptr<const Surface> surface, const GeometryContext& geoCtx,
      const Vector4& pos4, const Vector3& dir, Scalar qOverP,
      std::optional<CovarianceMatrix> cov,
      ParticleHypothesis particleHypothesis) {
    Result<BoundVector> bound = detail::transformFreeToBoundParameters(
        pos4.segment<3>(ePos0), pos4[eTime], dir, qOverP, *surface, geoCtx);

    if (!bound.ok()) {
      return bound.error();
    }

    return GenericBoundTrackParameters{std::move(surface), std::move(*bound),
                                       std::move(cov),
                                       std::move(particleHypothesis)};
  }

  /// Parameters are not default constructible due to the charge type.
  GenericBoundTrackParameters() = delete;
  GenericBoundTrackParameters(const GenericBoundTrackParameters&) = default;
  GenericBoundTrackParameters(GenericBoundTrackParameters&&) = default;
  ~GenericBoundTrackParameters() = default;
  GenericBoundTrackParameters& operator=(const GenericBoundTrackParameters&) =
      default;
  GenericBoundTrackParameters& operator=(GenericBoundTrackParameters&&) =
      default;

  /// Parameters vector.
  ParametersVector& parameters() { return m_params; }
  /// Parameters vector.
  const ParametersVector& parameters() const { return m_params; }
  /// Optional covariance matrix.
  std::optional<CovarianceMatrix>& covariance() { return m_cov; }
  /// Optional covariance matrix.
  const std::optional<CovarianceMatrix>& covariance() const { return m_cov; }

  /// Access a single parameter value identified by its index.
  ///
  /// @tparam kIndex Track parameter index
  template <BoundIndices kIndex>
  Scalar get() const {
    return m_params[kIndex];
  }

  /// Local spatial position two-vector.
  Vector2 localPosition() const { return m_params.segment<2>(eBoundLoc0); }
  /// Space-time position four-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
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
  /// This uses the associated surface to transform the local position on
  /// the surface to globalcoordinates. This requires a geometry context to
  /// select the appropriate transformation and might be a computationally
  /// expensive operation.
  Vector3 position(const GeometryContext& geoCtx) const {
    return m_surface->localToGlobal(geoCtx, localPosition(), direction());
  }
  /// Time coordinate.
  Scalar time() const { return m_params[eBoundTime]; }

  /// Phi direction.
  Scalar phi() const { return m_params[eBoundPhi]; }
  /// Theta direction.
  Scalar theta() const { return m_params[eBoundTheta]; }
  /// Charge over momentum.
  Scalar qOverP() const { return m_params[eBoundQOverP]; }

  /// Unit direction three-vector, i.e. the normalized momentum
  /// three-vector.
  Vector3 direction() const {
    return makeDirectionFromPhiTheta(m_params[eBoundPhi],
                                     m_params[eBoundTheta]);
  }
  /// Absolute momentum.
  Scalar absoluteMomentum() const {
    return m_particleHypothesis.extractMomentum(m_params[eBoundQOverP]);
  }
  /// Transverse momentum.
  Scalar transverseMomentum() const {
    return std::sin(m_params[eBoundTheta]) * absoluteMomentum();
  }
  /// Momentum three-vector.
  Vector3 momentum() const { return absoluteMomentum() * direction(); }

  /// Particle electric charge.
  Scalar charge() const {
    return m_particleHypothesis.extractCharge(get<eBoundQOverP>());
  }

  /// Particle hypothesis.
  const ParticleHypothesis& particleHypothesis() const {
    return m_particleHypothesis;
  }

  /// Reference surface onto which the parameters are bound.
  const Surface& referenceSurface() const { return *m_surface; }
  /// Reference frame in which the local error is defined.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
  ///
  /// For planar surfaces, this is the transformation local-to-global
  /// rotation matrix. For non-planar surfaces, it is the local-to-global
  /// rotation matrix of the tangential plane at the track position.
  RotationMatrix3 referenceFrame(const GeometryContext& geoCtx) const {
    return m_surface->referenceFrame(geoCtx, position(geoCtx), momentum());
  }

 private:
  BoundVector m_params;
  std::optional<BoundSquareMatrix> m_cov;
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
    return (lhs.m_params == rhs.m_params) and (lhs.m_cov == rhs.m_cov) and
           (lhs.m_surface == rhs.m_surface) and
           (lhs.m_particleHypothesis == rhs.m_particleHypothesis);
  }
  /// Compare two bound track parameters for bitwise in-equality.
  friend bool operator!=(const GenericBoundTrackParameters& lhs,
                         const GenericBoundTrackParameters& rhs) {
    return not(lhs == rhs);
  }
  /// Print information to the output stream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const GenericBoundTrackParameters& tp) {
    detail::printBoundParameters(
        os, tp.referenceSurface(), tp.parameters(),
        tp.covariance().has_value() ? &tp.covariance().value() : nullptr);
    return os;
  }
};

}  // namespace Acts
