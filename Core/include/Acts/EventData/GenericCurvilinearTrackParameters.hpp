// This file is part of the Acts project.
//
// Copyright (C) 2016-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParametersConcept.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

namespace Acts {

/// Curvilinear track parameters for a single track.
///
/// @tparam particle_hypothesis_t Helper type to interpret the particle charge/momentum
///
/// This is intended as a user-facing data class that adds additional accessors
/// and charge/momentum interpretation on top of the pure parameters vector. All
/// parameters and their corresponding covariance matrix are stored in
/// curvilinear parametrization.
///
/// @see GenericBoundTrackParameters
template <typename particle_hypothesis_t>
class GenericCurvilinearTrackParameters {
 public:
  using Scalar = ActsScalar;
  using ParametersVector = BoundVector;
  using CovarianceMatrix = BoundSquareMatrix;
  using ParticleHypothesis = particle_hypothesis_t;

  /// Construct from four-position, direction, and qOverP.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored.
  /// @param qOverP Charge over momentum
  /// @param cov Curvilinear bound parameters covariance matrix
  /// @param particleHypothesis Particle hypothesis
  GenericCurvilinearTrackParameters(const Vector4& pos4, const Vector3& dir,
                                    Scalar qOverP,
                                    std::optional<CovarianceMatrix> cov,
                                    ParticleHypothesis particleHypothesis)
      : m_surface(pos4.segment<3>(ePos0), dir),
        m_params{detail::transformFreeToCurvilinearParameters(pos4[eTime], dir,
                                                              qOverP)},
        m_cov{std::move(cov)},
        m_particleHypothesis{std::move(particleHypothesis)} {}

  /// Construct from four-position, angles, and qOverP.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param phi Transverse track direction angle
  /// @param theta Longitudinal track direction angle
  /// @param qOverP Charge over momentum
  /// @param cov Curvilinear bound parameters covariance matrix
  /// @param particleHypothesis Particle hypothesis
  GenericCurvilinearTrackParameters(const Vector4& pos4, Scalar phi,
                                    Scalar theta, Scalar qOverP,
                                    std::optional<CovarianceMatrix> cov,
                                    ParticleHypothesis particleHypothesis)
      : m_surface(pos4.segment<3>(ePos0),
                  makeDirectionFromPhiTheta(phi, theta)),
        m_params{detail::transformFreeToCurvilinearParameters(pos4[eTime], phi,
                                                              theta, qOverP)},
        m_cov{std::move(cov)},
        m_particleHypothesis{std::move(particleHypothesis)} {}

  /// Converts a bound track parameter with a different hypothesis.
  template <typename other_particle_hypothesis_t>
  GenericCurvilinearTrackParameters(
      const GenericCurvilinearTrackParameters<other_particle_hypothesis_t>&
          other)
      : GenericCurvilinearTrackParameters(other.fourPosition(),
                                          other.particleHypothesis(),
                                          other.covariance()) {}

  /// Converts an unknown bound track parameter.
  template <typename other_track_parameter_t>
  static GenericCurvilinearTrackParameters create(
      const other_track_parameter_t& other) {
    static_assert(
        Concepts::BoundTrackParametersConcept<other_track_parameter_t>);

    return GenericCurvilinearTrackParameters(
        other.fourPosition(), other.particleHypothesis(), other.covariance());
  }

  /// Parameters are not default constructible due to the charge type.
  GenericCurvilinearTrackParameters() = delete;
  GenericCurvilinearTrackParameters(const GenericCurvilinearTrackParameters&) =
      default;
  GenericCurvilinearTrackParameters(GenericCurvilinearTrackParameters&&) =
      default;
  ~GenericCurvilinearTrackParameters() = default;
  GenericCurvilinearTrackParameters& operator=(
      const GenericCurvilinearTrackParameters&) = default;
  GenericCurvilinearTrackParameters& operator=(
      GenericCurvilinearTrackParameters&&) = default;

  /// Parameters vector.
  ParametersVector& parameters() { return m_params; }
  /// Parameters vector.
  const ParametersVector& parameters() const { return m_params; }
  /// Vector of spatial impact parameters (i.e., d0 and z0)
  ActsVector<2> spatialImpactParameters() const { return m_params.head<2>(); }
  /// Vector of spatial and temporal impact parameters (i.e., d0, z0, and t)
  ActsVector<3> impactParameters() const {
    ActsVector<3> ip;
    ip.template head<2>() = m_params.template head<2>();
    ip(2) = m_params(eBoundTime);
    return ip;
  }

  /// Optional covariance matrix.
  std::optional<CovarianceMatrix>& covariance() { return m_cov; }
  /// Optional covariance matrix.
  const std::optional<CovarianceMatrix>& covariance() const { return m_cov; }
  /// Covariance matrix of the spatial impact parameters (i.e., of d0 and z0)
  std::optional<ActsSquareMatrix<2>> spatialImpactParameterCovariance() const {
    if (!m_cov.has_value()) {
      return std::nullopt;
    }

    return m_cov.value().template topLeftCorner<2, 2>();
  }

  /// Covariance matrix of the spatial and temporal impact parameters (i.e., of
  /// d0, z0, and t)
  std::optional<ActsSquareMatrix<3>> impactParameterCovariance() const {
    if (!m_cov.has_value()) {
      return std::nullopt;
    }

    ActsSquareMatrix<3> ipCov;
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
  template <BoundIndices kIndex>
  Scalar get() const {
    return m_params[kIndex];
  }

  /// Local spatial position two-vector.
  Vector2 localPosition() const { return m_params.segment<2>(eBoundLoc0); }
  /// Space-time position four-vector.
  ///
  /// This uses the associated surface to transform the local position on
  /// the surface to globalcoordinates.
  Vector4 fourPosition() const {
    Vector4 pos4;
    pos4.segment<3>(ePos0) = m_surface.center();
    pos4[eTime] = m_params[eBoundTime];
    return pos4;
  }
  /// Spatial position three-vector.
  ///
  /// This uses the associated surface to transform the local position on
  /// the surface to globalcoordinates.
  Vector3 position() const { return m_surface.center(); }
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

  GenericBoundTrackParameters<ParticleHypothesis> boundParameters() const {
    return GenericBoundTrackParameters<ParticleHypothesis>(
        m_surface.planeSurface(), m_params, m_cov, m_particleHypothesis);
  }

 private:
  /// reference surface
  CurvilinearSurface m_surface;
  BoundVector m_params;
  std::optional<BoundSquareMatrix> m_cov;
  // TODO use [[no_unique_address]] once we switch to C++20
  ParticleHypothesis m_particleHypothesis;
};

}  // namespace Acts
