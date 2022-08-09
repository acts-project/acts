// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/detail/PrintParameters.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <cassert>
#include <cmath>
#include <optional>
#include <type_traits>

namespace Acts {

/// Track parameters not bound to a surface for a single track.
///
/// @tparam charge_t Helper type to interpret the particle charge
///
/// Parameters and covariance matrix are stored using the free parametrization
/// defined in `enum FreeIndices`.
template <class charge_t>
class SingleFreeTrackParameters {
 public:
  using Scalar = ActsScalar;
  using ParametersVector = FreeVector;
  using CovarianceMatrix = FreeSymMatrix;

  /// Construct from a parameters vector and particle charge.
  ///
  /// @param params Free parameters vector
  /// @param q Particle charge
  /// @param cov Free parameters covariance matrix
  ///
  /// In principle, only the charge magnitude is needed her to allow unambigous
  /// extraction of the absolute momentum. The particle charge is required as
  /// an input here to be consistent with the other constructors below that
  /// that also take the charge as an input. The charge sign is only used in
  /// debug builds to check for consistency with the q/p parameter.
  SingleFreeTrackParameters(const ParametersVector& params, Scalar q,
                            std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_params(params),
        m_cov(std::move(cov)),
        m_chargeInterpreter(std::abs(q)) {
    assert((0 <= (params[eFreeQOverP] * q)) and "Inconsistent q/p and q signs");
  }

  /// Construct from a parameters vector.
  ///
  /// @param params Free parameters vector
  /// @param cov Free parameters covariance matrix
  /// @tparam T Internal helper template be able to check charge type
  ///
  /// This constructor is only available if there are no potential charge
  /// ambiguities, i.e. the charge interpretation type is default-constructible.
  template <typename T = charge_t,
            std::enable_if_t<std::is_default_constructible_v<T>, int> = 0>
  SingleFreeTrackParameters(const ParametersVector& params,
                            std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_params(params), m_cov(std::move(cov)) {}

  /// Construct from four-position, angles, absolute momentum, and charge.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param phi Transverse track direction angle
  /// @param theta Longitudinal track direction angle
  /// @param p Absolute momentum
  /// @param q Particle charge
  /// @param cov Free parameters covariance matrix
  SingleFreeTrackParameters(const Vector4& pos4, Scalar phi, Scalar theta,
                            Scalar p, Scalar q,
                            std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_params(FreeVector::Zero()),
        m_cov(std::move(cov)),
        m_chargeInterpreter(std::abs(q)) {
    assert((0 <= p) and "Absolute momentum must be positive");

    auto dir = makeDirectionUnitFromPhiTheta(phi, theta);
    m_params[eFreePos0] = pos4[ePos0];
    m_params[eFreePos1] = pos4[ePos1];
    m_params[eFreePos2] = pos4[ePos2];
    m_params[eFreeTime] = pos4[eTime];
    m_params[eFreeDir0] = dir[eMom0];
    m_params[eFreeDir1] = dir[eMom1];
    m_params[eFreeDir2] = dir[eMom2];
    m_params[eFreeQOverP] = (q != Scalar(0)) ? (q / p) : (1 / p);
  }

  /// Construct from four-position, angles, and charge-over-momentum.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param phi Transverse track direction angle
  /// @param theta Longitudinal track direction angle
  /// @param qOverP Charge-over-momentum-like parameter
  /// @param cov Free parameters covariance matrix
  ///
  /// This constructor is only available if there are no potential charge
  /// ambiguities, i.e. the charge interpretation type is default-constructible.
  template <typename T = charge_t,
            std::enable_if_t<std::is_default_constructible_v<T>, int> = 0>
  SingleFreeTrackParameters(const Vector4& pos4, Scalar phi, Scalar theta,
                            Scalar qOverP,
                            std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_params(FreeVector::Zero()), m_cov(std::move(cov)) {
    auto dir = makeDirectionUnitFromPhiTheta(phi, theta);
    m_params[eFreePos0] = pos4[ePos0];
    m_params[eFreePos1] = pos4[ePos1];
    m_params[eFreePos2] = pos4[ePos2];
    m_params[eFreeTime] = pos4[eTime];
    m_params[eFreeDir0] = dir[eMom0];
    m_params[eFreeDir1] = dir[eMom1];
    m_params[eFreeDir2] = dir[eMom2];
    m_params[eFreeQOverP] = qOverP;
  }

  /// Parameters are not default constructible due to the charge type.
  SingleFreeTrackParameters() = delete;
  SingleFreeTrackParameters(const SingleFreeTrackParameters&) = default;
  SingleFreeTrackParameters(SingleFreeTrackParameters&&) = default;
  ~SingleFreeTrackParameters() = default;
  SingleFreeTrackParameters& operator=(const SingleFreeTrackParameters&) =
      default;
  SingleFreeTrackParameters& operator=(SingleFreeTrackParameters&&) = default;

  /// Parameters vector.
  const ParametersVector& parameters() const { return m_params; }
  /// Optional covariance matrix.
  const std::optional<CovarianceMatrix>& covariance() const { return m_cov; }

  /// Access a single parameter value indentified by its index.
  ///
  /// @tparam kIndex Track parameter index
  template <FreeIndices kIndex>
  Scalar get() const {
    return m_params[kIndex];
  }

  /// Space-time position four-vector.
  Vector4 fourPosition() const {
    Vector4 pos4;
    pos4[ePos0] = m_params[eFreePos0];
    pos4[ePos1] = m_params[eFreePos1];
    pos4[ePos2] = m_params[eFreePos2];
    pos4[eTime] = m_params[eFreeTime];
    return pos4;
  }
  /// Spatial position three-vector.
  Vector3 position() const { return m_params.segment<3>(eFreePos0); }
  /// Time coordinate.
  Scalar time() const { return m_params[eFreeTime]; }

  /// Unit direction three-vector, i.e. the normalized momentum three-vector.
  Vector3 unitDirection() const {
    return m_params.segment<3>(eFreeDir0).normalized();
  }
  /// Absolute momentum.
  Scalar absoluteMomentum() const {
    return m_chargeInterpreter.extractMomentum(m_params[eFreeQOverP]);
  }
  /// Transverse momentum.
  Scalar transverseMomentum() const {
    // direction vector w/ arbitrary normalization can be parametrized as
    //   [f*sin(theta)*cos(phi), f*sin(theta)*sin(phi), f*cos(theta)]
    // w/ f,sin(theta) positive, the transverse magnitude is then
    //   sqrt(f^2*sin^2(theta)) = f*sin(theta)
    Scalar transverseMagnitude =
        std::hypot(m_params[eFreeDir0], m_params[eFreeDir1]);
    // absolute magnitude is f by construction
    Scalar magnitude = std::hypot(transverseMagnitude, m_params[eFreeDir2]);
    // such that we can extract sin(theta) = f*sin(theta) / f
    return (transverseMagnitude / magnitude) * absoluteMomentum();
  }
  /// Momentum three-vector.
  Vector3 momentum() const { return absoluteMomentum() * unitDirection(); }

  /// Particle electric charge.
  Scalar charge() const {
    return m_chargeInterpreter.extractCharge(get<eFreeQOverP>());
  }

 private:
  FreeVector m_params;
  std::optional<FreeSymMatrix> m_cov;
  // TODO use [[no_unique_address]] once we switch to C++20
  charge_t m_chargeInterpreter;

  /// Print information to the output stream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const SingleFreeTrackParameters& tp) {
    detail::printFreeParameters(
        os, tp.parameters(),
        tp.covariance().has_value() ? &tp.covariance().value() : nullptr);
    return os;
  }
};

}  // namespace Acts
