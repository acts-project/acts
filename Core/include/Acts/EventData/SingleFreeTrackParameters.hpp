// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/EventData/detail/PrintParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <cassert>
#include <cmath>
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
  using Scalar = FreeScalar;
  using ParametersVector = FreeVector;
  using CovarianceMatrix = FreeSymMatrix;

  /// Construct from a parameters vector and and particle charge.
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
      : m_paramSet(std::move(cov), params), m_chargeInterpreter(std::abs(q)) {
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
      : m_paramSet(std::move(cov), params), m_chargeInterpreter(T()) {}

  /// Construct from four-position, angles, absolute momentum, and charge.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param phi Transverse track direction angle
  /// @param theta Longitudinal track direction angle
  /// @param p Absolute momentum
  /// @param q Particle charge
  /// @param cov Free parameters covariance matrix
  SingleFreeTrackParameters(const Vector4D& pos4, Scalar phi, Scalar theta,
                            Scalar p, Scalar q,
                            std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_paramSet(std::move(cov), ParametersVector::Zero()),
        m_chargeInterpreter(std::abs(q)) {
    assert((0 <= p) and "Absolute momentum must be positive");

    m_paramSet.setParameter<eFreePos0>(pos4[ePos0]);
    m_paramSet.setParameter<eFreePos1>(pos4[ePos1]);
    m_paramSet.setParameter<eFreePos2>(pos4[ePos2]);
    m_paramSet.setParameter<eFreeTime>(pos4[eTime]);
    auto dir = makeDirectionUnitFromPhiTheta(phi, theta);
    m_paramSet.setParameter<eFreeDir0>(dir[eMom0]);
    m_paramSet.setParameter<eFreeDir1>(dir[eMom1]);
    m_paramSet.setParameter<eFreeDir2>(dir[eMom2]);
    m_paramSet.setParameter<eFreeQOverP>((q != Scalar(0)) ? (q / p) : (1 / p));
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
  SingleFreeTrackParameters(const Vector4D& pos4, Scalar phi, Scalar theta,
                            Scalar qOverP,
                            std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_paramSet(std::move(cov), ParametersVector::Zero()),
        m_chargeInterpreter(T()) {
    m_paramSet.setParameter<eFreePos0>(pos4[ePos0]);
    m_paramSet.setParameter<eFreePos1>(pos4[ePos1]);
    m_paramSet.setParameter<eFreePos2>(pos4[ePos2]);
    m_paramSet.setParameter<eFreeTime>(pos4[eTime]);
    auto dir = makeDirectionUnitFromPhiTheta(phi, theta);
    m_paramSet.setParameter<eFreeDir0>(dir[eMom0]);
    m_paramSet.setParameter<eFreeDir1>(dir[eMom1]);
    m_paramSet.setParameter<eFreeDir2>(dir[eMom2]);
    m_paramSet.setParameter<eFreeQOverP>(qOverP);
  }

  // this class does not have a custom default constructor and thus should not
  // provide any custom default cstors, dstor, or assignment. see ISOCPP C.20.

  /// Parameters vector.
  ParametersVector parameters() const { return m_paramSet.getParameters(); }
  /// Optional covariance matrix.
  const std::optional<CovarianceMatrix>& covariance() const {
    return m_paramSet.getCovariance();
  }

  /// Access a single parameter value indentified by its index.
  ///
  /// @tparam kIndex Track parameter index
  template <FreeIndices kIndex>
  Scalar get() const {
    return m_paramSet.template getParameter<kIndex>();
  }
  /// Access a single parameter uncertainty identified by its index.
  ///
  /// @tparam kIndex Track parameter index
  /// @retval zero if the track parameters have no associated covariance
  /// @retval parameter standard deviation if the covariance is available
  template <FreeIndices kIndex>
  Scalar uncertainty() const {
    return m_paramSet.template getUncertainty<kIndex>();
  }

  /// Space-time position four-vector.
  Vector4D fourPosition() const {
    Vector4D pos4;
    pos4[ePos0] = get<eFreePos0>();
    pos4[ePos1] = get<eFreePos1>();
    pos4[ePos2] = get<eFreePos2>();
    pos4[eTime] = get<eFreeTime>();
    return pos4;
  }
  /// Spatial position three-vector.
  Vector3D position() const {
    return parameters().template segment<3>(eFreePos0);
  }
  /// Time coordinate.
  Scalar time() const { return get<eFreeTime>(); }

  /// Unit direction three-vector, i.e. the normalized momentum three-vector.
  Vector3D unitDirection() const {
    return parameters().template segment<3>(eFreeDir0).normalized();
  }
  /// Absolute momentum.
  Scalar absoluteMomentum() const {
    return m_chargeInterpreter.extractMomentum(get<eFreeQOverP>());
  }
  /// Transverse momentum.
  Scalar transverseMomentum() const {
    // direction vector w/ arbitrary normalization can be parametrized as
    //   [f*sin(theta)*cos(phi), f*sin(theta)*sin(phi), f*cos(theta)]
    // w/ f,sin(theta) positive, the transverse magnitude is then
    //   sqrt(f^2*sin^2(theta)) = f*sin(theta)
    Scalar transverseMagnitude = std::hypot(get<eFreeDir0>(), get<eFreeDir1>());
    // absolute magnitude is f by construction
    Scalar magnitude = std::hypot(transverseMagnitude, get<eFreeDir2>());
    // such that we can extract sin(theta) = f*sin(theta) / f
    return (transverseMagnitude / magnitude) * absoluteMomentum();
  }
  /// Momentum three-vector.
  Vector3D momentum() const { return absoluteMomentum() * unitDirection(); }

  /// Particle electric charge.
  constexpr Scalar charge() const {
    return m_chargeInterpreter.extractCharge(get<eFreeQOverP>());
  }

 private:
  FullFreeParameterSet m_paramSet;
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
