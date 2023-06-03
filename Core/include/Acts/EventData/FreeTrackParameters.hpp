// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
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
class FreeTrackParameters {
 public:
  using Scalar = ActsScalar;
  using ParametersVector = FreeVector;
  using CovarianceMatrix = FreeSymMatrix;

  /// Construct from a parameters vector.
  ///
  /// @param params Free parameters vector
  /// @param cov Free parameters covariance matrix
  FreeTrackParameters(const ParametersVector& params,
                      std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_params(params), m_cov(std::move(cov)) {}

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
  FreeTrackParameters(const Vector4& pos4, Scalar phi, Scalar theta,
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
  FreeTrackParameters() = delete;
  FreeTrackParameters(const FreeTrackParameters&) = default;
  FreeTrackParameters(FreeTrackParameters&&) = default;
  ~FreeTrackParameters() = default;
  FreeTrackParameters& operator=(const FreeTrackParameters&) = default;
  FreeTrackParameters& operator=(FreeTrackParameters&&) = default;

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
  Vector3 direction() const {
    return m_params.segment<3>(eFreeDir0).normalized();
  }

 private:
  FreeVector m_params;
  std::optional<FreeSymMatrix> m_cov;

  /// Print information to the output stream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const FreeTrackParameters& tp) {
    detail::printFreeParameters(
        os, tp.parameters(),
        tp.covariance().has_value() ? &tp.covariance().value() : nullptr);
    return os;
  }
};

}  // namespace Acts
