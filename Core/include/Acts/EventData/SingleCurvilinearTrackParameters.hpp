// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SingleBoundTrackParameters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

namespace Acts {

/// Curvilinear track parameters for a single track.
///
/// @tparam charge_t Helper type to interpret the particle charge/momentum
///
/// This is intended as a user-facing data class that adds additional accessors
/// and charge/momentum interpretation on top of the pure parameters vector. All
/// parameters and their corresponding covariance matrix are stored in
/// curvilinear parametrization.
///
/// @see SingleBoundTrackParameters
template <typename charge_t>
class SingleCurvilinearTrackParameters
    : public SingleBoundTrackParameters<charge_t> {
  using Base = SingleBoundTrackParameters<charge_t>;

 public:
  using Scalar = BoundScalar;
  using ParametersVector = BoundVector;
  using CovarianceMatrix = BoundSymMatrix;

  /// Construct from four-position, direction, absolute momentum, and charge.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored.
  /// @param p Absolute momentum
  /// @param q Particle charge
  /// @param cov Curvilinear bound parameters covariance matrix
  SingleCurvilinearTrackParameters(
      const Vector4D& pos4, const Vector3D& dir, Scalar p, Scalar q,
      std::optional<CovarianceMatrix> cov = std::nullopt)
      : Base(Surface::makeShared<PlaneSurface>(pos4.segment<3>(ePos0), dir),
             detail::transformFreeToCurvilinearParameters(
                 pos4[eTime], dir, (q != Scalar(0)) ? (q / p) : (1 / p)),
             q, std::move(cov)) {
    assert((0 <= p) and "Absolute momentum must be positive");
  }

  /// Construct from four-position, direction, and charge-over-momentum.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored.
  /// @param qOverP Charge-over-momentum-like parameter
  /// @param cov Curvilinear bound parameters covariance matrix
  ///
  /// This constructor is only available if there are no potential charge
  /// ambiguities, i.e. the charge interpretation type is default-constructible.
  template <typename T = charge_t,
            std::enable_if_t<std::is_default_constructible_v<T>, int> = 0>
  SingleCurvilinearTrackParameters(
      const Vector4D& pos4, const Vector3D& dir, Scalar qOverP,
      std::optional<CovarianceMatrix> cov = std::nullopt)
      : Base(Surface::makeShared<PlaneSurface>(pos4.segment<3>(ePos0), dir),
             detail::transformFreeToCurvilinearParameters(pos4[eTime], dir,
                                                          qOverP),
             std::move(cov)) {}

  /// Construct from four-position, angles, absolute momentum, and charge.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param phi Transverse track direction angle
  /// @param theta Longitudinal track direction angle
  /// @param p Absolute momentum
  /// @param q Particle charge
  /// @param cov Curvilinear bound parameters covariance matrix
  SingleCurvilinearTrackParameters(
      const Vector4D& pos4, Scalar phi, Scalar theta, Scalar p, Scalar q,
      std::optional<CovarianceMatrix> cov = std::nullopt)
      : Base(Surface::makeShared<PlaneSurface>(
                 pos4.segment<3>(ePos0),
                 makeDirectionUnitFromPhiTheta(phi, theta)),
             detail::transformFreeToCurvilinearParameters(
                 pos4[eTime], phi, theta, (q != Scalar(0)) ? (q / p) : (1 / p)),
             q, std::move(cov)) {
    assert((0 <= p) and "Absolute momentum must be positive");
  }

  /// Construct from four-position, angles, and charge-over-momentum.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param phi Transverse track direction angle
  /// @param theta Longitudinal track direction angle
  /// @param qOverP Charge-over-momentum-like parameter
  /// @param cov Curvilinear bound parameters covariance matrix
  ///
  /// This constructor is only available if there are no potential charge
  /// ambiguities, i.e. the charge interpretation type is default-constructible.
  template <typename T = charge_t,
            std::enable_if_t<std::is_default_constructible_v<T>, int> = 0>
  SingleCurvilinearTrackParameters(
      const Vector4D& pos4, Scalar phi, Scalar theta, Scalar qOverP,
      std::optional<CovarianceMatrix> cov = std::nullopt)
      : Base(Surface::makeShared<PlaneSurface>(
                 pos4.segment<3>(ePos0),
                 makeDirectionUnitFromPhiTheta(phi, theta)),
             detail::transformFreeToCurvilinearParameters(pos4[eTime], phi,
                                                          theta, qOverP),
             std::move(cov)) {}

  // this class does not have a custom default constructor and thus should not
  // provide any custom default cstors, dstor, or assignment. see ISOCPP C.20.
};

}  // namespace Acts
