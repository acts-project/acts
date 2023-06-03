// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

namespace Acts {

/// Curvilinear track parameters for a single track.
///
/// This is intended as a user-facing data class that adds additional accessors
/// and charge/momentum interpretation on top of the pure parameters vector. All
/// parameters and their corresponding covariance matrix are stored in
/// curvilinear parametrization.
///
/// @see BoundTrackParameters
class CurvilinearTrackParameters : public BoundTrackParameters {
  using Base = BoundTrackParameters;

 public:
  using Scalar = ActsScalar;
  using ParametersVector = BoundVector;
  using CovarianceMatrix = BoundSymMatrix;

  /// Construct from four-position, direction, and charge-over-momentum.
  ///
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored.
  /// @param qOverP Charge-over-momentum-like parameter
  /// @param cov Curvilinear bound parameters covariance matrix
  CurvilinearTrackParameters(const Vector4& pos4, const Vector3& dir,
                             Scalar qOverP,
                             std::optional<CovarianceMatrix> cov = std::nullopt)
      : Base(Surface::makeShared<PlaneSurface>(pos4.segment<3>(ePos0), dir),
             detail::transformFreeToCurvilinearParameters(pos4[eTime], dir,
                                                          qOverP),
             std::move(cov)) {}

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
  CurvilinearTrackParameters(const Vector4& pos4, Scalar phi, Scalar theta,
                             Scalar qOverP,
                             std::optional<CovarianceMatrix> cov = std::nullopt)
      : Base(Surface::makeShared<PlaneSurface>(
                 pos4.segment<3>(ePos0),
                 makeDirectionUnitFromPhiTheta(phi, theta)),
             detail::transformFreeToCurvilinearParameters(pos4[eTime], phi,
                                                          theta, qOverP),
             std::move(cov)) {}

  /// Parameters are not default constructible due to the charge type.
  CurvilinearTrackParameters() = delete;
  CurvilinearTrackParameters(const CurvilinearTrackParameters&) = default;
  CurvilinearTrackParameters(CurvilinearTrackParameters&&) = default;
  ~CurvilinearTrackParameters() = default;
  CurvilinearTrackParameters& operator=(const CurvilinearTrackParameters&) =
      default;
  CurvilinearTrackParameters& operator=(CurvilinearTrackParameters&&) = default;
};

}  // namespace Acts
