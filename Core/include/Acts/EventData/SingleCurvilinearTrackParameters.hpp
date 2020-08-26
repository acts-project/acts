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

#include <memory>

namespace Acts {

/// Single track parameters bound to a curvilinear reference surface.
///
/// @tparam charge_policy_t Selection type for the particle charge
///
/// These are just regular bound parameters with a fixed surface type.
///
/// @see SingleBoundTrackParameters
template <typename charge_policy_t>
class SingleCurvilinearTrackParameters
    : public SingleBoundTrackParameters<charge_policy_t> {
  using Base = SingleBoundTrackParameters<charge_policy_t>;

 public:
  using Scalar = BoundParametersScalar;
  using ParametersVector = BoundVector;
  using CovarianceMatrix = BoundSymMatrix;

  /// Construct charged curvilinear parameters from global position, momentum.
  ///
  /// @param[in] cov Optional covariance matrix in the curvilinear frame
  /// @param[in] pos The global track three-position vector
  /// @param[in] mom The global track three-momentum vector
  /// @param[in] charge The particle charge
  /// @param[in] time The time coordinate
  template <typename T = charge_policy_t,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleCurvilinearTrackParameters(std::optional<CovarianceMatrix> cov,
                                   const Vector3D& pos, const Vector3D& mom,
                                   Scalar charge, Scalar time)
      : Base(Surface::makeShared<PlaneSurface>(pos, mom),
             detail::transformFreeToCurvilinearParameters(time, mom,
                                                          charge / mom.norm()),
             std::move(cov)) {}

  /// Construct neutral curvilinear parameters from global position, momentum.
  ///
  /// @param[in] cov Optional covariance matrix in the curvilinear frame
  /// @param[in] pos The global track three-position vector
  /// @param[in] mom The global track three-momentum vector
  /// @param[in] time The time coordinate
  template <typename T = charge_policy_t,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleCurvilinearTrackParameters(std::optional<CovarianceMatrix> cov,
                                   const Vector3D& pos, const Vector3D& mom,
                                   Scalar time)
      : Base(Surface::makeShared<PlaneSurface>(pos, mom),
             detail::transformFreeToCurvilinearParameters(time, mom,
                                                          1 / mom.norm()),
             std::move(cov)) {}

  // this class does not have a custom default constructor and thus should not
  // provide any custom default cstors, dstor, or assignment. see ISOCPP C.20.

  /// Access the spatial position vector.
  ///
  /// The surface is owned by the parameters object and thus is independent
  /// from the geometry context.
  Vector3D position() const { return Base::position(GeometryContext()); }
  // Make sure that the position access via geometry context is also available
  // so that bound and curvilinear parameters can be used interchangeably.
  using Base::position;
};

}  // namespace Acts
