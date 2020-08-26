// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/ChargePolicy.hpp"
#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/EventData/detail/PrintParameters.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <cassert>
#include <cmath>

namespace Acts {

/// Single track parameters bound to a reference surface.
///
/// @tparam charge_policy_t Selection type for the particle charge
///
/// This is a base class for bound parameters. All parameters and their
/// corresponding covariance matrix are stored in bound parametrization. The
/// specific definition of the local spatial parameters is defined by the
/// associated surface.
///
/// @note This class holds shared ownership on its reference surface.
template <class charge_policy_t>
class SingleBoundTrackParameters {
 public:
  using Scalar = BoundParametersScalar;
  using ParametersVector = BoundVector;
  using CovarianceMatrix = BoundSymMatrix;

  /// Construct charged bound parameters from parameters vector on the surface.
  ///
  /// @param[in] surface The reference surface the parameters are defined on
  /// @param[in] params The parameter vector
  /// @param[in] cov Optional covariance in the reference frame
  template <typename T = charge_policy_t,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleBoundTrackParameters(std::shared_ptr<const Surface> surface,
                             const ParametersVector& params,
                             std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_paramSet(std::move(cov), params),
        m_chargePolicy(std::copysign(1., params[eBoundQOverP])),
        m_surface(std::move(surface)) {
    assert(m_surface);
  }

  /// Construct neutral bound parameters from parameters vector on the surface.
  ///
  /// @param[in] surface The reference surface the parameters are defined on
  /// @param[in] params The parameter vector
  /// @param[in] cov Optional covariance in the reference frame
  template <typename T = charge_policy_t,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleBoundTrackParameters(std::shared_ptr<const Surface> surface,
                             const ParametersVector& params,
                             std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_paramSet(std::move(cov), params), m_surface(std::move(surface)) {
    assert(m_surface);
  }

  /// Construct charged bound parameters from global position and momentum.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global transformation
  /// @param[in] cov Optional covariance in the reference frame
  /// @param[in] pos The global track three-position vector
  /// @param[in] mom The global track three-momentum vector
  /// @param[in] charge The particle charge
  /// @param[in] time The time coordinate
  /// @param[in] surface The reference surface the parameters are bound to
  template <typename T = charge_policy_t,
            std::enable_if_t<std::is_same<T, ChargedPolicy>::value, int> = 0>
  SingleBoundTrackParameters(const GeometryContext& geoCtx,
                             std::optional<CovarianceMatrix> cov,
                             const Vector3D& pos, const Vector3D& mom,
                             Scalar charge, Scalar time,
                             std::shared_ptr<const Surface> surface)
      : m_paramSet(std::move(cov),
                   detail::transformFreeToBoundParameters(
                       pos, time, mom, charge / mom.norm(), *surface, geoCtx)),
        m_chargePolicy(charge),
        m_surface(std::move(surface)) {
    assert(m_surface);
  }

  /// Construct neutral bound parameters from global position and momentum.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global transformation
  /// @param[in] cov Optional covariance in the reference frame
  /// @param[in] pos The global track three-position vector
  /// @param[in] mom The global track three-momentum vector
  /// @param[in] time The time coordinate
  /// @param[in] surface The reference surface the parameters are bound to
  template <typename T = charge_policy_t,
            std::enable_if_t<std::is_same<T, NeutralPolicy>::value, int> = 0>
  SingleBoundTrackParameters(const GeometryContext& geoCtx,
                             std::optional<CovarianceMatrix> cov,
                             const Vector3D& pos, const Vector3D& mom,
                             Scalar time,
                             std::shared_ptr<const Surface> surface)
      : m_paramSet(std::move(cov),
                   detail::transformFreeToBoundParameters(
                       pos, time, mom, 1 / mom.norm(), *surface, geoCtx)),
        m_surface(std::move(surface)) {
    assert(m_surface);
  }

  // this class does not have a custom default constructor and thus should not
  // provide any custom default cstors, dstor, or assignment. see ISOCPP C.20.

  /// Access the parameter set holding the parameters vector and covariance.
  const FullParameterSet& getParameterSet() const { return m_paramSet; }
  /// Access the bound parameters vector.
  ParametersVector parameters() const { return m_paramSet.getParameters(); }
  /// Access the optional covariance matrix.
  const std::optional<CovarianceMatrix>& covariance() const {
    return m_paramSet.getCovariance();
  }

  /// Access a single parameter value indentified by its index.
  ///
  /// @tparam kIndex Track parameter index
  template <BoundParametersIndices kIndex>
  Scalar get() const {
    return m_paramSet.template getParameter<kIndex>();
  }
  /// Access a single parameter uncertainty identified by its index.
  ///
  /// @tparam kIndex Track parameter index
  /// @retval zero if the track parameters have no associated covariance
  /// @retval parameter standard deviation if the covariance is available
  template <BoundParametersIndices kIndex>
  Scalar uncertainty() const {
    return m_paramSet.getUncertainty<kIndex>();
  }

  /// Access the spatial position vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global transformation
  ///
  /// This uses the associated surface to transform the local position to global
  /// coordinates. This requires a geometry context to select the appropriate
  /// transformation and might be a computationally expensive operation.
  Vector3D position(const GeometryContext& geoCtx) const {
    const Vector2D loc(get<eBoundLoc0>(), get<eBoundLoc1>());
    const Vector3D dir =
        makeDirectionUnitFromPhiTheta(get<eBoundPhi>(), get<eBoundTheta>());
    Vector3D pos;
    m_surface->localToGlobal(geoCtx, loc, dir, pos);
    return pos;
  }
  /// Access the time coordinate.
  Scalar time() const { return get<eBoundTime>(); }

  /// Access the direction pseudo-rapidity.
  Scalar eta() const { return -std::log(std::tan(get<eBoundTheta>() / 2)); }
  /// Access the absolute transverse momentum.
  Scalar pT() const {
    return std::sin(get<eBoundTheta>()) / std::abs(get<eBoundQOverP>());
  }
  /// Access the momentum three-vector.
  Vector3D momentum() const {
    auto mom =
        makeDirectionUnitFromPhiTheta(get<eBoundPhi>(), get<eBoundTheta>());
    mom *= std::abs(1 / get<eBoundQOverP>());
    return mom;
  }

  /// Access the particle electric charge.
  Scalar charge() const { return m_chargePolicy.getCharge(); }

  /// Access the reference surface onto which the parameters are bound.
  const Surface& referenceSurface() const { return *m_surface; }
  /// Access the reference frame in which the local error is defined.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global transformation
  ///
  /// For planar surfaces, this is the transformation local-to-global rotation
  /// matrix. For non-planar surfaces, it is the local-to-global rotation matrix
  /// of the tangential plane at the track position.
  RotationMatrix3D referenceFrame(const GeometryContext& geoCtx) const {
    return m_surface->referenceFrame(geoCtx, position(geoCtx), momentum());
  }

 private:
  /// parameter set holding parameters vector and covariance.
  FullParameterSet m_paramSet;
  /// charge policy to distinguish differently charged particles
  charge_policy_t m_chargePolicy;
  /// reference surface
  std::shared_ptr<const Surface> m_surface;

  /// Compare two bound parameters for equality.
  friend bool operator==(const SingleBoundTrackParameters& lhs,
                         const SingleBoundTrackParameters& rhs) {
    return (lhs.m_paramSet == rhs.m_paramSet) and
           (lhs.m_chargePolicy == rhs.m_chargePolicy) and
           (lhs.m_surface == rhs.m_surface);
  }
  /// Compare two bound parameters for inequality.
  friend bool operator!=(const SingleBoundTrackParameters& lhs,
                         const SingleBoundTrackParameters& rhs) {
    return !(lhs == rhs);
  }
  /// Print information to the output stream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const SingleBoundTrackParameters& tp) {
    detail::printBoundParameters(
        os, tp.referenceSurface(), tp.parameters(),
        tp.covariance().has_value() ? &tp.covariance().value() : nullptr);
    return os;
  }
};

}  // namespace Acts
