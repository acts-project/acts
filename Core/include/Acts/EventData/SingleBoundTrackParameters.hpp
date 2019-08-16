// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/EventData/detail/PrintParameters.hpp"
#include "Acts/EventData/detail/TransformationFreeToBound.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <type_traits>

namespace Acts {

/// Track parameters bound to a reference surface for a single track.
///
/// @tparam charge_t Helper type to interpret the particle charge/momentum
///
/// This is intended as a user-facing data class that adds additional accessors
/// and charge/momentum interpretation on top of the pure parameters vector. All
/// parameters and their corresponding covariance matrix are stored in bound
/// parametrization. The specific definition of the local spatial parameters is
/// defined by the associated surface.
///
/// @note This class holds shared ownership on its reference surface.
template <class charge_t>
class SingleBoundTrackParameters {
 public:
  using Scalar = BoundScalar;
  using ParametersVector = BoundVector;
  using CovarianceMatrix = BoundSymMatrix;

  /// Construct from a parameters vector on the surface and particle charge.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param params Bound parameters vector
  /// @param q Particle charge
  /// @param cov Bound parameters covariance matrix
  ///
  /// In principle, only the charge magnitude is needed her to allow unambigous
  /// extraction of the absolute momentum. The particle charge is required as
  /// an input here to be consistent with the other constructors below that
  /// that also take the charge as an input. The charge sign is only used in
  /// debug builds to check for consistency with the q/p parameter.
  SingleBoundTrackParameters(std::shared_ptr<const Surface> surface,
                             const ParametersVector& params, Scalar q,
                             std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_paramSet(std::move(cov), params),
        m_surface(std::move(surface)),
        m_chargeInterpreter(std::abs(q)) {
    assert((0 <= (params[eBoundQOverP] * q)) and
           "Inconsistent q/p and q signs");
    assert(m_surface);
  }

  /// Construct from a parameters vector on the surface.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param params Bound parameters vector
  /// @param cov Bound parameters covariance matrix
  ///
  /// This constructor is only available if there are no potential charge
  /// ambiguities, i.e. the charge type is default-constructible.
  template <typename T = charge_t,
            std::enable_if_t<std::is_default_constructible_v<T>, int> = 0>
  SingleBoundTrackParameters(std::shared_ptr<const Surface> surface,
                             const ParametersVector& params,
                             std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_paramSet(std::move(cov), params),
        m_surface(std::move(surface)),
        m_chargeInterpreter(T()) {
    assert(m_surface);
  }

  /// Construct from four-position, direction, absolute momentum, and charge.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param geoCtx Geometry context for the local-to-global transformation
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored.
  /// @param p Absolute momentum
  /// @param q Particle charge
  /// @param cov Bound parameters covariance matrix
  SingleBoundTrackParameters(std::shared_ptr<const Surface> surface,
                             const GeometryContext& geoCtx,
                             const Vector4D& pos4, const Vector3D& dir,
                             Scalar p, Scalar q,
                             std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_paramSet(std::move(cov),
                   detail::transformFreeToBoundParameters(
                       pos4.segment<3>(ePos0), pos4[eTime], dir,
                       (q != Scalar(0)) ? (q / p) : (1 / p), *surface, geoCtx)),
        m_surface(std::move(surface)),
        m_chargeInterpreter(std::abs(q)) {
    assert((0 <= p) and "Absolute momentum must be positive");
    assert(m_surface);
  }

  /// Construct from four-position, direction, and charge-over-momentum.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param geoCtx Geometry context for the local-to-global transformation
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored.
  /// @param qOverP Charge-over-momentum-like parameter
  /// @param cov Bound parameters covariance matrix
  ///
  /// This constructor is only available if there are no potential charge
  /// ambiguities, i.e. the charge type is default-constructible. The position
  /// must be located on the surface.
  template <typename T = charge_t,
            std::enable_if_t<std::is_default_constructible_v<T>, int> = 0>
  SingleBoundTrackParameters(std::shared_ptr<const Surface> surface,
                             const GeometryContext& geoCtx,
                             const Vector4D& pos4, const Vector3D& dir,
                             Scalar qOverP,
                             std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_paramSet(std::move(cov), detail::transformFreeToBoundParameters(
                                       pos4.segment<3>(ePos0), pos4[eTime], dir,
                                       qOverP, *surface, geoCtx)),
        m_surface(std::move(surface)),
        m_chargeInterpreter(T()) {
    assert(m_surface);
  }

  // this class does not have a custom default constructor and thus should not
  // provide any custom default cstors, dstor, or assignment. see ISOCPP C.20.

  /// Access the parameter set holding the parameters vector and covariance.
  const FullBoundParameterSet& getParameterSet() const { return m_paramSet; }
  /// Parameters vector.
  ParametersVector parameters() const { return m_paramSet.getParameters(); }
  /// Optional covariance matrix.
  const std::optional<CovarianceMatrix>& covariance() const {
    return m_paramSet.getCovariance();
  }

  /// Access a single parameter value indentified by its index.
  ///
  /// @tparam kIndex Track parameter index
  template <BoundIndices kIndex>
  Scalar get() const {
    return m_paramSet.template getParameter<kIndex>();
  }
  /// Access a single parameter uncertainty identified by its index.
  ///
  /// @tparam kIndex Track parameter index
  /// @retval zero if the track parameters have no associated covariance
  /// @retval parameter standard deviation if the covariance is available
  template <BoundIndices kIndex>
  Scalar uncertainty() const {
    return m_paramSet.template getUncertainty<kIndex>();
  }

  /// Space-time position four-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global transformation
  ///
  /// This uses the associated surface to transform the local position on the
  /// surface to globalcoordinates. This requires a geometry context to select
  /// the appropriate transformation and might be a computationally expensive
  /// operation.
  Vector4D fourPosition(const GeometryContext& geoCtx) const {
    const Vector2D loc(get<eBoundLoc0>(), get<eBoundLoc1>());
    const Vector3D dir =
        makeDirectionUnitFromPhiTheta(get<eBoundPhi>(), get<eBoundTheta>());
    Vector4D pos4;
    pos4.segment<3>(ePos0) = m_surface->localToGlobal(geoCtx, loc, dir);
    pos4[eTime] = get<eBoundTime>();
    return pos4;
  }
  /// Spatial position three-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global transformation
  ///
  /// This uses the associated surface to transform the local position on the
  /// surface to globalcoordinates. This requires a geometry context to select
  /// the appropriate transformation and might be a computationally expensive
  /// operation.
  Vector3D position(const GeometryContext& geoCtx) const {
    const Vector2D loc(get<eBoundLoc0>(), get<eBoundLoc1>());
    const Vector3D dir =
        makeDirectionUnitFromPhiTheta(get<eBoundPhi>(), get<eBoundTheta>());
    return m_surface->localToGlobal(geoCtx, loc, dir);
  }
  /// Time coordinate.
  Scalar time() const { return get<eBoundTime>(); }

  /// Unit direction three-vector, i.e. the normalized momentum three-vector.
  Vector3D unitDirection() const {
    return makeDirectionUnitFromPhiTheta(get<eBoundPhi>(), get<eBoundTheta>());
  }
  /// Absolute momentum.
  Scalar absoluteMomentum() const {
    return m_chargeInterpreter.extractMomentum(get<eBoundQOverP>());
  }
  /// Transverse momentum.
  Scalar transverseMomentum() const {
    return std::sin(get<eBoundTheta>()) * absoluteMomentum();
  }
  /// Momentum three-vector.
  Vector3D momentum() const { return absoluteMomentum() * unitDirection(); }

  /// Particle electric charge.
  constexpr Scalar charge() const {
    return m_chargeInterpreter.extractCharge(get<eBoundQOverP>());
  }

  /// Reference surface onto which the parameters are bound.
  const Surface& referenceSurface() const { return *m_surface; }
  /// Reference frame in which the local error is defined.
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
  FullBoundParameterSet m_paramSet;
  /// reference surface
  std::shared_ptr<const Surface> m_surface;
  // TODO use [[no_unique_address]] once we switch to C++20
  charge_t m_chargeInterpreter;

  /// Compare two bound parameters for equality.
  friend bool operator==(const SingleBoundTrackParameters& lhs,
                         const SingleBoundTrackParameters& rhs) {
    return (lhs.m_paramSet == rhs.m_paramSet) and
           (lhs.m_surface == rhs.m_surface) and
           (lhs.m_chargeInterpreter == rhs.m_chargeInterpreter);
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
