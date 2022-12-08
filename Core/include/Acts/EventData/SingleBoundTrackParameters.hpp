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
  using Scalar = ActsScalar;
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
      : m_params(params),
        m_cov(std::move(cov)),
        m_surface(std::move(surface)),
        m_chargeInterpreter(std::abs(q)) {
    assert((0 <= (params[eBoundQOverP] * q)) and
           "Inconsistent q/p and q signs");
    assert(m_surface);
    normalizePhiTheta();
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
      : m_params(params), m_cov(std::move(cov)), m_surface(std::move(surface)) {
    assert(m_surface);
    normalizePhiTheta();
  }

  /// Factory to construct from four-position, direction, absolute momentum, and
  /// charge.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param geoCtx Geometry context for the local-to-global transformation
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored.
  /// @param p Absolute momentum
  /// @param q Particle charge
  /// @param cov Bound parameters covariance matrix
  ///
  /// @note The returned result indicates whether the free parameters could
  /// successfully be converted to on-surface parameters.
  static Result<SingleBoundTrackParameters<charge_t>> create(
      std::shared_ptr<const Surface> surface, const GeometryContext& geoCtx,
      const Vector4& pos4, const Vector3& dir, Scalar p, Scalar q,
      std::optional<CovarianceMatrix> cov = std::nullopt) {
    Result<BoundVector> bound = detail::transformFreeToBoundParameters(
        pos4.segment<3>(ePos0), pos4[eTime], dir,
        (q != Scalar(0)) ? (q / p) : (1 / p), *surface, geoCtx);

    if (!bound.ok()) {
      return bound.error();
    }

    return SingleBoundTrackParameters<charge_t>{std::move(surface), *bound, q,
                                                std::move(cov)};
  }

  /// Factory to construct from four-position, direction, and
  /// charge-over-momentum.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param geoCtx Geometry context for the local-to-global transformation
  /// @param pos4 Track position/time four-vector
  /// @param dir Track direction three-vector; normalization is ignored.
  /// @param qOverP Charge-over-momentum-like parameter
  /// @param cov Bound parameters covariance matrix
  ///
  /// @note This factory is only available if there are no potential charge
  /// ambiguities, i.e. the charge type is default-constructible. The
  /// position must be located on the surface.
  /// @note The returned result indicates whether the free parameters could
  /// successfully be converted to on-surface parameters.
  template <typename T = charge_t,
            std::enable_if_t<std::is_default_constructible_v<T>, int> = 0>
  static Result<SingleBoundTrackParameters<charge_t>> create(
      std::shared_ptr<const Surface> surface, const GeometryContext& geoCtx,
      const Vector4& pos4, const Vector3& dir, Scalar qOverP,
      std::optional<CovarianceMatrix> cov = std::nullopt) {
    Result<BoundVector> bound = detail::transformFreeToBoundParameters(
        pos4.segment<3>(ePos0), pos4[eTime], dir, qOverP, *surface, geoCtx);
    if (!bound.ok()) {
      return bound.error();
    }

    return SingleBoundTrackParameters<charge_t>{std::move(surface), *bound,
                                                std::move(cov)};
  }

  /// Parameters are not default constructible due to the charge type.
  SingleBoundTrackParameters() = delete;
  SingleBoundTrackParameters(const SingleBoundTrackParameters&) = default;
  SingleBoundTrackParameters(SingleBoundTrackParameters&&) = default;
  ~SingleBoundTrackParameters() = default;
  SingleBoundTrackParameters& operator=(const SingleBoundTrackParameters&) =
      default;
  SingleBoundTrackParameters& operator=(SingleBoundTrackParameters&&) = default;

  /// Parameters vector.
  ParametersVector& parameters() { return m_params; }
  /// Parameters vector.
  const ParametersVector& parameters() const { return m_params; }
  /// Optional covariance matrix.
  std::optional<CovarianceMatrix>& covariance() { return m_cov; }
  /// Optional covariance matrix.
  const std::optional<CovarianceMatrix>& covariance() const { return m_cov; }

  /// Access a single parameter value indentified by its index.
  ///
  /// @tparam kIndex Track parameter index
  template <BoundIndices kIndex>
  Scalar get() const {
    return m_params[kIndex];
  }

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
    const Vector2 loc(m_params[eBoundLoc0], m_params[eBoundLoc1]);
    const Vector3 dir = makeDirectionUnitFromPhiTheta(m_params[eBoundPhi],
                                                      m_params[eBoundTheta]);
    Vector4 pos4;
    pos4.segment<3>(ePos0) = m_surface->localToGlobal(geoCtx, loc, dir);
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
    const Vector2 loc(m_params[eBoundLoc0], m_params[eBoundLoc1]);
    const Vector3 dir = makeDirectionUnitFromPhiTheta(m_params[eBoundPhi],
                                                      m_params[eBoundTheta]);
    return m_surface->localToGlobal(geoCtx, loc, dir);
  }
  /// Time coordinate.
  Scalar time() const { return m_params[eBoundTime]; }

  /// Unit direction three-vector, i.e. the normalized momentum
  /// three-vector.
  Vector3 unitDirection() const {
    return makeDirectionUnitFromPhiTheta(m_params[eBoundPhi],
                                         m_params[eBoundTheta]);
  }
  /// Absolute momentum.
  Scalar absoluteMomentum() const {
    return m_chargeInterpreter.extractMomentum(m_params[eBoundQOverP]);
  }
  /// Transverse momentum.
  Scalar transverseMomentum() const {
    return std::sin(m_params[eBoundTheta]) * absoluteMomentum();
  }
  /// Momentum three-vector.
  Vector3 momentum() const { return absoluteMomentum() * unitDirection(); }

  /// Particle electric charge.
  Scalar charge() const {
    return m_chargeInterpreter.extractCharge(get<eBoundQOverP>());
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
  std::optional<BoundSymMatrix> m_cov;
  /// reference surface
  std::shared_ptr<const Surface> m_surface;
  // TODO use [[no_unique_address]] once we switch to C++20
  charge_t m_chargeInterpreter;

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
  ///   parameter values, or compare them for compability instead of
  ///   equality; you might also have different (floating point) thresholds
  ///   of equality in different contexts. None of that can be handled by
  ///   this operator. Users should think really hard if this is what they
  ///   want and we might decided that we will remove this in the future.
  friend bool operator==(const SingleBoundTrackParameters<charge_t>& lhs,
                         const SingleBoundTrackParameters<charge_t>& rhs) {
    return (lhs.m_params == rhs.m_params) and (lhs.m_cov == rhs.m_cov) and
           (lhs.m_surface == rhs.m_surface) and
           (lhs.m_chargeInterpreter == rhs.m_chargeInterpreter);
  }
  /// Compare two bound track parameters for bitwise in-equality.
  friend bool operator!=(const SingleBoundTrackParameters<charge_t>& lhs,
                         const SingleBoundTrackParameters<charge_t>& rhs) {
    return not(lhs == rhs);
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
