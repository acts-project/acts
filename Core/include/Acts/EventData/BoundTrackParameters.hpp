// This file is part of the Acts project.
//
// Copyright (C) 2016-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
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

/// Track parameters bound to a reference surface.
///
/// This is intended as a user-facing data class that adds additional accessors
/// and charge/momentum interpretation on top of the pure parameters vector. All
/// parameters and their corresponding covariance matrix are stored in bound
/// parametrization. The specific definition of the local spatial parameters is
/// defined by the associated surface.
///
/// @note This class holds shared ownership on its reference surface.
class BoundTrackParameters {
 public:
  using Scalar = ActsScalar;
  using ParametersVector = BoundVector;
  using CovarianceMatrix = BoundSymMatrix;

  /// Construct from a parameters vector on the surface.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param params Bound parameters vector
  /// @param cov Bound parameters covariance matrix
  BoundTrackParameters(std::shared_ptr<const Surface> surface,
                       const ParametersVector& params,
                       std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_params(params), m_cov(std::move(cov)), m_surface(std::move(surface)) {
    assert(m_surface);
    normalizePhiTheta();
  }

  /// Factory to construct from four-position, direction and qop.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param geoCtx Geometry context for the local-to-global transformation
  /// @param pos4 Track position/time four-vectorl
  /// @param dir Track direction three-vector; normalization is ignored
  /// @param qOverP Charge over absolute momentum
  /// @param cov Bound parameters covariance matrix
  ///
  /// @note The returned result indicates whether the free parameters could
  /// successfully be converted to on-surface parameters.
  static Result<BoundTrackParameters> create(
      std::shared_ptr<const Surface> surface, const GeometryContext& geoCtx,
      const Vector4& pos4, const Vector3& dir, Scalar qOverP,
      std::optional<CovarianceMatrix> cov = std::nullopt) {
    Result<BoundVector> bound = detail::transformFreeToBoundParameters(
        pos4.segment<3>(ePos0), pos4[eTime], dir, qOverP, *surface, geoCtx);

    if (!bound.ok()) {
      return bound.error();
    }

    return BoundTrackParameters{std::move(surface), *bound, std::move(cov)};
  }

  /// Parameters are not default constructible due to the charge type.
  BoundTrackParameters() = delete;
  BoundTrackParameters(const BoundTrackParameters&) = default;
  BoundTrackParameters(BoundTrackParameters&&) = default;
  ~BoundTrackParameters() = default;
  BoundTrackParameters& operator=(const BoundTrackParameters&) = default;
  BoundTrackParameters& operator=(BoundTrackParameters&&) = default;

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

  /// Local spatial position two-vector.
  Vector2 localPosition() const { return m_params.segment<2>(eBoundLoc0); }
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
    Vector4 pos4;
    pos4.segment<3>(ePos0) =
        m_surface->localToGlobal(geoCtx, localPosition(), direction());
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
    return m_surface->localToGlobal(geoCtx, localPosition(), direction());
  }
  /// Time coordinate.
  Scalar time() const { return m_params[eBoundTime]; }

  /// Unit direction three-vector, i.e. the normalized momentum
  /// three-vector.
  Vector3 direction() const {
    return makeDirectionUnitFromPhiTheta(m_params[eBoundPhi],
                                         m_params[eBoundTheta]);
  }

  Scalar phi() const { return get<eBoundPhi>(); }
  Scalar theta() const { return get<eBoundTheta>(); }
  Scalar qop() const { return get<eBoundQOverP>(); }

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
    return m_surface->referenceFrame(geoCtx, position(geoCtx), direction());
  }

 private:
  BoundVector m_params;
  std::optional<BoundSymMatrix> m_cov;
  /// reference surface
  std::shared_ptr<const Surface> m_surface;

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
  friend bool operator==(const BoundTrackParameters& lhs,
                         const BoundTrackParameters& rhs) {
    return (lhs.m_params == rhs.m_params) and (lhs.m_cov == rhs.m_cov) and
           (lhs.m_surface == rhs.m_surface);
  }
  /// Compare two bound track parameters for bitwise in-equality.
  friend bool operator!=(const BoundTrackParameters& lhs,
                         const BoundTrackParameters& rhs) {
    return not(lhs == rhs);
  }
  /// Print information to the output stream.
  friend std::ostream& operator<<(std::ostream& os,
                                  const BoundTrackParameters& tp) {
    detail::printBoundParameters(
        os, tp.referenceSurface(), tp.parameters(),
        tp.covariance().has_value() ? &tp.covariance().value() : nullptr);
    return os;
  }
};

}  // namespace Acts
