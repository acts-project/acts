// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SingleBoundTrackParameters.hpp"
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

/// TODO as a quick start, the accessors only return the first components, Also
/// internaly just a vector of SingleBoundTrackParameters. Is this performance
/// critical code?

/// Track parameters bound to a reference surface for a single track.
///
/// @tparam charge_t Helper type to interpret the particle charge/momentum
///
/// TODO Adjust description for Multi
/// This is intended as a user-facing data class that adds additional accessors
/// and charge/momentum interpretation on top of the pure parameters vector. All
/// parameters and their corresponding covariance matrix are stored in bound
/// parametrization. The specific definition of the local spatial parameters is
/// defined by the associated surface.
///
/// @note This class holds shared ownership on its reference surface.
template <class charge_t>
class MultiComponentBoundTrackParameters {
 public:
  using Scalar = ActsScalar;
  using ParametersVector = BoundVector;
  using CovarianceMatrix = BoundSymMatrix;
  using SingleParameters = SingleBoundTrackParameters<charge_t>;

 private:
  /// Components
  std::vector<std::tuple<double, SingleParameters>> m_components;

  // TODO Actually redundant
  std::shared_ptr<const Surface> m_surface;

  // TODO use [[no_unique_address]] once we switch to C++20
  charge_t m_chargeInterpreter;

 public:
  /// Construct from multiple components
  template <typename covariance_t>
  MultiComponentBoundTrackParameters(
      std::shared_ptr<const Surface> surface,
      const std::vector<std::tuple<double, BoundVector, covariance_t>>& cmps,
      Scalar q)
      : m_surface(surface), m_chargeInterpreter(std::abs(q)) {
    static_assert(
        std::is_same_v<Acts::BoundSymMatrix, covariance_t> ||
        std::is_same_v<std::optional<Acts::BoundSymMatrix>, covariance_t>);
    for (const auto& [weight, params, cov] : cmps) {
      m_components.push_back({weight, SingleBoundTrackParameters<charge_t>(
                                          surface, params, cov, q)});
    }
  }

  /// Construct from multiple components
  template <typename covariance_t, typename T = charge_t,
            std::enable_if_t<std::is_default_constructible_v<T>, int> = 0>
  MultiComponentBoundTrackParameters(
      std::shared_ptr<const Surface> surface,
      const std::vector<std::tuple<double, BoundVector, covariance_t>>& cmps)
      : m_surface(surface) {
    static_assert(
        std::is_same_v<Acts::BoundSymMatrix, covariance_t> ||
        std::is_same_v<std::optional<Acts::BoundSymMatrix>, covariance_t>);
    for (const auto& [weight, params, cov] : cmps) {
      m_components.push_back({weight, SingleParameters(surface, params, cov)});
    }
  }

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
  MultiComponentBoundTrackParameters(
      std::shared_ptr<const Surface> surface, const ParametersVector& params,
      Scalar q, std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_surface(surface), m_chargeInterpreter(std::abs(q)) {
    m_components.push_back({1., SingleParameters(surface, params, cov, q)});
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
  MultiComponentBoundTrackParameters(
      std::shared_ptr<const Surface> surface, const ParametersVector& params,
      std::optional<CovarianceMatrix> cov = std::nullopt)
      : m_surface(surface) {
    m_components.push_back({1., SingleParameters(surface, params, cov)});
  }

  /// Parameters are not default constructible due to the charge type.
  MultiComponentBoundTrackParameters() = delete;
  MultiComponentBoundTrackParameters(
      const MultiComponentBoundTrackParameters&) = default;
  MultiComponentBoundTrackParameters(MultiComponentBoundTrackParameters&&) =
      default;
  ~MultiComponentBoundTrackParameters() = default;
  MultiComponentBoundTrackParameters& operator=(
      const MultiComponentBoundTrackParameters&) = default;
  MultiComponentBoundTrackParameters& operator=(
      MultiComponentBoundTrackParameters&&) = default;

  /// Parameters vector.
  const ParametersVector& parameters() const {
    return std::get<SingleParameters>(m_components.front()).parameters();
  }

  /// Return the components vector
  const auto& components() const { return m_components; }

  /// Optional covariance matrix.
  const std::optional<CovarianceMatrix>& covariance() const {
    return std::get<SingleParameters>(m_components.front()).covariance();
  }

  /// Access a single parameter value indentified by its index.
  ///
  /// @tparam kIndex Track parameter index
  template <BoundIndices kIndex>
  Scalar get() const {
    return std::get<SingleParameters>(m_components.front())
        .template get<kIndex>();
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
    return std::get<SingleParameters>(m_components.front())
        .fourPosition(geoCtx);
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
    return std::get<SingleParameters>(m_components.front()).position(geoCtx);
  }
  /// Time coordinate.
  Scalar time() const {
    return std::get<SingleParameters>(m_components.front()).time();
  }

  /// Unit direction three-vector, i.e. the normalized momentum
  /// three-vector.
  Vector3 unitDirection() const {
    return std::get<SingleParameters>(m_components.front()).unitDirection();
  }
  /// Absolute momentum.
  Scalar absoluteMomentum() const {
    return std::get<SingleParameters>(m_components.front()).absoluteMomentum();
  }
  /// Transverse momentum.
  Scalar transverseMomentum() const {
    return std::get<SingleParameters>(m_components.front())
        .transverseMomentum();
  }
  /// Momentum three-vector.
  Vector3 momentum() const { return absoluteMomentum() * unitDirection(); }

  /// Particle electric charge.
  Scalar charge() const {
    return std::get<SingleParameters>(m_components.front()).charge();
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
};

}  // namespace Acts
