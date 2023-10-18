// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <cmath>
#include <memory>
#include <type_traits>
#include <utility>

namespace Acts {

/// This class is only a light wrapper around a surface and a vector of
/// parameters. Its main purpose is to provide many constructors for the
/// underlying vector. Most accessors are generated from the
/// BoundTrackParameters equivalent and thus may be expensive
/// @note This class holds shared ownership on its reference surface.
/// @note The accessors for parameters, covariance, position, etc.
/// are the weighted means of the components.
/// @note If all covariances are zero, the accessor for the total
/// covariance does return std::nullopt;
/// TODO Add constructor from range and projector maybe?
class MultiComponentBoundTrackParameters {
  using Parameters = BoundTrackParameters;
  using ParticleHypothesis = Parameters::ParticleHypothesis;

  std::vector<std::tuple<double, BoundVector, std::optional<BoundSquareMatrix>>>
      m_components;
  std::shared_ptr<const Surface> m_surface;

  // TODO use [[no_unique_address]] once we switch to C++20
  ParticleHypothesis m_particleHypothesis;

  /// Helper function to reduce the component vector to a single representation
  template <typename projector_t>
  auto reduce(projector_t&& proj) const {
    using Ret = std::decay_t<decltype(proj(std::declval<Parameters>()))>;

    Ret ret;

    if constexpr (std::is_floating_point_v<Ret>) {
      ret = 0.0;
    } else {
      ret = Ret::Zero();
    }

    for (auto i = 0ul; i < m_components.size(); ++i) {
      const auto [weight, single_pars] = (*this)[i];
      ret += weight * proj(single_pars);
    }

    return ret;
  }

 public:
  using Scalar = typename Parameters::Scalar;
  using ParametersVector = typename Parameters::ParametersVector;
  using CovarianceMatrix = typename Parameters::CovarianceMatrix;

  /// Construct from multiple components
  template <typename covariance_t>
  MultiComponentBoundTrackParameters(
      std::shared_ptr<const Surface> surface,
      const std::vector<std::tuple<double, BoundVector, covariance_t>>& cmps,
      ParticleHypothesis particleHypothesis)
      : m_surface(std::move(surface)),
        m_particleHypothesis(particleHypothesis) {
    static_assert(
        std::is_same_v<BoundSquareMatrix, covariance_t> ||
        std::is_same_v<std::optional<BoundSquareMatrix>, covariance_t>);
    if constexpr (std::is_same_v<BoundSquareMatrix, covariance_t>) {
      for (const auto& [weight, params, cov] : cmps) {
        m_components.push_back({weight, params, cov});
      }
    } else {
      m_components = cmps;
    }
  }

  /// Construct from a parameters vector on the surface and particle charge.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param params Bound parameters vector
  /// @param particleHypothesis Particle hypothesis for these parameters
  /// @param cov Bound parameters covariance matrix
  ///
  /// In principle, only the charge magnitude is needed her to allow
  /// unambiguous extraction of the absolute momentum. The particle charge is
  /// required as an input here to be consistent with the other constructors
  /// below that that also take the charge as an input. The charge sign is
  /// only used in debug builds to check for consistency with the q/p
  /// parameter.
  MultiComponentBoundTrackParameters(std::shared_ptr<const Surface> surface,
                                     const BoundVector& params,
                                     std::optional<BoundSquareMatrix> cov,
                                     ParticleHypothesis particleHypothesis)
      : m_surface(std::move(surface)),
        m_particleHypothesis(particleHypothesis) {
    m_components.push_back({1., params, std::move(cov)});
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

  /// Access the parameters
  const auto& components() const { return m_components; }

  /// Reference surface onto which the parameters are bound.
  const Surface& referenceSurface() const { return *m_surface; }

  /// Get the weight and a GenericBoundTrackParameters object for one component
  std::pair<double, Parameters> operator[](std::size_t i) const {
    return std::make_pair(
        std::get<double>(m_components[i]),
        Parameters(m_surface, std::get<BoundVector>(m_components[i]),
                   std::get<std::optional<BoundSquareMatrix>>(m_components[i]),
                   m_particleHypothesis));
  }

  /// Parameters vector.
  ParametersVector parameters() const {
    return reduce([](const Parameters& p) { return p.parameters(); });
  }

  /// Optional covariance matrix.
  std::optional<CovarianceMatrix> covariance() const {
    const auto ret = reduce([](const Parameters& p) {
      return p.covariance() ? *p.covariance() : CovarianceMatrix::Zero();
    });

    if (ret == CovarianceMatrix::Zero()) {
      return std::nullopt;
    } else {
      return ret;
    }
  }

  /// Access a single parameter value identified by its index.
  ///
  /// @tparam kIndex Track parameter index
  template <BoundIndices kIndex>
  Scalar get() const {
    return reduce([&](const Parameters& p) { return p.get<kIndex>(); });
  }

  /// Space-time position four-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
  Vector4 fourPosition(const GeometryContext& geoCtx) const {
    return reduce([&](const Parameters& p) { return p.fourPosition(geoCtx); });
  }

  /// Spatial position three-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
  Vector3 position(const GeometryContext& geoCtx) const {
    return reduce([&](const Parameters& p) { return p.position(geoCtx); });
  }

  /// Time coordinate.
  Scalar time() const {
    return reduce([](const Parameters& p) { return p.time(); });
  }

  /// Unit direction three-vector, i.e. the normalized momentum
  /// three-vector.
  Vector3 direction() const {
    return reduce([](const Parameters& p) { return p.direction(); })
        .normalized();
  }

  /// Phi direction.
  Scalar phi() const { return VectorHelpers::phi(direction()); }

  /// Theta direction.
  Scalar theta() const { return VectorHelpers::theta(direction()); }

  /// Charge over momentum.
  Scalar qOverP() const { return get<eBoundQOverP>(); }

  /// Absolute momentum.
  Scalar absoluteMomentum() const {
    return reduce([](const Parameters& p) { return p.absoluteMomentum(); });
  }

  /// Transverse momentum.
  Scalar transverseMomentum() const {
    return reduce([](const Parameters& p) { return p.transverseMomentum(); });
  }

  /// Momentum three-vector.
  Vector3 momentum() const {
    return reduce([](const Parameters& p) { return p.momentum(); });
  }

  /// Particle electric charge.
  Scalar charge() const {
    return reduce([](const Parameters& p) { return p.charge(); });
  }

  /// Particle hypothesis.
  const ParticleHypothesis& particleHypothesis() const {
    return m_particleHypothesis;
  }
};

}  // namespace Acts
