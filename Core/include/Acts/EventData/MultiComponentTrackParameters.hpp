// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

#include <cmath>
#include <memory>
#include <type_traits>
#include <utility>

namespace Acts {

/// This class is only a light wrapper around a surface and a vector of
/// parameters. Its main purpose is to provide many constructors for the
/// underlying vector. Most accessors are generated from the
/// GenericBoundTrackParameters equivalent and thus may be expensive
/// @tparam charge_t Helper type to interpret the particle charge/momentum
/// @note This class holds shared ownership on its reference surface.
/// @note The accessors for parameters, covariance, position, etc.
/// are the weighted means of the components.
/// TODO Add constructor from range and projector maybe?
template <typename charge_t>
class MultiComponentBoundTrackParameters {
  using SingleParameters = GenericBoundTrackParameters<charge_t>;

  std::vector<std::tuple<double, BoundVector, BoundSquareMatrix>>
      m_components;
  std::shared_ptr<const Surface> m_surface;

  // TODO use [[no_unique_address]] once we switch to C++20
  charge_t m_chargeInterpreter;

  /// Helper function to reduce the component vector to a single representation
  template <typename projector_t>
  auto reduce(projector_t&& proj) const {
    using Ret = std::decay_t<decltype(proj(std::declval<SingleParameters>()))>;

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
  using Scalar = typename SingleParameters::Scalar;
  using ParametersVector = typename SingleParameters::ParametersVector;
  using CovarianceMatrix = typename SingleParameters::CovarianceMatrix;

  /// Construct from multiple components with charge scalar q
  MultiComponentBoundTrackParameters(
      std::shared_ptr<const Surface> surface,
      const std::vector<std::tuple<double, BoundVector, CovarianceMatrix>>& cmps,
      ActsScalar q)
      : m_components(cmps), m_surface(std::move(surface)), m_chargeInterpreter(std::abs(q)) {
  }

  /// Construct from multiple components with charge_t
  template <typename T = charge_t,
            std::enable_if_t<std::is_default_constructible_v<T>, int> = 0>
  MultiComponentBoundTrackParameters(
      std::shared_ptr<const Surface> surface,
      const std::vector<std::tuple<double, BoundVector, CovarianceMatrix>>& cmps)
      : m_components(cmps), m_surface(std::move(surface)) {
  }

  /// Construct from a parameters vector on the surface and particle charge.
  ///
  /// @param surface Reference surface the parameters are defined on
  /// @param params Bound parameters vector
  /// @param q Particle charge
  /// @param cov Bound parameters covariance matrix
  ///
  /// In principle, only the charge magnitude is needed her to allow
  /// unambiguous extraction of the absolute momentum. The particle charge is
  /// required as an input here to be consistent with the other constructors
  /// below that that also take the charge as an input. The charge sign is
  /// only used in debug builds to check for consistency with the q/p
  /// parameter.
  MultiComponentBoundTrackParameters(
      std::shared_ptr<const Surface> surface, const BoundVector& params,
      ActsScalar q, CovarianceMatrix cov)
      : m_surface(std::move(surface)), m_chargeInterpreter(std::abs(q)) {
    m_components.push_back({1., params, cov});
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
      std::shared_ptr<const Surface> surface, const BoundVector& params, CovarianceMatrix cov)
      : m_surface(std::move(surface)) {
    m_components.push_back({1., params, cov});
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
  std::pair<double, SingleParameters> operator[](std::size_t i) const {
    return std::make_pair(
        std::get<double>(m_components[i]),
        SingleParameters(
            m_surface, std::get<BoundVector>(m_components[i]),
            std::get<BoundSquareMatrix>(m_components[i])));
  }

  /// Parameters vector.
  ParametersVector parameters() const {
    return reduce([](const SingleParameters& p) { return p.parameters(); });
  }

  /// Optional covariance matrix.
  CovarianceMatrix covariance() const {
    return reduce([](const SingleParameters& p) {
      return *p.covariance();
    });
  }

  /// Space-time position four-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
  Vector4 fourPosition(const GeometryContext& geoCtx) const {
    return reduce(
        [&](const SingleParameters& p) { return p.fourPosition(geoCtx); });
  }

  /// Spatial position three-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
  Vector3 position(const GeometryContext& geoCtx) const {
    return reduce(
        [&](const SingleParameters& p) { return p.position(geoCtx); });
  }

  /// Time coordinate.
  Scalar time() const {
    return reduce([](const SingleParameters& p) { return p.time(); });
  }

  /// Unit direction three-vector, i.e. the normalized momentum
  /// three-vector.
  Vector3 direction() const {
    return reduce([](const SingleParameters& p) { return p.direction(); })
        .normalized();
  }

  /// Absolute momentum.
  Scalar absoluteMomentum() const {
    return reduce(
        [](const SingleParameters& p) { return p.absoluteMomentum(); });
  }

  /// Transverse momentum.
  Scalar transverseMomentum() const {
    return reduce(
        [](const SingleParameters& p) { return p.transverseMomentum(); });
  }

  /// Momentum three-vector.
  Vector3 momentum() const {
    return reduce([](const SingleParameters& p) { return p.momentum(); });
  }

  /// Particle electric charge.
  Scalar charge() const {
    return reduce([](const SingleParameters& p) { return p.charge(); });
  }
};

/// This class mimics the behaviour of the curvilinear parameters for ordinary
/// track parameters. To adopt this concept, a "common surface" is constructed,
/// and all parameters are projected onto this surface. The use of this is
/// questionable, and if the result is reasonable depends largely on the initial
/// multi component state. However, the propagator infrastructure forces the
/// existance of this type
/// @tparam charge_t Helper type to interpret the particle charge/momentum
template <typename charge_t>
class MultiComponentCurvilinearTrackParameters
    : public MultiComponentBoundTrackParameters<charge_t> {
  using covariance_t = BoundSquareMatrix;
 public:
  // template <typename covariance_t>
  using ConstructionTuple = std::tuple<double, Acts::Vector4, Acts::Vector3,
                                       ActsScalar, covariance_t>;

 private:
  using Base = MultiComponentBoundTrackParameters<charge_t>;

  // template <typename covariance_t>
  using BaseConstructionTuple =
      std::tuple<std::shared_ptr<Acts::Surface>,
                 std::vector<std::tuple<double, BoundVector, covariance_t>>,
                 ActsScalar>;

  /// We need this helper function in order to construct the base class properly
  // template <typename covariance_t>
  static BaseConstructionTuple construct(
      const std::vector<ConstructionTuple>& curvi) {
    // TODO where to get a geometry context here
    Acts::GeometryContext gctx{};

    // Construct and average surface
    Acts::Vector3 avgPos = Acts::Vector3::Zero();
    Acts::Vector3 avgDir = Acts::Vector3::Zero();
    for (const auto& [w, pos4, dir, qop, cov] : curvi) {
      avgPos += w * pos4.template segment<3>(0);
      avgDir += w * dir;
    }

    auto s = Surface::makeShared<PlaneSurface>(avgPos, avgDir);

    std::vector<std::tuple<double, BoundVector, covariance_t>> bound;
    bound.reserve(curvi.size());

    // Project the position onto the surface, keep everything else as is
    for (const auto& [w, pos4, dir, qop, cov] : curvi) {
      Vector3 newPos =
          s->intersect(gctx, pos4.template segment<3>(eFreePos0), dir, false)
              .intersection.position;
      
      BoundVector bv =
          detail::transformFreeToCurvilinearParameters(pos4[eTime], dir, qop);
          
      // Because of the projection this should never fail
      bv.template segment<2>(eBoundLoc0) = *(s->globalToLocal(gctx, newPos, dir));
      
      bound.emplace_back(w, bv, cov);
    }

    // TODO bring this first to compile then check qop
    return {s, bound, {}};
  }

  // Private constructor from a tuple
  // template <typename covariance_t>
  MultiComponentCurvilinearTrackParameters(
      const BaseConstructionTuple& t)
      : Base(std::get<0>(t), std::get<1>(t), std::get<2>(t)) {}

 public:
  // template <typename covariance_t>
  MultiComponentCurvilinearTrackParameters(
      const std::vector<ConstructionTuple>& cmps)
      : MultiComponentCurvilinearTrackParameters(construct(cmps)) {}
};

}  // namespace Acts
