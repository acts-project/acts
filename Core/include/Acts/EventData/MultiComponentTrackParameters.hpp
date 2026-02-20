// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/GenericBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

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
 public:
  /// Type alias for bound track parameters
  using Parameters = BoundTrackParameters;
  /// Type alias for particle hypothesis
  using ParticleHypothesis = Parameters::ParticleHypothesis;
  /// Type alias for bound parameters vector
  using ParametersVector = typename Parameters::ParametersVector;
  /// Type alias for covariance matrix
  using CovarianceMatrix = typename Parameters::CovarianceMatrix;

 private:
  std::vector<std::tuple<double, BoundVector, std::optional<BoundMatrix>>>
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
  /// Type alias for construction tuple containing weight, position, direction,
  /// q/p, and covariance
  using ConstructionTuple = std::tuple<double, Acts::Vector4, Acts::Vector3,
                                       double, CovarianceMatrix>;

  /// We need this helper function in order to construct the base class properly
  /// @param geoCtx Geometry context for construction
  /// @param curvi Vector of construction tuples containing component data
  /// @param particleHypothesis Particle hypothesis for the parameters
  /// @return Multi-component bound track parameters in curvilinear representation
  static MultiComponentBoundTrackParameters createCurvilinear(
      const GeometryContext& geoCtx,
      const std::vector<ConstructionTuple>& curvi,
      ParticleHypothesis particleHypothesis) {
    // Construct and average surface
    Acts::Vector3 avgPos = Acts::Vector3::Zero();
    Acts::Vector3 avgDir = Acts::Vector3::Zero();
    for (const auto& [w, pos4, dir, qop, cov] : curvi) {
      avgPos += w * pos4.template segment<3>(0);
      avgDir += w * dir;
    }

    std::shared_ptr<PlaneSurface> s =
        CurvilinearSurface(avgPos, avgDir).planeSurface();

    std::vector<std::tuple<double, ParametersVector, CovarianceMatrix>> bound;
    bound.reserve(curvi.size());

    // Project the position onto the surface, keep everything else as is
    for (const auto& [w, pos4, dir, qop, cov] : curvi) {
      Intersection3D closestIntersection =
          s->intersect(geoCtx, pos4.template segment<3>(eFreePos0), dir,
                       BoundaryTolerance::Infinite())
              .closest();
      const Vector3& newPos = closestIntersection.position();

      ParametersVector bv =
          transformFreeToCurvilinearParameters(pos4[eTime], dir, qop);

      // Because of the projection this should never fail
      bv.template segment<2>(eBoundLoc0) =
          *(s->globalToLocal(geoCtx, newPos, dir));

      bound.emplace_back(w, bv, cov);
    }

    return MultiComponentBoundTrackParameters(s, bound, particleHypothesis);
  }

  /// Construct from multiple components
  /// @param surface Surface on which the parameters are bound
  /// @param cmps Vector of weight, parameters vector, and covariance components
  /// @param particleHypothesis Particle hypothesis for the parameters
  template <typename covariance_t>
  MultiComponentBoundTrackParameters(
      std::shared_ptr<const Surface> surface,
      const std::vector<std::tuple<double, ParametersVector, covariance_t>>&
          cmps,
      ParticleHypothesis particleHypothesis)
      : m_surface(std::move(surface)),
        m_particleHypothesis(particleHypothesis) {
    static_assert(std::is_same_v<BoundMatrix, covariance_t> ||
                  std::is_same_v<std::optional<BoundMatrix>, covariance_t>);
    if constexpr (std::is_same_v<BoundMatrix, covariance_t>) {
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
                                     const ParametersVector& params,
                                     std::optional<BoundMatrix> cov,
                                     ParticleHypothesis particleHypothesis)
      : m_surface(std::move(surface)),
        m_particleHypothesis(particleHypothesis) {
    m_components.push_back({1., params, std::move(cov)});
  }

  /// Parameters are not default constructible due to the charge type.
  MultiComponentBoundTrackParameters() = delete;
  /// Copy constructor
  MultiComponentBoundTrackParameters(
      const MultiComponentBoundTrackParameters&) = default;
  /// Move constructor
  MultiComponentBoundTrackParameters(MultiComponentBoundTrackParameters&&) =
      default;
  ~MultiComponentBoundTrackParameters() = default;
  /// Copy assignment operator
  /// @return Reference to this object after copying
  MultiComponentBoundTrackParameters& operator=(
      const MultiComponentBoundTrackParameters&) = default;
  /// Move assignment operator
  /// @return Reference to this object after moving
  MultiComponentBoundTrackParameters& operator=(
      MultiComponentBoundTrackParameters&&) = default;

  /// Comply with bound convertible, in this case return a copy
  /// @return Copy of this multi-component track parameters
  MultiComponentBoundTrackParameters toBound() const { return *this; }

  /// Access the parameters
  /// @return Reference to the vector of parameter components
  const auto& components() const { return m_components; }

  /// Reference surface onto which the parameters are bound.
  /// @return Reference to the bound reference surface
  const Surface& referenceSurface() const { return *m_surface; }

  /// Get the weight and a GenericBoundTrackParameters object for one component
  /// @param i Index of the component to access
  /// @return Pair of weight and bound track parameters for the component
  std::pair<double, Parameters> operator[](std::size_t i) const {
    return {
        std::get<double>(m_components[i]),
        Parameters(m_surface, std::get<ParametersVector>(m_components[i]),
                   std::get<std::optional<CovarianceMatrix>>(m_components[i]),
                   m_particleHypothesis)};
  }

  /// Convert to single component parameters by merging all components into one
  /// @return Single component bound track parameters representing the weighted average of all components
  BoundTrackParameters toSingleComponent() const;

  /// Parameters vector.
  /// @return Weighted average of parameters from all components
  ParametersVector parameters() const {
    return reduce([](const Parameters& p) { return p.parameters(); });
  }

  /// Optional covariance matrix.
  /// @return Optional weighted average covariance matrix, nullopt if all components have zero covariance
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
  /// @return Weighted average of the parameter at the specified index
  template <BoundIndices kIndex>
  double get() const {
    return reduce([&](const Parameters& p) { return p.get<kIndex>(); });
  }

  /// Space-time position four-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
  /// @return Weighted average four-dimensional position vector
  Vector4 fourPosition(const GeometryContext& geoCtx) const {
    return reduce([&](const Parameters& p) { return p.fourPosition(geoCtx); });
  }

  /// Spatial position three-vector.
  ///
  /// @param[in] geoCtx Geometry context for the local-to-global
  /// transformation
  /// @return Weighted average three-dimensional position vector
  Vector3 position(const GeometryContext& geoCtx) const {
    return reduce([&](const Parameters& p) { return p.position(geoCtx); });
  }

  /// Time coordinate.
  /// @return Weighted average time coordinate
  double time() const {
    return reduce([](const Parameters& p) { return p.time(); });
  }

  /// Unit direction three-vector, i.e. the normalized momentum
  /// three-vector.
  /// @return Weighted average normalized direction vector
  Vector3 direction() const {
    return reduce([](const Parameters& p) { return p.direction(); })
        .normalized();
  }

  /// Phi direction.
  /// @return Azimuthal angle phi derived from average direction
  double phi() const { return VectorHelpers::phi(direction()); }

  /// Theta direction.
  /// @return Polar angle theta derived from average direction
  double theta() const { return VectorHelpers::theta(direction()); }

  /// Charge over momentum.
  /// @return Weighted average charge over momentum ratio
  double qOverP() const { return get<eBoundQOverP>(); }

  /// Absolute momentum.
  /// @return Weighted average absolute momentum magnitude
  double absoluteMomentum() const {
    return reduce([](const Parameters& p) { return p.absoluteMomentum(); });
  }

  /// Transverse momentum.
  /// @return Weighted average transverse momentum magnitude
  double transverseMomentum() const {
    return reduce([](const Parameters& p) { return p.transverseMomentum(); });
  }

  /// Momentum three-vector.
  /// @return Weighted average three-dimensional momentum vector
  Vector3 momentum() const {
    return reduce([](const Parameters& p) { return p.momentum(); });
  }

  /// Particle electric charge.
  /// @return Weighted average particle electric charge
  double charge() const {
    return reduce([](const Parameters& p) { return p.charge(); });
  }

  /// Particle hypothesis.
  /// @return Reference to the particle hypothesis
  const ParticleHypothesis& particleHypothesis() const {
    return m_particleHypothesis;
  }
};

}  // namespace Acts
