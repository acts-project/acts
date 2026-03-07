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

enum class ComponentMergeMethod;

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
  /// Type alias for the component tuple
  using Component =
      std::tuple<double, ParametersVector, std::optional<CovarianceMatrix>>;

 private:
  std::vector<double> m_weights;
  std::vector<ParametersVector> m_parameters;
  std::vector<CovarianceMatrix> m_covariances;
  std::shared_ptr<const Surface> m_surface;
  ParticleHypothesis m_particleHypothesis;

 public:
  /// Type alias for construction tuple containing weight, position, direction,
  /// q/p, and covariance
  using ConstructionTuple =
      std::tuple<double, Vector4, Vector3, double, CovarianceMatrix>;

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
    Vector3 avgPos = Vector3::Zero();
    Vector3 avgDir = Vector3::Zero();
    for (const auto& [w, pos4, dir, qop, cov] : curvi) {
      avgPos += w * pos4.template segment<3>(0);
      avgDir += w * dir;
    }

    std::shared_ptr<PlaneSurface> s =
        CurvilinearSurface(avgPos, avgDir).planeSurface();

    std::vector<Component> bound;
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
    static_assert(
        std::is_same_v<CovarianceMatrix, covariance_t> ||
        std::is_same_v<std::optional<CovarianceMatrix>, covariance_t>);

    if (cmps.empty()) {
      throw std::invalid_argument(
          "Cannot construct MultiComponentBoundTrackParameters with empty "
          "components");
    }

    for (const auto& [weight, params, cov] : cmps) {
      m_weights.push_back(weight);
      m_parameters.push_back(params);
      if constexpr (std::is_same_v<CovarianceMatrix, covariance_t>) {
        m_covariances.push_back(cov);
      } else if (std::is_same_v<std::optional<CovarianceMatrix>,
                                covariance_t>) {
        if (cov.has_value()) {
          m_covariances.push_back(*cov);
        }
      }
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
                                     std::optional<CovarianceMatrix> cov,
                                     ParticleHypothesis particleHypothesis)
      : m_surface(std::move(surface)),
        m_particleHypothesis(particleHypothesis) {
    m_weights.push_back(1.);
    m_parameters.push_back(params);
    if (cov) {
      m_covariances.push_back(*cov);
    }
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

  /// Size of the multi-component parameters
  /// @return Number of components in the multi-component parameters
  std::size_t size() const { return m_weights.size(); }

  /// Access to the weights of the components
  /// @return Span of the weights for all components
  std::span<const double> weights() const { return m_weights; }

  /// Access to the parameters of the components
  /// @return Span of the parameters vectors for all components
  std::span<const ParametersVector> parameters() const { return m_parameters; }

  /// Access to the covariances of the components
  /// @return Span of the covariance matrices for all components
  std::span<const CovarianceMatrix> covariances() const {
    return m_covariances;
  }

  /// Reference surface onto which the parameters are bound.
  /// @return Reference to the bound reference surface
  const Surface& referenceSurface() const { return *m_surface; }

  /// Particle hypothesis.
  /// @return Reference to the particle hypothesis
  const ParticleHypothesis& particleHypothesis() const {
    return m_particleHypothesis;
  }

  /// Check if covariance matrices are available for the components
  /// @return True if covariance matrices are available, false otherwise
  bool hasCovariance() const { return !m_covariances.empty(); }

  /// Convert the multi-component parameters to a vector of components
  /// @return Vector of components, where each component is a tuple of weight, parameters vector, and optional covariance matrix
  std::vector<Component> toComponents() const {
    std::vector<Component> cmps;
    cmps.reserve(size());
    for (std::size_t i = 0; i < size(); ++i) {
      cmps.emplace_back(
          m_weights[i], m_parameters[i],
          hasCovariance() ? std::optional(m_covariances[i]) : std::nullopt);
    }
    return cmps;
  }

  /// Get the weight and a GenericBoundTrackParameters object for one component
  /// @param i Index of the component to access
  /// @return Pair of weight and bound track parameters for the component
  std::pair<double, Parameters> operator[](std::size_t i) const {
    return {m_weights[i],
            Parameters(m_surface, m_parameters[i],
                       hasCovariance() ? std::optional(m_covariances[i])
                                       : std::nullopt,
                       m_particleHypothesis)};
  }

  /// Merge to single component parameters
  /// @return Single component bound track parameters representing the weighted average of all components
  BoundTrackParameters merge(ComponentMergeMethod method) const;
};

}  // namespace Acts
