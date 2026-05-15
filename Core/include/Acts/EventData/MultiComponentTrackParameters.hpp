// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/detail/MultiComponentTrackParametersConcept.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <memory>
#include <ranges>
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
  bool m_hasCovariance{false};

 public:
  /// Type alias for construction tuple containing weight, position, direction,
  /// q/p, and covariance
  using ConstructionTuple =
      std::tuple<double, Vector4, Vector3, double, BoundMatrix>;

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
      const Intersection3D closestIntersection =
          s->intersect(geoCtx, pos4.template segment<3>(eFreePos0), dir,
                       BoundaryTolerance::Infinite())
              .closest();
      const Vector3& newPos = closestIntersection.position();

      BoundVector bv =
          transformFreeToCurvilinearParameters(pos4[eTime], dir, qop);

      // Because of the projection this should never fail
      bv.template segment<2>(eBoundLoc0) =
          *(s->globalToLocal(geoCtx, newPos, dir));

      bound.emplace_back(w, bv, cov);
    }

    return MultiComponentBoundTrackParameters(
        s, bound,
        [](const Component& cmp) {
          return std::tie(std::get<0>(cmp), std::get<1>(cmp),
                          *std::get<2>(cmp));
        },
        particleHypothesis);
  }

  /// Construct from a surface and particle hypothesis, without components
  /// @param surface Reference surface the parameters are defined on
  /// @param hasCovariance Flag indicating if covariance matrices will be provided for the
  /// @param particleHypothesis Particle hypothesis for the parameters
  MultiComponentBoundTrackParameters(std::shared_ptr<const Surface> surface,
                                     bool hasCovariance,
                                     ParticleHypothesis particleHypothesis)
      : m_surface(std::move(surface)),
        m_particleHypothesis(particleHypothesis),
        m_hasCovariance(hasCovariance) {}

  /// Construct from a single component
  /// @param surface Reference surface the parameters are defined on
  /// @param params Bound parameters vector
  /// @param cov Bound parameters covariance matrix
  /// @param particleHypothesis Particle hypothesis for these parameters
  MultiComponentBoundTrackParameters(std::shared_ptr<const Surface> surface,
                                     const BoundVector& params,
                                     std::optional<BoundMatrix> cov,
                                     ParticleHypothesis particleHypothesis)
      : m_surface(std::move(surface)),
        m_particleHypothesis(particleHypothesis) {
    m_weights.push_back(1.);
    m_parameters.push_back(params);
    if (cov) {
      m_covariances.push_back(*cov);
    }
    m_hasCovariance = cov.has_value();
  }

  /// Construct from multiple components without covariance
  /// @param surface Surface on which the parameters are bound
  /// @param cmps Vector of weight, parameters vector, and covariance components
  /// @param proj Projector to use for the parameters
  /// @param particleHypothesis Particle hypothesis for the parameters
  template <std::ranges::range component_range_t, typename projector_t>
  MultiComponentBoundTrackParameters(std::shared_ptr<const Surface> surface,
                                     const component_range_t& cmps,
                                     const projector_t& proj,
                                     ParticleHypothesis particleHypothesis)
    requires detail::ComponentRangeAndProjectorWithoutCovarianceConcept<
                 component_range_t, projector_t>
      : m_surface(std::move(surface)),
        m_particleHypothesis(particleHypothesis) {
    if (cmps.begin() == cmps.end()) {
      throw std::invalid_argument(
          "Cannot construct MultiComponentBoundTrackParameters without "
          "components");
    }

    for (const auto& cmp : cmps) {
      const auto& [weight, parameters] = proj(cmp);
      m_weights.push_back(weight);
      m_parameters.push_back(parameters);
    }
  }

  /// Construct from multiple components with covariance
  /// @param surface Surface on which the parameters are bound
  /// @param cmps Vector of weight, parameters vector, and covariance components
  /// @param proj Projector to use for the parameters
  /// @param particleHypothesis Particle hypothesis for the parameters
  template <std::ranges::range component_range_t, typename projector_t>
  MultiComponentBoundTrackParameters(std::shared_ptr<const Surface> surface,
                                     const component_range_t& cmps,
                                     const projector_t& proj,
                                     ParticleHypothesis particleHypothesis)
    requires detail::ComponentRangeAndProjectorWithCovarianceConcept<
                 component_range_t, projector_t>
      : m_surface(std::move(surface)),
        m_particleHypothesis(particleHypothesis),
        m_hasCovariance(true) {
    if (cmps.begin() == cmps.end()) {
      throw std::invalid_argument(
          "Cannot construct MultiComponentBoundTrackParameters without "
          "components");
    }

    for (const auto& cmp : cmps) {
      const auto& [weight, parameters, covariance] = proj(cmp);
      m_weights.push_back(weight);
      m_parameters.push_back(parameters);
      m_covariances.push_back(covariance);
    }
  }

  /// No default constructor, because we always need at least a surface and
  /// particle hypothesis
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
  [[deprecated(
      "You already have a universal bound track parameter at hand. You can "
      "drop `toBound()`.")]]
  MultiComponentBoundTrackParameters toBound() const {
    return *this;
  }

  /// Size of the multi-component parameters
  /// @return Number of components in the multi-component parameters
  std::size_t size() const { return m_weights.size(); }

  /// Check if the multi-component parameters are empty
  /// @return True if there are no components, false otherwise
  bool empty() const { return m_weights.empty(); }

  /// Access to the weights of the components
  /// @return The weights for all components
  const std::vector<double>& weights() const { return m_weights; }

  /// Access to the parameters of the components
  /// @return The parameters vectors for all components
  const std::vector<ParametersVector>& parameters() const {
    return m_parameters;
  }

  /// Access to the covariances of the components
  /// @return The covariance matrices for all components
  /// @throws std::runtime_error if covariance matrices are not available for this multi-component parameters
  const std::vector<CovarianceMatrix>& covariances() const {
    if (!hasCovariance()) {
      throw std::runtime_error(
          "Covariance matrices are not available for this "
          "MultiComponentBoundTrackParameters");
    }
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
  bool hasCovariance() const { return m_hasCovariance; }

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

  /// Merge component mixture into a single set of parameters using the
  /// specified method.
  /// @param method Method to use for merging the components into a single set of parameters
  /// @return Single component bound track parameters representing the mixture
  BoundTrackParameters merge(ComponentMergeMethod method) const;

  /// Clear all components from the multi-component parameters
  void clear();

  /// Reserve space for a number of components in the multi-component parameters
  /// @param n Number of components to reserve space for
  void reserve(std::size_t n);

  /// Add a component to the multi-component parameters
  /// @param weight Weight of the new component
  /// @param params Parameters vector of the new component
  void pushComponent(double weight, const ParametersVector& params);

  /// Add a component with covariance to the multi-component parameters
  /// @param weight Weight of the new component
  /// @param params Parameters vector of the new component
  /// @param cov Covariance matrix of the new component
  void pushComponent(double weight, const ParametersVector& params,
                     const CovarianceMatrix& cov);

  /// Add a component with optional covariance to the multi-component parameters
  /// @param weight Weight of the new component
  /// @param params Parameters vector of the new component
  /// @param cov Optional covariance matrix of the new component
  void pushComponent(double weight, const ParametersVector& params,
                     const std::optional<CovarianceMatrix>& cov);

  /// Normalize the weights of the components so that they sum up to 1
  void normalizeWeights();
};

}  // namespace Acts
