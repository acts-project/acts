// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Direction.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/BoundTrackParameters.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/detail/periodic.hpp"

#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace Acts {
class IVolumeMaterial;
class MagneticFieldProvider;
}  // namespace Acts

namespace Acts::Experimental {

/// Collection of bound parameter variations, each entry is a pair of the index
/// of the bound parameter
using BoundParameterVariation = std::vector<std::pair<std::size_t, double>>;

/// Base class for generating bound parameter variations
struct BoundParameterVariationGenerator {
  virtual ~BoundParameterVariationGenerator() = default;
  /// Generate a variation map based on the input bound parameters and
  /// covariance
  /// @param input the input bound parameters
  /// @param covariance the covariance of the input bound parameters
  /// @return the corresponding parameter variation
  virtual BoundParameterVariation variationMap(
      const BoundVector& input, const BoundMatrix& covariance) const = 0;
};

/// Generator for bound parameter variations based on fixed deltas
struct DeltaBoundParameterVariationGenerator
    : public BoundParameterVariationGenerator {
  /// The deltas
  std::vector<BoundVector> deltas;

  /// Construct a generator with a single delta applied to all parameters
  /// @param delta the delta
  explicit DeltaBoundParameterVariationGenerator(const double delta)
      : deltas(1, BoundVector::Constant(delta)) {}

  /// Construct a generator with multiple deltas applied to all parameters
  /// @param deltas_ the delta
  explicit DeltaBoundParameterVariationGenerator(
      const std::vector<double>& deltas_) {
    this->deltas.reserve(deltas_.size());
    for (const auto& delta : deltas_) {
      this->deltas.emplace_back(BoundVector::Constant(delta));
    }
  }

  /// Construct a generator with a single delta vector
  /// @param delta the delta vector
  explicit DeltaBoundParameterVariationGenerator(const BoundVector& delta)
      : deltas(1, delta) {}

  /// Construct a generator with multiple delta vectors
  /// @param deltas_ the delta vectors
  explicit DeltaBoundParameterVariationGenerator(
      std::vector<BoundVector> deltas_)
      : deltas(std::move(deltas_)) {}

  /// Generate a variation map based on the fixed deltas
  /// @return the corresponding parameter variation
  BoundParameterVariation variationMap(
      const BoundVector& /*input*/,
      const BoundMatrix& /*covariance*/) const override {
    BoundParameterVariation variationMap;
    variationMap.reserve(deltas.size() * eBoundSize);
    for (std::size_t i = 0; i < deltas.size(); ++i) {
      for (std::size_t j = 0; j < eBoundSize; ++j) {
        variationMap.emplace_back(j, deltas[i][j]);
      }
    }
    return variationMap;
  }
};

/// Generator for bound parameter variations based on the covariance
struct CovarianceBoundParameterVariationGenerator
    : public BoundParameterVariationGenerator {
  /// The sigma factors
  std::vector<BoundVector> sigmaFactors;

  /// Construct a generator with a single sigma factor applied to all parameters
  /// @param sigmaFactor the sigma factor
  explicit CovarianceBoundParameterVariationGenerator(const double sigmaFactor)
      : sigmaFactors(1, BoundVector::Constant(sigmaFactor)) {}

  /// Construct a generator with multiple sigma factors applied to all
  /// parameters
  /// @param sigmaFactors_ the sigma factors
  explicit CovarianceBoundParameterVariationGenerator(
      const std::vector<double>& sigmaFactors_) {
    sigmaFactors.reserve(sigmaFactors_.size());
    for (const auto& factor : sigmaFactors_) {
      sigmaFactors.emplace_back(BoundVector::Constant(factor));
    }
  }

  /// Construct a generator with a single sigma factor vector
  /// @param sigmaFactor the sigma factor vector
  explicit CovarianceBoundParameterVariationGenerator(
      const BoundVector& sigmaFactor)
      : sigmaFactors(1, sigmaFactor) {}

  /// Construct a generator with multiple sigma factor vectors
  /// @param sigmaFactors_ the sigma factor vectors
  explicit CovarianceBoundParameterVariationGenerator(
      std::vector<BoundVector> sigmaFactors_)
      : sigmaFactors(std::move(sigmaFactors_)) {}

  /// Generate a variation map based on the covariance and sigma factors
  /// @param covariance the covariance of the input bound parameters
  /// @return the corresponding parameter variation
  BoundParameterVariation variationMap(
      const BoundVector& /*input*/,
      const BoundMatrix& covariance) const override {
    BoundParameterVariation variationMap;
    variationMap.reserve(sigmaFactors.size() * eBoundSize);
    for (std::size_t i = 0; i < sigmaFactors.size(); ++i) {
      for (std::size_t j = 0; j < eBoundSize; ++j) {
        variationMap.emplace_back(
            j, sigmaFactors[i][j] * std::sqrt(covariance(j, j)));
      }
    }
    return variationMap;
  }
};

/// RiddersStepper implements the Ridders method for numerical differentiation
/// to compute the Jacobian of the bound to bound transformation. It uses a
/// primary stepper to perform the nominal step and multiple secondary steppers
/// with varied initial conditions to compute the Jacobian columns. The
/// variations are generated based on the input covariance and a configurable
/// variation generator.
/// @tparam stepper_impl_t the type of the underlying stepper implementation
template <Concepts::SingleStepper stepper_impl_t>
class RiddersStepper final {
 public:
  /// The type of the underlying stepper implementation
  using StepperImpl = stepper_impl_t;

  /// The bound parameters type used by this stepper
  using BoundParameters = BoundTrackParameters;
  /// The Jacobian type for the bound to bound transformation
  using Jacobian = BoundMatrix;
  /// The covariance type for the bound parameters
  using Covariance = BoundMatrix;
  /// The type of the bound state returned by the `boundState` method
  using BoundState = std::tuple<BoundParameters, Jacobian, double>;

  /// Configuration struct for the RiddersStepper
  struct Config {
    /// Configuration for the underlying stepper implementation
    StepperImpl::Config stepperConfig;

    /// The generator for the bound parameter variations
    std::shared_ptr<const BoundParameterVariationGenerator> parameterVariation{
        std::make_shared<CovarianceBoundParameterVariationGenerator>(
            std::vector<double>{-1e-3, 1e-3})};

    /// The maximum number of steps to attempt when stepping towards the
    /// curvilinear surface for the bound state transformation
    std::size_t maxStepsToCurvilinearSurface = 10;
  };

  /// The options type for the RiddersStepper
  using Options = typename StepperImpl::Options;

  /// The state struct for the RiddersStepper
  struct State {
    /// @param primaryStepperStateIn the state of the primary stepper
    explicit State(const StepperImpl::State& primaryStepperStateIn)
        : primaryStepperState(primaryStepperStateIn) {}

    /// The state of the primary stepper, which performs the nominal step
    StepperImpl::State primaryStepperState;
    /// The states of the secondary steppers, which perform steps with varied
    /// initial conditions to compute the Jacobian columns
    std::vector<typename StepperImpl::State> secondaryStepperStates;

    /// The variation map, which contains the indices of the varied parameters
    /// and the corresponding deltas applied to them
    BoundParameterVariation variationMap;

    /// Flag indicating whether covariance transport is enabled
    bool covTransport = false;
    /// The covariance of the last known bound parameters
    Covariance cov = Covariance::Zero();
    /// The Jacobian of the last bound to bound transformation
    Jacobian jacobian = Jacobian::Identity();

    /// The accumulated path length during the stepping process
    double pathAccumulated = 0;

    /// Statistics for the stepper, such as the number of steps taken and the
    /// number of successful steps
    StepperStatistics statistics;

    // Parameters to reuse for curvilinear state
    /// The last propagation direction
    Direction lastStepPropagationDirection = Direction::Forward();
    /// The last material
    const IVolumeMaterial* lastStepMaterial = nullptr;
    /// The last surface tolerance
    double lastSurfaceTolerance = 0.;
    /// The last constrained step type
    ConstrainedStep::Type lastStepConstraintType =
        ConstrainedStep::Type::Navigator;
  };

  /// Construct a RiddersStepper with the given stepper implementation
  /// @param stepperImpl the underlying stepper implementation
  explicit RiddersStepper(StepperImpl stepperImpl)
      : m_stepperImpl(std::move(stepperImpl)) {}

  /// Construct a RiddersStepper with the given magnetic field provider
  /// @param bField the magnetic field provider
  explicit RiddersStepper(std::shared_ptr<const MagneticFieldProvider> bField)
      : m_stepperImpl(std::move(bField)) {}

  /// Construct a RiddersStepper with the given configuration
  /// @param config the configuration for the RiddersStepper
  explicit RiddersStepper(const Config& config)
      : m_config(config), m_stepperImpl(config.stepperConfig) {}

  /// Create a new state for the RiddersStepper based on the given options
  /// @param options the options
  /// @return a new state for the RiddersStepper
  State makeState(const Options& options) const {
    State state(m_stepperImpl.makeState(options));
    return state;
  }

  /// Initialize the state of the RiddersStepper based on the given bound
  /// parameters
  /// @param state the state
  /// @param boundParameters the bound parameters to initialize the state with
  void initialize(State& state, const BoundParameters& boundParameters) const {
    initialize(state, boundParameters.parameters(),
               boundParameters.covariance(),
               boundParameters.particleHypothesis(),
               boundParameters.referenceSurface());
  }

  /// Initialize the state of the RiddersStepper based on the given bound vector
  /// and covariance
  /// @param state the state
  /// @param boundVector the bound vector to initialize the state with
  /// @param covariance the covariance of the bound parameters, used to generate the variations for the secondary steppers
  /// @param particleHypothesis the particle hypothesis
  /// @param surface the surface on which the bound parameters are defined
  void initialize(State& state, const BoundVector& boundVector,
                  const std::optional<BoundMatrix>& covariance,
                  ParticleHypothesis particleHypothesis,
                  const Surface& surface) const {
    state.variationMap.clear();
    state.secondaryStepperStates.clear();

    m_stepperImpl.initialize(state.primaryStepperState, boundVector,
                             std::nullopt, particleHypothesis, surface);

    if (!covariance.has_value()) {
      state.covTransport = false;
      state.cov = BoundMatrix::Zero();
      // intentionally not resetting `jacobian` here as otherwise `update`
      // break. this should be revisitted at some point - might be best to
      // rework the stepper interface and how the jacobian is accessed

      return;
    }

    state.variationMap =
        m_config.parameterVariation->variationMap(boundVector, *covariance);

    for (const auto& [index, delta] : state.variationMap) {
      BoundVector nudgedParams = boundVector;
      nudgedParams[index] += delta;

      state.secondaryStepperStates.push_back(
          m_stepperImpl.makeState(state.primaryStepperState.options));
      m_stepperImpl.initialize(state.secondaryStepperStates.back(),
                               nudgedParams, std::nullopt, particleHypothesis,
                               surface);
    }

    state.covTransport = true;
    state.cov = *covariance;
    // same as above:
    // intentionally not resetting `jacobian` here as otherwise `update`
    // break. this should be revisitted at some point - might be best to
    // rework the stepper interface and how the jacobian is accessed
  }

  /// Get the magnetic field at the given position for the primary stepper state
  /// @param state the state of the RiddersStepper
  /// @param position the position at which to get the magnetic field
  /// @return the magnetic field at the given position for the primary stepper state
  Result<Vector3> getField(State& state, const Vector3& position) const {
    return m_stepperImpl.getField(state.primaryStepperState, position);
  }

  /// Get the position of the primary stepper state
  /// @param state the state of the RiddersStepper
  /// @return the position of the primary stepper state
  Vector3 position(const State& state) const {
    return m_stepperImpl.position(state.primaryStepperState);
  }

  /// Get the direction of the primary stepper state
  /// @param state the state of the RiddersStepper
  /// @return the direction of the primary stepper state
  Vector3 direction(const State& state) const {
    return m_stepperImpl.direction(state.primaryStepperState);
  }

  /// Get the charge over momentum of the primary stepper state
  /// @param state the state of the RiddersStepper
  /// @return the charge over momentum of the primary stepper state
  double qOverP(const State& state) const {
    return m_stepperImpl.qOverP(state.primaryStepperState);
  }

  /// Get the absolute momentum of the primary stepper state
  /// @param state the state of the RiddersStepper
  /// @return the absolute momentum of the primary stepper state
  double absoluteMomentum(const State& state) const {
    return m_stepperImpl.absoluteMomentum(state.primaryStepperState);
  }

  /// Get the momentum vector of the primary stepper state
  /// @param state the state of the RiddersStepper
  /// @return the momentum vector of the primary stepper state
  Vector3 momentum(const State& state) const {
    return m_stepperImpl.momentum(state.primaryStepperState);
  }

  /// Get the charge of the primary stepper state
  /// @param state the state of the RiddersStepper
  /// @return the charge of the primary stepper state
  double charge(const State& state) const {
    return m_stepperImpl.charge(state.primaryStepperState);
  }

  /// Get the particle hypothesis of the primary stepper state
  /// @param state the state of the RiddersStepper
  /// @return the particle hypothesis of the primary stepper state
  const ParticleHypothesis& particleHypothesis(const State& state) const {
    return m_stepperImpl.particleHypothesis(state.primaryStepperState);
  }

  /// Get the time of the primary stepper state
  /// @param state the state of the RiddersStepper
  /// @return the time of the primary stepper state
  double time(const State& state) const {
    return m_stepperImpl.time(state.primaryStepperState);
  }

  /// Update the surface status of the primary and secondary stepper states
  /// based on the given surface and navigation direction
  /// @param state the state of the RiddersStepper
  /// @param surface the surface
  /// @param index the index of the surface
  /// @param navDir the navigation direction
  /// @param boundaryTolerance the boundary tolerance
  /// @param surfaceTolerance the surface tolerance
  /// @param stype the type of the constrained step
  /// @param logger the logger
  /// @return the overall intersection status after updating the surface
  ///    status for the primary and secondary stepper states
  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction navDir, const BoundaryTolerance& boundaryTolerance,
      double surfaceTolerance, ConstrainedStep::Type stype,
      const Logger& logger = getDummyLogger()) const {
    state.lastSurfaceTolerance = surfaceTolerance;
    state.lastStepConstraintType = stype;

    IntersectionStatus overallStatus = m_stepperImpl.updateSurfaceStatus(
        state.primaryStepperState, surface, index, navDir, boundaryTolerance,
        surfaceTolerance, stype, logger);
    if (overallStatus == IntersectionStatus::unreachable) {
      return overallStatus;
    }

    const auto combine = [](const IntersectionStatus a,
                            const IntersectionStatus b) {
      using enum IntersectionStatus;
      if ((a == unreachable) || (b == unreachable)) {
        return unreachable;
      }
      if ((a == reachable) || (b == reachable)) {
        return reachable;
      }
      return onSurface;
    };

    for (auto& secondaryState : state.secondaryStepperStates) {
      const IntersectionStatus status = m_stepperImpl.updateSurfaceStatus(
          secondaryState, surface, index, navDir, BoundaryTolerance::Infinite(),
          surfaceTolerance, stype, logger);
      overallStatus = combine(overallStatus, status);
    }

    return overallStatus;
  }

  /// Update the step size based on the given navigation target, direction, and
  /// constrained step type
  /// @param state the state of the RiddersStepper
  /// @param target the navigation target
  /// @param direction the direction
  /// @param stype the type of the constrained step
  void updateStepSize(State& state, const NavigationTarget& target,
                      Direction direction, ConstrainedStep::Type stype) const {
    m_stepperImpl.updateStepSize(state.primaryStepperState, target, direction,
                                 stype);

    for (auto& secondaryState : state.secondaryStepperStates) {
      m_stepperImpl.updateStepSize(secondaryState, target, direction, stype);
    }
  }

  /// Update the step size based on the given step size and constrained step
  /// type
  /// @param state the state of the RiddersStepper
  /// @param stepSize the step size
  /// @param stype the type of the constrained step
  void updateStepSize(State& state, double stepSize,
                      ConstrainedStep::Type stype) const {
    m_stepperImpl.updateStepSize(state.primaryStepperState, stepSize, stype);

    for (auto& secondaryState : state.secondaryStepperStates) {
      m_stepperImpl.updateStepSize(secondaryState, stepSize, stype);
    }
  }

  /// Get the step size for the given constrained step type
  /// @param state the state of the RiddersStepper
  /// @param stype the type of the constrained step for which to get the step size
  /// @return the step size of the primary stepper state for the given constrained step type
  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return m_stepperImpl.getStepSize(state.primaryStepperState, stype);
  }

  /// Release the step size for the given constrained step type
  /// @param state the state of the RiddersStepper
  /// @param stype the type of the constrained step for which to release the step size
  void releaseStepSize(State& state, ConstrainedStep::Type stype) const {
    m_stepperImpl.releaseStepSize(state.primaryStepperState, stype);

    for (auto& secondaryState : state.secondaryStepperStates) {
      m_stepperImpl.releaseStepSize(secondaryState, stype);
    }
  }

  /// Get a string representation of the step size constraints
  /// @param state the state of the RiddersStepper
  /// @return a string representation of the step size constraints
  std::string outputStepSize(const State& state) const {
    std::stringstream ss;
    ss << state.primaryStepperState.stepSize.toString();
    ss << " || ";
    for (const auto& secondaryState : state.secondaryStepperStates) {
      ss << secondaryState.stepSize.toString() << " | ";
    }
    return ss.str();
  }

  /// Compute the bound state
  /// @param state the state of the RiddersStepper
  /// @param surface the surface
  /// @param transportCovariance flag indicating whether to transport the covariance to the bound parameters
  /// @param freeToBoundCorrection the correction
  /// @return the bound state
  Result<BoundState> boundState(
      State& state, const Surface& surface, bool transportCovariance = true,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const {
    Result<BoundState> primaryResult = m_stepperImpl.boundState(
        state.primaryStepperState, surface, false, freeToBoundCorrection);

    if (!transportCovariance) {
      std::get<1>(*primaryResult) = state.jacobian;
      return primaryResult;
    }

    std::vector<BoundVector> boundVectors;
    std::vector<double> pathLenghts;
    boundVectors.reserve(state.secondaryStepperStates.size());
    pathLenghts.reserve(state.secondaryStepperStates.size());
    for (auto& secondaryState : state.secondaryStepperStates) {
      const Result<BoundState> result = m_stepperImpl.boundState(
          secondaryState, surface, false, freeToBoundCorrection);
      if (!result.ok()) {
        return result.error();
      }
      boundVectors.push_back(std::get<0>(*result).parameters());
      pathLenghts.push_back(std::get<2>(*result));
    }
    const auto& referenceVector = std::get<0>(*primaryResult).parameters();

    state.jacobian = BoundMatrix::Zero();
    std::array<std::size_t, eBoundSize> numberOfVariationsPerParameter{};

    for (std::size_t i = 0; i < state.variationMap.size(); ++i) {
      const auto& [index, delta] = state.variationMap[i];
      const auto& nudgedVector = boundVectors[i];

      const auto diff = [&]() {
        BoundVector d = nudgedVector - referenceVector;
        d[eBoundPhi] = detail::difference_periodic(nudgedVector[eBoundPhi],
                                                   referenceVector[eBoundPhi],
                                                   2 * std::numbers::pi);
        return d;
      };

      state.jacobian.col(index) += diff() / delta;
      ++numberOfVariationsPerParameter[index];
    }
    for (std::size_t i = 0; i < eBoundSize; ++i) {
      state.jacobian.col(i) /= numberOfVariationsPerParameter[i];
    }

    state.cov = state.jacobian * state.cov * state.jacobian.transpose();

    BoundTrackParameters newParams(surface.getSharedPtr(), referenceVector,
                                   state.cov, particleHypothesis(state));

    state.variationMap =
        m_config.parameterVariation->variationMap(referenceVector, state.cov);

    return Result<BoundState>::success(
        BoundState{std::move(newParams), state.jacobian, pathLenghts.front()});
  }

  /// Prepare the stepper for the curvilinear transformation
  /// @param state the state of the RiddersStepper
  /// @return true if the preparation was successful for all stepper states, false otherwise
  bool prepareCurvilinearState(State& state) const {
    bool result =
        m_stepperImpl.prepareCurvilinearState(state.primaryStepperState);
    for (auto& secondaryState : state.secondaryStepperStates) {
      result = result && m_stepperImpl.prepareCurvilinearState(secondaryState);
    }
    return result;
  }

  /// Compute the curvilinear state
  /// @param state the state of the RiddersStepper
  /// @param transportCovariance flag indicating whether to transport the covariance to the curvilinear parameters
  /// @return the curvilinear state
  BoundState curvilinearState(State& state,
                              bool transportCovariance = true) const {
    if (!transportCovariance) {
      return m_stepperImpl.curvilinearState(state.primaryStepperState, false);
    }

    const auto curvilinearSurface =
        CurvilinearSurface(position(state), direction(state)).surface();

    for (std::size_t i = 0; i < m_config.maxStepsToCurvilinearSurface; ++i) {
      const auto surfaceStatus = updateSurfaceStatus(
          state, *curvilinearSurface, 0, state.lastStepPropagationDirection,
          BoundaryTolerance::Infinite(), state.lastSurfaceTolerance,
          state.lastStepConstraintType);
      if (surfaceStatus == IntersectionStatus::onSurface) {
        break;
      }
      if (surfaceStatus == IntersectionStatus::unreachable) {
        throw std::runtime_error(
            "RiddersStepper: Curvilinear surface is unreachable");
      }
      const auto stepResult = step(state, state.lastStepPropagationDirection,
                                   state.lastStepMaterial);
      if (!stepResult.ok()) {
        throw std::runtime_error(
            "RiddersStepper: Failed to step towards curvilinear surface: " +
            stepResult.error().message());
      }
    }

    return boundState(state, *curvilinearSurface, true).value();
  }

  /// Update the state of the RiddersStepper
  /// @param state the state of the RiddersStepper
  /// @param boundParams the bound vector
  /// @param covariance the covariance of the bound parameters
  /// @param surface the surface on which the bound parameters are defined
  void update(State& state, const FreeVector& /*freeParams*/,
              const BoundVector& boundParams, const BoundMatrix& covariance,
              const Surface& surface) const {
    initialize(state, boundParams, covariance, particleHypothesis(state),
               surface);
  }

  /// Update the state of the RiddersStepper
  /// @param state the state of the RiddersStepper
  /// @param position the position
  /// @param direction the direction
  /// @param qOverP the charge over momentum
  /// @param time the time
  void update(State& state, const Vector3& position, const Vector3& direction,
              double qOverP, double time) const {
    const auto curvilinearSurface =
        CurvilinearSurface(position, direction).surface();
    const std::optional<BoundMatrix> cov =
        state.covTransport ? std::optional<BoundMatrix>(state.cov)
                           : std::nullopt;
    initialize(state,
               transformFreeToBoundParameters(
                   position, time, direction, qOverP, *curvilinearSurface,
                   state.primaryStepperState.options.geoContext)
                   .value(),
               cov, particleHypothesis(state), *curvilinearSurface);
  }

  /// Transport the covariance to curvilinear parameters
  /// @param state the state of the RiddersStepper
  void transportCovarianceToCurvilinear(State& state) const {
    curvilinearState(state, true);
  }

  /// Transport the covariance to bound parameters
  /// @param state the state of the RiddersStepper
  /// @param surface the surface
  /// @param freeToBoundCorrection the correction
  Result<void> transportCovarianceToBound(
      State& state, const Surface& surface,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const {
    Result<BoundState> result =
        boundState(state, surface, true, freeToBoundCorrection);
    if (!result.ok()) {
      return result.error();
    }
    return Result<void>::success();
  }

  /// Perform a step
  /// @param state the state of the RiddersStepper
  /// @param propagationDirection the direction
  /// @param material the material
  /// @return a result containing the step length or an error if the step failed
  Result<double> step(State& state, Direction propagationDirection,
                      const IVolumeMaterial* material) const {
    state.lastStepPropagationDirection = propagationDirection;
    state.lastStepMaterial = material;

    const Result<double> stepResult = m_stepperImpl.step(
        state.primaryStepperState, propagationDirection, material);
    if (!stepResult.ok()) {
      return stepResult.error();
    }
    for (auto& secondaryState : state.secondaryStepperStates) {
      const Result<double> secondaryStepResult =
          m_stepperImpl.step(secondaryState, propagationDirection, material);
      if (!secondaryStepResult.ok()) {
        return secondaryStepResult.error();
      }
    }
    state.pathAccumulated += *stepResult;
    return *stepResult;
  }

 private:
  /// The stepper configuration
  Config m_config;

  /// The underlying stepper implementation
  StepperImpl m_stepperImpl;
};

}  // namespace Acts::Experimental
