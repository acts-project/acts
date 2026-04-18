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
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"
#include "Acts/Utilities/Intersection.hpp"

#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace Acts {
class IVolumeMaterial;
class MagneticFieldProvider;
}  // namespace Acts

namespace Acts::Experimental {

using BoundParameterVariation = std::vector<std::pair<std::size_t, double>>;

struct BoundParameterVariationGenerator {
  virtual ~BoundParameterVariationGenerator() = default;
  virtual BoundParameterVariation variationMap(
      const BoundVector& input, const BoundMatrix& covariance) const = 0;
};

struct DeltaBoundParameterVariationGenerator
    : public BoundParameterVariationGenerator {
  std::vector<BoundVector> deltas;

  explicit DeltaBoundParameterVariationGenerator(double delta)
      : deltas(1, BoundVector::Constant(delta)) {}
  explicit DeltaBoundParameterVariationGenerator(
      const std::vector<double>& deltas_) {
    this->deltas.reserve(deltas_.size());
    for (const auto& delta : deltas_) {
      this->deltas.emplace_back(BoundVector::Constant(delta));
    }
  }
  explicit DeltaBoundParameterVariationGenerator(const BoundVector& delta)
      : deltas(1, delta) {}
  explicit DeltaBoundParameterVariationGenerator(
      std::vector<BoundVector> deltas_)
      : deltas(std::move(deltas_)) {}

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

struct CovarianceBoundParameterVariationGenerator
    : public BoundParameterVariationGenerator {
  std::vector<BoundVector> sigmaFactors;

  explicit CovarianceBoundParameterVariationGenerator(double sigmaFactor)
      : sigmaFactors(1, BoundVector::Constant(sigmaFactor)) {}
  explicit CovarianceBoundParameterVariationGenerator(
      const std::vector<double>& sigmaFactors_) {
    sigmaFactors.reserve(sigmaFactors_.size());
    for (const auto& factor : sigmaFactors_) {
      sigmaFactors.emplace_back(BoundVector::Constant(factor));
    }
  }
  explicit CovarianceBoundParameterVariationGenerator(
      const BoundVector& sigmaFactor)
      : sigmaFactors(1, sigmaFactor) {}
  explicit CovarianceBoundParameterVariationGenerator(
      std::vector<BoundVector> sigmaFactors_)
      : sigmaFactors(std::move(sigmaFactors_)) {}

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

template <Concepts::SingleStepper stepper_impl_t>
class RiddersStepper final {
 public:
  using StepperImpl = stepper_impl_t;

  using BoundParameters = BoundTrackParameters;
  using Jacobian = BoundMatrix;
  using Covariance = BoundMatrix;
  using BoundState = std::tuple<BoundParameters, Jacobian, double>;

  struct Config {
    StepperImpl::Config stepperConfig;

    std::shared_ptr<const BoundParameterVariationGenerator> parameterVariation{
        std::make_shared<CovarianceBoundParameterVariationGenerator>(
            std::vector<double>{-1e-1, 1e-1})};
  };

  using Options = typename StepperImpl::Options;

  struct State {
    explicit State(const StepperImpl::State& primaryStepperStateIn)
        : primaryStepperState(primaryStepperStateIn) {}

    StepperImpl::State primaryStepperState;
    std::vector<typename StepperImpl::State> secondaryStepperStates;

    BoundParameterVariation variationMap;

    bool covTransport = false;
    Covariance cov = Covariance::Zero();
    Jacobian jacobian = Jacobian::Identity();

    double pathAccumulated = 0.;

    StepperStatistics statistics;
  };

  explicit RiddersStepper(std::shared_ptr<const MagneticFieldProvider> bField)
      : m_stepperImpl(std::move(bField)) {}

  explicit RiddersStepper(const Config& config)
      : m_config(config), m_stepperImpl(config.stepperConfig) {}

  State makeState(const Options& options) const {
    State state(m_stepperImpl.makeState(options));
    return state;
  }

  void initialize(State& state, const BoundParameters& boundParameters) const {
    initialize(state, boundParameters.parameters(),
               boundParameters.covariance(),
               boundParameters.particleHypothesis(),
               boundParameters.referenceSurface());
  }

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
      // TODO otherwise `update` gets screwed
      // state.jacobian = BoundMatrix::Identity();

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
    // TODO otherwise `update` gets screwed
    // state.jacobian = BoundMatrix::Identity();
  }

  Result<Vector3> getField(State& state, const Vector3& position) const {
    return m_stepperImpl.getField(state.primaryStepperState, position);
  }

  Vector3 position(const State& state) const {
    return m_stepperImpl.position(state.primaryStepperState);
  }

  Vector3 direction(const State& state) const {
    return m_stepperImpl.direction(state.primaryStepperState);
  }

  double qOverP(const State& state) const {
    return m_stepperImpl.qOverP(state.primaryStepperState);
  }

  double absoluteMomentum(const State& state) const {
    return m_stepperImpl.absoluteMomentum(state.primaryStepperState);
  }

  Vector3 momentum(const State& state) const {
    return m_stepperImpl.momentum(state.primaryStepperState);
  }

  double charge(const State& state) const {
    return m_stepperImpl.charge(state.primaryStepperState);
  }

  const ParticleHypothesis& particleHypothesis(const State& state) const {
    return m_stepperImpl.particleHypothesis(state.primaryStepperState);
  }

  double time(const State& state) const {
    return m_stepperImpl.time(state.primaryStepperState);
  }

  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction navDir, const BoundaryTolerance& boundaryTolerance,
      double surfaceTolerance, ConstrainedStep::Type stype,
      const Logger& logger = getDummyLogger()) const {
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

  void updateStepSize(State& state, const NavigationTarget& target,
                      Direction direction, ConstrainedStep::Type stype) const {
    m_stepperImpl.updateStepSize(state.primaryStepperState, target, direction,
                                 stype);

    for (auto& secondaryState : state.secondaryStepperStates) {
      m_stepperImpl.updateStepSize(secondaryState, target, direction, stype);
    }
  }

  void updateStepSize(State& state, double stepSize,
                      ConstrainedStep::Type stype) const {
    m_stepperImpl.updateStepSize(state.primaryStepperState, stepSize, stype);

    for (auto& secondaryState : state.secondaryStepperStates) {
      m_stepperImpl.updateStepSize(secondaryState, stepSize, stype);
    }
  }

  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return m_stepperImpl.getStepSize(state.primaryStepperState, stype);
  }

  void releaseStepSize(State& state, ConstrainedStep::Type stype) const {
    m_stepperImpl.releaseStepSize(state.primaryStepperState, stype);

    for (auto& secondaryState : state.secondaryStepperStates) {
      m_stepperImpl.releaseStepSize(secondaryState, stype);
    }
  }

  std::string outputStepSize(const State& state) const {
    std::stringstream ss;
    ss << state.primaryStepperState.stepSize.toString();
    for (const auto& secondaryState : state.secondaryStepperStates) {
      ss << " || " << secondaryState.stepSize.toString();
    }
    return ss.str();
  }

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

      state.jacobian.col(index) += (nudgedVector - referenceVector) / delta;
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

  bool prepareCurvilinearState(State& state) const {
    bool result =
        m_stepperImpl.prepareCurvilinearState(state.primaryStepperState);
    for (auto& secondaryState : state.secondaryStepperStates) {
      result = result && m_stepperImpl.prepareCurvilinearState(secondaryState);
    }
    return result;
  }

  BoundState curvilinearState(State& state,
                              bool transportCovariance = true) const {
    if (!transportCovariance) {
      return m_stepperImpl.curvilinearState(state.primaryStepperState, false);
    }
    return boundState(
               state,
               *CurvilinearSurface(position(state), direction(state)).surface(),
               true)
        .value();
  }

  void update(State& state, const FreeVector& /*freeParams*/,
              const BoundVector& boundParams, const BoundMatrix& covariance,
              const Surface& surface) const {
    initialize(state, boundParams, covariance, particleHypothesis(state),
               surface);
  }

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

  void transportCovarianceToCurvilinear(State& state) const {
    curvilinearState(state, true);
  }

  void transportCovarianceToBound(
      State& state, const Surface& surface,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const {
    boundState(state, surface, true, freeToBoundCorrection);
  }

  Result<double> step(State& state, Direction propagationDirection,
                      const IVolumeMaterial* material) const {
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

  void setIdentityJacobian(State& state) const {
    m_stepperImpl.setIdentityJacobian(state.primaryStepperState);
    for (auto& secondaryState : state.secondaryStepperStates) {
      m_stepperImpl.setIdentityJacobian(secondaryState);
    }
  }

 private:
  Config m_config;

  StepperImpl m_stepperImpl;
};

}  // namespace Acts::Experimental
