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
#include "Acts/Propagator/MultiStepperLoop.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Propagator/StepperConcept.hpp"
#include "Acts/Propagator/StepperStatistics.hpp"
#include "Acts/Surfaces/CurvilinearSurface.hpp"

#include <iostream>
#include <tuple>

namespace Acts {
class IVolumeMaterial;
}

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

template <Concepts::SingleStepper stepper_impl_t,
          Concepts::MultiStepper multi_stepper_impl_t =
              MultiStepperLoop<stepper_impl_t>>
class RiddersStepper final {
 public:
  using StepperImpl = stepper_impl_t;
  using MultiStepperImpl = multi_stepper_impl_t;

  using BoundParameters = BoundTrackParameters;
  using Jacobian = BoundMatrix;
  using Covariance = BoundMatrix;
  using BoundState = std::tuple<BoundParameters, Jacobian, double>;

  using MultiBoundParameters = MultiStepperImpl::BoundParameters;
  using MultiBoundState = std::tuple<MultiBoundParameters, Jacobian, double>;

  struct Config {
    MultiStepperImpl::Config multiStepperConfig;

    std::shared_ptr<const BoundParameterVariationGenerator> parameterVariation{
        std::make_shared<DeltaBoundParameterVariationGenerator>(1e-4)};
  };

  using Options = typename MultiStepperImpl::Options;

  struct State {
    explicit State(const MultiStepperImpl::State& multiStepperStateIn)
        : multiStepperState(multiStepperStateIn) {}

    MultiStepperImpl::State multiStepperState;

    BoundParameterVariation variationMap;

    bool covTransport = false;
    Covariance cov = Covariance::Zero();
    Jacobian jacobian = Jacobian::Identity();

    double pathAccumulated = 0.;

    StepperStatistics statistics;
  };

  explicit RiddersStepper(std::shared_ptr<const MagneticFieldProvider> bField)
      : m_multiStepperImpl(std::move(bField)) {}

  explicit RiddersStepper(const Config& config)
      : m_config(config), m_multiStepperImpl(config.multiStepperConfig) {}

  State makeState(const Options& options) const {
    State state(m_multiStepperImpl.makeState(options));
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

    if (!covariance.has_value()) {
      state.covTransport = false;
      state.cov = BoundMatrix::Zero();
      // TODO otherwise `update` gets screwed
      // state.jacobian = BoundMatrix::Identity();

      m_multiStepperImpl.initialize(
          state.multiStepperState,
          MultiBoundParameters(surface.getSharedPtr(), boundVector,
                               std::nullopt, particleHypothesis));
      return;
    }

    state.variationMap =
        m_config.parameterVariation->variationMap(boundVector, *covariance);

    MultiBoundParameters multiBoundParameters(surface.getSharedPtr(), false,
                                              particleHypothesis);
    multiBoundParameters.reserve(1 + state.variationMap.size());

    multiBoundParameters.pushComponent(1., boundVector);

    for (const auto& [index, delta] : state.variationMap) {
      BoundVector nudgedParams = boundVector;
      nudgedParams[index] += delta;
      multiBoundParameters.pushComponent(0., nudgedParams);
    }

    state.covTransport = true;
    state.cov = *covariance;
    // TODO otherwise `update` gets screwed
    // state.jacobian = BoundMatrix::Identity();

    m_multiStepperImpl.initialize(state.multiStepperState,
                                  multiBoundParameters);
  }

  Result<Vector3> getField(State& state, const Vector3& position) const {
    return singleStepper().getField(primaryState(state), position);
  }

  Vector3 position(const State& state) const {
    return singleStepper().position(primaryState(state));
  }

  Vector3 direction(const State& state) const {
    return singleStepper().direction(primaryState(state));
  }

  double qOverP(const State& state) const {
    return singleStepper().qOverP(primaryState(state));
  }

  double absoluteMomentum(const State& state) const {
    return singleStepper().absoluteMomentum(primaryState(state));
  }

  Vector3 momentum(const State& state) const {
    return singleStepper().momentum(primaryState(state));
  }

  double charge(const State& state) const {
    return singleStepper().charge(primaryState(state));
  }

  const ParticleHypothesis& particleHypothesis(const State& state) const {
    return singleStepper().particleHypothesis(primaryState(state));
  }

  double time(const State& state) const {
    return singleStepper().time(primaryState(state));
  }

  IntersectionStatus updateSurfaceStatus(
      State& state, const Surface& surface, std::uint8_t index,
      Direction navDir, const BoundaryTolerance& boundaryTolerance,
      double surfaceTolerance, ConstrainedStep::Type stype,
      const Logger& logger = getDummyLogger()) const {
    return m_multiStepperImpl.updateSurfaceStatus(
        state.multiStepperState, surface, index, navDir, boundaryTolerance,
        surfaceTolerance, stype, logger);
  }

  void updateStepSize(State& state, const NavigationTarget& target,
                      Direction direction, ConstrainedStep::Type stype) const {
    m_multiStepperImpl.updateStepSize(state.multiStepperState, target,
                                      direction, stype);
  }

  void updateStepSize(State& state, double stepSize,
                      ConstrainedStep::Type stype) const {
    m_multiStepperImpl.updateStepSize(state.multiStepperState, stepSize, stype);
  }

  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return m_multiStepperImpl.getStepSize(state.multiStepperState, stype);
  }

  void releaseStepSize(State& state, ConstrainedStep::Type stype) const {
    m_multiStepperImpl.releaseStepSize(state.multiStepperState, stype);
  }

  std::string outputStepSize(const State& state) const {
    return m_multiStepperImpl.outputStepSize(state.multiStepperState);
  }

  Result<BoundState> boundState(
      State& state, const Surface& surface, bool transportCovariance = true,
      const FreeToBoundCorrection& freeToBoundCorrection =
          FreeToBoundCorrection(false)) const {
    std::cout << "RiddersStepper::boundState called with transportCovariance = "
              << transportCovariance << "\n";

    if (!transportCovariance) {
      Result<BoundState> result = singleStepper().boundState(
          primaryState(state), surface, false, freeToBoundCorrection);
      std::get<1>(*result) = state.jacobian;
      return result;
    }

    std::vector<BoundVector> boundVectors;
    std::vector<double> pathLenghts;
    boundVectors.reserve(state.multiStepperState.components.size());
    pathLenghts.reserve(state.multiStepperState.components.size());
    for (auto& component : state.multiStepperState.components) {
      const Result<BoundState> result = singleStepper().boundState(
          component.state, surface, false, freeToBoundCorrection);
      std::cout << "Bound state result: " << (result.ok() ? "success" : "error")
                << "\n";
      if (!result.ok()) {
        return result.error();
      }
      boundVectors.push_back(std::get<0>(*result).parameters());
      pathLenghts.push_back(std::get<2>(*result));
    }
    const auto& referenceVector = boundVectors.front();

    state.jacobian = BoundMatrix::Zero();
    std::array<std::size_t, eBoundSize> numberOfVariationsPerParameter{};

    for (std::size_t i = 0; i < state.variationMap.size(); ++i) {
      const auto& [index, delta] = state.variationMap[i];
      const auto& nudgedVector = boundVectors[1 + i];

      state.jacobian.col(index) += (nudgedVector - referenceVector) / delta;
      ++numberOfVariationsPerParameter[index];
    }
    for (std::size_t i = 0; i < eBoundSize; ++i) {
      if (numberOfVariationsPerParameter[i] > 0) {
        state.jacobian.col(i) /= numberOfVariationsPerParameter[i];
      } else {
        // TODO
        state.jacobian.col(i) = BoundVector::Unit(i);
      }
    }

    std::cout << "Jacobian:\n" << state.jacobian << "\n";

    state.cov = state.jacobian * state.cov * state.jacobian.transpose();

    std::cout << "Transported covariance:\n" << state.cov << "\n";

    BoundTrackParameters newParams(surface.getSharedPtr(), referenceVector,
                                   state.cov, particleHypothesis(state));

    state.variationMap =
        m_config.parameterVariation->variationMap(referenceVector, state.cov);

    return Result<BoundState>::success(
        BoundState{std::move(newParams), state.jacobian, pathLenghts.front()});
  }

  bool prepareCurvilinearState(State& state) const {
    return m_multiStepperImpl.prepareCurvilinearState(state.multiStepperState);
  }

  BoundState curvilinearState(State& state,
                              bool transportCovariance = true) const {
    if (!transportCovariance) {
      return singleStepper().curvilinearState(primaryState(state), false);
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
                   state.multiStepperState.options.geoContext)
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
    const Result<double> stepResult = m_multiStepperImpl.step(
        state.multiStepperState, propagationDirection, material);
    if (!stepResult.ok()) {
      return stepResult.error();
    }
    state.pathAccumulated += *stepResult;
    return *stepResult;
  }

  void setIdentityJacobian(State& state) const {
    m_multiStepperImpl.setIdentityJacobian(state.multiStepperState);
  }

 private:
  Config m_config;

  MultiStepperImpl m_multiStepperImpl;

  const auto& singleStepper() const {
    return m_multiStepperImpl.singleStepper();
  }

  auto& primaryState(State& state) const {
    return state.multiStepperState.components.front().state;
  }
  const auto& primaryState(const State& state) const {
    return state.multiStepperState.components.front().state;
  }
};

}  // namespace Acts::Experimental
