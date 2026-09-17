// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/MultiStepperLoop.hpp"

#include "Acts/Propagator/MultiStepperError.hpp"

#include <algorithm>
#include <optional>
#include <vector>

namespace Acts {

template <Concepts::SingleStepper S, typename R>
auto MultiStepperLoop<S, R>::boundParameters(const State& state,
                                             const Surface& surface) const
    -> Result<BoundParameters> {
  assert(!state.components.empty());

  MultiComponentBoundTrackParameters params(
      surface.getSharedPtr(), state.covTransport, state.particleHypothesis);
  params.reserve(numberComponents(state));

  for (const auto& cmp : state.components) {
    auto bp = m_singleStepper.boundParameters(cmp.state, surface);
    if (bp.ok()) {
      params.pushComponent(cmp.weight, bp->parameters(), bp->covariance());
    }
  }

  if (params.empty()) {
    return MultiStepperError::AllComponentsConversionToBoundFailed;
  }

  params.normalizeWeights();

  return params;
}

template <Concepts::SingleStepper S, typename R>
auto MultiStepperLoop<S, R>::curvilinearParameters(const State& state) const
    -> BoundParameters {
  assert(!state.components.empty());

  std::vector<std::tuple<double, Vector4, Vector3, double, BoundMatrix>> cmps;
  cmps.reserve(numberComponents(state));

  for (const auto& cmp : state.components) {
    const auto cp = m_singleStepper.curvilinearParameters(cmp.state);

    cmps.emplace_back(cmp.weight, cp.fourPosition(state.options.geoContext),
                      cp.direction(), cp.qOverP(),
                      cp.covariance().value_or(BoundMatrix::Zero()));
  }

  return MultiComponentBoundTrackParameters::createCurvilinear(
      state.options.geoContext, cmps, state.particleHypothesis);
}

template <Concepts::SingleStepper S, typename R>
auto MultiStepperLoop<S, R>::transportToBound(
    State& state, const Surface& surface,
    const FreeToBoundCorrection& freeToBoundCorrection) const
    -> Result<Jacobian> {
  assert(!state.components.empty());

  const std::size_t selected = Reducer::index(state);
  std::optional<Result<Jacobian>> jacobian;

  for (std::size_t i = 0; i < state.components.size(); ++i) {
    auto& cmpState = state.components[i].state;

    // Force the component to be on the surface
    // This needs to be done because of the `averageOnSurface`-option of the
    // `MultiStepperSurfaceReached`-Aborter, which can be configured to end the
    // propagation when the mean of all components reached the destination
    // surface. Thus, it is not garantueed that all states are actually
    // onSurface.
    cmpState.pars.template segment<3>(eFreePos0) =
        surface
            .intersect(state.options.geoContext,
                       cmpState.pars.template segment<3>(eFreePos0),
                       cmpState.pars.template segment<3>(eFreeDir0),
                       BoundaryTolerance::Infinite())
            .closest()
            .position();

    auto res = m_singleStepper.transportToBound(cmpState, surface,
                                                freeToBoundCorrection);
    if (i == selected) {
      jacobian = std::move(res);
    }
  }

  return std::move(*jacobian);
}

template <Concepts::SingleStepper S, typename R>
Result<double> MultiStepperLoop<S, R>::step(
    State& state, Direction propDir, const IVolumeMaterial* material) const {
  using Status = Acts::IntersectionStatus;

  auto& components = state.components;
  const Logger& logger = *m_logger;

  // Update step count
  state.steps++;

  // Check if we abort because of m_stepLimitAfterFirstComponentOnSurface
  if (state.stepCounterAfterFirstComponentOnSurface) {
    (*state.stepCounterAfterFirstComponentOnSurface)++;

    // If the limit is reached, remove all components which are not on a
    // surface, reweight the components, perform no step and return 0
    if (*state.stepCounterAfterFirstComponentOnSurface >=
        m_stepLimitAfterFirstComponentOnSurface) {
      //! [step limit]
      for (auto& cmp : components) {
        if (cmp.status != Status::onSurface) {
          cmp.status = Status::unreachable;
        }
      }

      ACTS_VERBOSE("Stepper performed "
                   << m_stepLimitAfterFirstComponentOnSurface
                   << " steps after the first component hit a surface.");
      ACTS_VERBOSE(
          "-> remove all components not on a surface, perform no step");

      removeMissedComponents(state);
      reweightComponents(state);
      //! [step limit]

      ACTS_VERBOSE(components.size()
                   << " components left after removing missed components");

      state.stepCounterAfterFirstComponentOnSurface.reset();

      return 0.0;
    }
  }

  // Flag indicating if we need to reweight in the end
  bool reweightNecessary = false;

  // If at least one component is on a surface, we can remove all missed
  // components before the step. If not, we must keep them for the case that all
  // components miss and we need to retarget
  const auto cmpsOnSurface = std::ranges::count_if(components, [&](auto& cmp) {
    return cmp.status == IntersectionStatus::onSurface;
  });

  if (cmpsOnSurface > 0) {
    removeMissedComponents(state);
    reweightNecessary = true;
  }

  // Loop over all components and collect results in vector, write some
  // summary information to a stringstream
  SmallVector<std::optional<Result<double>>> results;
  double accumulatedPathLength = 0.0;
  std::size_t errorSteps = 0;

  // Lambda that performs the step for a component and returns false if the step
  // went ok and true if there was an error
  auto errorInStep = [this, &results, propDir, material, &accumulatedPathLength,
                      &errorSteps, &reweightNecessary](auto& component) {
    if (component.status == Status::onSurface) {
      // We need to add these, so the propagation does not fail if we have only
      // components on surfaces and failing states
      results.emplace_back(std::nullopt);
      return false;
    }

    results.emplace_back(
        m_singleStepper.step(component.state, propDir, material));

    if (results.back()->ok()) {
      accumulatedPathLength += component.weight * results.back()->value();
      return false;
    } else {
      ++errorSteps;
      reweightNecessary = true;
      return true;
    }
  };

  // Loop over components and remove errorous components
  const auto removedTail = std::ranges::remove_if(components, errorInStep);
  components.erase(removedTail.begin(), removedTail.end());

  // Reweight if necessary
  if (reweightNecessary) {
    reweightComponents(state);
  }

  // Print the result vector to a string so we can log it
  auto summary = [](auto& result_vec) {
    std::stringstream ss;
    for (auto& optRes : result_vec) {
      if (!optRes) {
        ss << "on surface | ";
      } else if (optRes->ok()) {
        ss << optRes->value() << " | ";
      } else {
        ss << optRes->error() << " | ";
      }
    }
    auto str = ss.str();
    str.resize(str.size() - 3);
    return str;
  };

  // Print the summary
  if (errorSteps == 0) {
    ACTS_VERBOSE("Performed steps: " << summary(results));
  } else {
    ACTS_WARNING("Performed steps with errors: " << summary(results));
  }

  // Return error if there is no ok result
  if (components.empty()) {
    return MultiStepperError::AllComponentsSteppingError;
  }

  // Invalidate the component status after each step
  for (auto& cmp : components) {
    cmp.status = Status::unreachable;
  }

  // Return the weighted accumulated path length of all successful steps
  state.pathAccumulated += accumulatedPathLength;
  return accumulatedPathLength;
}

}  // namespace Acts
