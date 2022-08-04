// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

template <typename E, typename R, typename A>
auto MultiEigenStepperLoop<E, R, A>::boundState(
    State& state, const Surface& surface, bool transportCov,
    const FreeToBoundCorrection& freeToBoundCorrection) const
    -> Result<BoundState> {
  if (numberComponents(state) == 1) {
    return SingleStepper::boundState(state.components.front().state, surface,
                                     transportCov, freeToBoundCorrection);
  } else {  // Do the combinatio
    SmallVector<std::pair<double, BoundTrackParameters>> states;
    double accumulatedPathLength = 0.0;
    int failedBoundTransforms = 0;

    for (auto i = 0ul; i < numberComponents(state); ++i) {
      auto bs = SingleStepper::boundState(state.components[i].state, surface,
                                          transportCov, freeToBoundCorrection);

      if (bs.ok()) {
        states.push_back(
            {state.components[i].weight, std::get<BoundTrackParameters>(*bs)});
        accumulatedPathLength +=
            std::get<double>(*bs) * state.components[i].weight;
      } else {
        failedBoundTransforms++;
      }
    }

    if (states.size() == 0) {
      return MultiStepperError::AllComponentsConversionToBoundFailed;
    }

    if (failedBoundTransforms > 0) {
      return MultiStepperError::SomeComponentsConversionToBoundFailed;
    }

    // TODO At ATLAS, the final parameters seem to be computed with the mode of
    // the mixture. At the moment, we use the mean of the mixture here, but
    // there should be done a comparison sometimes in the future. This could
    // also be configurable maybe...
    const auto proj = [&](const auto& wbs) {
      return std::tie(wbs.first, wbs.second.parameters(),
                      wbs.second.covariance());
    };

    const auto [params, cov] =
        detail::angleDescriptionSwitch(surface, [&](const auto& desc) {
          return detail::combineGaussianMixture(states, proj, desc);
        });

    return BoundState{BoundTrackParameters(surface.getSharedPtr(), params, cov),
                      Jacobian::Zero(), accumulatedPathLength};
  }
}

template <typename E, typename R, typename A>
auto MultiEigenStepperLoop<E, R, A>::curvilinearState(State& state,
                                                      bool transportCov) const
    -> CurvilinearState {
  if (numberComponents(state) == 1) {
    return SingleStepper::curvilinearState(state.components.front().state,
                                           transportCov);
  } else {
    Vector4 pos4 = Vector4::Zero();
    Vector3 dir = Vector3::Zero();
    ActsScalar qop = 0.0;
    BoundSymMatrix cov = BoundSymMatrix::Zero();
    ActsScalar pathLenth = 0.0;

    // TODO At ATLAS, the final parameters seem to be computed with the mode of
    // the mixture. At the moment, we use the mean of the mixture here, but
    // there should be done a comparison sometimes in the future. This could
    // also be configurable maybe...
    for (auto i = 0ul; i < numberComponents(state); ++i) {
      const auto [cp, jac, pl] = SingleStepper::curvilinearState(
          state.components[i].state, transportCov);

      pos4 += state.components[i].weight * cp.fourPosition(state.geoContext);
      dir += state.components[i].weight * cp.unitDirection();
      qop += state.components[i].weight * (cp.charge() / cp.absoluteMomentum());
      if (cp.covariance()) {
        cov += state.components[i].weight * *cp.covariance();
      }
      pathLenth += state.components[i].weight * pathLenth;
    }

    return CurvilinearState{CurvilinearTrackParameters(pos4, dir, qop, cov),
                            Jacobian::Zero(), pathLenth};
  }
}

template <typename E, typename R, typename A>
template <typename propagator_state_t>
Result<double> MultiEigenStepperLoop<E, R, A>::step(
    propagator_state_t& state) const {
  const auto& logger = state.options.logger;
  State& stepping = state.stepping;

  // It is not possible to remove components from the vector, since the
  // GSF actor relies on the fact that the ordering and number of
  // components does not change
  auto invalidateComponent = [](auto& cmp) {
    cmp.status = Intersection3D::Status::missed;
    cmp.weight = 0.0;
    cmp.state.pars.template segment<3>(eFreeDir0) = Vector3::Zero();
  };

  // Lambda for reweighting the components
  auto reweight = [](auto& cmps) {
    ActsScalar sumOfWeights = 0.0;
    for (const auto& cmp : cmps) {
      sumOfWeights += cmp.weight;
    }
    for (auto& cmp : cmps) {
      cmp.weight /= sumOfWeights;
    }
  };

  // Update step count
  stepping.steps++;

  // Check if we abort because of m_stepLimitAfterFirstComponentOnSurface
  if (stepping.stepCounterAfterFirstComponentOnSurface) {
    (*stepping.stepCounterAfterFirstComponentOnSurface)++;

    // If the limit is reached, remove all components which are not on a
    // surface, reweight the components, perform no step and return 0
    if (*stepping.stepCounterAfterFirstComponentOnSurface >=
        m_stepLimitAfterFirstComponentOnSurface) {
      for (auto& cmp : stepping.components) {
        if (cmp.status != Intersection3D::Status::onSurface) {
          invalidateComponent(cmp);
        }
      }

      reweight(stepping.components);

      ACTS_VERBOSE("Stepper performed "
                   << m_stepLimitAfterFirstComponentOnSurface
                   << " after the first component hit a surface.");
      ACTS_VERBOSE(
          "-> remove all components not on a surface, perform no step");

      stepping.stepCounterAfterFirstComponentOnSurface.reset();

      return 0.0;
    }
  }

  // Loop over all components and collect results in vector, write some
  // summary information to a stringstream
  SmallVector<std::optional<Result<double>>> results;
  double accumulatedPathLength = 0.0;
  std::size_t errorSteps = 0;

  for (auto& component : stepping.components) {
    // We must also propagate missed components for the case that all
    // components miss the target and we need to re-target
    if (component.status == Intersection3D::Status::onSurface) {
      // We need to add these, so the propagation does not fail if we have only
      // components on surfaces and failing states
      results.emplace_back(std::nullopt);
      continue;
    }

    using ThisSinglePropState =
        SinglePropState<SingleState, decltype(state.navigation),
                        decltype(state.options), decltype(state.geoContext)>;

    ThisSinglePropState single_state(component.state, state.navigation,
                                     state.options, state.geoContext);

    results.emplace_back(SingleStepper::step(single_state));

    if (results.back()->ok()) {
      accumulatedPathLength += component.weight * results.back()->value();
    } else {
      ++errorSteps;
      invalidateComponent(component);
    }
  }

  // Since we have invalidated some components, we need to reweight
  if (errorSteps > 0) {
    reweight(stepping.components);
  }

  // Print the result vector to a string so we can log it
  auto summary = [](auto& result_vec) {
    std::stringstream ss;
    for (auto& optRes : result_vec) {
      if (not optRes) {
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
  if (errorSteps == results.size()) {
    return MultiStepperError::AllComponentsSteppingError;
  }

  // Return the weighted accumulated path length of all successful steps
  stepping.pathAccumulated += accumulatedPathLength;
  return accumulatedPathLength;
}
}  // namespace Acts
