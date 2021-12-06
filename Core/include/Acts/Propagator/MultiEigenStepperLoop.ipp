// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

namespace Acts {

template <typename E, typename R, typename A>
template <typename state_type_t>
state_type_t MultiEigenStepperLoop<E, R, A>::combineComponents(
    State& state, const Surface *surface, bool transportCov) const {
  // Check if the state_type_t matches the requirements
  static_assert(std::is_same_v<state_type_t, BoundState> ||
                std::is_same_v<state_type_t, CurvilinearState>);
  constexpr bool is_bound_state =
      std::is_same_v<state_type_t, EigenStepper<>::BoundState>;

  // Do the combination
  std::vector<std::pair<double, state_type_t>> states;
  double accumulatedPathLength = 0.0;
  int failedBoundTransforms = 0;

  for (auto i = 0ul; i < numberComponents(state); ++i) {
    if constexpr (is_bound_state) {
      auto bs = SingleStepper::boundState(state.components[i].state, *surface,
                                          transportCov);

      if (bs.ok()) {
        states.push_back({state.components[i].weight, *bs});
        accumulatedPathLength +=
            std::get<double>(*bs) * state.components[i].weight;
      } else {
        failedBoundTransforms++;
      }
    } else {
      auto cs = SingleStepper::curvilinearState(state.componenents[i].state,
                                                transportCov);

      states.push_back({state.components[i].weight, cs});
      accumulatedPathLength +=
          std::get<double>(cs) * state.components[i].weight;
    }
  }

  if constexpr (is_bound_state) {
    if (failedBoundTransforms > 0) {
      ACTS_ERROR("Multi component bound state: "
                 << failedBoundTransforms << " of " << numberComponents(state)
                 << " transforms failed");
    }
  }

  // TODO also implement a method of using the mode
  const auto [params, cov] = detail::combineBoundGaussianMixture(
      states.begin(), states.end(), [&](const auto& wbs) {
        const auto& bp = std::get<BoundTrackParameters>(wbs.second);
        return std::tie(wbs.first, bp.parameters(), bp.covariance());
      });

  if constexpr (is_bound_state) {
    return BoundState{BoundTrackParameters(surface->getSharedPtr(), params, cov),
                      Jacobian::Zero(), accumulatedPathLength};
  } else {
    return CurvilinearState{CurvilinearTrackParameters(params, cov),
                            Jacobian::Zero(), accumulatedPathLength};
  }
}

template <typename E, typename R, typename A>
auto MultiEigenStepperLoop<E, R, A>::boundState(State& state,
                                                const Surface& surface,
                                                bool transportCov) const
    -> Result<BoundState> {
  if (numberComponents(state) == 1) {
    return SingleStepper::boundState(state.components.front().state, surface,
                                     transportCov);
  } else {
    return combineComponents<BoundState>(state, &surface, transportCov);
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
    return combineComponents<CurvilinearState>(state, nullptr, transportCov);
  }
}

template <typename E, typename R, typename A>
template <typename propagator_state_t>
Result<double> MultiEigenStepperLoop<E, R, A>::step(
    propagator_state_t& state) const {
  State& stepping = state.stepping;

  // Emit warning if charge is not the same for all componenents
  {
    std::stringstream ss;
    bool charge_ambigous = false;
    for (const auto& cmp : stepping.components) {
      ss << cmp.state.q << " ";
      if (cmp.state.q != stepping.components.front().state.q) {
        charge_ambigous = true;
      }
    }

    if (charge_ambigous) {
      ACTS_VERBOSE(stepping.steps << "Charge of components is ambigous: "
                                  << ss.str());
    } else {
      ACTS_VERBOSE(stepping.steps << "Charge of components: " << ss.str());
    }
  }

  // Update step count
  stepping.steps++;

  // Check if we abort because of m_stepLimitAfterFirstComponentOnSurface
  if (stepping.stepCounterAfterFirstComponentOnSurface) {
    (*stepping.stepCounterAfterFirstComponentOnSurface)++;

    // If the limit is reached, remove all components which are not on a
    // surface, reweight the components, perform no step and return 0
    if (*stepping.stepCounterAfterFirstComponentOnSurface >=
        m_stepLimitAfterFirstComponentOnSurface) {
      auto& cmps = stepping.components;

      // It is not possible to remove components from the vector, since the
      // GSF actor relies on the fact that the ordering and number of
      // components does not change
      for (auto& cmp : cmps) {
        if (cmp.status != Intersection3D::Status::onSurface) {
          cmp.status = Intersection3D::Status::missed;
          cmp.weight = 0.0;
          cmp.state.pars.template segment<3>(eFreeDir0) = Vector3::Zero();
        }
      }

      // Reweight
      const auto sum_of_weights = std::accumulate(
          begin(cmps), end(cmps), ActsScalar{0},
          [](auto sum, const auto& cmp) { return sum + cmp.weight; });
      for (auto& cmp : cmps) {
        cmp.weight /= sum_of_weights;
      }

      ACTS_VERBOSE(
          "hit m_stepLimitAfterFirstComponentOnSurface, "
          "perform no step");

      stepping.stepCounterAfterFirstComponentOnSurface.reset();

      return 0.0;
    }
  }

  // Loop over all components and collect results in vector, write some
  // summary information to a stringstream
  std::vector<Result<double>> results;
  std::stringstream ss;

  for (auto& component : stepping.components) {
    // We must also propagate missed components for the case that all
    // components miss the target we need to retarget
    if (component.status == Intersection3D::Status::onSurface) {
      ss << "cmp skipped\t";
      continue;
    }

    SinglePropState single_state{component.state, state.navigation,
                                 state.options, state.geoContext};
    results.push_back(SingleStepper::step(single_state));

    if (results.back().ok()) {
      ss << *results.back() << "\t";
    } else {
      ss << "step error: " << results.back().error() << "\t";
    }
  }

  // Return no component was updated
  if (results.empty()) {
    return 0.0;
  }

  // Collect pointers to results which are ok, since Result is not copyable
  std::vector<Result<double>*> ok_results;
  for (auto& res : results) {
    if (res.ok()) {
      ok_results.push_back(&res);
    }
  }

  // Return error if there is no ok result
  if (ok_results.empty()) {
    return MultiStepperError::AllComponentsSteppingError;
  }

  // Print the summary
  if (ok_results.size() == results.size()) {
    ACTS_VERBOSE("Performed steps: " << ss.str());
  } else {
    ACTS_WARNING("Performed steps with errors: " << ss.str());
  }

  // Compute the average stepsize for the return value and the
  // pathAccumulated
  const auto avg_step =
      std::accumulate(begin(ok_results), end(ok_results), 0.,
                      [](auto sum, auto res) { return sum + res->value(); }) /
      static_cast<double>(ok_results.size());
  stepping.pathAccumulated += avg_step;

  return avg_step;
}
}  // namespace Acts
