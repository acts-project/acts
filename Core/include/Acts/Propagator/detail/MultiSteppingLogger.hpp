// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <memory>
#include <vector>

// TODO slightly modified for Multistepper
struct MultiSteppingLogger {
  struct this_result {
    std::vector<std::vector<Acts::detail::Step>> steps;
    std::vector<Acts::detail::Step> averaged_steps;
  };

  using result_type = this_result;
  bool sterile = false;

  template <typename propagator_state_t, typename stepper_t>
  void operator()(const propagator_state_t& state, const stepper_t& stepper,
                  result_type& result) const {
    if (sterile or state.navigation.targetReached) {
      return;
    }

    if (result.steps.empty())
      result.steps.resize(stepper.numberComponents(state.stepping));

    // Very poor workaround vor changing number of components... the problem is,
    // we do not have the information here in the logger which component is the
    // child of which older component...
    if (stepper.numberComponents(state.stepping) == result.steps.size()) {
      const auto N = std::min(stepper.numberComponents(state.stepping),
                              result.steps.size());

      // Individual components
      for (auto i = 0ul; i < N; ++i) {
        Acts::detail::Step step;
        //       step.stepSize = state.stepping.stepSize;
        step.position = stepper.position(i, state.stepping);
        step.momentum = stepper.momentum(i, state.stepping) *
                        stepper.direction(i, state.stepping);

        if (state.navigation.currentSurface != nullptr)
          step.surface = state.navigation.currentSurface->getSharedPtr();

        step.volume = state.navigation.currentVolume;
        result.steps[i].push_back(std::move(step));
      }
    }

    // Average
    Acts::detail::Step avg_step;
    //     avg_step.stepSize = state.stepping.stepSize;
    avg_step.position = stepper.position(state.stepping);
    avg_step.momentum =
        stepper.momentum(state.stepping) * stepper.direction(state.stepping);

    if (state.navigation.currentSurface != nullptr)
      avg_step.surface = state.navigation.currentSurface->getSharedPtr();

    avg_step.volume = state.navigation.currentVolume;
    result.averaged_steps.push_back(std::move(avg_step));
  }

  template <typename propagator_state_t, typename stepper_t>
  void operator()(propagator_state_t& /*unused*/,
                  const stepper_t& /*unused*/) const {}
};
