// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/Result.hpp"

#include <cassert>

namespace Acts {

class Surface;

// TODO what about ridders propagator?
// TODO should there also be Stepping and Navigation?
// TODO should the Propagation be composed of Stepping and Navigation objects?
// TODO make this a detail for now?
template <typename propagator_t, typename propagator_state_t>
class Propagation {
 public:
  using Propagator = propagator_t;
  using Stepper = typename propagator_t::Stepper;
  using Navigator = typename propagator_t::Navigator;
  using State = propagator_state_t;

  Propagation(const Propagator &propagator, State &state)
      : m_propagator(&propagator), m_state(&state) {
    assert(m_propagator != nullptr && "Propagator pointer is null");
    assert(m_state != nullptr && "State pointer is null");
  }

  const Propagator &propagator() const { return *m_propagator; }
  const Stepper &stepper() const { return propagator().stepper(); }
  const Navigator &navigator() const { return propagator().navigator(); }
  State &state() const { return *m_state; }

  template <typename parameters_t>
  Result<void> initialize(const parameters_t &parameters) {
    return propagator().initialize(*m_state, parameters);
  }
  Result<void> performStep() { return propagator().performStep(*m_state); }
  Result<void> reachNextSurface() {
    return propagator().reachNextSurface(*m_state);
  }

  Vector3 position() const { return stepper().position(m_state->stepping); }
  Vector3 direction() const { return stepper().direction(m_state->stepping); }
  const Surface *currentSurface() const {
    return navigator().currentSurface(m_state->navigation);
  }

  void transportCovarianceToBound() {
    stepper().transportCovarianceToBound(m_state->stepping, *currentSurface());
  }
  void transportCovarianceToCurvilinear() {
    stepper().transportCovarianceToCurvilinear(m_state->stepping);
  }

  auto boundState(bool transportCov) {
    return stepper().boundState(m_state->stepping, *currentSurface(),
                                transportCov);
  }
  auto curvilinearState(bool transportCov) {
    return stepper().curvilinearState(m_state->stepping, transportCov);
  }

 private:
  const Propagator *m_propagator{nullptr};
  State *m_state{nullptr};
};

}  // namespace Acts
