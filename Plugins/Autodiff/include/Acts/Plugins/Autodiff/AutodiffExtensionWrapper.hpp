// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

#include <autodiff/forward.hpp>
#include <autodiff/forward/eigen.hpp>

namespace Acts {

/// @brief Default RKN4 evaluator for autodiff
template <template <typename> typename basic_extension_t>
struct AutodiffExtensionWrapper {
  /// @brief Default constructor
  AutodiffExtensionWrapper() = default;

  // Some typedefs
  using AutodiffScalar = autodiff::dual;

  using AutodiffVector3 = ActsMatrix<AutodiffScalar, 3, 1>;
  using AutodiffFreeMatrix = ActsMatrix<AutodiffScalar, 8, 8>;
  using AutodiffFreeVector = ActsMatrix<AutodiffScalar, 8, 1>;

  // The double-extension is needed to communicate with the "outer world" (the
  // stepper) and ensures it behaves exactly as the underlying extension, with
  // the exception of the computation of the transport-matrix. The corresponding
  // autodiff-extension can be found in the RKN4step-member-function (since it
  // is only needed locally). Another advantage of this approach is, that we do
  // not differentiate through the adaptive stepsize estimation in the stepper.
  basic_extension_t<double> m_doubleExtension;

  // Just call underlying extension
  template <typename propagator_state_t, typename stepper_t>
  int bid(const propagator_state_t& ps, const stepper_t& st) const {
    return m_doubleExtension.bid(ps, st);
  }

  // Just call underlying extension
  template <typename propagator_state_t, typename stepper_t>
  bool k(const propagator_state_t& state, const stepper_t& stepper,
         Vector3D& knew, const Vector3D& bField, std::array<double, 4>& kQoP,
         const int i = 0, const double h = 0.,
         const Vector3D& kprev = Vector3D()) {
    return m_doubleExtension.k(state, stepper, knew, bField, kQoP, i, h, kprev);
  }

  // Just call underlying extension
  template <typename propagator_state_t, typename stepper_t>
  bool finalize(propagator_state_t& state, const stepper_t& stepper,
                const double h) const {
    return m_doubleExtension.finalize(state, stepper, h);
  }

  // Here we call a custom implementation to compute the transport matrix
  template <typename propagator_state_t, typename stepper_t>
  bool finalize(propagator_state_t& state, const stepper_t& stepper,
                const double h, FreeMatrix& D) const {
    m_doubleExtension.finalize(state, stepper, h);
    return transportMatrix(state, stepper, h, D);
  }

 private:
  // A fake stepper-state
  struct FakeStepperState {
    AutodiffScalar t, p;
    AutodiffVector3 pos, dir;
    AutodiffFreeVector derivative;
    double q;
    bool covTransport = false;
  };

  // A fake propagator state
  template <class options_t, class navigation_t>
  struct FakePropState {
    FakeStepperState stepping;
    const options_t& options;
    const navigation_t& navigation;
  };

  // A fake stepper
  struct FakeStepper {
    auto charge(const FakeStepperState& s) const { return s.q; }
    auto momentum(const FakeStepperState& s) const { return s.p; }
    auto direction(const FakeStepperState& s) const { return s.dir; }
    auto position(const FakeStepperState& s) const { return s.pos; }
  };

  // Here the autodiff jacobian is computed
  template <typename propagator_state_t, typename stepper_t>
  bool transportMatrix(propagator_state_t& state, const stepper_t& stepper,
                       const double h, FreeMatrix& D) const {
    // Initialize fake stepper
    using ThisFakePropState =
        FakePropState<decltype(state.options), decltype(state.navigation)>;

    ThisFakePropState fstate{FakeStepperState(), state.options,
                             state.navigation};

    fstate.stepping.q = stepper.charge(state.stepping);

    // Init dependent values for autodiff
    AutodiffFreeVector initial_params;
    initial_params.segment<3>(eFreePos0) = stepper.position(state.stepping);
    initial_params(eFreeTime) = stepper.time(state.stepping);
    initial_params.segment<3>(eFreeDir0) = stepper.direction(state.stepping);
    initial_params(eFreeQOverP) =
        (fstate.stepping.q != 0. ? fstate.stepping.q : 1.) /
        stepper.momentum(state.stepping);

    const auto& sd = state.stepping.stepData;

    // Compute jacobian
    using namespace autodiff::forward;

    D = jacobian([&](const auto& in) { return RKN4step(in, sd, fstate, h); },
                 wrt(initial_params), at(initial_params))
            .template cast<double>();

    return true;
  }

  template <typename step_data_t, typename fake_state_t>
  auto RKN4step(const AutodiffFreeVector& in, const step_data_t& sd,
                fake_state_t state, const double h) const {
    // Initialize fake stepper
    FakeStepper stepper;

    // Set dependent variables
    state.stepping.pos = in.segment<3>(eFreePos0);
    state.stepping.t = in(eFreeTime);
    state.stepping.dir = in.segment<3>(eFreeDir0);
    state.stepping.p =
        (state.stepping.q != 0. ? state.stepping.q : 1.) / in(eFreeQOverP);

    std::array<AutodiffScalar, 4> kQoP;
    std::array<AutodiffVector3, 4> k;

    // Autodiff instance of the extension
    basic_extension_t<AutodiffScalar> ext;

    // Compute k values. Assume all return true, since these parameters
    // are already validated by the "outer RKN4"
    ext.k(state, stepper, k[0], sd.B_first, kQoP);
    ext.k(state, stepper, k[1], sd.B_middle, kQoP, 1, h * 0.5, k[0]);
    ext.k(state, stepper, k[2], sd.B_middle, kQoP, 2, h * 0.5, k[1]);
    ext.k(state, stepper, k[3], sd.B_last, kQoP, 3, h, k[2]);

    // finalize
    ext.finalize(state, stepper, h);

    // Compute RKN4 integration
    AutodiffFreeVector out;

    // position
    out.segment<3>(eFreePos0) = in.segment<3>(eFreePos0) +
                                h * in.segment<3>(eFreeDir0) +
                                h * h / 6. * (k[0] + k[1] + k[2]);

    // direction
    auto final_dir =
        in.segment<3>(eFreeDir0) + h / 6. * (k[0] + 2. * (k[1] + k[2]) + k[3]);

    out.segment<3>(eFreeDir0) = final_dir / final_dir.norm();

    // qop
    out(eFreeQOverP) =
        (state.stepping.q != 0. ? state.stepping.q : 1.) / state.stepping.p;

    // time
    out(eFreeTime) = state.stepping.t;

    return out;
  }
};
}  // namespace Acts
