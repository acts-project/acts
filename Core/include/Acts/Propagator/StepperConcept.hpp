// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {
class Surface;

namespace concept {
  namespace Stepper {

    template <typename T>
    using state_t = typename T::State;

    template <typename T>
    using return_t = typename T::template return_parameter_type<void, void>;

    template <typename T>
    using bound_state_t = typename T::BoundState;
    template <typename T>
    using curvilinear_state_t = typename T::CurvilinearState;

    METHOD_TRAIT(get_field_t, getField);
    METHOD_TRAIT(position_t, position);
    METHOD_TRAIT(direction_t, direction);
    METHOD_TRAIT(momentum_t, momentum);
    METHOD_TRAIT(charge_t, charge);
    METHOD_TRAIT(surface_reached_t, surfaceReached);
    METHOD_TRAIT(bound_state_method_t, boundState);
    METHOD_TRAIT(curvilinear_state_method_t, curvilinearState);
    METHOD_TRAIT(update_t, update);
    METHOD_TRAIT(covariance_transport_t, covarianceTransport);
    METHOD_TRAIT(step_t, step);

    template <typename T>
    using cov_transport_t = decltype(std::declval<T>().covTransport);
    template <typename T>
    using cov_t = decltype(std::declval<T>().cov);
    template <typename T>
    using nav_dir_t = decltype(std::declval<T>().navDir);
    template <typename T>
    using path_accumulated_t = decltype(std::declval<T>().pathAccumulated);
    template <typename T>
    using step_size_t = decltype(std::declval<T>().stepSize);

    // clang-format off
    template <typename S>
    constexpr bool StepperStateConcept
      = require<has_member<S, cov_transport_t, bool>,
                has_member<S, cov_t, ActsSymMatrixD<5>>,
                has_member<S, nav_dir_t, NavigationDirection>,
                has_member<S, path_accumulated_t, double>,
                has_member<S, step_size_t, detail::ConstrainedStep>
               >;
    // clang-format on

    // clang-format off
    template <typename S, typename state = typename S::State>
    constexpr bool StepperConcept
      = require<exists<state_t, S>,
                exists<bound_state_t, S>,
                exists<curvilinear_state_t, S>,
                exists<return_t, S>,
                has_method<const S, Vector3D, get_field_t, state&, const Vector3D&>,
                has_method<const S, Vector3D, position_t, const state&>,
                has_method<const S, Vector3D, direction_t, const state&>,
                has_method<const S, double, momentum_t, const state&>,
                has_method<const S, double, charge_t, const state&>,
                has_method<const S, bool, surface_reached_t, const state&, const Surface*>,
                has_method<const S, typename S::BoundState, bound_state_method_t, state&, const Surface&, bool>,
                has_method<const S, typename S::CurvilinearState, curvilinear_state_method_t, state&, bool>,
                has_method<const S, void, update_t, state&, const BoundParameters&>,
                has_method<const S, void, update_t, state&, const Vector3D&, const Vector3D&, double>,
                has_method<const S, void, covariance_transport_t, state&, bool>,
                has_method<const S, void, covariance_transport_t, state&, const Surface&, bool>
               >;
    // clang-format on
  }
}

template <typename stepper, typename state = typename stepper::State>
constexpr bool StepperConcept
    = Acts::concept::Stepper::StepperConcept<stepper, state>;
template <typename stepper>
constexpr bool StepperStateConcept
    = Acts::concept::Stepper::StepperStateConcept<stepper>;
}
