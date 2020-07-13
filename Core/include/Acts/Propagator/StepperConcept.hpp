// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {
class Surface;

namespace concept {
  namespace Stepper {

  template <typename T>
  using state_t = typename T::State;

  template <typename T>
  using jacobian_t = typename T::Jacobian;
  template <typename T>
  using covariance_t = typename T::Covariance;
  template <typename T>
  using bound_state_t = typename T::BoundState;
  template <typename T>
  using curvilinear_state_t = typename T::CurvilinearState;
  template <typename T>
  using bfield_t = typename T::BField;

  METHOD_TRAIT(get_field_t, getField);
  METHOD_TRAIT(position_t, position);
  METHOD_TRAIT(direction_t, direction);
  METHOD_TRAIT(momentum_t, momentum);
  METHOD_TRAIT(charge_t, charge);
  METHOD_TRAIT(time_t, time);
  METHOD_TRAIT(overstep_t, overstepLimit);
  METHOD_TRAIT(bound_state_method_t, boundState);
  METHOD_TRAIT(curvilinear_state_method_t, curvilinearState);
  METHOD_TRAIT(update_t, update);
  METHOD_TRAIT(covariance_transport_t, covarianceTransport);
  METHOD_TRAIT(step_t, step);
  METHOD_TRAIT(update_surface_status_t, updateSurfaceStatus);
  METHOD_TRAIT(set_step_size_t, setStepSize);
  METHOD_TRAIT(release_step_size_t, releaseStepSize);
  METHOD_TRAIT(output_step_size_t, outputStepSize);

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
                has_member<S, cov_t, BoundSymMatrix>,
                has_member<S, nav_dir_t, NavigationDirection>,
                has_member<S, path_accumulated_t, double>,
                has_member<S, step_size_t, ConstrainedStep>
               >;
  // clang-format on

  // clang-format off
    template <typename S, typename state = typename S::State>
      struct StepperConcept {
        constexpr static bool state_exists = exists<state_t, S>;
        static_assert(state_exists, "State type not found");
        constexpr static bool jacobian_exists = exists<jacobian_t, S>;
        static_assert(jacobian_exists, "Jacobian type not found");
        constexpr static bool covariance_exists = exists<covariance_t, S>;
        static_assert(covariance_exists, "Covariance type not found");
        constexpr static bool bound_state_exists = exists<bound_state_t, S>;
        static_assert(bound_state_exists, "BoundState type not found");
        constexpr static bool curvilinear_state_exists = exists<curvilinear_state_t, S>;
        static_assert(curvilinear_state_exists, "CurvilinearState type not found");
        constexpr static bool bfield_exists = exists<bfield_t, S>;
        static_assert(bfield_exists, "BField type not found");
        constexpr static bool get_field_exists = has_method<const S, Vector3D, get_field_t, state&, const Vector3D&>;
        static_assert(get_field_exists, "getField method not found");
        constexpr static bool position_exists = has_method<const S, Vector3D, position_t, const state&>;
        static_assert(position_exists, "position method not found");
        constexpr static bool direction_exists = has_method<const S, Vector3D, direction_t, const state&>;
        static_assert(direction_exists, "direction method not found");
        constexpr static bool momentum_exists = has_method<const S, double, momentum_t, const state&>;
        static_assert(momentum_exists, "momentum method not found");
        constexpr static bool charge_exists = has_method<const S, double, charge_t, const state&>;
        static_assert(charge_exists, "charge method not found");
        constexpr static bool time_exists = has_method<const S, double, time_t, const state&>;
        static_assert(time_exists, "time method not found");
        constexpr static bool overstep_exists = has_method<const S, double, overstep_t, const state&>;
        static_assert(overstep_exists, "overstepLimit method not found");
        constexpr static bool bound_state_method_exists= has_method<const S, typename S::BoundState, bound_state_method_t, state&, const Surface&>;
        static_assert(bound_state_method_exists, "boundState method not found");
        constexpr static bool curvilinear_state_method_exists = has_method<const S, typename S::CurvilinearState, curvilinear_state_method_t, state&>;
        static_assert(curvilinear_state_method_exists, "curvilinearState method not found");
        constexpr static bool update_method_exists = require<has_method<const S, void, update_t, state&, const FreeVector&, const BoundSymMatrix&>,
                                                             has_method<const S, void, update_t, state&, const Vector3D&, const Vector3D&, double, double>>;
        static_assert(update_method_exists, "update method not found");
        constexpr static bool covariance_transport_exists = require<has_method<const S, void, covariance_transport_t, state&>,
                                                                    has_method<const S, void, covariance_transport_t, state&, const Surface&>>;
        static_assert(covariance_transport_exists, "covarianceTransport method not found");
        constexpr static bool update_surface_exists = has_method<const S, Intersection::Status, update_surface_status_t, state&, const Surface&, const BoundaryCheck&>;
        static_assert(update_surface_exists, "updateSurfaceStatus method not found");
        constexpr static bool set_step_size_exists = has_method<const S, void, set_step_size_t, state&, double, ConstrainedStep::Type>;
        static_assert(set_step_size_exists, "setStepSize method not found");
        constexpr static bool release_step_size_exists = has_method<const S, void, release_step_size_t, state&>;
        static_assert(release_step_size_exists, "releaseStepSize method not found");
        constexpr static bool output_step_size_exists = has_method<const S, std::string, output_step_size_t, const state&>;
        static_assert(output_step_size_exists, "outputStepSize method not found");

        constexpr static bool value = require<state_exists,
                                              jacobian_exists,
                                              covariance_exists,
                                              bound_state_exists,
                                              curvilinear_state_exists,
                                              bfield_exists,
                                              get_field_exists,
                                              position_exists,
                                              direction_exists,
                                              momentum_exists,
                                              charge_exists,
                                              time_exists,
                                              bound_state_method_exists,
                                              curvilinear_state_method_exists,
                                              update_method_exists,
                                              covariance_transport_exists,
                                              update_surface_exists,
                                              set_step_size_exists,
                                              release_step_size_exists,
                                              output_step_size_exists>;
      };
  // clang-format on
  }  // namespace Stepper
}  // namespace concept

template <typename stepper, typename state = typename stepper::State>
constexpr bool StepperConcept =
    Acts::concept ::Stepper::StepperConcept<stepper, state>::value;
template <typename stepper>
constexpr bool StepperStateConcept =
    Acts::concept ::Stepper::StepperStateConcept<stepper>;
}  // namespace Acts