// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Propagator/ConstrainedStep.hpp"
#include "Acts/Surfaces/BoundaryCheck.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {
class Surface;

namespace Concepts::Stepper {

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

METHOD_TRAIT(reset_state_t, resetState);
METHOD_TRAIT(get_field_t, getField);
METHOD_TRAIT(position_t, position);
METHOD_TRAIT(direction_t, direction);
METHOD_TRAIT(qop_t, qOverP);
METHOD_TRAIT(absolute_momentum_t, absoluteMomentum);
METHOD_TRAIT(momentum_t, momentum);
METHOD_TRAIT(charge_t, charge);
METHOD_TRAIT(time_t, time);
METHOD_TRAIT(overstep_t, overstepLimit);
METHOD_TRAIT(bound_state_method_t, boundState);
METHOD_TRAIT(curvilinear_state_method_t, curvilinearState);
METHOD_TRAIT(update_t, update);
METHOD_TRAIT(covariance_transport_bound_t, transportCovarianceToBound);
METHOD_TRAIT(covariance_transport_curvilinear_t,
             transportCovarianceToCurvilinear);
METHOD_TRAIT(step_t, step);
METHOD_TRAIT(update_surface_status_t, updateSurfaceStatus);
METHOD_TRAIT(update_step_size_t, updateStepSize);
METHOD_TRAIT(get_step_size_t, getStepSize);
METHOD_TRAIT(release_step_size_t, releaseStepSize);
METHOD_TRAIT(output_step_size_t, outputStepSize);

template <typename T>
using cov_transport_t = decltype(std::declval<T>().covTransport);
template <typename T>
using cov_t = decltype(std::declval<T>().cov);
template <typename T>
using path_accumulated_t = decltype(std::declval<T>().pathAccumulated);
template <typename T>
using step_size_t = decltype(std::declval<T>().stepSize);

// clang-format off
    template <typename S>
    constexpr bool StepperStateConcept
      = require<has_member<S, cov_transport_t, bool>,
                has_member<S, cov_t, BoundSquareMatrix>,
                has_member<S, path_accumulated_t, double>//,
//                 has_member<S, step_size_t, ConstrainedStep>
               >;
// clang-format on

// clang-format off
template <typename S>
constexpr bool MultiStepperStateConcept= require<
  has_member<S, cov_transport_t, bool>,
  has_member<S, path_accumulated_t, double>
>;
// clang-format on

// clang-format off
    template <typename S, typename state = typename S::State>
      struct CommonStepperConcept {
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
        constexpr static bool reset_state_exists = has_method<const S, void, reset_state_t, state&, const BoundVector&, const BoundSquareMatrix&, const Surface&, const double>;
        static_assert(reset_state_exists, "resetState method not found");
        constexpr static bool position_exists = has_method<const S, Vector3, position_t, const state&>;
        static_assert(position_exists, "position method not found");
        constexpr static bool direction_exists = has_method<const S, Vector3, direction_t, const state&>;
        static_assert(direction_exists, "direction method not found");
        constexpr static bool qop_exists = has_method<const S, double, qop_t, const state&>;
        static_assert(qop_exists, "qOverP method not found");
        constexpr static bool absolute_momentum_exists = has_method<const S, double, absolute_momentum_t, const state&>;
        static_assert(absolute_momentum_exists, "absoluteMomentum method not found");
        constexpr static bool momentum_exists = has_method<const S, Vector3, momentum_t, const state&>;
        static_assert(momentum_exists, "momentum method not found");
        constexpr static bool charge_exists = has_method<const S, double, charge_t, const state&>;
        static_assert(charge_exists, "charge method not found");
        constexpr static bool time_exists = has_method<const S, double, time_t, const state&>;
        static_assert(time_exists, "time method not found");
        constexpr static bool overstep_exists = has_method<const S, double, overstep_t, const state&>;
        static_assert(overstep_exists, "overstepLimit method not found");
        constexpr static bool bound_state_method_exists= has_method<const S, Result<typename S::BoundState>, bound_state_method_t, state&, const Surface&, bool, const FreeToBoundCorrection&>;
        static_assert(bound_state_method_exists, "boundState method not found");
        constexpr static bool curvilinear_state_method_exists = has_method<const S, typename S::CurvilinearState, curvilinear_state_method_t, state&, bool>;
        static_assert(curvilinear_state_method_exists, "curvilinearState method not found");
        constexpr static bool covariance_transport_exists = require<has_method<const S, void, covariance_transport_curvilinear_t, state&>,
                                                                    has_method<const S, void, covariance_transport_bound_t, state&, const Surface&, const FreeToBoundCorrection&>>;
        static_assert(covariance_transport_exists, "covarianceTransport method not found");
        constexpr static bool update_surface_exists = has_method<const S, Intersection3D::Status, update_surface_status_t, state&, const Surface&, std::uint8_t, Direction, const BoundaryCheck&, ActsScalar, const Logger&>;
        static_assert(update_surface_exists, "updateSurfaceStatus method not found");
        constexpr static bool update_step_size_exists = has_method<const S, void, update_step_size_t, state&, double, ConstrainedStep::Type, bool>;
        static_assert(update_step_size_exists, "updateStepSize method not found");
        constexpr static bool get_step_size_exists = has_method<const S, double, get_step_size_t, const state &, ConstrainedStep::Type>;
        static_assert(get_step_size_exists, "getStepSize method not found");
        constexpr static bool release_step_size_exists = has_method<const S, void, release_step_size_t, state&, ConstrainedStep::Type>;
        static_assert(release_step_size_exists, "releaseStepSize method not found");
        constexpr static bool output_step_size_exists = has_method<const S, std::string, output_step_size_t, const state&>;
        static_assert(output_step_size_exists, "outputStepSize method not found");

        constexpr static bool value = require<state_exists,
                                              jacobian_exists,
                                              covariance_exists,
                                              bound_state_exists,
                                              curvilinear_state_exists,
                                              position_exists,
                                              direction_exists,
                                              qop_exists,
                                              absolute_momentum_exists,
                                              momentum_exists,
                                              charge_exists,
                                              time_exists,
                                              bound_state_method_exists,
                                              curvilinear_state_method_exists,
                                              covariance_transport_exists,
                                              update_surface_exists,
                                              update_step_size_exists,
                                              release_step_size_exists,
                                              output_step_size_exists>;

      };
// clang-format on

// clang-format off
    // NOTE This static_asserts in here must be commented out, since it would break the compilation for the MultiStepper
    template <typename S, typename state = typename S::State>
      struct SingleStepperConcept {
        constexpr static bool common_stepper_concept_fullfilled = CommonStepperConcept<S, state>::value;
        static_assert(common_stepper_concept_fullfilled, "Stepper does not fulfill common stepper concept");
        constexpr static bool update_method_exists = require<has_method<const S, void, update_t, state&, const FreeVector&, const BoundVector&, const BoundSquareMatrix&, const Surface&>, has_method<const S, void, update_t, state&, const Vector3&, const Vector3&, double, double>>;
        // static_assert(update_method_exists, "update method not found");
        constexpr static bool get_field_exists = has_method<const S, Result<Vector3>, get_field_t, state&, const Vector3&>;
        // static_assert(get_field_exists, "getField method not found");

        constexpr static bool value = require<common_stepper_concept_fullfilled,
                                              update_method_exists,
                                              get_field_exists>;
      };
// clang-format on

// clang-format off
    template <typename S, typename state = typename S::State>
      struct MultiStepperConcept {
        constexpr static bool common_stepper_concept_fullfilled = CommonStepperConcept<S, state>::value;
        static_assert(common_stepper_concept_fullfilled, "Common stepper concept not fulfilled");

        // TODO for now we do not check if the ComponentProxy does fulfill a concept
        template <typename T> using component_proxy_t = typename T::ComponentProxy;
        constexpr static bool component_proxy_exists = exists<component_proxy_t, S>;
        // static_assert(component_proxy_exists, "!component_proxy_exists");

        // TODO for now we do not check if the ConstComponentProxy does fulfill a concept
        template <typename T> using const_component_proxy_t = typename T::ConstComponentProxy;
        constexpr static bool const_component_proxy_exists = exists<const_component_proxy_t, S>;
        // static_assert(const_component_proxy_exists, "!const_component_proxy_exists");

        METHOD_TRAIT(number_components_t, numberComponents);
        constexpr static bool number_components_exists = has_method<const S, std::size_t, number_components_t, const state&>;
        // static_assert(number_components_exists, "!num_components_exists");

        // TODO We cannot check addComponents since it is a template member function

        METHOD_TRAIT(clear_components_t, clearComponents);
        constexpr static bool clear_components_exists = has_method<const S, void, clear_components_t, state&>;
        // static_assert(clear_components_exists, "!clear_components_exists");

        METHOD_TRAIT(remove_missed_components_t, removeMissedComponents);
        constexpr static bool remove_missed_components_exists = has_method<const S, void, remove_missed_components_t, state&>;
        // static_assert(remove_missed_components_exists, "!remove_missed_components_exists");

        constexpr static bool value = require<common_stepper_concept_fullfilled,
                                              component_proxy_exists,
                                              const_component_proxy_exists,
                                              number_components_exists,
                                              clear_components_exists,
                                              remove_missed_components_exists>;
        };
// clang-format on

}  // namespace Concepts::Stepper

template <typename stepper, typename state = typename stepper::State>
constexpr bool StepperConcept =
    Acts::Concepts ::Stepper::SingleStepperConcept<stepper, state>::value ||
    Acts::Concepts ::Stepper::MultiStepperConcept<stepper, state>::value;
template <typename stepper>
constexpr bool StepperStateConcept =
    Acts::Concepts ::Stepper::StepperStateConcept<stepper> ||
    Acts::Concepts ::Stepper::MultiStepperStateConcept<stepper>;
}  // namespace Acts
