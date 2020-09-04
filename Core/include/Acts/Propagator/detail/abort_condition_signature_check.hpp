// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/detail/MPL/type_collector.hpp"

#include <type_traits>

namespace Acts {

/// The following operators have to be implemented in order to satisfy
/// as an abort condition
///
/// clang-format off
///
/// IF the aborter declares an `action_type` upon whose result it will depend:
/// @code
///
/// template <typename propagator_state_t, typename stepper_t,
///           typename result_t>
/// bool
/// operator()(propagator_state_t& state,
///            const stepper_t, const result_t& r) const
/// {
///   return false;
/// }
///
/// @endcode
///
/// IF the aborter does NOT declare an `action_type`:
///
/// @code
///
/// template <typename propagator_state_t, typename stepper_t>
/// bool
/// operator()(propagator_state_t& state, const stepper_t& stepper) const
/// {
///   return false;
/// }
///
/// @endcode
///
/// clang-format off

namespace Concepts {
namespace detail_aborter {

/// Detection helper for call operator WITHOUT result
template <typename A, typename propagator_state_t, typename stepper_t>
using call_op_no_result_t = decltype(std::declval<const A>()(
    std::declval<propagator_state_t&>(), std::declval<const stepper_t&>()));

/// Detection helper for call operator WITH result
template <typename A, typename result_t, typename propagator_state_t,
          typename stepper_t>
using call_op_with_result_t = decltype(std::declval<const A>()(
    std::declval<propagator_state_t&>(), std::declval<const stepper_t&>(),
    std::declval<const result_t&>()));

// This is basically an if:
// if ( !Aborter.hasResult() ) { // has no result
template <typename T, typename propagator_state_t, typename stepper_t,
          bool has_result = false>
struct ConceptConditional {
  // check the existence of the correct call operator
  constexpr static bool value =
      Acts::Concepts ::exists<call_op_no_result_t, T, propagator_state_t,
                              stepper_t>;
};

// } else { // has a result
template <typename T, typename propagator_state_t, typename stepper_t>
struct ConceptConditional<T, propagator_state_t, stepper_t, true> {
  // unpack the result type from the action contained in the aborter type
  using result_type =
      Acts::detail::result_type_t<Acts::detail::action_type_t<T>>;
  // check the existence of the correct call operator
  constexpr static bool value =
      Acts::Concepts ::exists<call_op_with_result_t, T, result_type,
                              propagator_state_t, stepper_t>;
};
// } // endif

// Calls the 'if' above, depending on the value of `has_action_type_v`.
template <typename T, typename propagator_state_t, typename stepper_t>
struct Concept {
  constexpr static bool value =
      ConceptConditional<T, propagator_state_t, stepper_t,
                         Acts::detail::has_action_type_v<T>>::value;
};

}  // namespace detail_aborter

/// Meta function for checking if an aborter has a valid interface
template <typename T, typename propagator_state_t, typename stepper_t>
constexpr bool abort_condition_signature_check_v =
    detail_aborter::Concept<T, propagator_state_t, stepper_t>::value;
}  // namespace Concepts
}  // namespace Acts
