// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Include all core actors
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/tuple.hpp"
#include "detray/utils/tuple_helpers.hpp"

// System include(s)
#include <type_traits>

namespace detray::detail {

/// Extract the tuple of actor states from an actor type
/// @{
// Simple actor: No observers
template <typename actor_t>
struct get_state_tuple {
 private:
  using state_t = typename actor_t::state;

  // Remove empty default state of base actor type from tuple
  using principal = std::conditional_t<std::same_as<state_t, base_actor::state>,
                                       dtuple<>, dtuple<state_t>>;
  using principal_ref =
      std::conditional_t<std::same_as<state_t, base_actor::state>, dtuple<>,
                         dtuple<state_t &>>;

 public:
  using type = principal;
  using ref_type = principal_ref;
};

// Composite actor: Has observers
template <typename actor_t>
  requires(!std::same_as<typename std::remove_cvref_t<actor_t>::observer_states,
                         void>)
struct get_state_tuple<actor_t> {
 private:
  using principal_actor_t = typename actor_t::actor_type;

  using principal = typename get_state_tuple<principal_actor_t>::type;
  using principal_ref = typename get_state_tuple<principal_actor_t>::ref_type;

  using observers = typename actor_t::observer_states;
  using observer_refs = typename actor_t::observer_state_refs;

 public:
  using type = detail::tuple_cat_t<principal, observers>;
  using ref_type = detail::tuple_cat_t<principal_ref, observer_refs>;
};

/// Tuple of state types
template <typename actor_t>
using state_tuple_t = get_state_tuple<actor_t>::type;

/// Tuple of references
template <typename actor_t>
using state_ref_tuple_t = get_state_tuple<actor_t>::ref_type;
/// @}

}  // namespace detray::detail
