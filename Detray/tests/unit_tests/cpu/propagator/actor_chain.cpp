// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s).
#include "detray/propagator/actor_chain.hpp"

#include "detray/definitions/units.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/concepts.hpp"

// GTest include(s).
#include <gtest/gtest.h>

// System include(s).
#include <sstream>
#include <string>

using namespace detray;

/// Actor that prints its call chain and subject data
struct print_actor : public detray::base_actor {
  /// State keeps an internal string representation
  struct state {
    std::stringstream stream{};

    state() = default;
    state(state &&) noexcept = default;
    state(const state &) = delete;
    state &operator=(state const &) = delete;
    state &operator=(state &&) noexcept = default;

    std::string to_string() const { return stream.str(); }
  };

  // Actor result: status code
  using result = detray::actor::result;

  /// Actor implementation: append call notification to internal string
  template <typename propagator_state_t>
  result operator()(state &printer_state,
                    const propagator_state_t & /*p_state*/) const {
    printer_state.stream << "[print actor]:";

    return {detray::actor::status::e_notify};
  }

  /// Observing actor implementation: append call notification to internal
  /// string
  template <typename subj_result_t, typename propagator_state_t>
  result operator()(state &printer_state,
                    const propagator_state_t & /*p_state*/,
                    const subj_result_t &res) const {
    printer_state.stream << "[print actor obs " << res.buffer->back() << "]:";

    return {detray::actor::status::e_notify};
  }
};

/// Example actor that couts the number of elements in its buffer
template <template <typename...> class vector_t>
struct example_actor : public detray::base_actor {
  /// Actor state
  struct state {
    // Keep dynamic data per propagation stream
    vector_t<float> buffer = {};
  };

  // Actor result
  struct result : detray::actor::result {
    vector_t<float> *buffer{nullptr};
  };

  /// Actor implementation: Counts vector elements
  template <typename propagator_state_t>
  result operator()(state &example_state,
                    const propagator_state_t & /*p_state*/) const {
    example_state.buffer.push_back(
        static_cast<float>(example_state.buffer.size()));

    return {{detray::actor::status::e_notify}, &example_state.buffer};
  }

  /// Observing actor implementation: Counts vector elements (division)
  template <typename propagator_state_t>
  result operator()(state &example_state,
                    const propagator_state_t & /*p_state*/,
                    const result &res) const {
    example_state.buffer.push_back(static_cast<float>(res.buffer->size()) *
                                   0.1f);

    return {{detray::actor::status::e_notify}, &example_state.buffer};
  }

  /// Observing actor implementation to printer: do nothing
  template <typename subj_result_t, typename propagator_state_t>
    requires(!std::is_same_v<subj_result_t, state>)
  result operator()(state &example_state,
                    const propagator_state_t & /*p_state*/,
                    const subj_result_t & /*subject_result*/) const {
    return {{detray::actor::status::e_notify}, &example_state.buffer};
  }
};

using example_actor_t = example_actor<std::vector>;
// Implements example_actor with two print observers
using composite1 = composite_actor<example_actor_t, print_actor, print_actor>;
// Implements example_actor with one print observer
using composite2 = composite_actor<example_actor_t, print_actor>;
// Implements example_actor through composite2 and has composite1 as observer
using composite3 = composite_actor<example_actor_t, composite1>;
// Implements example_actor through composite2<-composite3 with composite1 obs.
using composite4 = composite_actor<example_actor_t, composite1>;

template <int I>
struct ordered_actor : public detray::base_actor {
  struct state {};
};

using ordered_actor_chain =
    actor_chain<ordered_actor<0>, ordered_actor<1>, ordered_actor<2>>;

/* Test chaining of multiple actors
 * The chain goes as follows (depth first):
 *                          example_actor1
 *                              1.|
 *                          observer_lvl1 (print)
 *                              2.|
 *                          observer_lvl2 (example_actor observing print actor)
 *                      3./     5.|     6.\
 *            observer_lvl3 example_actor2 print
 *          (example_actor3)
 *               4.|
 *               print
 */
using observer_lvl3 = composite_actor<example_actor_t, print_actor>;
using observer_lvl2 = composite_actor<example_actor_t, observer_lvl3,
                                      example_actor_t, print_actor>;
using observer_lvl1 = composite_actor<print_actor, observer_lvl2>;
using chain = composite_actor<example_actor_t, observer_lvl1>;

// Test the actor chain on some dummy actor types
GTEST_TEST(detray_propagator, actor_chain) {
  static_assert(detray::concepts::actor<print_actor>);
  static_assert(detray::concepts::actor<example_actor_t>);
  static_assert(detray::concepts::actor<composite1>);
  static_assert(detray::concepts::actor<composite2>);
  static_assert(detray::concepts::actor<composite3>);
  static_assert(detray::concepts::actor<composite4>);
  static_assert(detray::concepts::actor<observer_lvl3>);
  static_assert(detray::concepts::actor<observer_lvl2>);
  static_assert(detray::concepts::actor<observer_lvl1>);
  static_assert(detray::concepts::actor<ordered_actor<0>>);
  static_assert(detray::concepts::actor<ordered_actor<1>>);
  static_assert(detray::concepts::actor<ordered_actor<2>>);

  static_assert(!detray::concepts::composite_actor<print_actor>);
  static_assert(!detray::concepts::composite_actor<example_actor_t>);
  static_assert(detray::concepts::composite_actor<composite1>);
  static_assert(detray::concepts::composite_actor<composite2>);
  static_assert(detray::concepts::composite_actor<composite3>);
  static_assert(detray::concepts::composite_actor<composite4>);
  static_assert(detray::concepts::composite_actor<observer_lvl3>);
  static_assert(detray::concepts::composite_actor<observer_lvl2>);
  static_assert(detray::concepts::composite_actor<observer_lvl1>);
  static_assert(
      std::same_as<typename ordered_actor_chain::state_tuple,
                   dtuple<ordered_actor<0>::state, ordered_actor<1>::state,
                          ordered_actor<2>::state>>);
  static_assert(
      std::same_as<typename ordered_actor_chain::state_ref_tuple,
                   dtuple<ordered_actor<0>::state &, ordered_actor<1>::state &,
                          ordered_actor<2>::state &>>);

  // The actor states (can be reused between actors)
  example_actor_t::state example_state{};
  print_actor::state printer_state{};

  // Aggregate actor states to be able to pass them through the chain
  auto actor_states = detray::tie(example_state, printer_state);

  // Propagator state
  struct empty_prop_state {};
  empty_prop_state prop_state{};

  // Chain of actors
  using actor_chain_t = actor_chain<example_actor_t, composite1, composite2,
                                    composite3, composite4>;

  static_assert(detray::concepts::actor_chain<actor_chain<>>);
  static_assert(detray::concepts::actor_chain<actor_chain_t>);

  // Run
  actor_chain_t run_actors{};
  run_actors(actor_states, prop_state);

  ASSERT_TRUE(printer_state.to_string().compare(
                  "[print actor obs 1]:[print actor obs 1]:[print actor obs "
                  "2]:[print actor obs 0.4]:[print actor obs 0.4]:[print "
                  "actor obs 0.6]:[print actor obs 0.6]:") == 0)
      << "Printer call chain: " << printer_state.to_string() << std::endl;

  // Test chaining of multiple actors

  // Reset example actor state
  example_state.buffer.clear();
  printer_state.stream.str("");
  printer_state.stream.clear();

  // Run the chain
  actor_chain<chain> run_chain{};
  run_chain(actor_states, prop_state);

  ASSERT_TRUE(printer_state.to_string().compare(
                  "[print actor obs 0]:[print actor obs 0.1]:[print actor "
                  "obs 0.2]:") == 0)
      << "Printer call chain: " << printer_state.to_string() << std::endl;
}
