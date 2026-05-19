// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project includes(s)
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/base_actor.hpp"

// System include(s)
#include <iostream>
#include <sstream>
#include <string>

namespace {

// Propagator state: Does nothing
struct empty_prop_state {};
const empty_prop_state prop_state{};

/// Actor that prints its call chain and subject data
struct print_actor : public detray::base_actor {
  /// State keeps an internal string representation
  struct state {
    std::stringstream stream{};

    std::string to_string() const { return stream.str(); }
  };

  // Actor result: status code
  using result = detray::actor::result;

  /// Actor implementation: append call notification to internal string
  template <typename propagator_state_t>
  result operator()(state &printer_state,
                    const propagator_state_t & /*p_state*/,
                    const detray::actor::result /*res*/) const {
    printer_state.stream << "[print actor]:";

    return {detray::actor::status::e_notify};
  }

  /// Observing actor implementation: append call notification to internal
  /// string
  template <typename subj_result_t, typename propagator_state_t>
  result operator()(state &printer_state,
                    const propagator_state_t & /*p_state*/,
                    const subj_result_t res) const {
    printer_state.stream << "[" << res << ": print actor obs "
                         << res.buffer->back() << "]:";

    return {detray::actor::status::e_notify};
  }
};

/// Example actor that counts the number of elements in its buffer
template <template <typename...> class vector_t>
struct example_actor : public detray::base_actor {
  /// Actor state
  struct state {
    // Keep dynamic data per propagation stream
    vector_t<float> buffer = {};
  };

  /// Actor result
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
    example_state.buffer.push_back(static_cast<float>(res.buffer->size()) /
                                   10.f);

    return {{detray::actor::status::e_notify}, &example_state.buffer};
  }

  /// Observing actor implementation to printer: do nothing
  template <typename subj_result_t, typename propagator_state_t>
    requires(!std::is_same_v<subj_result_t, result>)
  result operator()(state &example_state,
                    const propagator_state_t & /*p_state*/,
                    const subj_result_t & /*res*/) const {
    return {{detray::actor::status::e_notify}, &example_state.buffer};
  }
};

}  // anonymous namespace

// Run the actor chain on some dummy actor types
int main() {
  using example_actor_t = example_actor<std::vector>;
  // Implements example_actor with two print observers
  using composite =
      detray::composite_actor<example_actor_t, print_actor, print_actor>;

  std::clog << "Actor Definition Tutorial\n=========================\n";

  example_actor_t::state example_state{};
  print_actor::state printer_state{};

  // Aggregate actor states to be able to pass them through the chain
  auto actor_states = detray::tie(example_state, printer_state);

  // Chain of actors
  using actor_chain_t = detray::actor_chain<example_actor_t, composite>;
  // Run
  actor_chain_t run_actors{};
  run_actors(actor_states, prop_state);

  std::clog << "\nactor list: " << printer_state.to_string() << std::endl;

  // Test chaining of multiple actors

  /* Test chaining of multiple actors
   * The chain goes as follows (depth first):
   *                          example_actor1
   *                              1.|
   *                          observer_lvl1 (print)
   *                              2.|
   *                          observer_lvl2 (example_actor observing print
   * actor)
   *                      3./     5.|     6.\
   *            observer_lvl3 example_actor2 print
   *          (example_actor3)
   *               4.|
   *               print
   */
  using observer_lvl3 = detray::composite_actor<example_actor_t, print_actor>;
  using observer_lvl2 = detray::composite_actor<example_actor_t, observer_lvl3,
                                                example_actor_t, print_actor>;
  using observer_lvl1 = detray::composite_actor<print_actor, observer_lvl2>;
  using chain = detray::composite_actor<example_actor_t, observer_lvl1>;

  // Reset example actor state
  example_state.buffer.clear();
  printer_state.stream.str("");
  printer_state.stream.clear();

  // Run the chain
  detray::actor_chain<chain> run_chain{};
  run_chain(actor_states, prop_state);

  std::clog << "actor chain: " << printer_state.to_string() << std::endl;
}
