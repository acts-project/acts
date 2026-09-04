// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/utils/logging.hpp"

// Vecmem include(s)
#include <vecmem/containers/device_vector.hpp>

// System includes
#include <utility>

namespace detray::actor {

template <typename sf_descriptor_t>
struct surface_sequencer : public base_actor {
  struct state {
    using surface_type = sf_descriptor_t;
    using sequence_type = vecmem::device_vector<sf_descriptor_t>;

    /// Constructor with the vector of track states
    DETRAY_HOST_DEVICE
    explicit state(sequence_type seq) : m_sequence(seq) {}

    /// Flags that the internal buffer has experienced an overflow
    DETRAY_HOST_DEVICE void set_overflow() { m_overflow = true; }

    /// @returns true if the internal buffer has reached its capacity
    DETRAY_HOST_DEVICE bool overflow() const { return m_overflow; }

    /// Configure whether to collect only material surfaces
    DETRAY_HOST_DEVICE void collect_only_material_surfaces(bool material_only) {
      m_collect_only_material_surfaces = material_only;
    }

    /// @returns whether to only collect material surfaces
    DETRAY_HOST_DEVICE bool collect_only_material_surfaces() const {
      return m_collect_only_material_surfaces;
    }

    /// @returns access to the surface sequence buffer
    DETRAY_HOST_DEVICE const sequence_type& sequence() const {
      return m_sequence;
    }

    /// @returns access to the surface sequence buffer
    DETRAY_HOST_DEVICE sequence_type& sequence() { return m_sequence; }

   private:
    /// Sequence of surfaces the navigation encountered
    sequence_type m_sequence;
    /// Internal buffer overflow
    bool m_overflow = false;
    /// Whether to collect only surfaces with material
    bool m_collect_only_material_surfaces = false;
  };

  /// Record detector surfaces in sequence if they are sensitive or
  /// have material
  template <typename propagator_state_t>
    requires std::same_as<
        sf_descriptor_t,
        typename propagator_state_t::detector_type::surface_type>
  DETRAY_HOST_DEVICE void operator()(state& actor_state,
                                     propagator_state_t& propagation) const {
    const auto& navigation = propagation.navigation();

    DETRAY_VERBOSE_HOST_DEVICE(
        "Actor: Checking for next surface in sequence...");

    if (actor_state.collect_only_material_surfaces() &&
        navigation.is_on_sensitive()) {
      return;
    } else if (!(navigation.is_on_sensitive() ||
                 navigation.encountered_sf_material())) {
      return;
    }

    if (actor_state.sequence().size() == actor_state.sequence().capacity()) {
      if (actor_state.sequence().capacity() == 0) {
        DETRAY_VERBOSE_HOST_DEVICE(
            "Actor: No capacity to store surface candidate - skipping");
      } else {
        DETRAY_VERBOSE_HOST_DEVICE("Actor: Sequence overflow!");
        actor_state.set_overflow();
        propagation.navigation().exit();
        propagation.heartbeat(false);
      }
      return;
    }

    const auto& sf_desc = std::as_const(navigation).current().surface();
    assert(!sf_desc.identifier().is_invalid());

    actor_state.sequence().push_back(sf_desc);
    DETRAY_VERBOSE_HOST("Actor: Added surface: " << sf_desc);
  }
};

}  // namespace detray::actor
