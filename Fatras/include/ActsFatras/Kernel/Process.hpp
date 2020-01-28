// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>

#include "Acts/Utilities/Definitions.hpp"

namespace ActsFatras {

/// @brief Process class that turns a parameterized or
/// tabularized fast simulation module into a process that
/// can be plugged into the PhysicsList
///
/// This is plugin for physics processes
///  - scattering
///  - energy loss
///  - pair production
///  - hadronic interaction
///  - decay
///
/// The type (and actual trigger) of the particle
/// and interaction is steered via the Selector list
/// for in and out.
template <typename physics_t, typename selector_in_t, typename selector_out_t,
          typename selector_child_t>

struct Process {
  /// The actual physics that is happening
  physics_t process;

  /// The selector list
  selector_in_t selectorIn;
  selector_out_t selectorOut;
  selector_child_t selectorChild;

  /// This is the scattering call operator
  template <typename generator_t, typename detector_t, typename particle_t>
  bool operator()(generator_t &gen, const detector_t &det, particle_t &in,
                  std::vector<particle_t> &out) const {
    // check if the process applies
    if (selectorIn(det, in)) {
      // apply energy loss and get eventual children
      auto children = process(gen, det, in);
      if (children.size()) {
        // copy the children that comply with the child selector
        std::copy_if(
            children.begin(), children.end(), std::back_inserter(out),
            [this, det](const particle_t &p) { return selectorChild(det, p); });
      }
    }
    // check if this killed the particle,
    // or pushed below threshold
    return (!selectorOut(det, in));
  }
};

}  // namespace ActsFatras
