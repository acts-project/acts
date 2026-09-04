// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/units.hpp"

// System include(s).
#include <limits>
#include <ostream>

namespace detray::stepping {

enum class id {
  // False for non-charged tracks
  e_linear = 0,
  // True for charged tracks
  e_rk = 1,
};

struct config {
  /// Minimum step size
  float min_stepsize{1e-4f * unit<float>::mm};
  /// Runge-Kutta numeric error tolerance
  float rk_error_tol{1e-4f * unit<float>::mm};
  /// Step size constraint
  float step_constraint{std::numeric_limits<float>::max()};
  /// Maximal path length of track
  float path_limit{5.f * unit<float>::m};
  /// Maximum number of Runge-Kutta step trials
  std::size_t max_rk_updates{10000u};
  /// Use mean energy loss (Bethe)
  /// if false, most probable energy loss (Landau) will be used
  bool use_mean_loss{true};
  /// Use eloss gradient in error propagation
  bool use_eloss_gradient{false};
  /// Use b field gradient in error propagation
  bool use_field_gradient{false};
  /// Do covariance transport
  bool do_covariance_transport{true};

  /// Print the stepping configuration
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& out, const config& cfg) {
    out << "  Min. Stepsize         : "
        << cfg.min_stepsize / detray::unit<float>::mm << " [mm]\n"
        << "  Runge-Kutta tolerance : "
        << cfg.rk_error_tol / detray::unit<float>::mm << " [mm]\n"
        << "  Max. step updates     : " << cfg.max_rk_updates << "\n"
        << "  Stepsize  constraint  : "
        << cfg.step_constraint / detray::unit<float>::mm << " [mm]\n"
        << "  Path limit            : "
        << cfg.path_limit / detray::unit<float>::m << " [m]\n"
        << std::boolalpha << "  Use Bethe energy loss : " << cfg.use_mean_loss
        << "\n"
        << "  Do cov. transport     : " << cfg.do_covariance_transport << "\n";

    if (cfg.do_covariance_transport) {
      out << std::boolalpha
          << "  Use eloss gradient    : " << cfg.use_eloss_gradient << "\n"
          << "  Use B-field gradient  : " << cfg.use_field_gradient << "\n";
    }
    // Reset state
    out << std::noboolalpha;

    return out;
  }
};

}  // namespace detray::stepping
