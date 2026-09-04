// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/propagator/propagation_config.hpp"

// System include(s)
#include <limits>
#include <ostream>

namespace detray::test {

/// Test configuration type
template <concepts::scalar scalar_t>
struct configuration {
  /// General testing
  /// @{
  /// Tolerance to compare two floating point values
  scalar_t m_tolerance{1e-5f};
  /// Shorthand for infinity
  scalar_t inf{std::numeric_limits<scalar_t>::infinity()};
  /// Shorthand for the floating point epsilon
  scalar_t epsilon{std::numeric_limits<scalar_t>::epsilon()};
  /// @}

  /// Propagation
  /// @{
  propagation::config m_prop_cfg{};
  /// @}

  /// Setters
  /// @{
  configuration& tol(scalar_t t) {
    m_tolerance = t;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  scalar_t tol() const { return m_tolerance; }
  propagation::config& propagation() { return m_prop_cfg; }
  const propagation::config& propagation() const { return m_prop_cfg; }
  /// @}

 private:
  /// Print configuration
  friend std::ostream& operator<<(std::ostream& out, const configuration& cfg) {
    out << "Test Configuration\n"
        << "----------------------------\n"
        << "  test tolerance        : " << cfg.tol() << "\n\n";
    out << cfg.propagation();

    return out;
  }
};

}  // namespace detray::test
