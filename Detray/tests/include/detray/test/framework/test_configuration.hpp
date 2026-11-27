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
#include <string>
#include <string_view>

namespace detray::test {

/// Test configuration type
template <concepts::scalar scalar_t, typename derived_config_t = void>
struct configuration {
 private:
  /// Cast to the derived type to return the correct reference from setters
  /// @{
  /// @returns reference to the derived type
  DETRAY_HOST_DEVICE
  constexpr auto cast_impl() {
    if constexpr (!std::same_as<derived_config_t, void>) {
      return static_cast<derived_config_t&>(*this);
    } else {
      return *this;
    }
  }
  /// @returns reference to the derived type - const
  DETRAY_HOST_DEVICE
  constexpr auto cast_impl() const {
    if constexpr (!std::same_as<derived_config_t, void>) {
      return static_cast<const derived_config_t&>(*this);
    } else {
      return *this;
    }
  }
  /// @}

 public:
  /// Default constructor
  configuration() = default;

  /// Default constructor
  /// @param test_name the name of the test
  configuration(const std::string_view test_name) : m_name{test_name} {};

  /// Propagation
  /// @{
  propagation::config m_prop_cfg{};
  /// @}

  /// General testing
  /// @{
  /// Description of the test
  std::string m_name{"detray_test"};
  /// Tolerance to compare two floating point values
  scalar_t m_tolerance{1e-5f};
  /// Shorthand for infinity
  scalar_t inf{std::numeric_limits<scalar_t>::infinity()};
  /// Shorthand for the floating point epsilon
  scalar_t epsilon{std::numeric_limits<scalar_t>::epsilon()};
  /// @}

  /// Setters
  /// @{
  /// Set the test name.
  /// @param n The configured test name.
  auto name(const std::string& n) {
    m_name = n;
    return cast_impl();
  }
  /// Set the tolerance for floating point comparisons.
  /// @param t New comparison tolerance.
  /// @return This configuration.
  auto tol(scalar_t t) {
    m_tolerance = t;
    return cast_impl();
  }
  /// @}

  /// Getters
  /// @{
  /// Get the test name.
  /// @return The configured test name.
  const std::string& name() const { return m_name; }
  /// Get the tolerance for floating point comparisons.
  /// @return The configured comparison tolerance.
  scalar_t tol() const { return m_tolerance; }
  /// Get the propagation configuration.
  /// @return Mutable propagation configuration.
  propagation::config& propagation() { return m_prop_cfg; }
  /// Get the propagation configuration.
  /// @return Const propagation configuration.
  const propagation::config& propagation() const { return m_prop_cfg; }
  /// @}

 private:
  /// Print configuration
  friend std::ostream& operator<<(std::ostream& out, const configuration& cfg) {
    out << "Test Configuration\n"
        << "----------------------------\n"
        << "  name                  : " << cfg.name() << "\n"
        << "  test tolerance        : " << cfg.tol() << "\n\n";
    out << cfg.propagation();

    return out;
  }
};

}  // namespace detray::test
