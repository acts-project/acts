// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/units.hpp"
#include "detray/propagator/propagator.hpp"

// Detray test include(s)
#include "detray/test/framework/test_configuration.hpp"
#include "detray/test/framework/types.hpp"

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>
#include <ostream>
#include <string>

namespace detray::test {

/// Base type for test fixtures with google test
template <typename scope = ::testing::Test>
class fixture_base : public scope {
 public:
  /// Linear algebra typedefs
  /// @{
  using algebra_t = test::algebra;
  using index = test::index;
  using scalar = test::scalar;
  using point2 = test::point2;
  using point3 = test::point3;
  using vector2 = test::vector2;
  using vector3 = test::vector3;
  using transform3 = test::transform3;
  template <index ROW, index COL>
  using matrix = test::matrix<ROW, COL>;
  /// @}

  using configuration = detray::test::configuration<scalar>;

  /// Constructor
  explicit fixture_base(const configuration& cfg = {})
      : m_tolerance{cfg.tol()},
        m_inf{cfg.inf},
        m_epsilon{cfg.epsilon},
        m_path_limit{cfg.propagation().stepping.path_limit},
        m_overstep_tolerance{
            cfg.propagation().navigation.intersection.overstep_tolerance},
        m_step_constraint{cfg.propagation().stepping.step_constraint} {}

  /// Default destructor
  ~fixture_base() override = default;

  /// @returns the benchmark name
  std::string name() const { return "detray_test"; };

 protected:
  static void SetUpTestSuite() { /* Do nothing */ }
  static void TearDownTestSuite() { /* Do nothing */ }
  void SetUp() override { /* Do nothing */ }
  void TearDown() override { /* Do nothing */ }

  /// Getters
  /// @{
  scalar tolerance() const { return m_tolerance; }
  scalar inf() const { return m_inf; }
  scalar epsilon() const { return m_epsilon; }
  scalar path_limit() const { return m_path_limit; }
  scalar overstep_tolerance() const { return m_overstep_tolerance; }
  scalar step_constraint() const { return m_step_constraint; }
  /// @}

 private:
  scalar m_tolerance{};
  scalar m_inf{};
  scalar m_epsilon{};
  scalar m_path_limit{};
  scalar m_overstep_tolerance{};
  scalar m_step_constraint{};
};

}  // namespace detray::test
