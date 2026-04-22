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
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <limits>

namespace detray::detail {

template <concepts::scalar scalar_t, typename = void>
class quadratic_equation {};

/// Class to solve a quadratic equation of type a * x^2 + b * x + c = 0
///
/// @note If there are no real solutions, the result is undefined
/// @note The solutions are sorted by default. If there is only one solution,
/// the larger value is undefined.
template <concepts::scalar scalar_t>
  requires std::is_arithmetic_v<scalar_t>
class quadratic_equation<scalar_t> {
 public:
  quadratic_equation() = delete;

  /// Solve the quadratic equation with the coefficients @param a, @param b
  /// and @param c
  ///
  /// @param tolerance threshold to compare the discrimant against to decide
  ///                  if we have two separate solutions.
  DETRAY_HOST_DEVICE
  constexpr quadratic_equation(
      const scalar_t a, const scalar_t b, const scalar_t c,
      const scalar_t tolerance = std::numeric_limits<scalar_t>::epsilon()) {
    // linear case
    if (math::fabs(a) <= tolerance) {
      if (math::fabs(b) <= tolerance) {
        m_solutions = 0;
      } else {
        m_solutions = 1;
        m_values[0] = -c / b;
      }
    } else {
      const scalar_t discriminant{b * b - 4.f * a * c};
      // If there is more than one solution, then a != 0 and q != 0
      if (discriminant > tolerance) {
        m_solutions = 2;
        const scalar_t q{-0.5f *
                         (b + detail::copysign(math::sqrt(discriminant), b))};
        m_values = {q / a, c / q};
        // Sort the two solutions
        if (m_values[0] > m_values[1]) {
          m_values = {m_values[1], m_values[0]};
        }
      }
      // Only one solution and a != 0
      else if (discriminant >= 0.f) {
        m_solutions = 1;
        m_values[0] = -0.5f * b / a;
      }
      // discriminant < 0 is not allowed, since all solutions should be
      // real
    }
  }

  /// Getters for the solution(s)
  /// @{
  constexpr int solutions() const { return m_solutions; }
  constexpr scalar_t smaller() const { return m_values[0]; }
  constexpr scalar_t larger() const { return m_values[1]; }
  /// @}

 private:
  /// Number of solutions of the equation
  int m_solutions{0};
  /// The solutions
  darray<scalar_t, 2> m_values{detail::invalid_value<scalar_t>(),
                               detail::invalid_value<scalar_t>()};
};

/// Class to solve a quadratic equation of type a * x^2 + b * x + c = 0
///
/// @note If there are no real solutions, the result is undefined
/// @note The solutions are sorted by default. If there is only one
/// solution, the larger value is undefined.
template <concepts::scalar scalar_t>
  requires(!std::is_arithmetic_v<scalar_t>)
class quadratic_equation<scalar_t> {
 public:
  quadratic_equation() = delete;

  DETRAY_HOST_DEVICE
  constexpr quadratic_equation(const scalar_t &a, const scalar_t &b,
                               const scalar_t &c,
                               const scalar_t &tolerance = 1e-6f) {
    // Linear case
    auto one_sol = (math::fabs(a) <= tolerance);
    m_solutions(one_sol) = 1.f;
    m_values[0] = -c / b;

    // Early exit
    if (detray::detail::all_of(one_sol)) {
      return;
    }

    const scalar_t discriminant = b * b - (4.f * a) * c;

    const auto two_sol = (discriminant > tolerance);
    one_sol = !two_sol && (discriminant >= 0.f);

    // If there is more than one solution, then a != 0 and q != 0
    if (detray::detail::any_of(two_sol)) {
      m_solutions = 2.f;
      m_solutions.setZeroInverted(two_sol);

      const scalar_t q =
          -0.5f * (b + math::copysign(math::sqrt(discriminant), b));

      scalar_t first = q / a;
      scalar_t second = c / q;
      first.setZeroInverted(two_sol);
      second.setZeroInverted(two_sol);

      // Sort the solutions
      const auto do_swap = (second < first);
      if (detray::detail::all_of(do_swap)) {
        m_values = {second, first};
      } else if (detray::detail::none_of(do_swap)) {
        m_values = {first, second};
      } else {
        const auto tmp = second;
        second(do_swap) = first;
        first(do_swap) = tmp;
        m_values = {first, second};
      }
    }

    // Only one solution and a != 0
    if (detray::detail::any_of(one_sol)) {
      scalar_t sol = 1.f;
      scalar_t result = -0.5f * b / a;
      sol.setZeroInverted(one_sol);
      result.setZeroInverted(one_sol);

      m_solutions += sol;
      m_values[0] += result;
    }
    // discriminant < 0 is not allowed, since all solutions should
    // be real
  }

  /// Getters for the solution(s)
  /// @{
  constexpr const auto &solutions() const { return m_solutions; }
  constexpr const scalar_t &smaller() const { return m_values[0]; }
  constexpr const scalar_t &larger() const { return m_values[1]; }
  /// @}

 private:
  /// Number of solutions of the equation (needs to be floating point to
  /// apply the masks correctly)
  scalar_t m_solutions = 0.f;
  /// The solutions
  darray<scalar_t, 2> m_values{static_cast<scalar_t>(0.f),
                               static_cast<scalar_t>(0.f)};
};

template <typename S>
quadratic_equation(const S a, const S &b, const S &c, const S &tolerance)
    -> quadratic_equation<S>;

}  // namespace detray::detail
