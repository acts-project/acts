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
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/navigation/intersection/intersection.hpp"
#include "detray/navigation/intersection/intersection_config.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/logging.hpp"

// System include(s).
#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace detray {

/// @brief Try to find a bracket around a root
///
/// @param [in] a lower initial boundary
/// @param [in] b upper initial boundary
/// @param [in] f function for which to find the root
/// @param [out] bracket bracket around the root
/// @param [in] k scale factor with which to widen the bracket at every step
///
/// @see Numerical Recepes pp. 445
///
/// @return whether a bracket was found
template <concepts::scalar scalar_t, typename function_t>
DETRAY_HOST_DEVICE inline bool expand_bracket(const scalar_t a,
                                              const scalar_t b, function_t &f,
                                              darray<scalar_t, 2> &bracket,
                                              const scalar_t k = 1.f) {
  if (a == b) {
    throw std::invalid_argument(
        "Root bracketing: Not a valid start interval [" + std::to_string(a) +
        ", " + std::to_string(b) + "]");
  }

  scalar_t lower{a > b ? b : a};
  scalar_t upper{a > b ? a : b};

  // Sample function points at interval
  scalar_t f_l{f(lower)};
  scalar_t f_u{f(upper)};
  std::size_t n_tries{0u};

  DETRAY_DEBUG_HOST("Initial bracket: [" << lower << ", " << upper << "]");

  /// Check if the bracket has become invalid
  const auto check_bracket = [a, b, &bracket](std::size_t n, scalar_t fl,
                                              scalar_t fu, scalar_t l,
                                              scalar_t u) {
    if ((n == 1000u) || !std::isfinite(fl) || !std::isfinite(fu) ||
        !std::isfinite(l) || !std::isfinite(u)) {
      DETRAY_VERBOSE_HOST("Could not bracket a root (a="
                          << l << ", b=" << u << ", f(a)=" << fl
                          << ", f(b)=" << fu
                          << ", root might not exist). Running Newton-Raphson "
                             "without bisection.");
      // Reset value
      bracket = {a, b};
      return false;
    }
    return true;
  };

  // If there is no sign change in interval, we don't know if there is a root
  while (!math::signbit(f_l * f_u)) {
    // No interval could be found to bracket the root
    // Might be correct, if there is not root
    if (!check_bracket(n_tries, f_l, f_u, lower, upper)) {
      return false;
    }
    scalar_t d{k * (upper - lower)};
    // Make interval larger in the direction where the function is smaller
    if (math::fabs(f_l) < math::fabs(f_u)) {
      lower -= d;
      f_l = f(lower);
    } else {
      upper += d;
      f_u = f(upper);
    }
    ++n_tries;
  }

  if (!check_bracket(n_tries, f_l, f_u, lower, upper)) {
    return false;
  } else {
    bracket = {lower, upper};
    return true;
  }
}

/// @brief Find a root using the Newton-Raphson algorithm
///
/// @param evaluate_func evaluate the function and its derivative
/// @param s initial guess for the root
/// @param convergence_tolerance max distance between from root before finished
/// @param max_n_tries max number of Newton-Bisection step to try
/// @param max_path don't consider root if it is too far away
///
/// @see Numerical Recepes pp. 445
///
/// @return pathlength to root and the last step size
template <typename scalar_t, typename function_t>
DETRAY_HOST_DEVICE inline std::pair<scalar_t, scalar_t> newton_raphson(
    function_t &evaluate_func, scalar_t s,
    const scalar_t convergence_tolerance = 1.f * unit<scalar_t>::um,
    const std::size_t max_n_tries = 1000u,
    const scalar_t max_path = 5.f * unit<scalar_t>::m) {
  constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
  constexpr scalar_t epsilon{std::numeric_limits<scalar_t>::epsilon()};

  DETRAY_DEBUG_HOST("Initial path estimate: s=" << s << "mm");

  if (math::fabs(s) >= max_path) {
    DETRAY_VERBOSE_HOST("Initial path estimate outside search area: s=" << s);
  }
  if (math::fabs(s) >= inv) {
    const std::string err_msg{"Initial path estimate invalid"};
    DETRAY_FATAL_HOST(err_msg);
    throw std::invalid_argument(err_msg);
  }

  // Run the iteration on s
  scalar_t s_prev{0.f};
  std::size_t n_tries{0u};
  auto [f_s, df_s] = evaluate_func(s);

  while (math::fabs(s - s_prev) > convergence_tolerance) {
    // Root already found?
    if (math::fabs(f_s) < convergence_tolerance) {
      DETRAY_DEBUG_HOST("Found: s = " << s << "mm, err = " << epsilon << "mm");
      return std::make_pair(s, epsilon);
    }

    // No intersection can be found if dividing by zero
    if (math::fabs(df_s) == 0.f) {
      DETRAY_ERROR_HOST(
          "Newton step encountered invalid derivative - skipping");
      return std::make_pair(inv, inv);
    }

    // Newton step
    s_prev = s;
    s -= f_s / df_s;

    // Update function evaluation
    std::tie(f_s, df_s) = evaluate_func(s);

    ++n_tries;

    // No intersection found within max number of trials
    if (n_tries >= max_n_tries) {
      DETRAY_VERBOSE_HOST("Helix intersector did not converge after "
                          << n_tries << " steps - skipping");
      return std::make_pair(inv, inv);
    }
  }
  // Final pathlengt to root and latest step size
  DETRAY_DEBUG_HOST("Found: s = " << s << "mm, err = " << math::fabs(s - s_prev)
                                  << "mm");
  return std::make_pair(s, math::fabs(s - s_prev));
}

/// @brief Find a root using the Newton-Raphson and Bisection algorithms
///
/// @param evaluate_func evaluate the function and its derivative
/// @param s initial guess for the root
/// @param convergence_tolerance max distance between from root before finished
/// @param max_n_tries max number of Newton-Bisection step to try
/// @param max_path don't consider root if it is too far away
///
/// @see Numerical Recepes pp. 445
///
/// @return pathlength to root and the last step size
template <concepts::scalar scalar_t, typename function_t>
DETRAY_HOST_DEVICE inline std::pair<scalar_t, scalar_t> newton_raphson_safe(
    function_t &evaluate_func, scalar_t s,
    const scalar_t convergence_tolerance = 1.f * unit<scalar_t>::um,
    const std::size_t max_n_tries = 1000u,
    const scalar_t max_path = 5.f * unit<scalar_t>::m) {
  constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
  constexpr scalar_t epsilon{std::numeric_limits<scalar_t>::epsilon()};

  DETRAY_DEBUG_HOST("Initial path estimate: s=" << s << "mm");

  // Evaluate the test function at point 'x'
  auto f = [&evaluate_func](const scalar_t x) {
    auto [f_x, df_x] = evaluate_func(x);

    return f_x;
  };

  if (math::fabs(s) >= max_path) {
    DETRAY_VERBOSE_HOST(
        "Initial path estimate outside search area: s=" << s << "mm");
  }
  if (math::fabs(s) >= inv) {
    const std::string err_msg{"Initial path estimate invalid"};
    DETRAY_FATAL_HOST(err_msg);
    throw std::invalid_argument(err_msg);
  }

  // Initial bracket (test a certain range around 's')
  scalar_t a{math::fabs(s) == 0.f ? -0.2f : 0.8f * s};
  scalar_t b{math::fabs(s) == 0.f ? 0.2f : 1.2f * s};
  darray<scalar_t, 2> br{};
  bool is_bracketed = expand_bracket(a, b, f, br);

  // Update initial guess on the root after bracketing
  s = is_bracketed ? 0.5f * (br[1] + br[0]) : s;

  if (!is_bracketed) {
    DETRAY_VERBOSE_HOST("Bracketing failed for initial path estimate: s=" << s);
  } else {
    // Check bracket
    [[maybe_unused]] auto [f_a, df_a] = evaluate_func(br[0]);
    [[maybe_unused]] auto [f_b, df_b] = evaluate_func(br[1]);

    // Bracket is not guaranteed to contain a root
    if (!math::signbit(f_a * f_b)) {
      throw std::runtime_error(
          "Incorrect bracket around root: No sign change!");
    }

    // No bisection algorithm possible if one bracket boundary is inf
    // (is already checked in bracketing alg)
    if ((math::fabs(br[0]) >= inv) || (math::fabs(br[1]) >= inv)) {
      throw std::runtime_error(
          "Incorrect bracket around root: Boundary reached inf!");
    }

    // Root is not within the maximal pathlength
    bool bracket_outside_tol{math::fabs(s) > max_path &&
                             math::fabs(br[0]) >= max_path &&
                             math::fabs(br[1]) >= max_path};
    if (bracket_outside_tol) {
      DETRAY_VERBOSE_HOST("Root outside maximum search area (s = "
                          << s << ", a: " << br[0] << ", b: " << br[1]
                          << ") - skipping");
      return std::make_pair(inv, inv);
    }

    // Root already found?
    if (math::fabs(f_a) < convergence_tolerance) {
      DETRAY_DEBUG_HOST("Found: s = " << a << "mm, err = " << epsilon << "mm");
      return std::make_pair(a, epsilon);
    }
    if (math::fabs(f_b) < convergence_tolerance) {
      DETRAY_DEBUG_HOST("Found: s = " << b << "mm, err = " << epsilon << "mm");
      return std::make_pair(b, epsilon);
    }

    // Make 'a' the boundary for the negative function value -> easier to
    // update
    bool is_lower_a{math::signbit(f_a)};
    a = br[is_lower_a ? 0u : 1u];
    b = br[is_lower_a ? 1u : 0u];
  }

  // Run the iteration on s
  scalar_t s_prev{0.f};
  std::size_t n_tries{0u};
  auto [f_s, df_s] = evaluate_func(s);
  if (math::fabs(f_s) < convergence_tolerance) {
    DETRAY_DEBUG_HOST(
        "Found: s = " << s << "mm, err = " << convergence_tolerance << "mm");
    return std::make_pair(s, epsilon);
  }
  if (math::signbit(f_s)) {
    a = s;
  } else {
    b = s;
  }

  while (math::fabs(s - s_prev) > convergence_tolerance) {
    // Does Newton step escape bracket?
    bool bracket_escape{true};
    scalar_t s_newton{0.f};
    if (math::fabs(df_s) != 0.f) {
      s_newton = s - f_s / df_s;
      bracket_escape = math::signbit((s_newton - a) * (b - s_newton));
    }

    // This criterion from Numerical Recipes seems to work, but why?
    /*const bool slow_convergence{math::fabs(2.f * f_s) >
                                math::fabs((s_prev - s) * df_s)};*/

    // Take a bisection step if it converges faster than Newton
    // |f(next_newton_s)| > |f(next_bisection_s)|
    bool slow_convergence{true};
    // The criterion is only well defined if the step lengths are small
    if (const scalar_t ds_bisection{0.5f * (a + b) - s};
        is_bracketed &&
        (math::fabs(ds_bisection) < 10.f * unit<scalar_t>::mm)) {
      slow_convergence =
          (2.f * math::fabs(f_s) > math::fabs(df_s * ds_bisection + f_s));
    }

    s_prev = s;

    // Run bisection if Newton-Raphson would be poor
    if (is_bracketed &&
        (bracket_escape || slow_convergence || math::fabs(df_s) == 0.f)) {
      // Test the function sign in the middle of the interval
      s = 0.5f * (a + b);
    } else {
      // No intersection can be found if dividing by zero
      if (!is_bracketed && math::fabs(df_s) == 0.f) {
        DETRAY_VERBOSE_HOST("Newton step encountered invalid derivative at s="
                            << s << " after " << n_tries
                            << " steps - skipping");

        return std::make_pair(inv, inv);
      }

      s = s_newton;
    }

    // Update function and bracket
    std::tie(f_s, df_s) = evaluate_func(s);
    if (is_bracketed && math::signbit(f_s)) {
      a = s;
    } else {
      b = s;
    }

    // Converges to a point outside the search space - early stop
    if (math::fabs(s) > max_path && math::fabs(s_prev) > max_path &&
        ((a < -max_path && b < -max_path) || (a > max_path && b > max_path))) {
      DETRAY_VERBOSE_HOST("WARNING: Root finding left the search space at (s = "
                          << s << ", a: " << a << ", b: " << b << ") after "
                          << n_tries << " steps - skipping");

      return std::make_pair(inv, inv);
    }

    ++n_tries;

    // No intersection found within max number of trials
    if (n_tries >= max_n_tries) {
      // Should have found the root
      if (is_bracketed) {
        std::stringstream err_str{};
        err_str << "Helix intersector did not find root for s=" << s << " in ["
                << a << ", " << b << "]";

        DETRAY_FATAL_HOST(err_str.str());
        throw std::runtime_error(err_str.str());
      } else {
        DETRAY_VERBOSE_HOST("Helix intersector did not converge after "
                            << n_tries
                            << " steps unbracketed search - skipping");
      }
      return std::make_pair(inv, inv);
    }
  }
  // Final pathlengt to root and latest step size
  DETRAY_DEBUG_HOST("Found: s = " << s << "mm, err = " << math::fabs(s - s_prev)
                                  << "mm");
  return std::make_pair(s, math::fabs(s - s_prev));
}

/// @brief Fill an intersection with the result of the root finding
///
/// @param [out] sfi the surface intersection
/// @param [in] traj the test trajectory that intersects the surface
/// @param [in] s path length to the root
/// @param [in] ds approximation error for the root
/// @param [in] mask the mask of the surface
/// @param [in] trf the transform of the surface
/// @param [in] mask_tolerance minimal and maximal mask tolerance
template <typename intersection_t, concepts::algebra algebra_t,
          typename surface_descr_t, typename mask_t, typename trajectory_t,
          concepts::transform3D transform3_t, concepts::scalar scalar_t>
DETRAY_HOST_DEVICE constexpr void resolve_mask(
    intersection_t &is, const trajectory_t &traj,
    const intersection_point_err<algebra_t> &ip, const surface_descr_t sf_desc,
    const mask_t &mask, const transform3_t &trf,
    const intersection::config &intr_cfg,
    const scalar_t /*external_mask_tol*/ = 0.f) {
  assert((intr_cfg.min_mask_tolerance == intr_cfg.max_mask_tolerance) &&
         "Helix intersectors use only one mask tolerance value");

  // Build intersection struct from test trajectory, if the distance is valid
  if (!detail::is_invalid_value(ip.path)) {
    is.set_path(ip.path);
    is.set_local(
        mask_t::to_local_frame3D(trf, traj.pos(ip.path), traj.dir(ip.path)));

    const scalar_t cos_incidence_angle = vector::dot(
        mask_t::get_local_frame().normal(trf, is.local()), traj.dir(ip.path));

    scalar_t tol{sf_desc.is_portal() ? 0.f : intr_cfg.min_mask_tolerance};
    // If tolerance is inf, use tolerance estimation (intr_cfg is 'float'!)
    if (tol >= detail::invalid_value<float>()) {
      // Due to floating point errors this can be negative if cos ~ 1
      const scalar_t sin_inc2{
          math::fabs(1.f - cos_incidence_angle * cos_incidence_angle)};

      tol = math::fabs(ip.path_err * math::sqrt(sin_inc2));
    }
    // Make sure the tol. has been estimated/configured in a sensible way
    assert(!math::signbit(tol));
    assert(tol < 1000.f * unit<scalar_t>::mm);

    is.set_status(mask.is_inside(is.local(), tol)
                      ? intersection::status::e_inside
                      : intersection::status::e_outside);
    is.set_surface(sf_desc);
    is.set_direction(!math::signbit(ip.path));
    is.set_volume_link(mask.volume_link());
  } else {
    // Not a valid intersection
    is.set_status(intersection::status::e_outside);
  }

  std::string intr_status{"unknown"};
  if (is.is_inside()) {
    intr_status = "inside";
  } else if (is.is_edge()) {
    intr_status = "edge";
  } else {
    intr_status = "outside";
  }

  DETRAY_DEBUG_HOST("Intersection status: " << intr_status);
}

}  // namespace detray
