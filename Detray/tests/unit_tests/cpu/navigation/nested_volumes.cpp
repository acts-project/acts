// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Project include(s)
#include "detray/definitions/units.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/navigation/caching_navigator.hpp"
#include "detray/navigation/navigation_state.hpp"
#include "detray/propagator/actor_chain.hpp"
#include "detray/propagator/base_actor.hpp"
#include "detray/propagator/line_stepper.hpp"
#include "detray/propagator/propagator.hpp"
#include "detray/tracks/tracks.hpp"

// Detray test include(s)
#include "detray/test/common/build_nested_detector.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/utils/inspectors.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>

// GTest include(s)
#include <gtest/gtest.h>

// System include(s)
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace detray;

namespace {

using algebra_t = test::algebra;
using scalar = test::scalar;
using point3 = test::point3;
using vector3 = test::vector3;

/// A surface that was reached during the navigation
struct surface_hit {
  /// Index of the surface in the detector
  dindex sf_idx{dindex_invalid};
  /// Volume that the navigation was in when the surface was reached
  dindex nav_volume{dindex_invalid};
  /// Whether the surface is a sensitive surface
  bool is_sensitive{false};
};

/// Navigation inspector that records the surfaces that are reached, together
/// with the current volume of the navigation
struct surface_hit_tracer {
  /// No device view needed for this host-only inspector
  struct void_view : public detray::detail::dbase_view {};
  using view_type = void_view;
  using const_view_type = const void_view;

  /// The recorded surface hits
  std::vector<surface_hit> hits{};

  /// Record a surface, if the navigation is on a surface
  template <typename state_t, typename point3_t, typename vector3_t,
            typename... Args>
  void operator()(const state_t &state, const navigation::config & /*cfg*/,
                  const point3_t & /*pos*/, const vector3_t & /*dir*/,
                  const char * /*message*/, Args &&.../*args*/) {
    if (!(state.status() == navigation::status::e_on_object ||
          state.status() == navigation::status::e_on_portal)) {
      return;
    }
    const auto &sf = state.current().surface();
    const surface_hit hit{sf.index(), static_cast<dindex>(state.volume()),
                          sf.is_sensitive()};

    // Do not record the same surface hit twice
    if (!hits.empty() && hits.back().sf_idx == hit.sf_idx &&
        hits.back().nav_volume == hit.nav_volume) {
      return;
    }
    hits.push_back(hit);
  }

  /// Ignore status messages
  template <typename state_t>
  void operator()(const state_t & /*state*/, const char * /*message*/) const {
    /*Do nothing*/
  }
};

/// Abort the propagation after a maximal number of steps (catches endless
/// loops between volumes, which do not advance the path length)
struct step_limit_aborter : public base_actor {
  /// Actor state: step counter
  struct state {
    /// Number of steps so far
    unsigned int n_steps{0u};
    /// Maximal number of steps
    unsigned int max_steps{1000u};
  };

  /// Count the steps and abort, if the limit is reached
  template <typename propagator_state_t>
  void operator()(state &abrt_state, propagator_state_t &prop_state) const {
    if (prop_state.navigation().finished()) {
      return;
    }
    if (++abrt_state.n_steps > abrt_state.max_steps) {
      prop_state.navigation().abort("Step limit reached");
      prop_state.heartbeat(false);
    }
  }
};

/// Truth: The sensitive surfaces that a ray crosses, sorted by distance
///
/// All sensitive surfaces of the nested detector are rectangles in the x-y
/// plane without rotation.
///
/// @param det the detector
/// @param ori the ray origin
/// @param dir the ray direction
/// @param edge_tol distance to a rectangle edge that is considered ambiguous
/// @param is_ambiguous set to true, if the ray passes close to an edge
///
/// @returns the sorted sensitive surface indices
template <typename detector_t>
std::vector<dindex> truth_sensitives(const detector_t &det, const point3 &ori,
                                     const vector3 &dir, const scalar edge_tol,
                                     bool &is_ambiguous) {
  typename detector_t::geometry_context ctx{};
  std::vector<std::pair<scalar, dindex>> crossings{};

  for (const auto &sf_desc : det.surfaces()) {
    if (!sf_desc.is_sensitive() || std::abs(dir[2]) < 1e-6f) {
      continue;
    }
    const auto sf = geometry::surface{det, sf_desc};
    const point3 t = sf.transform(ctx).translation();

    const scalar path{(t[2] - ori[2]) / dir[2]};
    if (path <= 0.f) {
      continue;
    }
    const point3 glob = ori + path * dir;
    const scalar loc_x{glob[0] - t[0]};
    const scalar loc_y{glob[1] - t[1]};

    // One mask per sensitive surface
    const scalar min_x{sf.boundary(rectangle2D::e_min_x).front()};
    const scalar min_y{sf.boundary(rectangle2D::e_min_y).front()};
    const scalar max_x{sf.boundary(rectangle2D::e_max_x).front()};
    const scalar max_y{sf.boundary(rectangle2D::e_max_y).front()};

    // Distance to the closest rectangle edge (negative: outside)
    const scalar dist{
        std::min({loc_x - min_x, max_x - loc_x, loc_y - min_y, max_y - loc_y})};
    if (std::abs(dist) < edge_tol) {
      is_ambiguous = true;
    }
    if (dist > 0.f) {
      crossings.emplace_back(path, sf_desc.index());
    }
  }

  std::ranges::sort(crossings);

  std::vector<dindex> sf_indices{};
  for (const auto &c : crossings) {
    sf_indices.push_back(c.second);
  }
  return sf_indices;
}

/// Propagate straight tracks through the nested detector and compare the
/// sensitive surfaces with the truth
///
/// @param cfg the nested detector configuration
/// @param origin the origin of the tracks
/// @param start_volume the volume that contains the track origin
/// @param nav_dir the navigation direction
void run_nested_navigation_test(const nested_detector_config &cfg,
                                const point3 &origin, const dindex start_volume,
                                const navigation::direction nav_dir) {
  vecmem::host_memory_resource host_mr;
  const auto [det, names] = build_nested_detector<algebra_t>(host_mr, cfg);

  using detector_t = std::remove_cvref_t<decltype(det)>;
  using inspector_t =
      aggregate_inspector<surface_hit_tracer, navigation::print_inspector>;
  using navigator_t =
      caching_navigator<detector_t, navigation::default_cache_size,
                        inspector_t>;
  using stepper_t = line_stepper<algebra_t>;
  using propagator_t =
      propagator<stepper_t, navigator_t, actor_chain<step_limit_aborter>>;

  propagation::config prop_cfg{};
  propagator_t prop{prop_cfg};

  // Ignore tracks that pass close to a rectangle edge
  constexpr scalar edge_tol{1e-2f * unit<scalar>::mm};

  constexpr unsigned int n_theta{60u};
  constexpr unsigned int n_phi{60u};

  unsigned int n_tracks{0u};
  unsigned int n_failed{0u};
  unsigned int n_sensitives{0u};
  for (unsigned int i_theta = 0u; i_theta < n_theta; ++i_theta) {
    // Avoid the poles and directions parallel to the sensitive surfaces
    const scalar theta{0.02f + (constant<scalar>::pi - 0.04f) *
                                   static_cast<scalar>(i_theta) /
                                   static_cast<scalar>(n_theta - 1u)};
    for (unsigned int i_phi = 0u; i_phi < n_phi; ++i_phi) {
      // Avoid symmetric directions
      const scalar phi{0.013f + 2.f * constant<scalar>::pi *
                                    static_cast<scalar>(i_phi) /
                                    static_cast<scalar>(n_phi)};

      const vector3 dir{std::sin(theta) * std::cos(phi),
                        std::sin(theta) * std::sin(phi), std::cos(theta)};

      // The direction in which the track is followed
      const vector3 nav_vec{
          nav_dir == navigation::direction::e_forward ? dir : -1.f * dir};

      bool is_ambiguous{false};
      const auto exp_sensitives{
          truth_sensitives(det, origin, nav_vec, edge_tol, is_ambiguous)};
      if (is_ambiguous) {
        continue;
      }

      const free_track_parameters<algebra_t> track{
          origin, 0.f, 10.f * unit<scalar>::GeV * dir, -1.f};

      typename propagator_t::state state(track, det, prop_cfg.context);
      state.navigation().set_volume(start_volume);
      state.navigation().set_direction(nav_dir);

      ++n_tracks;
      step_limit_aborter::state limit_state{};
      const bool success{prop.propagate(state, detray::tie(limit_state))};

      // Recorded sensitive surfaces and the navigation volume
      const auto &hits = state.navigation()
                             .inspector()
                             .template get<surface_hit_tracer>()
                             .hits;
      std::vector<dindex> rec_sensitives{};
      bool wrong_volume{false};
      std::stringstream hit_stream{};
      for (const auto &hit : hits) {
        hit_stream << "  sf " << hit.sf_idx << " (vol "
                   << det.surface(hit.sf_idx).volume() << ") in nav. volume "
                   << names.at(hit.nav_volume) << "\n";
        if (hit.is_sensitive) {
          rec_sensitives.push_back(hit.sf_idx);
          wrong_volume |= (det.surface(hit.sf_idx).volume() != hit.nav_volume);
        }
      }

      const bool is_ok{success && state.navigation().finished() &&
                       !wrong_volume && rec_sensitives == exp_sensitives};
      n_sensitives += static_cast<unsigned int>(exp_sensitives.size());

      if (!is_ok) {
        ++n_failed;
        // Only print the first failures
        if (n_failed <= 5u) {
          std::stringstream exp_stream{};
          for (const dindex idx : exp_sensitives) {
            exp_stream << idx << " ";
          }
          ADD_FAILURE() << "Track (theta=" << theta << ", phi=" << phi
                        << ") failed: success=" << success
                        << ", finished=" << state.navigation().finished()
                        << ", wrong volume=" << wrong_volume
                        << "\n expected sensitives: " << exp_stream.str()
                        << "\n recorded surfaces:\n"
                        << hit_stream.str();
          // Full navigation log of the first failed track
          if (n_failed == 1u) {
            std::cout << state.navigation()
                             .inspector()
                             .template get<navigation::print_inspector>()
                             .to_string()
                      << std::endl;
          }
        }
      }
    }
  }

  EXPECT_EQ(n_failed, 0u) << "of " << n_tracks << " tracks";
  // Make sure the test is not trivial
  EXPECT_GT(n_tracks, n_theta * n_phi / 2u);
  EXPECT_GT(n_sensitives, n_tracks);
}

}  // anonymous namespace

/// Start in the world volume, above box A
GTEST_TEST(detray_navigation, nested_volumes_from_world) {
  const point3 origin{0.37f, -0.61f, -20.3f};

  for (const float inset : {0.f, 0.5f * unit<float>::um, 1.f}) {
    SCOPED_TRACE("inset = " + std::to_string(inset));
    const auto cfg = nested_detector_config{}.inset(inset);

    run_nested_navigation_test(cfg, origin, 0u,
                               navigation::direction::e_forward);
    run_nested_navigation_test(cfg, origin, 0u,
                               navigation::direction::e_backward);
  }
}

/// Start in the world volume, below box A: Tracks enter box A through the
/// face that box G touches
GTEST_TEST(detray_navigation, nested_volumes_flush_entry) {
  const point3 origin{0.37f, -0.61f, -160.3f};

  for (const float inset : {0.f, 0.5f * unit<float>::um, 1.f}) {
    SCOPED_TRACE("inset = " + std::to_string(inset));
    const auto cfg = nested_detector_config{}.inset(inset);

    run_nested_navigation_test(cfg, origin, 0u,
                               navigation::direction::e_forward);
  }
}

/// Start in the daughter volumes
GTEST_TEST(detray_navigation, nested_volumes_from_daughters) {
  // Box G (second hierarchy level)
  for (const float inset : {0.f, 1.f}) {
    SCOPED_TRACE("box G, inset = " + std::to_string(inset));
    const auto cfg = nested_detector_config{}.inset(inset);

    run_nested_navigation_test(cfg, {0.37f, -0.61f, -127.3f + inset}, 2u,
                               navigation::direction::e_forward);
  }

  // Box D (touches box C)
  SCOPED_TRACE("box D");
  run_nested_navigation_test(nested_detector_config{}, {0.37f, 129.39f, 70.3f},
                             5u, navigation::direction::e_forward);
}
