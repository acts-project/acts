// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/coordinates/cartesian3D.hpp"
#include "detray/geometry/coordinates/polar2D.hpp"
#include "detray/geometry/detail/shape_utils.hpp"
#include "detray/geometry/detail/vertexer.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <string_view>

namespace detray {

/// @brief Geometrical shape of a stereo annulus that is used for the ITk
/// strip endcaps.
///
/// The stereo annulus is defined in two different(!) polar coordinate systems
/// that differ by an origin shift. The boundaries are the inner and outer
/// radius (bounds[0] and bounds[1]) in the polar coordinate system of an
/// endcap disc (called beam system in the following) that is centered on
/// the beam axis, as well as the two phi boundaries that are defined in the
/// system that is shifted into the focal point of the strips on a given sensor
/// (called focal system in the following). Note, that the local coordinate
/// system of the annulus surface is the same as the shifted disc (focal)
/// system! The mask phi boundaries (bounds[2] and bounds[3]) are defined
/// relative to the average phi position ( bounds[6]) of the strips in the
/// focal system.
/// Due to the focal polar coordinate system of the strips needing a different
/// origin from the beam polar system, two additional conversion parameters are
/// included (bounds[4], bounds[5]). These are the origin shift in x and y
/// respectively.
class annulus2D {
 public:
  /// The name for this shape
  static constexpr std::string_view name = "(stereo) annulus2D";

  /// Names for the mask boundary values
  enum boundaries : unsigned int {
    e_min_r = 0u,
    e_max_r = 1u,
    e_min_phi_rel = 2u,
    e_max_phi_rel = 3u,
    e_average_phi = 4u,
    e_shift_x = 5u,
    e_shift_y = 6u,
    e_size = 7u,
  };

  /// Container definition for the shape boundary values
  template <concepts::scalar scalar_t>
  using bounds_type = darray<scalar_t, boundaries::e_size>;

  /// Local coordinate frame ( focal system )
  template <concepts::algebra algebra_t>
  using local_frame_type = polar2D<algebra_t>;

  /// Result type of a boundary check
  template <typename bool_t>
  using result_type = detail::boundary_check_result<bool_t>;

  /// Dimension of the local coordinate system
  static constexpr std::size_t dim{2u};

  /// @returns the stereo angle calculated from the mask @param bounds .
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE darray<scalar_t, 8> stereo_angle(
      const bounds_type<scalar_t> &bounds) const {
    // Half stereo angle (phi_s / 2) (y points in the long strip direction)
    return 2.f * math::atan(bounds[e_shift_y] / bounds[e_shift_x]);
  }

  /// @returns The phi position in relative to the average phi of the annulus.
  template <concepts::scalar scalar_t, concepts::point point_t>
  DETRAY_HOST_DEVICE inline scalar_t get_phi_rel(
      const bounds_type<scalar_t> &bounds, const point_t &loc_p) const {
    // Rotate by avr phi in the focal system (this is usually zero)
    return loc_p[1] - bounds[e_average_phi];
  }

  /// @returns The squared radial position in the beam frame.
  template <concepts::scalar scalar_t, concepts::point point_t>
  DETRAY_HOST_DEVICE inline scalar_t get_r2_beam_frame(
      const bounds_type<scalar_t> &bounds, const point_t &loc_p) const {
    // Go to beam frame to check r boundaries. Use the origin
    // shift in polar coordinates for that
    // TODO: Put shift in r-phi into the bounds?
    point_t shift_xy;
    shift_xy[0u] = -bounds[e_shift_x];
    shift_xy[1u] = -bounds[e_shift_y];
    const scalar_t shift_r = vector::perp(shift_xy);
    const scalar_t shift_phi = vector::phi(shift_xy);

    return shift_r * shift_r + loc_p[0] * loc_p[0] +
           2.f * shift_r * loc_p[0] *
               math::cos(get_phi_rel(bounds, loc_p) - shift_phi);
  }

  /// @brief Find the minimum distance to any boundary.
  ///
  /// @note the point is expected to be given in local coordinates by the
  /// caller. For the annulus shape, the local coordinate system of the
  /// strips is used (focal system).
  ///
  /// @param bounds the boundary values for this shape
  /// @param loc_p the point to be checked in the local coordinate system
  ///
  /// @return the minimum distance.
  template <concepts::scalar scalar_t, concepts::point point_t>
  DETRAY_HOST_DEVICE inline scalar_t min_dist_to_boundary(
      const bounds_type<scalar_t> &bounds, const point_t &loc_p) const {
    // The two quantities to check: r^2 in beam system, phi in focal system:

    // Rotate by avr phi in the focal system (this is usually zero)
    const scalar_t phi_rel_focal = get_phi_rel(bounds, loc_p);

    // Check phi boundaries, which are well def. in focal frame
    const scalar_t min_phi_dist =
        math::min(math::fabs(phi_rel_focal - bounds[e_min_phi_rel]),
                  math::fabs(bounds[e_max_phi_rel] - phi_rel_focal));

    const auto r_beam = math::sqrt(get_r2_beam_frame(bounds, loc_p));

    const scalar_t min_r_dist = math::min(math::fabs(r_beam - bounds[e_min_r]),
                                          math::fabs(bounds[e_max_r] - r_beam));

    // Compare the radius with the chord
    return math::min(min_r_dist,
                     2.f * loc_p[0] * math::sin(0.5f * min_phi_dist));
  }

  /// @brief Check boundary values for a local point.
  /// @{
  /// @param bounds the boundary values for this shape
  /// @param trf the surface transform
  /// @param glob_p the point to be checked in the global coordinate system
  /// @param tol dynamic tolerance determined by caller
  ///
  /// @return true if the local point lies within the given boundaries.
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE constexpr result_type<dbool<algebra_t>> check_boundaries(
      const bounds_type<dscalar<algebra_t>> &bounds,
      const dtransform3D<algebra_t> &trf, const dpoint3D<algebra_t> &glob_p,
      const dscalar<algebra_t> tol =
          std::numeric_limits<dscalar<algebra_t>>::epsilon(),
      const dscalar<algebra_t> edge_tol = 0.f) const {
    using scalar_t = dscalar<algebra_t>;

    // Move point to local plane: Focal frame in cartesian coordinates
    const auto loc_3D{cartesian3D<algebra_t>::global_to_local(trf, glob_p, {})};

    // Shift local 3D position into beam frame to check the radius
    const scalar_t new_x{loc_3D[0] + bounds[e_shift_x]};
    const scalar_t new_y{loc_3D[1] + bounds[e_shift_y]};

    const scalar_t r_beam{math::sqrt(math::fma(new_x, new_x, new_y * new_y))};

    auto inside_mask = ((r_beam + tol) >= bounds[e_min_r]) &&
                       (r_beam <= (bounds[e_max_r] + tol));

    // Try to avoid the costly phi calculation
    auto phi_focal{detail::invalid_value<scalar_t>()};
    if (detail::any_of(inside_mask)) {
      // Get phi for phi-bounds check and rotate by average phi
      phi_focal = vector::phi(loc_3D) - bounds[e_average_phi];
      // Estimate angular tolerance along r
      const scalar_t phi_tol{detail::phi_tolerance(tol, r_beam)};

      inside_mask = (phi_focal >= (bounds[e_min_phi_rel] - phi_tol)) &&
                    (phi_focal <= (bounds[e_max_phi_rel] + phi_tol)) &&
                    inside_mask;
    }

    decltype(inside_mask) inside_edge{false};
    if (detail::any_of(edge_tol > 0.f)) {
      // Edge tolerance
      const scalar_t full_tol{tol + edge_tol};

      inside_edge = ((r_beam + full_tol) >= bounds[e_min_r]) &&
                    (r_beam <= (bounds[e_max_r] + full_tol));

      if (detail::any_of(inside_edge)) {
        // If phi had not been calculated before, do it now
        if (detail::is_invalid_value(phi_focal)) {
          phi_focal = vector::phi(loc_3D) - bounds[e_average_phi];
        }

        const scalar_t phi_tol_full{detail::phi_tolerance(full_tol, r_beam)};

        inside_edge = (phi_focal >= (bounds[e_min_phi_rel] - phi_tol_full)) &&
                      (phi_focal <= (bounds[e_max_phi_rel] + phi_tol_full)) &&
                      inside_edge;
      }
    }

    return result_type<decltype(inside_mask)>{inside_mask, inside_edge};
  }

  /// @note the point is expected to be given in local coordinates by the
  /// caller. For the annulus shape, the local coordinate system of the
  /// strips is used (focal system).
  ///
  /// @param bounds the boundary values for this shape
  /// @param loc_p the point to be checked in the local coordinate system
  /// @param tol dynamic tolerance determined by caller
  ///
  /// @return true if the local point lies within the given boundaries.
  template <concepts::scalar scalar_t, concepts::point point_t>
  DETRAY_HOST_DEVICE constexpr auto check_boundaries(
      const bounds_type<scalar_t> &bounds, const point_t &loc_p,
      const scalar_t tol = std::numeric_limits<scalar_t>::epsilon(),
      const scalar_t edge_tol = 0.f) const {
    // The two quantities to check: r^2 in beam system, phi in focal system:

    // Rotate by avr phi in the focal system (this is usually zero)
    const scalar_t phi_focal = get_phi_rel(bounds, loc_p);

    // Check phi boundaries, which are well def. in focal frame
    const scalar_t phi_tol = detail::phi_tolerance(tol, loc_p[0]);
    auto inside_mask = !((phi_focal < (bounds[e_min_phi_rel] - phi_tol)) ||
                         (phi_focal > (bounds[e_max_phi_rel] + phi_tol)));

    // Try to avoid the costly r_beam calculation
    auto r_beam2{detail::invalid_value<scalar_t>()};
    if (detail::any_of(inside_mask)) {
      r_beam2 = get_r2_beam_frame(bounds, loc_p);

      // Apply tolerances as squares: 0 <= a, 0 <= b: a^2 <= b^2 <=> a <=
      // b
      const scalar_t minR_tol = math::max(bounds[e_min_r] - tol, scalar_t(0.f));
      const scalar_t maxR_tol = bounds[e_max_r] + tol;

      assert(detail::all_of(minR_tol >= scalar_t(0.f)));

      inside_mask = (r_beam2 >= (minR_tol * minR_tol)) &&
                    (r_beam2 <= (maxR_tol * maxR_tol)) && inside_mask;
    }

    decltype(inside_mask) inside_edge{false};
    if (detail::any_of(edge_tol > 0.f)) {
      // Edge tolerance
      const scalar_t full_tol{tol + edge_tol};
      const scalar_t phi_tol_full = detail::phi_tolerance(full_tol, loc_p[0]);

      const auto phi_check_edge =
          (phi_focal >= (bounds[e_min_phi_rel] - phi_tol_full)) &&
          (phi_focal <= (bounds[e_max_phi_rel] + phi_tol_full));

      if (detail::any_of(inside_edge)) {
        // If phi had not been calculated before, do it now
        if (detail::is_invalid_value(r_beam2)) {
          r_beam2 = get_r2_beam_frame(bounds, loc_p);
        }

        const scalar_t minR_tol_edge =
            math::max(bounds[e_min_r] - full_tol, scalar_t(0.f));
        const scalar_t maxR_tol_edge = bounds[e_max_r] + full_tol;

        assert(detail::all_of(minR_tol_edge >= scalar_t(0.f)));

        inside_edge = (r_beam2 >= (minR_tol_edge * minR_tol_edge)) &&
                      (r_beam2 <= (maxR_tol_edge * maxR_tol_edge)) &&
                      phi_check_edge;
      }
    }

    return result_type<decltype(inside_mask)>{inside_mask, inside_edge};
  }
  /// @}

  /// @brief Measure of the shape: Area
  ///
  /// @note (not yet implemented!)
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the stereo annulus area on the plane.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t measure(
      const bounds_type<scalar_t> &bounds) const {
    return area(bounds);
  }

  /// @brief The area of a the shape
  ///
  /// @note (not yet implemented!)
  ///
  /// @param bounds the boundary values for this shape
  ///
  /// @returns the stereo annulus area.
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr scalar_t area(
      const bounds_type<scalar_t> &) const {
    return detail::invalid_value<scalar_t>();
  }

  /// @brief Merge two annulus shapes
  ///
  /// @param bounds the boundary values for this shape
  /// @param o_bounds the boundary values for the other shape
  ///
  /// @returns merged bound values
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE constexpr bounds_type<scalar_t> merge(
      const bounds_type<scalar_t> &bounds,
      const bounds_type<scalar_t> &o_bounds) const {
    assert(bounds[e_average_phi] == o_bounds[e_average_phi]);
    assert(bounds[e_shift_x] == o_bounds[e_shift_x]);
    assert(bounds[e_shift_y] == o_bounds[e_shift_y]);

    bounds_type<scalar_t> new_bounds{};

    new_bounds[e_min_r] = math::min(bounds[e_min_r], o_bounds[e_min_r]);
    new_bounds[e_max_r] = math::max(bounds[e_max_r], o_bounds[e_max_r]);
    new_bounds[e_min_phi_rel] =
        math::min(bounds[e_min_phi_rel], o_bounds[e_min_phi_rel]);
    new_bounds[e_max_phi_rel] =
        math::max(bounds[e_max_phi_rel], o_bounds[e_max_phi_rel]);
    new_bounds[e_average_phi] = bounds[e_average_phi];
    new_bounds[e_shift_x] = bounds[e_shift_x];
    new_bounds[e_shift_y] = bounds[e_shift_y];

    return new_bounds;
  }

  /// @brief Lower and upper point for minimal axis aligned bounding box.
  ///
  /// Computes the min and max vertices in a local cartesian frame.
  ///
  /// @param bounds the boundary values for this shape
  /// @param env dynamic envelope around the shape
  ///
  /// @returns and array of coordinates that contains the lower point (first
  /// three values) and the upper point (latter three values).
  // @TODO: this is a terrible approximation: restrict to annulus corners
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE darray<dscalar<algebra_t>, 6> local_min_bounds(
      const bounds_type<dscalar<algebra_t>> &bounds,
      const dscalar<algebra_t> env =
          std::numeric_limits<dscalar<algebra_t>>::epsilon()) const {
    using scalar_t = dscalar<algebra_t>;
    using point_t = dpoint2D<algebra_t>;

    assert(env > 0.f);

    const auto c_pos = corners(bounds);

    const scalar_t o_x{bounds[e_shift_x]};
    const scalar_t o_y{bounds[e_shift_y]};

    // Corner points 'b' and 'c' in local cartesian beam system
    const point_t b{c_pos[4] * math::cos(c_pos[5]) - o_x,
                    c_pos[4] * math::sin(c_pos[5]) - o_y};
    const point_t c{c_pos[6] * math::cos(c_pos[7]) - o_x,
                    c_pos[6] * math::sin(c_pos[7]) - o_y};

    // bisector = 0.5 * (c + b). Scale to the length of the full circle to
    // get the outermost point
    const point_t t = bounds[e_max_r] * vector::normalize(c + b);

    // Find the min/max positions in x and y
    darray<scalar_t, 5> x_pos{c_pos[2] * math::cos(c_pos[3]) - o_x, b[0], c[0],
                              c_pos[0] * math::cos(c_pos[1]) - o_x, t[0]};
    darray<scalar_t, 5> y_pos{c_pos[2] * math::sin(c_pos[3]) - o_y, b[1], c[1],
                              c_pos[0] * math::sin(c_pos[1]) - o_y, t[1]};

    constexpr scalar_t inv{detail::invalid_value<scalar_t>()};
    scalar_t min_x{inv};
    scalar_t min_y{inv};
    scalar_t max_x{-inv};
    scalar_t max_y{-inv};
    for (unsigned int i{0u}; i < 5u; ++i) {
      min_x = x_pos[i] < min_x ? x_pos[i] : min_x;
      max_x = x_pos[i] > max_x ? x_pos[i] : max_x;
      min_y = y_pos[i] < min_y ? y_pos[i] : min_y;
      max_y = y_pos[i] > max_y ? y_pos[i] : max_y;
    }

    return {min_x - env, min_y - env, -env, max_x + env, max_y + env, env};
  }

  /// @brief Stereo annulus corners in polar strip system.
  ///
  /// @param bounds the boundary values for the stereo annulus
  ///
  /// @note see calculation of strip lengths in
  ///       https://cds.cern.ch/record/1514636?ln=en p10-13
  ///
  /// @returns an array of coordinates that contains the lower point (first
  /// four values) and the upper point (latter four values).
  template <concepts::scalar scalar_t>
  DETRAY_HOST_DEVICE darray<scalar_t, 8> corners(
      const bounds_type<scalar_t> &bounds) const {
    // Calculate the r-coordinate of a point in the strip system from the
    // circle arc radius (e.g. min_r) and the phi position in the strip
    // system (e.g. for the corners these are average_phi + min_phi_rel and
    // average_phi + max_phi_rel).
    auto get_strips_pc_r = [&bounds](const scalar_t R,
                                     const scalar_t phi) -> scalar_t {
      // f: Shift distance between beamline and focal system origin
      const scalar_t f2{bounds[e_shift_x] * bounds[e_shift_x] +
                        bounds[e_shift_y] * bounds[e_shift_y]};

      // f * sin(phi_s / 2 + phi) using: f_y / f = sin(phi_s / 2) and
      // sin(a + b) = sin(a)*cos(b) + cos(a)*sin(b)
      const scalar_t f_sin_phi{bounds[e_shift_x] * math::cos(phi) +
                               bounds[e_shift_y] * math::sin(phi)};

      return f_sin_phi + math::sqrt(f_sin_phi * f_sin_phi - f2 + R * R);
    };

    // Calculate the polar coordinates for the corners
    const scalar_t min_phi{bounds[e_average_phi] + bounds[e_min_phi_rel]};
    const scalar_t max_phi{bounds[e_average_phi] + bounds[e_max_phi_rel]};
    darray<scalar_t, 8> corner_pos;
    // bottom left: min_r, min_phi_rel
    corner_pos[0] = get_strips_pc_r(bounds[e_min_r], min_phi);
    corner_pos[1] = min_phi;
    // bottom right: min_r, max_phi_rel
    corner_pos[2] = get_strips_pc_r(bounds[e_min_r], max_phi);
    corner_pos[3] = max_phi;
    // top right: max_r, max_phi_rel
    corner_pos[4] = get_strips_pc_r(bounds[e_max_r], max_phi);
    corner_pos[5] = max_phi;
    // top left: max_r, min_phi_rel
    corner_pos[6] = get_strips_pc_r(bounds[e_max_r], min_phi);
    corner_pos[7] = min_phi;

    return corner_pos;
  }

  /// @returns the shapes centroid in local cartesian coordinates
  /// @note the calculated centroid position is only an approximation
  /// (centroid of the four corner points)!
  template <concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE dpoint3D<algebra_t> centroid(
      const bounds_type<dscalar<algebra_t>> &bounds) const {
    using scalar_t = dscalar<algebra_t>;

    // Strip polar system
    const auto crns = corners(bounds);

    // Coordinates of the centroid position in strip system
    const scalar_t r{0.25f * (crns[0] + crns[2] + crns[4] + crns[6])};
    const scalar_t phi{bounds[e_average_phi]};

    return r * dpoint3D<algebra_t>{math::cos(phi), math::sin(phi), scalar_t(0)};
  }

  /// Generate vertices in local cartesian frame
  ///
  /// @param bounds the boundary values for the stereo annulus
  /// @param n_seg is the number of line segments
  ///
  /// @return a generated list of vertices
  template <concepts::algebra algebra_t>
  DETRAY_HOST dvector<dpoint3D<algebra_t>> vertices(
      const bounds_type<dscalar<algebra_t>> &bounds, dindex n_seg) const {
    using scalar_t = dscalar<algebra_t>;
    using point2_t = dpoint2D<algebra_t>;
    using point3_t = dpoint3D<algebra_t>;

    scalar_t min_r = bounds[e_min_r];
    scalar_t max_r = bounds[e_max_r];
    scalar_t min_phi = bounds[e_average_phi] + bounds[e_min_phi_rel];
    scalar_t max_phi = bounds[e_average_phi] + bounds[e_max_phi_rel];
    scalar_t origin_x = bounds[e_shift_x];
    scalar_t origin_y = bounds[e_shift_y];

    point2_t origin_m = {origin_x, origin_y};

    /// Helper method: find inner outer radius at edges in STRIP PC
    auto circIx = [](scalar_t O_x, scalar_t O_y, scalar_t r,
                     scalar_t phi) -> point2_t {
      //                      ____________________________________________
      //                     /      2  2                    2    2  2    2
      //     O_x + O_y*m - \/  - O_x *m  + 2*O_x*O_y*m - O_y  + m *r  + r
      // x =
      // --------------------------------------------------------------
      //                                  2
      //                                 m  + 1
      //
      // y = m*x
      //
      scalar_t m = math::tan(phi);
      point2_t dir = {math::cos(phi), math::sin(phi)};
      scalar_t x1 = (O_x + O_y * m -
                     math::sqrt(-math::pow(O_x, 2.f) * math::pow(m, 2.f) +
                                2.f * O_x * O_y * m - math::pow(O_y, 2.f) +
                                math::pow(m, 2.f) * math::pow(r, 2.f) +
                                math::pow(r, 2.f))) /
                    (math::pow(m, 2.f) + 1.f);
      scalar_t x2 = (O_x + O_y * m +
                     math::sqrt(-math::pow(O_x, 2.f) * math::pow(m, 2.f) +
                                2.f * O_x * O_y * m - math::pow(O_y, 2.f) +
                                math::pow(m, 2.f) * math::pow(r, 2.f) +
                                math::pow(r, 2.f))) /
                    (math::pow(m, 2.f) + 1.f);

      if (point2_t v1 = {x1, m * x1}; vector::dot(v1, dir) > 0.f) {
        return v1;
      }
      return {x2, m * x2};
    };

    // calculate corners in STRIP XY
    point2_t ul_xy = circIx(origin_x, origin_y, max_r, max_phi);
    point2_t ll_xy = circIx(origin_x, origin_y, min_r, max_phi);
    point2_t ur_xy = circIx(origin_x, origin_y, max_r, min_phi);
    point2_t lr_xy = circIx(origin_x, origin_y, min_r, min_phi);

    auto inner_phi = detail::phi_values(vector::phi(lr_xy - origin_m),
                                        vector::phi(ll_xy - origin_m), n_seg);
    auto outer_phi = detail::phi_values(vector::phi(ul_xy - origin_m),
                                        vector::phi(ur_xy - origin_m), n_seg);

    dvector<point3_t> annulus_vertices;
    annulus_vertices.reserve(inner_phi.size() + outer_phi.size());
    for (scalar_t iphi : inner_phi) {
      annulus_vertices.push_back(point3_t{min_r * math::cos(iphi) + origin_x,
                                          min_r * math::sin(iphi) + origin_y,
                                          0.f});
    }

    for (scalar_t ophi : outer_phi) {
      annulus_vertices.push_back(point3_t{max_r * math::cos(ophi) + origin_x,
                                          max_r * math::sin(ophi) + origin_y,
                                          0.f});
    }

    return annulus_vertices;
  }

  /// @brief Check consistency of boundary values.
  ///
  /// @param bounds the boundary values for this shape
  /// @param os output stream for error messages
  ///
  /// @return true if the bounds are consistent.
  template <concepts::scalar scalar_t>
  DETRAY_HOST constexpr bool check_consistency(
      const bounds_type<scalar_t> &bounds, std::ostream &os) const {
    constexpr auto tol{10.f * std::numeric_limits<scalar_t>::epsilon()};

    if (math::signbit(bounds[e_min_r]) || bounds[e_max_r] < tol) {
      os << "DETRAY ERROR (HOST): Radial bounds must be in the range [0, "
            "numeric_max)";
      return false;
    }
    if (bounds[e_min_r] >= bounds[e_max_r] ||
        math::fabs(bounds[e_min_r] - bounds[e_max_r]) < tol) {
      os << "DETRAY ERROR (HOST): Min radius must be smaller than max "
            "radius.";
      return false;
    }
    if ((bounds[e_min_phi_rel] < -constant<scalar_t>::pi ||
         bounds[e_min_phi_rel] > constant<scalar_t>::pi) ||
        (bounds[e_max_phi_rel] < -constant<scalar_t>::pi ||
         bounds[e_max_phi_rel] > constant<scalar_t>::pi) ||
        (bounds[e_average_phi] < -constant<scalar_t>::pi ||
         bounds[e_average_phi] > constant<scalar_t>::pi)) {
      os << "DETRAY ERROR (HOST): Angles must map onto [-pi, pi] range.";
      return false;
    }
    if (bounds[e_min_phi_rel] >= bounds[e_max_phi_rel] ||
        math::fabs(bounds[e_min_phi_rel] - bounds[e_max_phi_rel]) < tol) {
      os << "DETRAY ERROR (HOST): Min relative angle must be smaller "
            "than max relative "
            "angle.";
      return false;
    }

    return true;
  }
};

}  // namespace detray
