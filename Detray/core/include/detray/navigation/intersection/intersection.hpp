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
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/definitions/units.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <cstdint>
#include <limits>
#include <ostream>

namespace detray {

namespace intersection {

/// Whether the intersector result will contain the local position
inline constexpr bool contains_pos{true};

/// Intersection status
enum class status : std::uint_least8_t {
  e_inside = 0u,   ///< Inside the mask within numeric uncertainty
  e_edge = 1u,     ///< Inside the mask's external tolerance band
  e_outside = 2u,  ///< Outside of mask (including all tolerances)
};

}  // namespace intersection

/// Result of intersector: Point of intersection on a surface
template <concepts::algebra algebra_t, concepts::point point_t,
          bool has_pos = true>
struct intersection_point {};

/// Result of intersector: Only contains the distance (production)
template <concepts::algebra algebra_t, concepts::point point_t>
struct intersection_point<algebra_t, point_t, !intersection::contains_pos> {
  using scalar_type = dscalar<algebra_t>;

  /// @returns true if debug information needs to be filled
  static consteval bool contains_pos() { return !intersection::contains_pos; }

  /// Distance between track position and surface along test trajectory
  scalar_type path = detail::invalid_value<dvalue<algebra_t>>();

  /// @returns true if the data represents a valid intersection solution
  constexpr bool is_valid() const {
    constexpr auto inv_path{10.f * unit<dvalue<algebra_t>>::m};
    return detray::detail::any_of(math::fabs(path) < inv_path);
  }

  /// Transform to a string for output debugging
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &out_stream,
                                  const intersection_point &ip) {
    out_stream << "dist:" << ip.path;

    return out_stream;
  }
};

/// Result of intersector: Distance and intersection point (debug)
template <concepts::algebra algebra_t, concepts::point point_t>
struct intersection_point<algebra_t, point_t, intersection::contains_pos>
    : public intersection_point<algebra_t, point_t,
                                !intersection::contains_pos> {
 private:
  using base_type =
      intersection_point<algebra_t, point_t, !intersection::contains_pos>;

 public:
  using scalar_type = dscalar<algebra_t>;
  using point_type = point_t;

  constexpr intersection_point() = default;

  DETRAY_HOST_DEVICE
  explicit constexpr intersection_point(const base_type ip) : base_type{ip} {}

  DETRAY_HOST_DEVICE
  constexpr intersection_point(const scalar_type p, const point_type &pnt)
      : base_type{p}, point{pnt} {}

  /// @returns true if debug information needs to be filled
  static consteval bool contains_pos() { return intersection::contains_pos; }

  /// Local position of the intersection on the surface
  point_type point{init_point()};

  /// Transform to a string for output debugging
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &out_stream,
                                  const intersection_point &ip) {
    using base_t =
        intersection_point<algebra_t, point_t, !intersection::contains_pos>;

    out_stream << static_cast<base_t>(ip);

    out_stream << ", point [" << ip.point[0] << ", " << ip.point[1];

    if constexpr (std::same_as<point_t, dpoint3D<algebra_t>>) {
      out_stream << ", " << ip.point[2] << "]";
    } else {
      out_stream << "]";
    }

    return out_stream;
  }

 private:
  /// Initialize points of different dimensionality correctly
  DETRAY_HOST_DEVICE
  constexpr point_t init_point() const {
    constexpr auto inv{detail::invalid_value<dvalue<algebra_t>>()};
    if constexpr (std::same_as<point_t, dpoint2D<algebra_t>>) {
      return {inv, inv};
    } else if constexpr (std::same_as<point_t, dpoint3D<algebra_t>>) {
      return {inv, inv, inv};
    } else {
      assert(false);
      return {};
    }
  }
};

/// Result of intersector: Distance, intersection point (2D or 3D,
/// local/global) and error estimate on distance (for mask resolution)
template <concepts::algebra algebra_t>
struct intersection_point_err
    : public intersection_point<algebra_t, dpoint3D<algebra_t>,
                                intersection::contains_pos> {
 private:
  using base_type = intersection_point<algebra_t, dpoint3D<algebra_t>,
                                       intersection::contains_pos>;

  static constexpr auto inv{detail::invalid_value<dscalar<algebra_t>>()};

 public:
  using scalar_type = dscalar<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
  using point_type = point3_type;

  constexpr intersection_point_err() = default;

  DETRAY_HOST_DEVICE
  explicit constexpr intersection_point_err(const base_type &base_ip)
      : base_type{base_ip} {}

  DETRAY_HOST_DEVICE
  explicit constexpr intersection_point_err(
      const intersection_point<algebra_t, point3_type,
                               !intersection::contains_pos> &ip)
      : base_type{ip.path, point3_type{inv, inv, inv}} {}

  DETRAY_HOST_DEVICE
  explicit constexpr intersection_point_err(
      const intersection_point<algebra_t, point2_type,
                               !intersection::contains_pos> &ip)
      : base_type{ip.path, point3_type{inv, inv, inv}} {}

  DETRAY_HOST_DEVICE
  constexpr intersection_point_err(const scalar_type p, const point3_type &pnt,
                                   const scalar_type p_err)
      : base_type{p, pnt}, path_err{p_err} {}

  /// Error estimation on the path
  scalar_type path_err{inv};
};

/// @brief This class holds the intersection information.
///
/// @tparam surface_t is the type of surface descriptor
/// @tparam algebra_t linear algebra and memory layout
/// @tparam has_pos whether the local position is saved or not
template <typename surface_t, concepts::algebra algebra_t,
          bool has_pos = intersection::contains_pos>
class intersection2D {
  using T = dvalue<algebra_t>;
  using bool_type = dbool<algebra_t>;
  using scalar_type = dscalar<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;

 public:
  using algebra_type = algebra_t;
  using surface_type = surface_t;
  // This is needed for the SIMD implementation, where the boolean type that
  // results from the mask tolerance check has to match the SIMD vector type
  // on which it is used on for a masked assignment of the different status
  // codes
  using status_type =
      std::conditional_t<concepts::soa<algebra_t>, T, intersection::status>;
  using nav_link_type = typename surface_type::navigation_link;

  /// Default constructor
  constexpr intersection2D() = default;

  /// Fully parametrized constructor - contains local position
  template <bool B = has_pos>
    requires B
  DETRAY_HOST_DEVICE constexpr intersection2D(
      const surface_type &sf, const scalar_type path, const point3_type &local,
      const nav_link_type nl, const dsimd<algebra_t, status_type> status,
      const bool_type dir)
      : m_surface{sf},
        m_ip{path, local},
        m_volume_link{nl},
        m_status{status},
        m_direction{dir} {}

  /// Fully parametrized constructor - contains no local position
  template <bool B = has_pos>
    requires(!B)
  DETRAY_HOST_DEVICE constexpr intersection2D(
      const surface_type &sf, const scalar_type path, const nav_link_type nl,
      const dsimd<algebra_t, status_type> status, const bool_type dir)
      : m_surface{sf},
        m_ip{path},
        m_volume_link{nl},
        m_status{status},
        m_direction{dir} {}

  /// @returns true if debug information needs to be filled
  static consteval bool contains_pos() { return has_pos; }

  /// @returns the intersected surface
  DETRAY_HOST_DEVICE
  constexpr surface_type surface() const { return m_surface; }

  /// @returns the intersected surface - non-const
  DETRAY_HOST_DEVICE
  constexpr surface_type &surface() { return m_surface; }

  /// Set the surface for this intersection to @param sf
  DETRAY_HOST_DEVICE
  constexpr void set_surface(const surface_type sf) { m_surface = sf; }

  /// @returns the link to the target volume to the navigator
  DETRAY_HOST_DEVICE
  constexpr nav_link_type volume_link() const { return m_volume_link; }

  /// Set the volume link to @param nl
  DETRAY_HOST_DEVICE
  constexpr void set_volume_link(const nav_link_type nl) { m_volume_link = nl; }

  /// @returns the direction of the intersection (before or behind the track)
  DETRAY_HOST_DEVICE
  constexpr bool_type is_along() const { return m_direction; }

  /// Set the direction to @param dir
  DETRAY_HOST_DEVICE
  constexpr void set_direction(const bool_type dir) { m_direction = dir; }

  /// @returns the distance of to the intersection point along the test traj.
  DETRAY_HOST_DEVICE
  constexpr scalar_type path() const { return m_ip.path; }

  /// Set the path to @param p
  DETRAY_HOST_DEVICE
  constexpr void set_path(scalar_type p) { m_ip.path = p; }

  /// @returns the local 3D intersection point (only debug)
  template <bool B = has_pos>
    requires B
  DETRAY_HOST_DEVICE constexpr point3_type local() const {
    return m_ip.point;
  }

  /// Set the intersection position to @param p
  template <bool B = has_pos>
    requires B
  DETRAY_HOST_DEVICE constexpr void set_local(const point3_type p) {
    m_ip.point = p;
  }

  /// @returns the intersection status
  DETRAY_HOST_DEVICE
  constexpr dsimd<algebra_t, status_type> status() const { return m_status; }

  /// Set the intersection status according to enum value @param s
  DETRAY_HOST_DEVICE
  constexpr void set_status(intersection::status s) {
    m_status = static_cast<status_type>(s);
  }

  /// Set the intersection status according to enum value @param s
  DETRAY_HOST_DEVICE
  constexpr void set_status_if(intersection::status s,
                               dbool<algebra_t> result_mask) {
    // @TODO find a unified conditional assignment in algebra_plugins
    if constexpr (concepts::soa<algebra_t>) {
      m_status(result_mask) = static_cast<status_type>(s);
    } else {
      m_status = result_mask ? static_cast<status_type>(s) : m_status;
    }
  }

  /// @note: Three way comparison cannot be used easily with SoA boolean masks
  /// @{
  /// @param rhs is the right hand side intersection for comparison
  DETRAY_HOST_DEVICE
  friend constexpr bool_type operator<(const intersection2D &lhs,
                                       const intersection2D &rhs) noexcept {
    return (math::fabs(lhs.path()) < math::fabs(rhs.path()));
  }

  /// @param rhs is the right hand side intersection for comparison
  DETRAY_HOST_DEVICE
  friend constexpr bool_type operator<=(const intersection2D &lhs,
                                        const intersection2D &rhs) noexcept {
    return (math::fabs(lhs.path()) <= math::fabs(rhs.path()));
  }

  /// @param rhs is the left hand side intersection for comparison
  DETRAY_HOST_DEVICE
  friend constexpr bool_type operator>(const intersection2D &lhs,
                                       const intersection2D &rhs) noexcept {
    return (math::fabs(lhs.path()) > math::fabs(rhs.path()));
  }

  /// @param rhs is the left hand side intersection for comparison
  DETRAY_HOST_DEVICE
  friend constexpr bool_type operator>=(const intersection2D &lhs,
                                        const intersection2D &rhs) noexcept {
    return (math::fabs(lhs.path()) > math::fabs(rhs.path()));
  }
  /// @}

  /// @param rhs is the left hand side intersection for comparison
  DETRAY_HOST_DEVICE
  friend constexpr bool_type operator==(const intersection2D &lhs,
                                        const intersection2D &rhs) noexcept {
    return math::fabs(lhs.path() - rhs.path()) <
           std::numeric_limits<float>::epsilon();
  }

  /// @returns true if any of the intersection results is 'inside'
  DETRAY_HOST_DEVICE
  constexpr bool is_inside() const {
    const dsimd<algebra_t, status_type> comp(
        static_cast<status_type>(intersection::status::e_inside));
    return detail::any_of(this->m_status == comp);
  }

  /// @returns true if any of the intersection results is 'edge'
  DETRAY_HOST_DEVICE
  constexpr bool is_edge() const {
    const dsimd<algebra_t, status_type> comp(
        static_cast<status_type>(intersection::status::e_edge));
    return detail::any_of(this->m_status == comp);
  }

  /// @returns true if any of the intersection results is 'inside' or 'edge'
  DETRAY_HOST_DEVICE
  constexpr bool is_probably_inside() const {
    const dsimd<algebra_t, status_type> comp(
        static_cast<status_type>(intersection::status::e_edge));
    return detail::any_of(this->m_status <= comp);
  }

  /// @returns true if all of the intersection results are 'outside'
  DETRAY_HOST_DEVICE
  constexpr bool is_outside() const {
    const dsimd<algebra_t, status_type> comp(
        static_cast<status_type>(intersection::status::e_outside));
    return detail::all_of(this->m_status == comp);
  }

 private:
  /// Transform to a string for output debugging
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &out_stream,
                                  const intersection2D &is) {
    out_stream << is.m_ip << ", "
               << ", surface: " << is.surface().identifier()
               << ", type: " << static_cast<int>(is.surface().mask().id())
               << ", links to vol:" << is.volume_link() << ")";

    if (is.is_inside()) {
      out_stream << ", status: inside";
    } else if (is.is_edge()) {
      out_stream << ", status: edge";
    } else {
      out_stream << ", status: outside";
    }
    if constexpr (std::is_scalar_v<bool_type>) {
      out_stream << (is.is_along() ? ", direction: along"
                                   : ", direction: opposite");
    } else {
      out_stream << ", status: " << is.status();
      out_stream << ", direction: " << is.is_along();
    }

    return out_stream;
  }

  /// Descriptor of the surface this intersection belongs to
  surface_type m_surface{};

  /// The intersection point (only saves the local point in debug mode)
  intersection_point<algebra_type, dpoint3D<algebra_type>, has_pos> m_ip{};

  /// Navigation information (next volume to go to)
  nav_link_type m_volume_link{detail::invalid_value<nav_link_type>()};

  /// Result of the intersection
  dsimd<algebra_t, status_type> m_status =
      static_cast<status_type>(intersection::status::e_outside);

  /// Direction of the intersection with respect to the track (true = along,
  /// false = opposite)
  bool_type m_direction{true};
};

}  // namespace detray
