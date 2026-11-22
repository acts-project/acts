// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/core/name_map.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/geometry/detail/volume_kernels.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/material/material.hpp"
#include "detray/navigation/accelerators/search_window.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <algorithm>
#include <iostream>
#include <limits>
#include <sstream>

namespace detray {

/// @brief Facade for a detray detector volume.
///
/// Volume class that acts as a logical container in the detector for geometry
/// objects, i.e. surfaces. The volume boundary surfaces, the so-called portals,
/// carry index links that join adjacent volumes. The volume class itself does
/// not contain any data itself, but keeps a descriptor with index-based links
/// into the data containers that are managed by the detector.
/// Every type of surface that is known by the volume (determined by the type
/// and size of the ID enum in the descriptor) lives in it's own geometry
/// accelerator data structure, e.g. portals reside in a brute force
/// accelerator (a simple vector), while sensitive surfaces are usually sorted
/// into a spatial grid.
template <concepts::detector detector_t>  // @TODO: This needs a concept
class tracking_volume {
  /// Linear algebra types
  using algebra_type = typename detector_t::algebra_type;
  using scalar_type = dscalar<algebra_type>;
  using point3_type = dpoint3D<algebra_type>;

  /// Volume descriptor type
  using descr_t = typename detector_t::volume_type;
  using context_t = typename detector_t::geometry_context;

 public:
  /// In case the geometry needs to be printed
  using name_map = detray::name_map;
  /// Object identifier enum used by the volume descriptor's accelerator links.
  using object_id = descr_t::object_id;

  /// Not allowed: always needs a detector and a descriptor.
  tracking_volume() = delete;

  /// Construct a volume facade from its detector and descriptor.
  /// @param det Detector that owns the volume data.
  /// @param desc Descriptor containing links for this volume.
  constexpr tracking_volume(const detector_t &det, const descr_t &desc)
      : m_detector{det}, m_desc{desc} {
    assert(m_desc.index() < det.volumes().size());
    assert(m_desc.id() != volume_id::e_unknown);
  }

  /// Constructor from detector @param det and volume index @param vol_idx in
  /// that detector.
  constexpr tracking_volume(const detector_t &det, const dindex vol_idx)
      : tracking_volume(det, det.volume(vol_idx)) {}

  /// @returns access to the underlying detector
  DETRAY_HOST_DEVICE
  auto detector() const -> const detector_t & { return m_detector; }

  /// Equality operator
  ///
  /// @param rhs is the right hand side to be compared to
  /// @returns True when both facades refer to the same detector and descriptor.
  DETRAY_HOST_DEVICE
  constexpr auto operator==(const tracking_volume &rhs) const -> bool {
    return (&m_detector == &(rhs.m_detector) && m_desc == rhs.m_desc);
  }

  /// @returns the volume shape id, e.g. 'cylinder'.
  DETRAY_HOST_DEVICE
  constexpr auto id() const -> volume_id { return m_desc.id(); }

  /// @returns the index of the volume in the detector volume container.
  DETRAY_HOST_DEVICE
  constexpr auto index() const -> dindex { return m_desc.index(); }

  /// @returns the index of the mother volume in the detector volume
  /// container. The index is invalid for a top-level volume.
  DETRAY_HOST_DEVICE
  constexpr auto mother() const -> dindex { return m_desc.mother(); }

  /// @returns true if the volume is not placed inside a mother volume
  DETRAY_HOST_DEVICE
  constexpr auto is_top_level() const -> bool { return !m_desc.has_mother(); }

  /// @returns true if the volume contains daughter volumes
  DETRAY_HOST_DEVICE
  constexpr auto has_daughters() const -> bool {
    return !m_desc.template accel_link<object_id::e_volume>().is_invalid();
  }

  /// Check whether a point lies inside the volume
  ///
  /// The check uses only the portals of the volume. The volume of the
  /// daughter volumes is part of the volume. Depending on the volume shape:
  /// - cuboid: The point must lie on the inner side of every portal plane.
  ///   The inner side is the side of the mean portal position.
  /// - cylinder: The discs give the radial bounds and the cylinder mantles
  ///   give the bounds in z (in the volume frame).
  /// - other shapes: A ray from the point crosses the portals an odd number
  ///   of times, if the point is inside. If a ray passes close to a portal
  ///   edge, the check is repeated with another ray direction.
  ///
  /// @param glob_p the point in global coordinates
  /// @param ctx the geometry context
  /// @param tol tolerance for points on the boundary (and for portal edges)
  ///
  /// @returns whether the point is inside, outside or on the boundary
  DETRAY_HOST_DEVICE constexpr auto is_inside(
      const point3_type &glob_p, const context_t &ctx = {},
      const scalar_type tol = 1.f * unit<scalar_type>::um) const
      -> volume_containment {
    if (id() == volume_id::e_cuboid) {
      return is_inside_cuboid(glob_p, ctx, tol);
    }
    if (id() == volume_id::e_cylinder) {
      bool is_supported{true};
      const volume_containment result{
          is_inside_cylinder(glob_p, ctx, tol, is_supported)};
      if (is_supported) {
        return result;
      }
    }
    return is_inside_ray_test(glob_p, ctx, tol);
  }

  /// @returns the (non contextual) transform for the placement of the
  /// volume in the detector geometry.
  DETRAY_HOST_DEVICE
  constexpr auto transform() const -> const
      typename detector_t::transform3_type & {
    return m_detector.transform_store().at(m_desc.transform());
  }

  /// @returns the center point of the volume.
  DETRAY_HOST_DEVICE
  constexpr auto center() const -> point3_type {
    return transform().translation();
  }

  /// @returns a pointer to the material parameters at the local position
  /// @param loc_p
  DETRAY_HOST_DEVICE constexpr auto material_parameters(
      const point3_type &loc_p) const -> const detray::material<scalar_type> * {
    return visit_material<typename detail::get_material_params>(loc_p);
  }

  /// @returns true if the volume carries material.
  DETRAY_HOST_DEVICE
  constexpr bool has_material() const { return m_desc.has_material(); }

  /// @returns an iterator pair for the requested type of surfaces.
  template <surface_id sf_type = surface_id::e_all>
  DETRAY_HOST_DEVICE constexpr decltype(auto) surfaces() const {
    if constexpr (sf_type == surface_id::e_all) {
      return detray::ranges::subrange{m_detector.surfaces(),
                                      m_desc.full_sf_range()};
    } else {
      return detray::ranges::subrange{m_detector.surfaces(),
                                      m_desc.template sf_link<sf_type>()};
    }
  }

  /// @returns an iterator pair for the volume portals.
  DETRAY_HOST_DEVICE constexpr decltype(auto) portals() const {
    return surfaces<surface_id::e_portal>();
  }

  /// @returns the total number of portal surfaces contained in the volume
  DETRAY_HOST_DEVICE constexpr dindex n_portals() const {
    const auto pt_idx_range = m_desc.template sf_link<surface_id::e_portal>();

    assert(pt_idx_range[1] > pt_idx_range[0]);
    return pt_idx_range[1] - pt_idx_range[0];
  }

  /// @returns the total number of sensitive surfaces contained in the volume
  DETRAY_HOST_DEVICE constexpr dindex n_sensitives() const {
    const auto sens_idx_range =
        m_desc.template sf_link<surface_id::e_sensitive>();

    assert(sens_idx_range[1] >= sens_idx_range[0]);
    return sens_idx_range[1] - sens_idx_range[0];
  }

  /// @returns the total number of passive surfaces contained in the volume
  DETRAY_HOST_DEVICE constexpr dindex n_passives() const {
    const auto ps_idx_range = m_desc.template sf_link<surface_id::e_passive>();

    assert(ps_idx_range[1] >= ps_idx_range[0]);
    return ps_idx_range[1] - ps_idx_range[0];
  }

  /// Apply a functor to all surfaces in one of the volume's acceleration
  /// structures
  ///
  /// @tparam I type of object to retrieve (passive, portal, sensitive etc)
  /// @tparam functor_t the prescription to be applied to the surfaces
  /// @tparam Args      types of additional arguments to the functor
  /// @param args Arguments forwarded to the accelerator visitor.
  template <object_id I, typename functor_t, typename... Args>
  DETRAY_HOST_DEVICE constexpr void visit_accelerator(Args &&...args) const {
    static_assert(I < object_id::e_all);

    if (const auto &link{m_desc.template accel_link<I>()}; !link.is_invalid()) {
      // Run over the surfaces in a single acceleration data structure
      // and apply the functor to the resulting neighborhood
      m_detector.accelerator_store().template visit<functor_t>(
          link, std::forward<Args>(args)...);
    }
  }

  /// Apply a functor to all acceleration structures of this volume.
  ///
  /// @tparam I type of object to retrieve (surface types, daughters etc)
  /// @tparam functor_t the prescription to be applied to the acc structure
  /// @tparam Args      types of additional arguments to the functor
  /// @param args Arguments forwarded to each accelerator visitor.
  template <typename functor_t, int I = static_cast<int>(object_id::e_all) - 1,
            typename... Args>
  DETRAY_HOST_DEVICE constexpr void visit_accelerators(Args &&...args) const {
    // Get the acceleration data structures for this volume and only visit,
    // if object type is contained in volume
    for (std::size_t id = 0u; id < static_cast<std::size_t>(object_id::e_size);
         ++id) {
      if (const auto &link{m_desc.accel_link()[id]}; !link.is_invalid()) {
        // Run over the surfaces in a single acceleration data structure
        // and apply the functor to the resulting neighborhood
        m_detector.accelerator_store().template visit<functor_t>(
            link, std::forward<Args>(args)...);
      }
    }
  }

  /// Apply a functor to all surfaces of a given surface id (portal, passive,
  /// sensitive) in the volume
  ///
  /// Translates the detray surface type id to the volume geometry object id
  ///
  /// @tparam functor_t the prescription to be applied to the surfaces
  /// @tparam Args      types of additional arguments to the functor
  /// @param args Arguments forwarded to the surface visitor.
  template <surface_id I, typename functor_t, typename... Args>
  DETRAY_HOST_DEVICE constexpr void visit_surfaces(Args &&...args) const {
    using surface_getter_t = detail::apply_to_surfaces<functor_t>;

    // Dispatch to the correct acceleration structure
    if constexpr (I == surface_id::e_portal) {
      visit_accelerator<object_id::e_portal, surface_getter_t>(
          std::forward<Args>(args)...);
    } else if constexpr (I == surface_id::e_sensitive) {
      visit_accelerator<object_id::e_sensitive, surface_getter_t>(
          std::forward<Args>(args)...);
    } else if constexpr (I == surface_id::e_passive) {
      visit_accelerator<object_id::e_passive, surface_getter_t>(
          std::forward<Args>(args)...);
    } else {
      // Visit all surface types, but not other geometric objects
      // (e.g. daughter volumes)
      visit_accelerator<object_id::e_portal, surface_getter_t>(
          std::forward<Args>(args)...);
      if constexpr (object_id::e_portal != object_id::e_passive) {
        visit_accelerator<object_id::e_passive, surface_getter_t>(
            std::forward<Args>(args)...);
      }
      visit_accelerator<object_id::e_sensitive, surface_getter_t>(
          std::forward<Args>(args)...);
    }
  }

  /// Apply a functor to all daughter volumes (every daughter volume once)
  ///
  /// @tparam functor_t the prescription to be applied to the daughter volumes
  /// @tparam Args      types of additional arguments to the functor
  /// @param args Arguments forwarded to the daughter-volume visitor.
  template <typename functor_t, typename... Args>
  DETRAY_HOST constexpr void visit_daughter_volumes(Args &&...args) const {
    using volume_getter_t = detail::apply_to_volumes<functor_t>;

    visit_accelerator<object_id::e_volume, volume_getter_t>(
        std::forward<Args>(args)...);
  }

  /// Apply a functor to a neighborhood of geometric objects around a
  /// track position in the volume.
  ///
  /// @note: The acceleration structures in the volume might return different
  /// geometric objects (e.g. surfaces vs. volumes). The passed functor must
  /// provide corresponding overloads of the call operator.
  ///
  /// @tparam functor_t the prescription to be applied to the surfaces
  ///                   (customization point for the navigation)
  /// @tparam track_t   the track around which to build up the neighborhood
  /// @tparam Args      types of additional arguments to the functor
  /// @param track Track position and direction used for neighborhood lookup.
  /// @param win_size Search window around the track position.
  /// @param ctx Geometry context used to resolve transforms.
  /// @param args Additional arguments forwarded to the neighborhood visitor.
  template <object_id I, typename functor_t, typename track_t,
            concepts::arithmetic window_size_t, typename... Args>
  DETRAY_HOST_DEVICE constexpr void visit_neighborhood(
      const track_t &track, const search_window<window_size_t, 2> &win_size,
      const context_t &ctx, Args &&...args) const {
    if constexpr (I == object_id::e_all) {
      visit_accelerators<detail::apply_to_neighbourhood<functor_t>>(
          m_detector, m_desc, track, win_size, ctx,
          std::forward<Args>(args)...);
    } else {
      visit_accelerator<I, detail::apply_to_neighbourhood<functor_t>>(
          m_detector, m_desc, track, win_size, ctx,
          std::forward<Args>(args)...);
    }
  }

  /// Call a functor on the volume material with additional arguments.
  ///
  /// @tparam functor_t the prescription to be applied to the material
  /// @tparam Args      types of additional arguments to the functor
  /// @param args Arguments forwarded to the material visitor.
  /// @returns The result produced by the material visitor.
  template <typename functor_t, typename... Args>
  DETRAY_HOST_DEVICE constexpr auto visit_material(Args &&...args) const {
    assert(has_material());
    const auto &materials = m_detector.material_store();
    return materials.template visit<functor_t>(m_desc.material(),
                                               std::forward<Args>(args)...);
  }

  /// Do a consistency check on the volume after building the detector.
  ///
  /// @param os output stream for error messages.
  ///
  /// @returns true if the volume is consistent
  DETRAY_HOST bool self_check(std::ostream &os) const {
    if (id() == volume_id::e_unknown) {
      os << "DETRAY ERROR (HOST): Unknown volume shape type in volume:\n"
         << *this << std::endl;
      return false;
    }
    if (detail::is_invalid_value(index())) {
      os << "DETRAY ERROR (HOST): Volume index undefined in volume:\n"
         << *this << std::endl;
      return false;
    }
    if (index() >= m_detector.volumes().size()) {
      os << "DETRAY ERROR (HOST): Volume index out of bounds in volume:\n"
         << *this << std::endl;
      return false;
    }
    if (detail::is_invalid_value(m_desc.transform())) {
      os << "DETRAY ERROR (HOST): Volume transform undefined in volume:\n"
         << *this << std::endl;
      return false;
    }
    if (m_desc.transform() >= m_detector.transform_store().size()) {
      os << "DETRAY ERROR (HOST): Volume transform index out of bounds "
            "in volume:\n"
         << *this << std::endl;
      return false;
    }
    if (m_desc.has_mother()) {
      if (m_desc.mother() >= m_detector.volumes().size()) {
        os << "DETRAY ERROR (HOST): Mother volume index out of bounds in "
              "volume:\n"
           << *this << std::endl;
        return false;
      }
      if (m_desc.mother() == index()) {
        os << "DETRAY ERROR (HOST): Volume is its own mother:\n"
           << *this << std::endl;
        return false;
      }
    }
    // Only check, if there is material in the detector
    if (!m_detector.material_store().all_empty() && has_material() &&
        m_desc.material().is_invalid_index()) {
      os << "DETRAY ERROR (HOST): Volume does not have valid material "
            "link:\n"
         << *this << std::endl;
      return false;
    }
    const auto &acc_link = m_desc.accel_link();
    if (detail::is_invalid_value(acc_link[0])) {
      os << "DETRAY ERROR (HOST): Link to portal lookup broken: " << acc_link[0]
         << "\n in volume: " << *this << std::endl;
      return false;
    }
    if (const auto &pt_link = m_desc.template sf_link<surface_id::e_portal>();
        detail::is_invalid_value(pt_link)) {
      os << "DETRAY ERROR (HOST): Link to portal surfaces broken: " << pt_link
         << "\n in volume: " << *this << std::endl;
      return false;
    }
    // Check consistency of surface ranges
    std::vector<dindex_range> sf_ranges = {
        m_desc.template sf_link<surface_id::e_portal>()};

    // Only add the other ranges in case they are not empty
    if (const auto &sens_range =
            m_desc.template sf_link<surface_id::e_sensitive>();
        sens_range[0] != sens_range[1]) {
      sf_ranges.push_back(sens_range);
    }
    if (const auto &psv_range =
            m_desc.template sf_link<surface_id::e_passive>();
        psv_range[0] != psv_range[1]) {
      sf_ranges.push_back(psv_range);
    }

    // Sort and check that the ranges are contiguous
    auto compare_ranges = [](const dindex_range &rg1, const dindex_range &rg2) {
      return rg1[0] < rg2[0];
    };

    std::ranges::sort(sf_ranges, compare_ranges);

    if ((sf_ranges.size() > 1 && sf_ranges[0][1] != sf_ranges[1][0]) ||
        (sf_ranges.size() > 2 && sf_ranges[1][1] != sf_ranges[2][0])) {
      os << "DETRAY ERROR (HOST): Surface index ranges not contiguous: "
         << m_desc.sf_link() << "\n in volume: " << *this << std::endl;
      return false;
    }

    // Warnings
    bool suspicious_links = false;
    std::stringstream warnings{};
    for (std::size_t i = 1u; i < acc_link.size(); ++i) {
      // An acceleration data structure link was set, but is invalid
      if (!acc_link[i].is_invalid_id() && acc_link[i].is_invalid_index()) {
        suspicious_links = true;
        warnings << "Link to acceleration data structure "
                 << static_cast<int>(acc_link[i].id()) << " is invalid"
                 << std::endl;
      }
    }
    if (suspicious_links) {
      DETRAY_WARN_HOST(warnings.str() << " in volume: " << *this);
    }

    return true;
  }

  /// @param names Optional detector volume names indexed by volume index.
  /// @returns the volume name (add an offset for the detector name).
  DETRAY_HOST_DEVICE
  auto name(const name_map &names) const -> std::string {
    return names.empty() ? "volume " + std::to_string(m_desc.index())
                         : names.at(m_desc.index());
  }

  /// @returns a string stream that prints the volume details
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &os, const tracking_volume &v) {
    os << v.m_desc;
    return os;
  }

 private:
  /// @returns the containment of a point from its minimal distance @p dist
  /// to the boundary (negative: outside) with tolerance @p tol
  DETRAY_HOST_DEVICE static constexpr auto to_containment(
      const scalar_type dist, const scalar_type tol) -> volume_containment {
    if (dist < -tol) {
      return volume_containment::e_outside;
    }
    return (dist <= tol) ? volume_containment::e_on_boundary
                         : volume_containment::e_inside;
  }

  /// Inside check for a cuboid volume: The volume is the intersection of
  /// the half-spaces of its portal planes. The inner side of every plane is
  /// the side of the mean portal position, which lies inside of the (convex)
  /// volume.
  ///
  /// @param glob_p the point in global coordinates
  /// @param ctx the geometry context
  /// @param tol tolerance for points on the boundary
  ///
  /// @returns whether the point is inside, outside or on the boundary
  DETRAY_HOST_DEVICE constexpr auto is_inside_cuboid(
      const point3_type &glob_p, const context_t &ctx,
      const scalar_type tol) const -> volume_containment {
    // Reference point inside of the volume
    point3_type center{0.f, 0.f, 0.f};
    for (const auto &pt_desc : portals()) {
      center =
          center +
          geometry::surface{m_detector, pt_desc}.transform(ctx).translation();
    }
    center = (1.f / static_cast<scalar_type>(n_portals())) * center;

    // Minimal distance to the planes (negative: outside)
    scalar_type dist{std::numeric_limits<scalar_type>::max()};
    for (const auto &pt_desc : portals()) {
      const auto &trf = geometry::surface{m_detector, pt_desc}.transform(ctx);
      const auto &n = trf.z();
      const scalar_type d_center{vector::dot(n, center - trf.translation())};
      const scalar_type d_point{vector::dot(n, glob_p - trf.translation())};
      dist = math::min(dist, d_center < 0.f ? -d_point : d_point);
    }

    return to_containment(dist, tol);
  }

  /// Inside check for a cylinder volume, whose portals are discs and
  /// cylinder mantles
  ///
  /// @param glob_p the point in global coordinates
  /// @param ctx the geometry context
  /// @param tol tolerance for points on the boundary
  /// @param is_supported set to false, if the portals are not supported
  ///                     (e.g. planes for phi sectors) (output)
  ///
  /// @returns whether the point is inside, outside or on the boundary
  DETRAY_HOST_DEVICE constexpr auto is_inside_cylinder(
      const point3_type &glob_p, const context_t &ctx, const scalar_type tol,
      bool &is_supported) const -> volume_containment {
    constexpr scalar_type inf{std::numeric_limits<scalar_type>::max()};

    const auto &vol_trf =
        m_detector.transform_store().at(m_desc.transform(), ctx);

    // r_min, r_max, z_min, z_max
    darray<scalar_type, 4> bounds{inf, -inf, inf, -inf};
    for (const auto &pt_desc : portals()) {
      const auto pt = geometry::surface{m_detector, pt_desc};
      const auto &pt_trf = pt.transform(ctx);

      const scalar_type z_offset{
          vol_trf.point_to_local(pt_trf.translation())[2]};
      const scalar_type z_sign{
          vector::dot(vol_trf.z(), pt_trf.z()) < 0.f ? -1.f : 1.f};
      pt.template visit_mask<detail::cylinder_volume_bounds>(
          z_offset, z_sign, bounds, is_supported);
    }

    // At least one disc and one mantle are needed
    if (!is_supported || bounds[1] < bounds[0] || bounds[3] < bounds[2]) {
      is_supported = false;
      return volume_containment::e_ambiguous;
    }

    const point3_type loc{vol_trf.point_to_local(glob_p)};
    const scalar_type r{vector::perp(loc)};

    // Without an inner mantle, the volume reaches down to the z-axis
    const scalar_type dist_r_min{bounds[0] > tol ? r - bounds[0] : inf};
    const scalar_type dist{
        math::min(math::min(dist_r_min, bounds[1] - r),
                  math::min(loc[2] - bounds[2], bounds[3] - loc[2]))};

    return to_containment(dist, tol);
  }

  /// Inside check for volumes of any shape: A ray from the point crosses the
  /// portals an odd number of times, if the point is inside. If a ray passes
  /// close to a portal edge, the check is repeated with another ray direction.
  ///
  /// @param glob_p the point in global coordinates
  /// @param ctx the geometry context
  /// @param tol tolerance for points on a portal and for portal edges
  ///
  /// @returns whether the point is inside, outside, on the boundary, or
  /// whether the check is ambiguous
  DETRAY_HOST_DEVICE constexpr auto is_inside_ray_test(
      const point3_type &glob_p, const context_t &ctx,
      const scalar_type tol) const -> volume_containment {
    using vector3_t = dvector3D<algebra_type>;

    // Test directions that are not aligned with the coordinate axes
    const darray<vector3_t, 3> directions{
        vector::normalize(vector3_t{0.31f, 0.47f, 0.83f}),
        vector::normalize(vector3_t{-0.71f, 0.29f, -0.64f}),
        vector::normalize(vector3_t{0.53f, -0.81f, 0.25f})};

    for (const vector3_t &dir : directions) {
      const detail::ray<algebra_type> ray{glob_p, 0.f, dir, 0.f};

      unsigned int n_crossings{0u};
      bool is_ambiguous{false};
      bool is_on_portal{false};
      for (const auto &pt_desc : portals()) {
        const auto pt = geometry::surface{m_detector, pt_desc};
        pt.template visit_mask<detail::count_portal_crossings>(
            ray, pt.transform(ctx), tol, n_crossings, is_ambiguous,
            is_on_portal);
      }

      if (is_on_portal) {
        return volume_containment::e_on_boundary;
      }
      if (!is_ambiguous) {
        return (n_crossings % 2u == 1u) ? volume_containment::e_inside
                                        : volume_containment::e_outside;
      }
    }

    return volume_containment::e_ambiguous;
  }

  /// Access to the detector stores
  const detector_t &m_detector;
  /// Access to the descriptor
  const descr_t &m_desc;
};

}  // namespace detray
