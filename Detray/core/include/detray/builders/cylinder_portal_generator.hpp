// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/shapes/cuboid3D.hpp"
#include "detray/utils/bounding_volume.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <cassert>
#include <limits>

namespace detray {

/// @brief configuration for the cylinder portal generator
template <concepts::scalar scalar_t>
struct cylinder_portal_config {
  /// Build inner cylinder portal (will use the same distance to the layer
  /// that was found for the outer cylinder portal)
  bool m_build_inner{true};
  /// Autofit the lower/upper z extend and inner/outer radii
  bool m_do_autofit{true};
  /// Minimal envelope for the portals (used in autofitting)
  scalar_t m_envelope{100.f * unit<scalar_t>::um};
  /// Fixed inner radius during autofit
  scalar_t m_fixed_inner_r{0.f};
  /// Fixed outer radius during autofit
  scalar_t m_fixed_outer_r{0.f};
  /// Fixed length of the cylinder
  scalar_t m_fixed_z{0.f};
  /// The portal volumes links (north, south, east, west)
  std::vector<dindex> m_volume_links{dindex_invalid, dindex_invalid,
                                     dindex_invalid, dindex_invalid};

  /// Setters
  /// @{
  constexpr cylinder_portal_config &build_inner(const bool b) {
    m_build_inner = b;
    return *this;
  }
  constexpr cylinder_portal_config &do_autofit(const bool b) {
    m_do_autofit = b;
    return *this;
  }
  constexpr cylinder_portal_config &envelope(const scalar_t e) {
    m_envelope = e;
    return *this;
  }
  constexpr cylinder_portal_config &fixed_inner_radius(const scalar_t r) {
    m_fixed_inner_r = r;
    return *this;
  }
  constexpr cylinder_portal_config &fixed_outer_radius(const scalar_t r) {
    m_fixed_outer_r = r;
    return *this;
  }
  constexpr cylinder_portal_config &fixed_half_length(const scalar_t z) {
    m_fixed_z = z;
    return *this;
  }
  constexpr cylinder_portal_config &link_north(const dindex l) {
    m_volume_links[0] = l;
    return *this;
  }
  constexpr cylinder_portal_config &link_south(const dindex l) {
    m_volume_links[1] = l;
    return *this;
  }
  constexpr cylinder_portal_config &link_east(const dindex l) {
    m_volume_links[2] = l;
    return *this;
  }
  constexpr cylinder_portal_config &link_west(const dindex l) {
    m_volume_links[3] = l;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  constexpr bool build_inner() const { return m_build_inner; }
  constexpr bool do_autofit() const { return m_do_autofit; }
  constexpr scalar_t envelope() const { return m_envelope; }
  constexpr scalar_t fixed_inner_radius() const { return m_fixed_inner_r; }
  constexpr scalar_t fixed_outer_radius() const { return m_fixed_outer_r; }
  constexpr scalar_t fixed_half_length() const { return m_fixed_z; }
  constexpr const auto &volume_links() const { return m_volume_links; }
  constexpr const auto &link_north() const { return m_volume_links[0]; }
  constexpr const auto &link_south() const { return m_volume_links[1]; }
  constexpr const auto &link_east() const { return m_volume_links[2]; }
  constexpr const auto &link_west() const { return m_volume_links[3]; }
  /// @}
};

/// @brief Generates a portal box around a volume that already contains surfaces
///
/// @tparam detector_t the type of detector the volume belongs to.
template <typename detector_t>
class cylinder_portal_generator final
    : public surface_factory_interface<detector_t> {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;
  using transform3_t = dtransform3D<algebra_t>;

  using surface_type = typename detector_t::surface_type;
  using mask_id = typename detector_t::masks::id;
  using mask_link_t = typename surface_type::mask_link;
  using material_id = typename detector_t::material::id;
  using material_link_t = typename surface_type::material_link;

  /// A functor to construct global bounding boxes around masks
  struct bounding_box_creator {
    using aabb_t = axis_aligned_bounding_volume<cuboid3D, algebra_t>;

    template <typename mask_group_t, typename idx_range_t>
    DETRAY_HOST_DEVICE inline void operator()(
        const mask_group_t &mask_group, const idx_range_t &idx_range,
        const scalar_t envelope, const transform3_t &trf,
        std::vector<aabb_t> &boxes) const {
      for (const auto &mask : detray::ranges::subrange(mask_group, idx_range)) {
        // Local minimum bounding box
        aabb_t box{mask, boxes.size(), envelope};
        // Bounding box in global coordinates (might no longer be
        // minimum)
        boxes.push_back(box.transform(trf));
      }
    }
  };

 public:
  /// Save the boundaries of the cylinder after autofitting the portals
  struct boundaries {
    scalar_t inner_radius{0.f};
    scalar_t outer_radius{0.f};
    scalar_t lower_z{0.f};
    scalar_t upper_z{0.f};
  };

  /// Construct from configuration @param cfg
  DETRAY_HOST
  explicit cylinder_portal_generator(const cylinder_portal_config<scalar_t> cfg)
      : m_cfg{cfg} {}

  /// @returns the number of portals this factory will produce
  DETRAY_HOST
  auto size() const -> dindex override { return 4u; }

  DETRAY_HOST
  void clear() override { /*Do nothing*/ };

  DETRAY_HOST
  void push_back(surface_data<detector_t> &&) override { /*Do nothing*/ }
  DETRAY_HOST
  auto push_back(std::vector<surface_data<detector_t>> &&)
      -> void override { /*Do nothing*/ }

  /// @brief Access the volume boundaries after fitting
  const boundaries &volume_boundaries() const { return m_boundaries; }

  /// Create minimum aabbs around all surfaces that are passed and then
  /// construct cylinder and disc portals using the measures of the global
  /// aabb.
  ///
  /// @param volume the volume the portals need to be added to.
  /// @param surfaces the surface collection to wrap and to add the portals to
  /// @param transforms the transforms of the surfaces.
  /// @param masks the masks of the surfaces.
  /// @param ctx the geometry context (not needed for portals).
  DETRAY_HOST
  auto operator()(typename detector_t::volume_type &volume,
                  typename detector_t::surface_lookup_container &surfaces,
                  typename detector_t::transform_container &transforms,
                  typename detector_t::mask_container &masks,
                  typename detector_t::geometry_context ctx = {})
      -> dindex_range override {
    DETRAY_VERBOSE_HOST("Generate cylinder portals...");

    using aabb_t = axis_aligned_bounding_volume<cuboid3D, algebra_t>;

    // Only build portals for cylinder volumes
    assert(volume.id() == volume_id::e_cylinder);
    const dindex vol_idx{volume.index()};

    const std::size_t n_surfaces{surfaces.size()};

    if (!m_cfg.do_autofit()) {
      // Without autofit, the portal bounds have to be given explicitly
      assert(!(m_cfg.fixed_inner_radius() == 0.f &&
               m_cfg.fixed_outer_radius() == 0.f) ||
             m_cfg.fixed_half_length() != 0.f);
    } else {
      // Need surfaces in volume to do autofit
      assert(n_surfaces != 0u);
    }

    scalar_t inner_r{m_cfg.fixed_inner_radius()};
    scalar_t outer_r{m_cfg.fixed_outer_radius()};
    scalar_t lower_z{-m_cfg.fixed_half_length()};
    scalar_t upper_z{m_cfg.fixed_half_length()};

    if (m_cfg.do_autofit()) {
      // The bounding boxes around the module surfaces
      std::vector<aabb_t> boxes;
      boxes.reserve(n_surfaces);

      for (const auto &sf : surfaces) {
        masks.template visit<bounding_box_creator>(
            sf.mask(), m_cfg.envelope(), transforms.at(sf.transform(), ctx),
            boxes);
      }

      // Build an aabb in the global space around the surface aabbs
      aabb_t world_box{boxes, boxes.size(), m_cfg.envelope()};
      // translation
      const point3_t center = world_box.template center<point3_t>();

      // The world box local frame is the global coordinate frame
      const point3_t box_min = world_box.template loc_min<point3_t>();
      const point3_t box_max = world_box.template loc_max<point3_t>();

      // Get the half lengths for the cylinder height and disc translation
      const point3_t h_lengths = 0.5f * (box_max - box_min);
      const scalar_t h_x{math::fabs(h_lengths[0])};
      const scalar_t h_y{math::fabs(h_lengths[1])};
      const scalar_t h_z{math::fabs(h_lengths[2])};

      const scalar_t outer_r_min{math::max(h_x, h_y)};
      const scalar_t mean_radius{get_mean_radius(surfaces, transforms)};

      outer_r = outer_r_min;
      inner_r = mean_radius - (outer_r - mean_radius);
      lower_z = math::min(center[2] - h_z, center[2] + h_z);
      upper_z = math::max(center[2] - h_z, center[2] + h_z);
    }

    // Observe boundary conditions
    if (m_cfg.fixed_inner_radius() > 0.f) {
      inner_r = m_cfg.fixed_inner_radius();
    }
    if (m_cfg.fixed_outer_radius() > 0.f) {
      outer_r = m_cfg.fixed_outer_radius();
    }
    if (m_cfg.fixed_half_length() > 0.f) {
      lower_z = -m_cfg.fixed_half_length();
      upper_z = m_cfg.fixed_half_length();
    }

    m_boundaries = {inner_r, outer_r, lower_z, upper_z};

    // If inner radius is 0, skip adding the inner cylinder
    if (m_cfg.build_inner()) {
      add_cylinder_portal(surfaces, transforms, masks, ctx, vol_idx,
                          m_cfg.link_south(), inner_r, lower_z, upper_z);
    }
    add_cylinder_portal(surfaces, transforms, masks, ctx, vol_idx,
                        m_cfg.link_north(), outer_r, lower_z, upper_z);

    add_disc_portal(surfaces, transforms, masks, ctx, vol_idx,
                    m_cfg.link_west(), inner_r, outer_r, lower_z);
    add_disc_portal(surfaces, transforms, masks, ctx, vol_idx,
                    m_cfg.link_east(), inner_r, outer_r, upper_z);

    return {static_cast<dindex>(n_surfaces),
            static_cast<dindex>(surfaces.size())};
  }

 private:
  /// @brief Add a single cylinder portal to the internal data storage
  void add_cylinder_portal(
      typename detector_t::surface_lookup_container &surfaces,
      typename detector_t::transform_container &transforms,
      typename detector_t::mask_container &masks,
      typename detector_t::geometry_context ctx, dindex vol_idx,
      dindex vol_link, const scalar_t r, const scalar_t lower_z,
      const scalar_t upper_z) const {
    constexpr auto invalid_src_link{detail::invalid_value<std::uint64_t>()};

    const scalar_t min_z{math::min(lower_z, upper_z)};
    const scalar_t max_z{math::max(lower_z, upper_z)};

    // translation
    const point3_t tsl{0.f, 0.f, 0.f};

    // Add transform and mask data
    transforms.emplace_back(ctx, tsl);
    masks.template emplace_back<mask_id::e_concentric_cylinder2D>(
        empty_context{}, vol_link, r, min_z, max_z);

    // Add surface links
    mask_link_t mask_link{};
    const auto mask_idx{static_cast<dindex>(
        masks.template size<mask_id::e_concentric_cylinder2D>() - 1u)};
    set_mask_link(mask_link, mask_id::e_concentric_cylinder2D, mask_idx);

    material_link_t material_link{material_id::e_none, dindex_invalid};

    surfaces.push_back(
        {static_cast<dindex>(transforms.size(ctx) - 1u), mask_link,
         material_link, vol_idx, surface_id::e_portal},
        invalid_src_link);
  }

  /// @brief Add a single disc portal
  void add_disc_portal(typename detector_t::surface_lookup_container &surfaces,
                       typename detector_t::transform_container &transforms,
                       typename detector_t::mask_container &masks,
                       typename detector_t::geometry_context ctx,
                       dindex vol_idx, dindex vol_link, const scalar_t inner_r,
                       const scalar_t outer_r, const scalar_t z) const {
    constexpr auto invalid_src_link{detail::invalid_value<std::uint64_t>()};

    const scalar_t min_r{math::min(inner_r, outer_r)};
    const scalar_t max_r{math::max(inner_r, outer_r)};

    // translation
    point3_t tsl{scalar_t(0), scalar_t(0), z};

    // Add transform and mask data
    transforms.emplace_back(ctx, tsl);
    masks.template emplace_back<mask_id::e_ring2D>(empty_context{}, vol_link,
                                                   min_r, max_r);

    // Add surface links
    mask_link_t mask_link{};
    const auto mask_idx{
        static_cast<dindex>(masks.template size<mask_id::e_ring2D>() - 1u)};
    set_mask_link(mask_link, mask_id::e_ring2D, mask_idx);

    material_link_t material_link{material_id::e_none, dindex_invalid};

    surfaces.push_back(
        {static_cast<dindex>(transforms.size(ctx) - 1u), mask_link,
         material_link, vol_idx, surface_id::e_portal},
        invalid_src_link);
  }

  /// @returns Calculate the mean radius of all sensitive surfaces
  template <typename surface_container_t, typename transform_container_t>
  scalar_t get_mean_radius(const surface_container_t &surfaces,
                           const transform_container_t &transforms) const {
    double mean{0.};

    for (const auto &sf_desc : surfaces) {
      const auto &trf = transforms.at(sf_desc.transform());
      mean += vector::perp(trf.translation());
    }

    return static_cast<scalar_t>(mean / static_cast<double>(surfaces.size()));
  }

  /// Set the mask link for a single portal, either range or single index
  void set_mask_link(mask_link_t &mask_link, const mask_id id,
                     const dindex lower_idx) const {
    mask_link.set_id(id);
    if constexpr (concepts::interval<typename mask_link_t::index_type>) {
      mask_link.set_index({lower_idx, 1u});
    } else {
      mask_link.set_index(lower_idx);
    }
  }

  /// Portal generator configuration
  cylinder_portal_config<scalar_t> m_cfg;
  /// Save the dimensions of the volume after autofitting
  boundaries m_boundaries{};
};

}  // namespace detray
