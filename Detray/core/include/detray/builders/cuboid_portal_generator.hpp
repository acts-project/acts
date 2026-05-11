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
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/math.hpp"
#include "detray/geometry/shapes/cuboid3D.hpp"
#include "detray/utils/bounding_volume.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <cassert>
#include <limits>

namespace detray {

/// @brief Generates a portal box around a volume that already contains surfaces
///
/// @tparam detector_t the type of detector the volume belongs to.
template <typename detector_t>
class cuboid_portal_generator final
    : public surface_factory_interface<detector_t> {
  using algebra_type = typename detector_t::algebra_type;
  using scalar_type = dscalar<algebra_type>;
  using transform3_type = dtransform3D<algebra_type>;

  /// A functor to construct global bounding boxes around masks
  struct bounding_box_creator {
    using aabb_t = axis_aligned_bounding_volume<cuboid3D, algebra_type>;

    template <typename mask_group_t, typename idx_range_t>
    DETRAY_HOST_DEVICE inline void operator()(
        const mask_group_t &mask_group, const idx_range_t &idx_range,
        const scalar_type envelope, const transform3_type &trf,
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
  /// Use @param env as portal envelope
  DETRAY_HOST
  explicit cuboid_portal_generator(const scalar_type env) : m_envelope{env} {}

  /// @returns the number of rectangle portals this factory will produce
  DETRAY_HOST
  auto size() const -> dindex override { return 6u; }

  DETRAY_HOST
  void clear() override { /*Do nothing*/ };

  DETRAY_HOST
  void push_back(
      surface_data<detector_t> && /*unused*/) override { /*Do nothing*/ }
  DETRAY_HOST
  auto push_back(std::vector<surface_data<detector_t>> && /*unused*/)
      -> void override { /*Do nothing*/ }

  /// Create minimum aabbs around all surfaces that are passed and then
  /// construct world portals from the global bounding box.
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
    DETRAY_VERBOSE_HOST("Generate cuboid portals...");

    using point3_t = dpoint3D<algebra_type>;
    using vector3_t = dvector3D<algebra_type>;

    using surface_t = typename detector_t::surface_type;
    using nav_link_t = typename surface_t::navigation_link;
    using mask_link_t = typename surface_t::mask_link;
    using material_link_t = typename surface_t::material_link;

    using aabb_t = axis_aligned_bounding_volume<cuboid3D, algebra_type>;

    constexpr auto invalid_src_link{detail::invalid_value<std::uint64_t>()};

    // Only build box portals for cuboid volumes
    assert(volume.id() == volume_id::e_cuboid);
    // Need surfaces to wrap
    std::size_t n_surfaces{surfaces.size()};
    assert(n_surfaces != 0u);

    // The surfaces container is prefilled with other surfaces
    auto surfaces_offset{static_cast<dindex>(n_surfaces)};

    // Fetch the position in the mask tuple for the rectangle portals
    constexpr auto rectangle_id{detector_t::masks::id::e_rectangle2D};

    // The material will be added in a later step
    constexpr auto no_material = surface_t::material_id::e_none;
    material_link_t material_link{no_material, dindex_invalid};

    // Max distance in case of infinite bounds
    constexpr scalar_type max_shift{0.01f *
                                    std::numeric_limits<scalar_type>::max()};

    // The bounding boxes around the module surfaces
    std::vector<aabb_t> boxes;
    boxes.reserve(n_surfaces);

    for (const auto &sf : surfaces) {
      masks.template visit<bounding_box_creator>(
          sf.mask(), m_envelope, transforms.at(sf.transform(), ctx), boxes);
    }
    // Build an aabb in the global space around the surface aabbs
    aabb_t world_box{boxes, boxes.size(), m_envelope};

    // translation
    const point3_t center = world_box.template center<point3_t>();

    // The world box local frame is the global coordinate frame
    const point3_t box_min = world_box.template loc_min<point3_t>();
    const point3_t box_max = world_box.template loc_max<point3_t>();

    // Get the half lengths for the rectangle sides and translation
    const point3_t h_lengths = 0.5f * (box_max - box_min);
    const scalar_type h_x{math::fabs(h_lengths[0])};
    const scalar_type h_y{math::fabs(h_lengths[1])};
    const scalar_type h_z{math::fabs(h_lengths[2])};

    // Volume links for the portal descriptors and the masks
    const dindex volume_idx{volume.index()};
    const nav_link_t volume_link{detail::invalid_value<nav_link_t>()};

    // Construct portals in the...

    //
    // ... x-y plane
    //
    // Only one rectangle needed for both surfaces
    mask_link_t mask_link{};
    mask_link.set_id(rectangle_id);
    const auto mask_idx{
        static_cast<dindex>(masks.template size<rectangle_id>())};
    if constexpr (concepts::interval<typename mask_link_t::index_type>) {
      mask_link.set_index({mask_idx, 1u});
    } else {
      mask_link.set_index(mask_idx);
    }
    masks.template emplace_back<rectangle_id>(empty_context{}, volume_link, h_x,
                                              h_y);

    // No rotation, but shift in z for both faces
    vector3_t shift{static_cast<scalar_type>(0), static_cast<scalar_type>(0),
                    detail::is_invalid_value(h_z) ? max_shift : h_z};
    transforms.emplace_back(ctx, static_cast<vector3_t>(center + shift));
    transforms.emplace_back(ctx, static_cast<vector3_t>(center - shift));

    // Build the portal surfaces
    dindex trf_idx{transforms.size(ctx) - 2};
    surfaces.push_back(
        {trf_idx, mask_link, material_link, volume_idx, surface_id::e_portal},
        invalid_src_link);

    surfaces.push_back(
        {++trf_idx, mask_link, material_link, volume_idx, surface_id::e_portal},
        invalid_src_link);

    //
    // ... x-z plane
    //
    mask_link.shift(1u);
    masks.template emplace_back<rectangle_id>(empty_context{}, volume_link, h_x,
                                              h_z);

    // Rotate by 90deg around x-axis, plus shift in y
    shift = {static_cast<scalar_type>(0),
             detail::is_invalid_value(h_y) ? max_shift : h_y,
             static_cast<scalar_type>(0)};
    vector3_t new_x{1.f, 0.f, 0.f};
    vector3_t new_z{0.f, -1.f, 0.f};
    transforms.emplace_back(ctx, static_cast<vector3_t>(center + shift), new_z,
                            new_x);
    transforms.emplace_back(ctx, static_cast<vector3_t>(center - shift), new_z,
                            new_x);

    surfaces.push_back(
        {++trf_idx, mask_link, material_link, volume_idx, surface_id::e_portal},
        invalid_src_link);

    surfaces.push_back(
        {++trf_idx, mask_link, material_link, volume_idx, surface_id::e_portal},
        invalid_src_link);

    //
    // ... y-z plane
    //
    mask_link.shift(1u);
    masks.template emplace_back<rectangle_id>(empty_context{}, volume_link, h_z,
                                              h_y);

    // Rotate by 90deg around y-axis, plus shift in x
    shift = {detail::is_invalid_value(h_x) ? max_shift : h_x,
             static_cast<scalar_type>(0), static_cast<scalar_type>(0)};
    new_x = {0.f, 0.f, -1.f};
    new_z = {1.f, 0.f, 0.f};
    transforms.emplace_back(ctx, static_cast<vector3_t>(center + shift), new_z,
                            new_x);
    transforms.emplace_back(ctx, static_cast<vector3_t>(center - shift), new_z,
                            new_x);

    surfaces.push_back(
        {++trf_idx, mask_link, material_link, volume_idx, surface_id::e_portal},
        invalid_src_link);

    surfaces.push_back(
        {++trf_idx, mask_link, material_link, volume_idx, surface_id::e_portal},
        invalid_src_link);

    return {surfaces_offset, static_cast<dindex>(surfaces.size())};
  }

 private:
  /// Portal envelope (min distance between portals and volume surfaces)
  scalar_type m_envelope;
};

}  // namespace detray
