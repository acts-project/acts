// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/line.hpp"
#include "detray/utils/axis_rotation.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/unit_vectors.hpp"

// System include(s)
#include <cassert>
#include <limits>

namespace detray {

/// @brief configuration for the wire layer generator
template <concepts::scalar scalar_t>
struct wire_layer_generator_config {
  /// Half length of the layer (z)
  scalar_t m_half_z{1000.f * unit<scalar_t>::mm};
  /// Cell size/radius
  scalar_t m_cell_size{10.f * unit<scalar_t>::mm};
  /// Stereo angle between wires
  scalar_t m_stereo_angle{50.f * unit<scalar_t>::mrad};
  /// Layer radius
  scalar_t m_layer_radius{500.f * unit<scalar_t>::mm};

  /// Setters
  /// @{
  constexpr wire_layer_generator_config &half_length(const scalar_t hz) {
    m_half_z = hz;
    return *this;
  }
  constexpr wire_layer_generator_config &cell_size(const scalar_t c) {
    m_cell_size = c;
    return *this;
  }
  constexpr wire_layer_generator_config &stereo_angle(const scalar_t angle) {
    m_stereo_angle = angle;
    return *this;
  }
  constexpr wire_layer_generator_config &inner_layer_radius(const scalar_t r) {
    m_layer_radius = r;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  constexpr scalar_t half_length() const { return m_half_z; }
  constexpr scalar_t cell_size() const { return m_cell_size; }
  constexpr scalar_t stereo_angle() const { return m_stereo_angle; }
  constexpr scalar_t inner_layer_radius() const { return m_layer_radius; }
  /// @}
};

/// @brief Generates a number of drift cells or straw tubes for a barrel layer
///
/// @tparam detector_t the type of detector the layer should be added to
/// @tparam mask_shape_t the shape for the drift cells or straw tubes
template <typename detector_t, typename mask_shape_t = line_square>
class wire_layer_generator final
    : public surface_factory_interface<detector_t> {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using transform3_t = dtransform3D<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;

 public:
  /// Build a wire chamber/straw tube layer according to the parameters given
  /// in @param cfg
  DETRAY_HOST
  explicit wire_layer_generator(
      const wire_layer_generator_config<scalar_t> &cfg)
      : m_cfg{cfg},
        m_layer_central_rad{m_cfg.inner_layer_radius() + m_cfg.cell_size()} {}

  /// @returns the number of wires this factory will produce
  DETRAY_HOST
  auto size() const -> dindex override {
    const scalar_t circumference{2.f * constant<scalar_t>::pi *
                                 m_layer_central_rad};
    return static_cast<dindex>(
        std::floor(circumference / 2.f * m_cfg.cell_size()));
  }

  /// This is a surface generator, no external surface data needed
  /// @{
  DETRAY_HOST
  void clear() override { /*Do nothing*/ };
  DETRAY_HOST
  void push_back(surface_data<detector_t> &&) override { /*Do nothing*/ }
  DETRAY_HOST
  auto push_back(std::vector<surface_data<detector_t>> &&)
      -> void override { /*Do nothing*/ }
  /// @}

  /// Create a wire chamber barrel layer.
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
    DETRAY_VERBOSE_HOST("Generate modules for wire chamber layer...");
    DETRAY_VERBOSE_HOST("-> Generate " << size() << " surfaces");

    using surface_t = typename detector_t::surface_type;
    using nav_link_t = typename surface_t::navigation_link;
    using mask_link_t = typename surface_t::mask_link;
    using material_link_t = typename surface_t::material_link;

    // Volume that this layer will be added to
    auto volume_idx{volume.index()};

    const dindex surfaces_offset{static_cast<dindex>(surfaces.size())};
    constexpr auto invalid_src_link{detail::invalid_value<std::uint64_t>()};

    // The type id of the surface mask shape (drift cell or straw tube)
    constexpr auto mask_id{
        types::id<typename detector_t::masks, mask<mask_shape_t, algebra_t>>};
    // Modules link back to mother volume in navigation
    const auto mask_volume_link{static_cast<nav_link_t>(volume_idx)};

    // Cell center positions
    detray::dvector<point3_t> cell_centers{};
    // Distance between wires in phi
    const scalar_t delta = 2.f * m_cfg.cell_size() / m_layer_central_rad;
    scalar_t phi{0.f};
    while (phi <= 2.f * constant<scalar_t>::pi) {
      const scalar_t x = m_layer_central_rad * math::cos(phi);
      const scalar_t y = m_layer_central_rad * math::sin(phi);
      const scalar_t z = 0.f;

      cell_centers.push_back({x, y, z});
      phi += delta;
    }

    // Generate a wire surface at every position
    for (const point3_t &cell_center : cell_centers) {
      // Mask link for the cell
      mask_link_t mask_link{mask_id, {masks.template size<mask_id>(), 1u}};
      // The material will be added in a later step
      material_link_t material_link{surface_t::material_id::e_none,
                                    dindex_invalid};
      // Link the placement transform to the surface
      const auto trf_index = transforms.size(ctx);

      surfaces.push_back({trf_index, mask_link, material_link, volume_idx,
                          surface_id::e_sensitive},
                         invalid_src_link);

      // The wire bounds
      masks.template emplace_back<mask_id>(empty_context{}, mask_volume_link,
                                           m_cfg.cell_size(),
                                           m_cfg.half_length());

      // Build the transform
      vector3_t z_axis{0.f, 0.f, 1.f};
      const vector3_t r_axis = vector::normalize(cell_center);
      z_axis = axis_rotation<algebra_t>(r_axis, m_cfg.stereo_angle())(z_axis);
      const vector3_t x_axis =
          unit_vectors<vector3_t>().make_curvilinear_unit_u(z_axis);

      transforms.emplace_back(ctx, cell_center, z_axis, x_axis);
    }

    return {surfaces_offset, static_cast<dindex>(surfaces.size())};
  }

 private:
  /// The generator configuration
  wire_layer_generator_config<scalar_t> m_cfg{};
  /// Central radius of the layer under construction
  /// TODO: Add envelope to avoid overlaps for wire cells
  scalar_t m_layer_central_rad{};
};

}  // namespace detray
