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
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <algorithm>
#include <cassert>
#include <limits>

namespace detray {

/// @brief configuration for the barrel module generator
///
/// The default values correspond to the toy detector inner pixel layer
template <concepts::scalar scalar_t>
struct barrel_generator_config {
  /// Half length of the barrel (z)
  scalar_t m_half_z{500.f * unit<scalar_t>::mm};
  /// Boundary values for the module masks
  std::vector<scalar_t> m_mask_values{8.4f * unit<scalar_t>::mm,
                                      36.f * unit<scalar_t>::mm};
  /// Phi tilt for modules
  scalar_t m_tilt_phi{0.14f};
  /// Layer radius
  scalar_t m_radius{32.f * unit<scalar_t>::mm};
  /// Radial module stagger
  scalar_t m_radial_stagger{0.5f * unit<scalar_t>::mm};
  /// Module overlap in z
  scalar_t m_z_overlap{2.f * unit<scalar_t>::mm};
  /// Number of modules in phi and z
  std::pair<unsigned int, unsigned int> m_binning = {16u, 14u};

  /// Setters
  /// @{
  constexpr barrel_generator_config &half_length(const scalar_t hz) {
    m_half_z = hz;
    return *this;
  }
  barrel_generator_config &module_bounds(const std::vector<scalar_t> &bnds) {
    m_mask_values.clear();
    std::ranges::copy(bnds, std::back_inserter(m_mask_values));
    return *this;
  }
  constexpr barrel_generator_config &tilt_phi(const scalar_t tilt) {
    m_tilt_phi = tilt;
    return *this;
  }
  constexpr barrel_generator_config &radius(const scalar_t r) {
    m_radius = r;
    return *this;
  }
  constexpr barrel_generator_config &radial_stagger(const scalar_t r_st) {
    m_radial_stagger = r_st;
    return *this;
  }
  constexpr barrel_generator_config &z_overlap(const scalar_t zo) {
    m_z_overlap = zo;
    return *this;
  }
  constexpr barrel_generator_config &binning(const unsigned int n_phi,
                                             const unsigned int n_z) {
    m_binning = {n_phi, n_z};
    return *this;
  }
  constexpr barrel_generator_config &binning(
      const std::pair<unsigned int, const unsigned int> &binning) {
    m_binning = binning;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  constexpr scalar_t half_length() const { return m_half_z; }
  const std::vector<scalar_t> &module_bounds() const { return m_mask_values; }
  constexpr scalar_t tilt_phi() const { return m_tilt_phi; }
  constexpr scalar_t radius() const { return m_radius; }
  constexpr scalar_t radial_stagger() const { return m_radial_stagger; }
  constexpr scalar_t z_overlap() const { return m_z_overlap; }
  constexpr const auto &binning() const { return m_binning; }
  /// @}
};

/// @brief Generates a number of surfaces in a barrel shape
///
/// @tparam detector_t the type of detector the layer should be added to
template <typename detector_t, typename mask_shape_t = rectangle2D>
class barrel_generator final : public surface_factory_interface<detector_t> {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using transform3_t = dtransform3D<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;

 public:
  /// Build a barrel layer according to the parameters given in @param cfg
  DETRAY_HOST
  explicit barrel_generator(const barrel_generator_config<scalar_t> &cfg)
      : m_cfg{cfg} {}

  /// @returns the number of surfaces this factory will produce
  DETRAY_HOST
  auto size() const -> dindex override {
    return static_cast<dindex>(m_cfg.binning().first * m_cfg.binning().second);
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

  /// Create a pixel tracker barrel layer.
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
    DETRAY_VERBOSE_HOST("Generate silicon tracker barrel modules...");
    DETRAY_VERBOSE_HOST("-> Generate " << size() << " surfaces");

    using surface_t = typename detector_t::surface_type;
    using nav_link_t = typename surface_t::navigation_link;
    using mask_link_t = typename surface_t::mask_link;
    using material_link_t = typename surface_t::material_link;

    const dindex surfaces_offset{static_cast<dindex>(surfaces.size())};
    constexpr auto invalid_src_link{detail::invalid_value<std::uint64_t>()};

    // The type id of the surface mask shape
    constexpr auto mask_id{
        types::id<typename detector_t::masks, mask<mask_shape_t, algebra_t>>};

    // The material will be added in a later step
    constexpr auto no_material{surface_t::material_id::e_none};

    auto volume_idx{volume.index()};
    // Modules link back to mother volume in navigation
    const auto mask_volume_link{static_cast<nav_link_t>(volume_idx)};

    // surface grid bins
    const unsigned int n_phi_bins{m_cfg.binning().first};
    const unsigned int n_z_bins{m_cfg.binning().second};

    // module positions
    std::vector<point3_t> mod_centers;
    mod_centers.reserve(n_phi_bins * n_z_bins);

    // prep work
    const scalar_t phi_step{2.f * constant<scalar_t>::pi /
                            static_cast<scalar_t>(n_phi_bins)};
    const scalar_t min_phi{-constant<scalar_t>::pi + 0.5f * phi_step};

    // @TODO: Only work for rectangles
    const scalar_t z_start{
        -0.5f * static_cast<scalar_t>(n_z_bins - 1u) *
        (2.f * m_cfg.module_bounds().at(1) - m_cfg.z_overlap())};
    const scalar_t z_step{(math::fabs(z_start) - z_start) /
                          static_cast<scalar_t>(n_z_bins - 1)};

    // loop over the z bins
    for (unsigned int z_bin = 0u; z_bin < n_z_bins; ++z_bin) {
      // prepare z and r
      const scalar_t mod_z{z_start + static_cast<scalar_t>(z_bin) * z_step};
      const scalar_t mod_r{
          (z_bin % 2u) != 0u ? m_cfg.radius() - 0.5f * m_cfg.radial_stagger()
                             : m_cfg.radius() + 0.5f * m_cfg.radial_stagger()};

      for (unsigned int phiBin = 0u; phiBin < n_phi_bins; ++phiBin) {
        // calculate the current phi value
        const scalar_t mod_phi{min_phi +
                               static_cast<scalar_t>(phiBin) * phi_step};
        mod_centers.push_back(point3_t{mod_r * math::cos(mod_phi),
                                       mod_r * math::sin(mod_phi), mod_z});
      }
    }

    // Create geometry data
    for (auto &mod_center : mod_centers) {
      // Surfaces with the linking into the local containers
      mask_link_t mask_link = {mask_id, {masks.template size<mask_id>(), 1u}};
      material_link_t material_link{no_material, dindex_invalid};
      const auto trf_index = transforms.size(ctx);

      surfaces.push_back({trf_index, mask_link, material_link, volume_idx,
                          surface_id::e_sensitive},
                         invalid_src_link);

      // Build the transform
      // The local phi
      const scalar_t mod_phi{vector::phi(mod_center)};
      // Local z axis is the normal vector
      const scalar_t tilt_phi{m_cfg.tilt_phi()};
      const vector3_t mod_local_z{math::cos(mod_phi + tilt_phi),
                                  math::sin(mod_phi + tilt_phi), scalar_t(0)};
      // Local x axis the normal to local y,z
      const vector3_t mod_local_x{-math::sin(mod_phi + tilt_phi),
                                  math::cos(mod_phi + tilt_phi), scalar_t(0)};

      // Create the module transform
      transforms.emplace_back(ctx, mod_center, mod_local_z, mod_local_x);
    }

    // Add the mask
    masks.template emplace_back<mask_id>(empty_context{}, m_cfg.module_bounds(),
                                         mask_volume_link);

    return {surfaces_offset, static_cast<dindex>(surfaces.size())};
  }

 private:
  /// The generator configuration
  barrel_generator_config<scalar_t> m_cfg{};
};

}  // namespace detray
