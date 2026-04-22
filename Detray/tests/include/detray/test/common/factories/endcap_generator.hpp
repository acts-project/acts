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
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/material/material.hpp"
#include "detray/utils/logging.hpp"

// System include(s)
#include <cassert>
#include <iostream>
#include <limits>

namespace detray {

/// @brief configuration for the barrel module generator
///
/// The default values correspond to the toy detector inner pixel layer
template <concepts::scalar scalar_t>
struct endcap_generator_config {
  /// Construct positive or negative endcap
  int m_side{-1};
  /// Center position (|z|)
  scalar_t m_center_z{0.f};
  /// Inner layer radius
  scalar_t m_inner_radius{25.f * unit<scalar_t>::mm};
  /// Outer layer radius
  scalar_t m_outer_radius{180.f * unit<scalar_t>::mm};
  /// Boundary values for the module masks (per ring)
  std::vector<std::vector<scalar_t>> m_mask_values{
      {3.f * unit<scalar_t>::mm, 9.5f * unit<scalar_t>::mm,
       39.f * unit<scalar_t>::mm},
      {6.f * unit<scalar_t>::mm, 10.f * unit<scalar_t>::mm,
       39.f * unit<scalar_t>::mm}};
  /// Stagger between the two rings
  scalar_t m_ring_stagger{2.f * unit<scalar_t>::mm};
  /// Stagger in phi (per ring)
  std::vector<scalar_t> m_phi_stagger = {4.f * unit<scalar_t>::mm,
                                         4.f * unit<scalar_t>::mm};
  /// Substagger in phi (per ring)
  std::vector<scalar_t> m_phi_sub_stagger = {0.5f * unit<scalar_t>::mm,
                                             0.5f * unit<scalar_t>::mm};
  /// Module tilt (per ring)
  std::vector<scalar_t> m_tilt = {0.f, 0.f};
  /// Number of modules in phi (per ring)
  std::vector<unsigned int> m_binning = {40u, 68u};

  /// Setters
  /// @{
  constexpr endcap_generator_config &side(const int s) {
    m_side = s;
    return *this;
  }
  constexpr endcap_generator_config &center(const scalar_t z) {
    m_center_z = std::abs(z);
    return *this;
  }
  constexpr endcap_generator_config &inner_radius(const scalar_t r) {
    m_inner_radius = r;
    return *this;
  }
  constexpr endcap_generator_config &outer_radius(const scalar_t r) {
    m_outer_radius = r;
    return *this;
  }
  endcap_generator_config &module_bounds(
      const std::vector<std::vector<scalar_t>> &bnds) {
    m_mask_values.clear();
    m_mask_values.push_back(bnds.at(0));
    m_mask_values.push_back(bnds.at(1));
    return *this;
  }
  constexpr endcap_generator_config &ring_stagger(const scalar_t rs) {
    m_ring_stagger = rs;
    return *this;
  }
  constexpr endcap_generator_config &phi_stagger(
      const std::vector<scalar_t> &ps) {
    m_phi_stagger = ps;
    return *this;
  }
  constexpr endcap_generator_config &phi_sub_stagger(
      const std::vector<scalar_t> &ps) {
    m_phi_sub_stagger = ps;
    return *this;
  }
  constexpr endcap_generator_config &module_tilt(
      const std::vector<scalar_t> &t) {
    m_tilt = t;
    return *this;
  }
  constexpr endcap_generator_config &binning(
      const std::vector<unsigned int> &b) {
    m_binning = b;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  constexpr int side() const { return m_side; }
  constexpr scalar_t center() const { return m_center_z; }
  constexpr scalar_t inner_radius() const { return m_inner_radius; }
  constexpr scalar_t outer_radius() const { return m_outer_radius; }
  const auto &module_bounds() const { return m_mask_values; }
  constexpr scalar_t ring_stagger() const { return m_ring_stagger; }
  constexpr const auto &phi_stagger() const { return m_phi_stagger; }
  constexpr const auto &phi_sub_stagger() const { return m_phi_sub_stagger; }
  constexpr const auto &module_tilt() const { return m_tilt; }
  constexpr const auto &binning() const { return m_binning; }
  /// @}
};

/// @brief Generates surfaces for an endcap layer consisting of two rings
///
/// @tparam detector_t the type of detector the layer should be added to
template <typename detector_t, typename mask_shape_t = trapezoid2D>
class endcap_generator final : public surface_factory_interface<detector_t> {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using transform3_t = dtransform3D<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;

 public:
  /// Build an endcap layer according to the parameters given in @param cfg
  DETRAY_HOST
  explicit endcap_generator(const endcap_generator_config<scalar_t> &cfg)
      : m_cfg{cfg} {}

  /// @returns the number of surfaces this factory will produce
  DETRAY_HOST
  auto size() const -> dindex override {
    dindex n_modules{0u};
    for (const unsigned int n_bins : m_cfg.binning()) {
      n_modules += n_bins;
    }
    return n_modules;
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

  /// Create a pixel tracker endcap layer with two rings.
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
    DETRAY_VERBOSE_HOST("Generate silicon tracker endcap modules...");
    DETRAY_VERBOSE_HOST("-> Generate " << size() << " surfaces");

    using surface_t = typename detector_t::surface_type;
    using nav_link_t = typename surface_t::navigation_link;
    using mask_link_t = typename surface_t::mask_link;
    using material_link_t = typename surface_t::material_link;

    auto volume_idx{volume.index()};
    const dindex surfaces_offset{static_cast<dindex>(surfaces.size())};
    constexpr auto invalid_src_link{detail::invalid_value<std::uint64_t>()};

    // The type id of the surface mask shape
    constexpr auto mask_id{
        types::id<typename detector_t::masks, mask<mask_shape_t, algebra_t>>};
    // The material will be added in a later step
    constexpr auto no_material{surface_t::material_id::e_none};
    // Modules link back to mother volume in navigation
    const auto mask_volume_link{static_cast<nav_link_t>(volume_idx)};

    // The radii of the rings
    std::vector<scalar_t> radii{};
    // The radial span of the disc
    scalar_t delta_r{m_cfg.outer_radius() - m_cfg.inner_radius()};

    //
    // Calculate the radial positions of the rings
    //

    // Only one ring
    if (m_cfg.binning().size() == 1u) {
      radii.push_back(0.5f * (m_cfg.inner_radius() + m_cfg.outer_radius()));
    } else {
      // Sum up the total length of the modules along r
      scalar_t tot_length{0.f};
      for (const auto &bounds : m_cfg.module_bounds()) {
        // Add an extra buffer, so that the trapezoid corners don't poke
        // out of the cylinder portals
        tot_length += 2.f * bounds.at(trapezoid2D::e_half_length_2) +
                      1.f * unit<scalar_t>::mm;
      }

      // Calculate the overlap (equal pay)
      scalar_t r_overlap{(tot_length - delta_r) /
                         static_cast<scalar_t>(m_cfg.module_bounds().size())};

      // Fill the radii and gaps
      scalar_t prev_r{m_cfg.inner_radius() + r_overlap};
      scalar_t prev_hl{0.f};
      scalar_t prev_ol{0.f};

      for (const auto &bounds : m_cfg.module_bounds()) {
        const scalar_t mod_hlength{bounds.at(trapezoid2D::e_half_length_2)};
        // Calculate the radius
        radii.push_back(prev_r + prev_hl - prev_ol + mod_hlength);
        prev_r = radii.back();
        prev_ol = 2.f * r_overlap;
        prev_hl = mod_hlength;
      }
    }

    //
    // Build the modules in every ring
    //
    for (unsigned int ir = 0u; ir < radii.size(); ++ir) {
      // Generate the position in z (observe ring stagger)
      // Convention: inner ring is closer to origin
      const scalar_t center{m_cfg.center()};
      const scalar_t staggered_center{
          (ir % 2u ? center + 0.5f * m_cfg.ring_stagger()
                   : center - 0.5f * m_cfg.ring_stagger())};
      const scalar_t rz{radii.size() == 1u ? center : staggered_center};

      // Generate the ring module positions (observe phi stagger)
      const scalar_t sub_stagger{m_cfg.phi_sub_stagger().size() > 1u
                                     ? m_cfg.phi_sub_stagger().at(ir)
                                     : 0.f};

      std::vector<point3_t> module_positions =
          module_positions_ring(rz, radii.at(ir), m_cfg.phi_stagger().at(ir),
                                sub_stagger, m_cfg.binning().at(ir));

      // Build the modules
      for (const point3_t &mod_position : module_positions) {
        // Module mask
        mask_link_t mask_link{mask_id, {masks.template size<mask_id>(), 1u}};
        material_link_t material_link{no_material, dindex_invalid};

        // Surface descriptor
        surfaces.push_back({transforms.size(ctx), mask_link, material_link,
                            volume_idx, surface_id::e_sensitive},
                           invalid_src_link);

        // Build the module transform

        // The center position of the module
        point3_t mod_center{mod_position};
        mod_center[2] *= static_cast<scalar_t>(m_cfg.side());
        // The rotation matrix of the module
        const scalar_t mod_phi{vector::phi(mod_position)};
        const vector3_t mod_loc_y{math::cos(mod_phi), math::sin(mod_phi),
                                  scalar_t(0)};
        // Take different axis to have the same readout direction
        const vector3_t mod_loc_z{scalar_t(0), scalar_t(0),
                                  static_cast<scalar_t>(m_cfg.side())};
        const vector3_t mod_loc_x{vector::cross(mod_loc_y, mod_loc_z)};

        // Create the module transform object
        transforms.emplace_back(ctx, mod_center, mod_loc_z, mod_loc_x);
      }

      // Build the mask for this ring
      std::vector<scalar_t> mask_values{m_cfg.module_bounds().at(ir)};

      // Precompute trapezoid divisor
      if constexpr (std::is_same_v<mask_shape_t, trapezoid2D>) {
        const scalar_t div{
            1.f / (2.f * mask_values.at(trapezoid2D::e_half_length_2))};

        mask_values.insert(mask_values.begin() + trapezoid2D::e_divisor, div);
      }

      masks.template emplace_back<mask_id>(empty_context{}, mask_values,
                                           mask_volume_link);
    }

    return {surfaces_offset, static_cast<dindex>(surfaces.size())};
  }

 private:
  /// Helper method for positioning of modules in an endcap ring
  ///
  /// @param z is the z position of the ring
  /// @param radius is the ring radius
  /// @param phi_stagger is the radial staggering along phi
  /// @param phi_sub_stagger is the overlap of the modules
  /// @param n_phi_bins is the number of bins in phi
  ///
  /// @return a vector of the module positions in a ring
  auto module_positions_ring(const scalar_t z, const scalar_t radius,
                             const scalar_t phi_stagger,
                             const scalar_t phi_sub_stagger,
                             const unsigned int n_phi_bins) const {
    // create and fill the positions
    std::vector<point3_t> mod_positions;
    mod_positions.reserve(n_phi_bins);

    // prep work
    const scalar_t phi_step{2.f * constant<scalar_t>::pi /
                            static_cast<scalar_t>(n_phi_bins)};
    const scalar_t min_phi{-constant<scalar_t>::pi + 0.5f * phi_step};

    for (unsigned int iphi = 0u; iphi < n_phi_bins; ++iphi) {
      // If we have a phi sub stagger presents
      scalar_t rzs{0.f};
      // Phi stagger affects 0 vs 1, 2 vs 3 ... etc
      // -> only works if it is a %4
      // Phi sub stagger affects 2 vs 4, 1 vs 3 etc.
      if (phi_sub_stagger != 0.f && !(n_phi_bins % 4u)) {
        // switch sides
        if (!(iphi % 4u)) {
          rzs = phi_sub_stagger;
        } else if (!((iphi + 1u) % 4u)) {
          rzs = -phi_sub_stagger;
        }
      }
      // The module phi
      const scalar_t phi{min_phi + static_cast<scalar_t>(iphi) * phi_step};
      // Main z position depending on phi bin
      const scalar_t rz{iphi % 2u ? z - 0.5f * phi_stagger
                                  : z + 0.5f * phi_stagger};
      mod_positions.push_back(
          {radius * math::cos(phi), radius * math::sin(phi), rz + rzs});
    }
    return mod_positions;
  }

  /// The generator configuration
  endcap_generator_config<scalar_t> m_cfg{};
};

}  // namespace detray
