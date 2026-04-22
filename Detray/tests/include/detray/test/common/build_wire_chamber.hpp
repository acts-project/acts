// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/cylinder_portal_generator.hpp"
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/grid_builder.hpp"
#include "detray/builders/homogeneous_material_builder.hpp"
#include "detray/builders/homogeneous_material_generator.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/wire_chamber_metadata.hpp"
#include "detray/geometry/shapes/line.hpp"
#include "detray/material/predefined_materials.hpp"
#include "detray/utils/consistency_checker.hpp"
#include "detray/utils/print_detector.hpp"

// Detray test include(s)
#include "detray/test/common/factories/wire_layer_generator.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

namespace detray {

/// Configuration for building a wire chamber detector
template <concepts::scalar scalar_t, typename wire_shape_t = line_square>
struct wire_chamber_config {
  /// Number of layers
  unsigned int m_n_layers{10u};
  /// The inner radius of the first layer
  scalar_t m_first_layer_inner_rad{500.f * unit<scalar_t>::mm};
  /// Half z of cylinder chamber
  scalar_t m_half_z{1000.f * unit<scalar_t>::mm};
  /// Radius of the material rods
  scalar_t m_mat_radius{15.f * unit<scalar_t>::um};
  /// Type of material for the material rods
  material<scalar_t> m_wire_mat{tungsten<scalar_t>()};
  /// Config for the wire generation (barrel)
  wire_layer_generator_config<scalar_t> m_wire_factory_cfg{};
  /// Configuration for the homogeneous material generator
  hom_material_config<scalar_t> m_material_config{};
  /// Do a full detector consistency check after building
  bool m_do_check{true};

  /// Setters
  /// @{
  constexpr wire_chamber_config &n_layers(const unsigned int n) {
    m_n_layers = n;
    return *this;
  }
  constexpr wire_chamber_config &first_layer_inner_radius(const scalar_t r) {
    m_first_layer_inner_rad = r;
    return *this;
  }
  constexpr wire_chamber_config &half_z(const scalar_t hz) {
    m_half_z = hz;
    return *this;
  }
  constexpr wire_chamber_config &cell_size(const scalar_t c) {
    m_wire_factory_cfg.cell_size(c);
    return *this;
  }
  constexpr wire_chamber_config &stereo_angle(const scalar_t s) {
    m_wire_factory_cfg.stereo_angle(s);
    return *this;
  }
  constexpr wire_chamber_config &mat_radius(const scalar_t r) {
    m_mat_radius = r;
    return *this;
  }
  constexpr wire_chamber_config &wire_material(const material<scalar_t> &m) {
    m_wire_mat = m;
    return *this;
  }
  constexpr wire_chamber_config &do_check(const bool check) {
    m_do_check = check;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  constexpr unsigned int n_layers() const { return m_n_layers; }
  constexpr scalar_t first_layer_inner_radius() const {
    return m_first_layer_inner_rad;
  }
  constexpr scalar_t half_z() const { return m_half_z; }
  constexpr scalar_t cell_size() const {
    return m_wire_factory_cfg.cell_size();
  }
  constexpr scalar_t stereo_angle() const {
    return m_wire_factory_cfg.stereo_angle();
  }
  constexpr scalar_t mat_radius() const { return m_mat_radius; }
  constexpr const material<scalar_t> &wire_material() const {
    return m_wire_mat;
  }
  constexpr wire_layer_generator_config<scalar_t> &layer_config() {
    return m_wire_factory_cfg;
  }
  constexpr const wire_layer_generator_config<scalar_t> &layer_config() const {
    return m_wire_factory_cfg;
  }
  constexpr auto &material_config() { return m_material_config; }
  constexpr const auto &material_config() const { return m_material_config; }
  constexpr bool do_check() const { return m_do_check; }
  /// @}

 private:
  /// Print the wire chamber configuration
  friend inline std::ostream &operator<<(std::ostream &out,
                                         const wire_chamber_config &cfg) {
    out << "\nWire Chamber\n"
        << "----------------------------\n"
        << "  No. layers            : " << cfg.n_layers() << "\n"
        << "  First layer inner rad.: " << cfg.first_layer_inner_radius()
        << " [mm]\n"
        << "  Half length z         : " << cfg.half_z() << " [mm]\n";

    if constexpr (std::same_as<wire_shape_t, line_square>) {
      out << "  Shape                 : wire cell\n";
    } else {
      out << "  Shape                 : straw tube\n";
    }

    out << "  Cell size             : " << cfg.cell_size() << " [mm]\n"
        << "  Stereo angle          : " << cfg.stereo_angle() << " [rad]\n"
        << "  Wire material         : " << cfg.wire_material() << "\n"
        << "  Material rad.         : " << cfg.mat_radius() << " [mm]\n";

    return out;
  }

};  // wire chamber config

template <concepts::algebra algebra_t, typename wire_shape_t>
inline auto build_wire_chamber(
    vecmem::memory_resource &resource,
    wire_chamber_config<dscalar<algebra_t>, wire_shape_t> &cfg) {
  using builder_t =
      detector_builder<wire_chamber_metadata<algebra_t>, volume_builder>;
  using detector_t = typename builder_t::detector_type;
  using scalar_t = dscalar<typename detector_t::algebra_type>;

  // Wire chamber detector builder
  builder_t det_builder;
  det_builder.set_name("wire_chamber");

  // Geometry context object
  typename detector_t::geometry_context gctx{};

  // Navigation link when leaving the detector
  using nav_link_t = typename detector_t::surface_type::navigation_link;
  constexpr auto leaving_world{detail::invalid_value<nav_link_t>()};
  const scalar_t inner_rad{cfg.first_layer_inner_radius()};
  const scalar_t cell_size{cfg.cell_size()};

  // Prepare grid building
  constexpr auto grid_id{
      detector_t::accel::id::e_surface_concentric_cylinder2D_grid};
  using cyl_grid_t = types::get<typename detector_t::accel, grid_id>;
  using loc_bin_idx_t = typename cyl_grid_t::loc_bin_index;
  static_assert(cyl_grid_t::dim == 2);
  using grid_builder_t =
      grid_builder<detector_t, cyl_grid_t, detray::fill_by_pos>;

  // Binning of the grid
  axis::multi_bin_range<cyl_grid_t::dim> bin_range{};
  // Min, max bin indices per axis
  bin_range[static_cast<std::size_t>(axis::label::e_rphi)] = {0, 100};
  bin_range[static_cast<std::size_t>(axis::label::e_cyl_z)] = {0, 1};

  // Spans of the grid axes
  std::vector<scalar_t> sf_grid_spans{-constant<scalar_t>::pi,
                                      constant<scalar_t>::pi, -cfg.half_z(),
                                      cfg.half_z()};

  constexpr unsigned int bin_capacity{3u};
  std::vector<std::pair<loc_bin_idx_t, dindex>> capacities{};
  capacities.reserve(
      static_cast<std::size_t>(bin_range[0][1] * bin_range[1][1]));

  //
  // Build empty inner volume, where silicon subdetectors would sit
  //
  auto inner_v_builder = det_builder.new_volume(volume_id::e_cylinder);
  inner_v_builder->add_volume_placement(/*identity*/);
  const dindex inner_vol_idx{inner_v_builder->vol_index()};
  inner_v_builder->set_name("inner_vol_" + std::to_string(inner_vol_idx));

  // Configure the portal factory
  // TODO: Add material maps that model the silicon detector budget
  cylinder_portal_config<scalar_t> inner_pt_cfg{};
  inner_pt_cfg.do_autofit(false)
      .fixed_half_length(cfg.half_z())
      .fixed_inner_radius(0.f)
      .fixed_outer_radius(inner_rad)
      // No inner subdetectors present -> don't build inner portal
      .build_inner(false)
      .link_north(inner_vol_idx + 1u)
      .link_south(leaving_world)
      .link_east(leaving_world)
      .link_west(leaving_world);

  auto inner_pt_factory =
      std::make_shared<cylinder_portal_generator<detector_t>>(inner_pt_cfg);
  inner_v_builder->add_surfaces(inner_pt_factory, gctx);

  //
  // Build layer volumes
  //

  // Configure the homogeneous material generation
  auto &mat_cfg = cfg.material_config();
  // Only build material on sensitive surfaces
  constexpr auto vac{vacuum<scalar_t>{}};
  mat_cfg.portal_material(vac).passive_material(vac);
  mat_cfg.sensitive_material(tungsten<scalar_t>{}).thickness(cfg.mat_radius());

  for (unsigned int i_lay = 0; i_lay < cfg.n_layers(); i_lay++) {
    // New volume for layer
    auto v_builder = det_builder.new_volume(volume_id::e_cylinder);
    auto vm_builder =
        det_builder.template decorate<homogeneous_material_builder<detector_t>>(
            v_builder);

    const dindex vol_idx{vm_builder->vol_index()};
    vm_builder->set_name("layer_vol_" + std::to_string(vol_idx));

    // The barrel volumes are centered at the origin
    vm_builder->add_volume_placement(/*identity*/);

    // The maximal inner and outer radius of the volume
    const scalar_t inner_layer_rad =
        inner_rad + static_cast<scalar_t>(i_lay) * 2.f * cell_size;
    const scalar_t outer_layer_rad =
        inner_rad + static_cast<scalar_t>(i_lay + 1u) * 2.f * cell_size;

    // Configure the wire layer factory for this layer
    auto &layer_cfg = cfg.layer_config();
    layer_cfg.inner_layer_radius(inner_layer_rad).half_length(cfg.half_z());

    const scalar_t sign = (i_lay % 2 == 0) ? 1 : -1;
    layer_cfg.stereo_angle(sign * math::fabs(layer_cfg.stereo_angle()));

    // Configure the portal factory
    cylinder_portal_config<scalar_t> layer_portal_cfg{};
    // Limit to maximum valid link
    auto link_north{i_lay == cfg.n_layers() - 1u ? leaving_world
                                                 : vol_idx + 1u};

    layer_portal_cfg.do_autofit(false)
        .fixed_half_length(cfg.half_z())
        .fixed_inner_radius(inner_layer_rad)
        .fixed_outer_radius(outer_layer_rad)
        // Link the volume portals to its neighbors
        .link_north(link_north)
        .link_south(vol_idx - 1u)
        .link_east(leaving_world)
        .link_west(leaving_world);

    // Register sensitive wires together with material
    auto wire_factory =
        std::make_unique<wire_layer_generator<detector_t, wire_shape_t>>(
            layer_cfg);
    auto wire_mat_factory =
        std::make_shared<homogeneous_material_generator<detector_t>>(
            std::move(wire_factory), mat_cfg);

    // Register portals
    auto portal_mat_factory =
        std::make_shared<cylinder_portal_generator<detector_t>>(
            layer_portal_cfg);

    // Add surfaces and material to volume builder
    vm_builder->add_surfaces(portal_mat_factory);
    vm_builder->add_surfaces(wire_mat_factory, gctx);

    // Add a cylinder grid to every barrel layer
    auto vgr_builder =
        det_builder.template decorate<grid_builder_t>(vm_builder);

    // Determine bin capacities
    capacities.clear();
    auto bin_indexer2D = axis::detail::get_bin_indexer(
        bin_range, std::make_integer_sequence<std::size_t, cyl_grid_t::dim>{});
    for (const auto [bin_idx0, bin_idx1] : bin_indexer2D) {
      // @Todo: fine-tune capacity
      loc_bin_idx_t loc_bin{static_cast<unsigned int>(bin_idx0),
                            static_cast<unsigned int>(bin_idx1)};
      capacities.emplace_back(loc_bin, bin_capacity);
    }

    vgr_builder->set_type(detector_t::geo_obj_ids::e_sensitive);
    vgr_builder->init_grid(sf_grid_spans,
                           {static_cast<std::size_t>(bin_range[0][1]),
                            static_cast<std::size_t>(bin_range[1][1])},
                           capacities);
  }

  // Build and return the detector and fill name map
  typename detector_t::name_map name_map{};
  auto det = det_builder.build(resource, name_map);

  if (cfg.do_check()) {
    const bool verbose_check{false};
    detray::detail::check_consistency(det, verbose_check, name_map);
  }

  DETRAY_DEBUG_HOST("\n" << detray::utils::print_detector(det, name_map));

  return std::make_pair(std::move(det), std::move(name_map));
}

}  // namespace detray
