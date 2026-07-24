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
#include "detray/builders/material_map_builder.hpp"
#include "detray/builders/material_map_generator.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/units.hpp"
#include "detray/detectors/toy_metadata.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/material/mixture.hpp"
#include "detray/material/predefined_materials.hpp"
#include "detray/utils/consistency_checker.hpp"
#include "detray/utils/print_detector.hpp"
#include "detray/utils/ranges.hpp"

// Detray test include(s)
#include "detray/test/common/factories/barrel_generator.hpp"
#include "detray/test/common/factories/endcap_generator.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace detray {

/// Configure the toy detector
template <concepts::scalar scalar_t>
struct toy_det_config {
  /// Default toy detector configuration
  toy_det_config() {
    // Barrel module creator
    m_barrel_factory_cfg.half_length(500.f * unit<scalar_t>::mm)
        .module_bounds({8.4f * unit<scalar_t>::mm, 36.f * unit<scalar_t>::mm})
        .tilt_phi(0.14f /*0.145*/)
        .radial_stagger(0.5f * unit<scalar_t>::mm /*2.f*/)
        .z_overlap(2.f * unit<scalar_t>::mm /*5.f*/);

    // Endcap module creator
    m_endcap_factory_cfg.inner_radius(m_beampipe_volume_radius)
        .outer_radius(m_outer_radius)
        .module_bounds({{3.f * unit<scalar_t>::mm, 9.5f * unit<scalar_t>::mm,
                         39.f * unit<scalar_t>::mm},
                        {6.f * unit<scalar_t>::mm, 10.f * unit<scalar_t>::mm,
                         39.f * unit<scalar_t>::mm}})
        .ring_stagger(2.f * unit<scalar_t>::mm)
        .phi_stagger({4.f * unit<scalar_t>::mm, 4.f * unit<scalar_t>::mm})
        .phi_sub_stagger({0.5f * unit<scalar_t>::mm, 0.5f * unit<scalar_t>::mm})
        .module_tilt({0.f, 0.f})
        .binning({40u, 68u});

    // Configure the material generation
    m_material_config.sensitive_material(silicon_tml<scalar_t>())
        .passive_material(beryllium_tml<scalar_t>())  // < beampipe
        .portal_material(vacuum<scalar_t>())
        .thickness(1.5f * unit<scalar_t>::mm);

    // Configure the material map generation
    m_beampipe_map_cfg.n_bins = {20u, 20u};
    m_beampipe_map_cfg.axis_index = 1u;
    m_beampipe_map_cfg.mapped_material = beryllium_tml<scalar_t>();
    m_beampipe_map_cfg.thickness = 0.8f * unit<scalar_t>::mm;
    // Don't scale the generation of the material thickness
    m_beampipe_map_cfg.scalor = 0.f;
    m_beampipe_map_cfg.mat_generator =
        detray::detail::generate_cyl_mat<scalar_t>;

    m_disc_map_cfg.n_bins = {5u, 20u};
    m_disc_map_cfg.axis_index = 0u;
    m_disc_map_cfg.mapped_material =
        mixture<scalar_t, silicon_tml<scalar_t, std::ratio<9, 10>>,
                aluminium<scalar_t, std::ratio<1, 10>>>{};
    m_disc_map_cfg.thickness = 1.f * unit<scalar_t>::mm;
    m_disc_map_cfg.scalor = 1e-4f;
    m_disc_map_cfg.mat_generator = detray::detail::generate_disc_mat<scalar_t>;

    m_cyl_map_cfg.n_bins = {20u, 20u};
    m_cyl_map_cfg.axis_index = 1u;
    m_cyl_map_cfg.mapped_material =
        mixture<scalar_t, silicon_tml<scalar_t, std::ratio<9, 10>>,
                aluminium<scalar_t, std::ratio<1, 10>>>{};
    m_cyl_map_cfg.thickness = 5.f * unit<scalar_t>::mm;
    m_cyl_map_cfg.scalor = 1e-8f;
    m_cyl_map_cfg.mat_generator = detray::detail::generate_cyl_mat<scalar_t>;
  }

  /// No. of barrel layers the detector should be built with
  unsigned int m_n_brl_layers{4u};
  /// No. of endcap layers (on either side) the detector should be built with
  unsigned int m_n_edc_layers{3u};
  /// Total outer radius of the pixel subdetector
  scalar_t m_outer_radius{180.f * unit<scalar_t>::mm};
  // Radius of the innermost volume that contains the beampipe
  scalar_t m_beampipe_volume_radius{25.f * unit<scalar_t>::mm};
  // Envelope around the modules used by the cylinder portal generator
  scalar_t m_portal_envelope{2.f * unit<scalar_t>::mm};
  /// Configuration for the homogeneous material generator
  hom_material_config<scalar_t> m_material_config{};
  /// Build spatial grid acceleration structures (otherwise brute force search)
  bool m_use_grids{true};
  /// Put homogeneneous material on the sensitive surfaces
  bool m_use_homogeneous_material{true};
  /// Put material maps on portals or use homogeneous material on modules
  bool m_use_material_maps{false};
  /// Configuration for the material map generator (beampipe)
  typename material_map_config<scalar_t>::map_config m_beampipe_map_cfg{};
  /// Configuration for the material map generator (disc)
  typename material_map_config<scalar_t>::map_config m_disc_map_cfg{};
  /// Configuration for the material map generator (cylinder)
  typename material_map_config<scalar_t>::map_config m_cyl_map_cfg{};
  /// Thickness of the beampipe material
  scalar_t m_beampipe_mat_thickness{0.8f * unit<scalar_t>::mm};
  /// Thickness of the material slabs in the homogeneous material description
  scalar_t m_module_mat_thickness{1.5f * unit<scalar_t>::mm};
  /// Radii at which to place the barrel module layers (including beampipe)
  std::vector<scalar_t> m_barrel_layer_radii = {
      19.f * unit<scalar_t>::mm, 32.f * unit<scalar_t>::mm,
      72.f * unit<scalar_t>::mm, 116.f * unit<scalar_t>::mm,
      172.f * unit<scalar_t>::mm};
  /// Number of modules in phi and z for the barrel
  std::vector<std::pair<unsigned int, unsigned int>> m_barrel_binning = {
      {0u, 0u}, {16u, 14u}, {32u, 14u}, {52u, 14u}, {78u, 14u}};
  /// Positions at which to place the endcap module layers on either side
  std::vector<scalar_t> m_endcap_layer_positions = {
      600.f * unit<scalar_t>::mm,  700.f * unit<scalar_t>::mm,
      820.f * unit<scalar_t>::mm,  960.f * unit<scalar_t>::mm,
      1100.f * unit<scalar_t>::mm, 1300.f * unit<scalar_t>::mm,
      1500.f * unit<scalar_t>::mm};
  /// Config for the module generation (barrel)
  barrel_generator_config<scalar_t> m_barrel_factory_cfg{};
  /// Config for the module generation (endcaps)
  endcap_generator_config<scalar_t> m_endcap_factory_cfg{};
  /// Run detector consistency check after reading
  bool m_do_check{true};

  /// Setters
  /// @{
  constexpr toy_det_config &n_brl_layers(const unsigned int n) {
    m_n_brl_layers = n;
    return *this;
  }
  constexpr toy_det_config &n_edc_layers(const unsigned int n) {
    m_n_edc_layers = n;
    return *this;
  }
  constexpr toy_det_config &envelope(const scalar_t env) {
    m_portal_envelope = env;
    return *this;
  }
  constexpr toy_det_config &use_grids(const bool b) {
    m_use_grids = b;
    return *this;
  }
  constexpr toy_det_config &use_homogeneous_material(const bool b) {
    m_use_homogeneous_material = b;
    return *this;
  }
  constexpr toy_det_config &use_material_maps(const bool b) {
    m_use_material_maps = b;
    return *this;
  }
  constexpr toy_det_config &cyl_map_bins(const std::size_t n_phi,
                                         const std::size_t n_z) {
    m_cyl_map_cfg.n_bins = {n_phi, n_z};
    return *this;
  }
  constexpr toy_det_config &disc_map_bins(const std::size_t n_r,
                                          const std::size_t n_phi) {
    m_disc_map_cfg.n_bins = {n_r, n_phi};
    return *this;
  }
  constexpr toy_det_config &material_map_min_thickness(const scalar_t t) {
    assert(t > 0.f);
    m_cyl_map_cfg.thickness = t;
    m_disc_map_cfg.thickness = t;
    return *this;
  }
  constexpr toy_det_config &beampipe_mat_thickness(const scalar_t t) {
    assert(t > 0.f);
    m_beampipe_mat_thickness = t;
    return *this;
  }
  constexpr toy_det_config &module_mat_thickness(const scalar_t t) {
    assert(t > 0.f);
    m_module_mat_thickness = t;
    return *this;
  }
  constexpr toy_det_config &mapped_material(const material<scalar_t> &mat) {
    m_cyl_map_cfg.mapped_material = mat;
    m_disc_map_cfg.mapped_material = mat;
    return *this;
  }
  constexpr toy_det_config &do_check(const bool check) {
    m_do_check = check;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  constexpr unsigned int n_brl_layers() const { return m_n_brl_layers; }
  constexpr unsigned int n_edc_layers() const { return m_n_edc_layers; }
  constexpr const auto &outer_radius() const { return m_outer_radius; }
  constexpr scalar_t envelope() const { return m_portal_envelope; }
  constexpr scalar_t beampipe_vol_radius() const {
    return m_beampipe_volume_radius;
  }
  constexpr auto &material_config() { return m_material_config; }
  constexpr const auto &material_config() const { return m_material_config; }
  constexpr bool use_grids() const { return m_use_grids; }
  constexpr bool use_homogeneous_material() const {
    return m_use_homogeneous_material;
  }
  constexpr bool use_material_maps() const { return m_use_material_maps; }
  constexpr auto &beampipe_material_map() { return m_beampipe_map_cfg; }
  constexpr const auto &beampipe_material_map() const {
    return m_beampipe_map_cfg;
  }
  constexpr const auto &cyl_material_map() const { return m_cyl_map_cfg; }
  constexpr auto &cyl_material_map() { return m_cyl_map_cfg; }
  constexpr const auto &disc_material_map() const { return m_disc_map_cfg; }
  constexpr auto &disc_material_map() { return m_disc_map_cfg; }
  constexpr const darray<std::size_t, 2> &cyl_map_bins() const {
    return m_cyl_map_cfg.n_bins;
  }
  constexpr const darray<std::size_t, 2> &disc_map_bins() const {
    return m_disc_map_cfg.n_bins;
  }
  constexpr scalar_t material_map_min_thickness() const {
    assert(m_cyl_map_cfg.thickness == m_disc_map_cfg.thickness);
    return m_cyl_map_cfg.thickness;
  }
  constexpr scalar_t beampipe_mat_thickness() const {
    return m_beampipe_mat_thickness;
  }
  constexpr scalar_t module_mat_thickness() const {
    return m_module_mat_thickness;
  }
  auto barrel_mat_generator() const { return m_cyl_map_cfg.mat_generator; }
  auto edc_mat_generator() const { return m_disc_map_cfg.mat_generator; }
  constexpr material<scalar_t> mapped_material() const {
    assert(m_cyl_map_cfg.mapped_material == m_disc_map_cfg.mapped_material);
    return m_cyl_map_cfg.mapped_material;
  }
  constexpr const auto &barrel_layer_radii() const {
    return m_barrel_layer_radii;
  }
  constexpr const auto &endcap_layer_positions() const {
    return m_endcap_layer_positions;
  }
  constexpr const auto &barrel_layer_binning() const {
    return m_barrel_binning;
  }
  constexpr barrel_generator_config<scalar_t> &barrel_config() {
    return m_barrel_factory_cfg;
  }
  constexpr endcap_generator_config<scalar_t> &endcap_config() {
    return m_endcap_factory_cfg;
  }
  constexpr bool do_check() const { return m_do_check; }
  /// @}

  /// Print the toy detector configuration
  friend std::ostream &operator<<(std::ostream &out,
                                  const toy_det_config &cfg) {
    out << "\nToy Detector\n"
        << "----------------------------\n"
        << "  No. barrel layers     : " << cfg.n_brl_layers() << "\n"
        << "  No. endcap layers     : " << cfg.n_edc_layers() << "\n"
        << "  Portal envelope       : " << cfg.envelope() << " [mm]\n";

    if (cfg.use_material_maps()) {
      const auto &cyl_map_bins = cfg.cyl_map_bins();
      const auto &disc_map_bins = cfg.disc_map_bins();

      out << "  Material maps \n"
          << "    -> cyl. map bins    : (phi: " << cyl_map_bins[0]
          << ", z: " << cyl_map_bins[1] << ")\n"
          << "    -> disc map bins    : (r: " << disc_map_bins[0]
          << ", phi: " << disc_map_bins[1] << ")\n"
          << "    -> cyl. min. thickness: "
          << cfg.cyl_material_map().thickness / detray::unit<scalar_t>::mm
          << " [mm]\n"
          << "    -> disc min. thickness: "
          << cfg.disc_material_map().thickness / detray::unit<scalar_t>::mm
          << " [mm]\n"
          << "    -> Material         : " << cfg.mapped_material() << "\n";
    } else if (cfg.use_homogeneous_material()) {
      out << "  Homogeneous material \n"
          << "    -> Thickness        : "
          << cfg.module_mat_thickness() / detray::unit<scalar_t>::mm
          << " [mm]\n"
          << "    -> Material         : " << silicon_tml<scalar_t>() << "\n";
    } else {
      out << "  No material\n";
    }

    return out;
  }
};

namespace detail {

// Helper type, used to define the r- or z-extent of detector volumes
template <concepts::scalar scalar_t>
struct extent2D {
  scalar_t lower;
  scalar_t upper;
};

/// Helper method to decorate a volume builder with material
///
/// @param cfg config for the toy detector
/// @param det_builder detector builder the volume belongs to
/// @param v_builder the builder of the volume that should be decorated
///
/// @returns the decorated volume builder and surface factory
template <typename detector_builder_t, typename detector_t, typename config_t>
volume_builder_interface<detector_t> *decorate_material(
    config_t &cfg, detector_builder_t &det_builder,
    volume_builder_interface<detector_t> *v_builder) {
  static_assert(
      std::is_same_v<detector_t, typename detector_builder_t::detector_type>,
      "Detector builder and volume builder/surface factory have different "
      "detector type");

  volume_builder_interface<detector_t> *vm_builder{v_builder};

  // Decorate the builder with material
  if (cfg.use_material_maps()) {
    // Build the volume with material maps
    vm_builder =
        det_builder.template decorate<material_map_builder<detector_t>>(
            v_builder);
  } else if (cfg.use_homogeneous_material()) {
    // Build the volume with a homogeneous material description
    vm_builder =
        det_builder.template decorate<homogeneous_material_builder<detector_t>>(
            v_builder);
  }

  if (!vm_builder) {
    throw std::runtime_error("Material decoration failed");
  }

  return vm_builder;
}

/// Helper method to decorate a surface factory with material generators
///
/// @param cfg config for the toy detector
/// @param sf_factory surface factory that should be decorated with material
/// @param is_module_factory should homogeneous material be added to surfaces
///                          of this factory (assumes no material maps are used)
///
/// @returns the decorated volume builder and surface factory
template <typename detector_t>
std::shared_ptr<surface_factory_interface<detector_t>> decorate_material(
    toy_det_config<typename detector_t::scalar_type> &cfg,
    std::unique_ptr<surface_factory_interface<detector_t>> sf_factory,
    bool is_module_factory = false) {
  using scalar_t = dscalar<typename detector_t::algebra_type>;
  using mask_id = typename detector_t::masks::id;
  using material_id = typename detector_t::material::id;

  // Decorate the surfaces with material
  if (cfg.use_material_maps()) {
    // How to generate the specific material map for every surface
    material_map_config<scalar_t> material_map_config{};

    // Add the complete configuration for the beampipe map generation
    cfg.beampipe_material_map().map_id =
        static_cast<dindex>(material_id::e_concentric_cylinder2D_map);
    material_map_config.set_map_config(mask_id::e_concentric_cylinder2D,
                                       surface_id::e_passive,
                                       cfg.beampipe_material_map());

    // Add the complete configuration for disc map generation
    cfg.disc_material_map().map_id =
        static_cast<dindex>(material_id::e_ring2D_map);
    material_map_config.set_map_config(mask_id::e_ring2D, surface_id::e_portal,
                                       cfg.disc_material_map());

    // Add the complete configuration for cylinder map generation
    cfg.cyl_material_map().map_id =
        static_cast<dindex>(material_id::e_concentric_cylinder2D_map);
    material_map_config.set_map_config(mask_id::e_concentric_cylinder2D,
                                       surface_id::e_portal,
                                       cfg.cyl_material_map());

    auto mat_generator = std::make_shared<material_map_generator<detector_t>>(
        std::move(sf_factory), material_map_config);

    return mat_generator;

  } else if (cfg.use_homogeneous_material() && is_module_factory) {
    // How to generate the specific material for every surface
    auto mat_generator =
        std::make_shared<homogeneous_material_generator<detector_t>>(
            std::move(sf_factory), cfg.material_config());

    return mat_generator;
  } else {
    // If no material should be added, return the factory as is
    return sf_factory;
  }
}

/// Add the portals for a cylinder volume from explicit parameters
///
/// @param v_builder the volume builder to add the portals to
/// @param cfg config for the toy detector
/// @param vol_z half length of the cylinder
/// @param h_z half length of the cylinder
/// @param inner_r inner volume radius
/// @param outer_r outer volume radius
/// @param link_north portal volume link of the outer cylinder
/// @param link_south portal volume link of the inner cylinder
/// @param link_east portal volume link of the left disc
/// @param link_west portal volume link of the right disc
/// @param add_material decorate material maps to portals
template <typename detector_t>
void add_cylinder_portals(volume_builder_interface<detector_t> *v_builder,
                          toy_det_config<typename detector_t::scalar_type> &cfg,
                          const typename detector_t::scalar_type lower_z,
                          const typename detector_t::scalar_type upper_z,
                          const typename detector_t::scalar_type inner_r,
                          const typename detector_t::scalar_type outer_r,
                          const dindex link_north, const dindex link_south,
                          const dindex link_east, const dindex link_west) {
  using algebra_t = typename detector_t::algebra_type;
  using transform3_t = dtransform3D<algebra_t>;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  const transform3_t identity{};
  const dindex vol_idx{v_builder->vol_index()};

  scalar_t min_r{math::min(inner_r, outer_r)};
  scalar_t max_r{math::max(inner_r, outer_r)};
  scalar_t min_z{math::min(lower_z, upper_z)};
  scalar_t max_z{math::max(lower_z, upper_z)};

  using factory_interface_t = surface_factory_interface<detector_t>;
  using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
  using disc_factory_t = surface_factory<detector_t, ring2D>;

  std::shared_ptr<factory_interface_t> pt_cyl_factory{nullptr};
  std::shared_ptr<factory_interface_t> pt_disc_factory{nullptr};

  if (cfg.use_material_maps()) {
    pt_cyl_factory =
        decorate_material<detector_t>(cfg, std::make_unique<cyl_factory_t>());
    pt_disc_factory =
        decorate_material<detector_t>(cfg, std::make_unique<disc_factory_t>());
  } else {
    pt_cyl_factory = std::make_shared<cyl_factory_t>();
    pt_disc_factory = std::make_shared<disc_factory_t>();
  }

  // Inner cylinder portal
  pt_cyl_factory->push_back({surface_id::e_portal, identity,
                             static_cast<nav_link_t>(link_south),
                             std::vector<scalar_t>{min_r, min_z, max_z}});
  // Outer cylinder portal
  pt_cyl_factory->push_back({surface_id::e_portal, identity,
                             static_cast<nav_link_t>(link_north),
                             std::vector<scalar_t>{max_r, min_z, max_z}});

  // Left disc portal
  pt_disc_factory->push_back(
      {surface_id::e_portal,
       transform3_t{
           point3_t{static_cast<scalar_t>(0), static_cast<scalar_t>(0), min_z}},
       static_cast<nav_link_t>(link_west),
       std::vector<scalar_t>{min_r, max_r}});
  // Right disc portal
  pt_disc_factory->push_back(
      {surface_id::e_portal,
       transform3_t{
           point3_t{static_cast<scalar_t>(0), static_cast<scalar_t>(0), max_z}},
       static_cast<nav_link_t>(link_east),
       std::vector<scalar_t>{min_r, max_r}});

  v_builder->add_surfaces(pt_cyl_factory);
  v_builder->add_surfaces(pt_disc_factory);

  v_builder->set_name("gap_" + std::to_string(vol_idx));
}

/// Helper method for creating the barrel surface grids.
///
/// @param det_builder detector builder the barrel section should be added to
/// @param cfg config for the toy detector
/// @param vol_index index of the volume to which the grid should be added
template <typename detector_builder_t>
inline void add_cylinder_grid(
    detector_builder_t &det_builder,
    toy_det_config<typename detector_builder_t::detector_type::scalar_type>
        &cfg,
    const dindex vol_index) {
  using detector_t = typename detector_builder_t::detector_type;
  using scalar_t = dscalar<typename detector_t::algebra_type>;

  constexpr auto grid_id =
      detector_t::accel::id::e_surface_concentric_cylinder2D_grid;
  using cyl_grid_t = types::get<typename detector_t::accel, grid_id>;

  using grid_builder_t =
      grid_builder<detector_t, cyl_grid_t, detray::fill_by_pos>;

  const auto &barrel_cfg{cfg.barrel_config()};
  const scalar_t h_z{barrel_cfg.half_length()};

  auto vgr_builder = det_builder.template decorate<grid_builder_t>(vol_index);

  if (!vgr_builder) {
    throw std::runtime_error("Grid decoration failed");
  }

  vgr_builder->set_type(detector_t::geo_obj_ids::e_sensitive);
  vgr_builder->init_grid(
      {-constant<scalar_t>::pi, constant<scalar_t>::pi, -h_z, h_z},
      {barrel_cfg.binning().first, barrel_cfg.binning().second});
}

/// Helper method for creating the endcap surface grids.
///
/// @param det_builder detector builder the barrel section should be added to
/// @param cfg config for the toy detector
/// @param vol_index index of the volume to which the grid should be added
template <typename detector_builder_t>
inline void add_disc_grid(
    detector_builder_t &det_builder,
    toy_det_config<typename detector_builder_t::detector_type::scalar_type>
        &cfg,
    const dindex vol_index) {
  using detector_t = typename detector_builder_t::detector_type;
  using scalar_t = dscalar<typename detector_t::algebra_type>;

  constexpr auto grid_id = detector_t::accel::id::e_surface_ring2D_grid;
  using disc_grid_t = types::get<typename detector_t::accel, grid_id>;

  using grid_builder_t =
      grid_builder<detector_t, disc_grid_t, detray::fill_by_pos>;

  const auto &endcap_cfg{cfg.endcap_config()};
  const scalar_t inner_r{cfg.beampipe_vol_radius()};
  const scalar_t outer_r{cfg.outer_radius()};

  auto vgr_builder = det_builder.template decorate<grid_builder_t>(vol_index);

  if (!vgr_builder) {
    throw std::runtime_error("Grid decoration failed");
  }

  vgr_builder->set_type(detector_t::geo_obj_ids::e_sensitive);
  vgr_builder->init_grid(
      {inner_r, outer_r, -constant<scalar_t>::pi, constant<scalar_t>::pi},
      {endcap_cfg.binning().size(), endcap_cfg.binning().back()});
}

/// Helper method to retrieve the volume extent from a potentially decorated
/// cylinder portal factory
///
/// @param[in] cfg config for the toy detector
/// @param[in] sf_factory (material) factory to get the cylinder extent from
/// @param[out] vol_bounds boundary struct
template <typename detector_t>
inline void get_volume_extent(
    toy_det_config<typename detector_t::scalar_type> &cfg,
    const std::shared_ptr<surface_factory_interface<detector_t>> &sf_factory,
    typename cylinder_portal_generator<detector_t>::boundaries &vol_bounds) {
  // Get the volume extent from the cylinder factory
  const cylinder_portal_generator<detector_t> *cyl_factory{nullptr};
  if (cfg.use_material_maps()) {
    // Retrieve the underlying portal factory
    auto decorator =
        std::dynamic_pointer_cast<const factory_decorator<detector_t>>(
            sf_factory);
    if (!decorator) {
      throw std::bad_cast();
    }
    cyl_factory = dynamic_cast<const cylinder_portal_generator<detector_t> *>(
        decorator->get_factory());
  } else {
    // No material was decorated
    cyl_factory = dynamic_cast<const cylinder_portal_generator<detector_t> *>(
        sf_factory.get());
  }
  if (!cyl_factory) {
    throw std::bad_cast();
  }

  // Set the new current boundaries, to construct the next gap
  vol_bounds = cyl_factory->volume_boundaries();
}

/// Helper method for creating the barrel section.
///
/// @param det_builder detector builder the barrel section should be added to
/// @param gctx geometry context
/// @param cfg config for the toy detector
/// @param beampipe_idx index of the beampipe outermost volume
///
/// @returns the radial extents of the barrel module layers and gap volumes
template <typename detector_builder_t>
inline auto add_barrel_detector(
    detector_builder_t &det_builder,
    typename detector_builder_t::detector_type::geometry_context &gctx,
    toy_det_config<typename detector_builder_t::detector_type::scalar_type>
        &cfg,
    dindex beampipe_idx) {
  using detector_t = typename detector_builder_t::detector_type;
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using transform3_t = dtransform3D<algebra_t>;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  // Register the sizes in z per volume index
  std::vector<std::pair<dindex, extent2D<scalar_t>>> volume_sizes{};

  // Mask volume link for portals that exit the detector
  constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};
  const transform3_t identity{};

  // Links to the connector gaps. Their volume index depends on the
  // number of endcap layer that are being constructed (before and after
  // the barrel volumes are constructed)
  auto link_east{end_of_world};
  auto link_west{end_of_world};

  if (cfg.n_edc_layers() > 0) {
    link_east = static_cast<nav_link_t>(det_builder.n_volumes() +
                                        2u * cfg.n_brl_layers() + 2u);
    link_west = static_cast<nav_link_t>(det_builder.n_volumes() -
                                        2u * cfg.n_edc_layers() + 1u);
  }

  const scalar_t h_z{cfg.barrel_config().half_length()};
  // Set the inner radius of the first gap to the radius of the beampipe vol.
  scalar_t gap_inner_r{cfg.beampipe_vol_radius()};

  typename cylinder_portal_generator<detector_t>::boundaries vol_bounds{};

  // Alternate barrel module layers and gap volumes
  bool is_gap = true;
  for (unsigned int i = 0u; i < 2u * cfg.n_brl_layers(); ++i) {
    // New volume
    auto v_builder = det_builder.new_volume(volume_id::e_cylinder);
    auto vm_builder = decorate_material(cfg, det_builder, v_builder);
    const dindex vol_idx{vm_builder->vol_index()};

    // The barrel volumes are centered at the origin
    vm_builder->add_volume_placement(identity);

    // Every second layer is a gap volume
    is_gap = !is_gap;
    if (is_gap) {
      auto link_north{vol_idx - 1};
      auto link_south{vol_idx - 3};

      // The first time a gap is built, it needs to link to the beampipe
      link_south = (i == 1u) ? beampipe_idx : link_south;

      detail::add_cylinder_portals(vm_builder, cfg, -h_z, h_z, gap_inner_r,
                                   vol_bounds.inner_radius, link_north,
                                   link_south, link_east, link_west);
      volume_sizes.push_back({vol_idx, {gap_inner_r, vol_bounds.inner_radius}});

      // Set the inner gap radius for the next gap volume
      gap_inner_r = vol_bounds.outer_radius;

      vm_builder->set_name("gap_" + std::to_string(vol_idx));

    } else {
      // Limit to maximum valid link
      auto link_north{std::min(
          vol_idx + 3, 2 * cfg.n_edc_layers() + 2 * cfg.n_brl_layers() + 1)};
      auto link_south{vol_idx + 1};

      // Configure the module factory for this layer
      auto &barrel_cfg = cfg.barrel_config();

      const unsigned int j{(i + 2u) / 2u};
      barrel_cfg.binning(cfg.barrel_layer_binning().at(j))
          .radius(cfg.barrel_layer_radii().at(j));

      // Configure the portal factory
      cylinder_portal_config<scalar_t> portal_cfg{};

      portal_cfg.envelope(cfg.envelope())
          .fixed_half_length(h_z)
          // Link the volume portals to its neighbors
          .link_north(link_north)
          .link_south(link_south)
          .link_east(link_east)
          .link_west(link_west);

      // Configure the material
      cfg.material_config().thickness(cfg.module_mat_thickness());

      // Add a layer of module surfaces (may have material)
      auto module_mat_factory = decorate_material<detector_t>(
          cfg,
          std::make_unique<barrel_generator<detector_t, rectangle2D>>(
              barrel_cfg),
          true);

      // Add cylinder and disc portals (may have material, depending on
      // whether material maps are being used or not)
      auto portal_mat_factory = decorate_material<detector_t>(
          cfg,
          std::make_unique<cylinder_portal_generator<detector_t>>(portal_cfg));

      vm_builder->add_surfaces(module_mat_factory, gctx);
      vm_builder->add_surfaces(portal_mat_factory);

      // Fill the volume extent into 'vol_bounds'
      get_volume_extent(cfg, portal_mat_factory, vol_bounds);

      volume_sizes.push_back(
          {vol_idx, {vol_bounds.inner_radius, vol_bounds.outer_radius}});

      vm_builder->set_name("barrel_" + std::to_string(vol_idx));

      if (cfg.use_grids()) {
        // Add a cylinder grid to every barrel module layer
        add_cylinder_grid(det_builder, cfg, vol_idx);
      }
    }
  }

  // Add a final gap volume to get to the full barrel radius
  auto v_builder = det_builder.new_volume(volume_id::e_cylinder);
  const dindex vol_idx{v_builder->vol_index()};
  v_builder->add_volume_placement(identity);

  detail::add_cylinder_portals(
      v_builder, cfg, -h_z, h_z, vol_bounds.outer_radius, cfg.outer_radius(),
      end_of_world, vol_idx - 2u, link_east, link_west);
  volume_sizes.push_back(
      {vol_idx, {vol_bounds.outer_radius, cfg.outer_radius()}});

  v_builder->set_name("gap_" + std::to_string(vol_idx));

  return volume_sizes;
}

/// Helper method for creating one of the two endcaps.
///
/// @param det_builder detector builder the barrel section should be added to
/// @param gctx geometry context
/// @param cfg config for the toy detector
/// @param beampipe_idx index of the beampipe outermost volume
///
/// @returns the z extents of the endcap module layers and gap volumes
template <typename detector_builder_t>
inline auto add_endcap_detector(
    detector_builder_t &det_builder,
    typename detector_builder_t::detector_type::geometry_context &gctx,
    toy_det_config<typename detector_builder_t::detector_type::scalar_type>
        &cfg,
    dindex beampipe_idx) {
  using detector_t = typename detector_builder_t::detector_type;
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  std::vector<std::pair<dindex, extent2D<scalar_t>>> volume_sizes{};

  // Mask volume link for portals that exit the detector
  constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};
  // Disc portal link to the barrel section (via connector gap volume, which
  // is the second volume initialized by this function, but built in a
  // later step, when all necessaty information is available)
  const dindex connector_link{det_builder.n_volumes() + 1u};

  const auto sign{static_cast<scalar_t>(cfg.endcap_config().side())};

  // Mask volume links
  auto link_north{end_of_world};
  auto link_south{beampipe_idx};

  // Set the inner radius of the first gap to the radius of the beampipe vol.
  const scalar_t inner_radius{cfg.beampipe_vol_radius()};
  const scalar_t outer_radius{cfg.outer_radius()};

  scalar_t gap_east_z{sign * cfg.barrel_config().half_length()};

  typename cylinder_portal_generator<detector_t>::boundaries vol_bounds{};

  // Alternate endcap module layers and gap volumes
  bool is_gap = true;
  for (dindex i = 0u; i < 2u * cfg.n_edc_layers(); ++i) {
    // New volume
    auto v_builder = det_builder.new_volume(volume_id::e_cylinder);
    auto vm_builder = decorate_material(cfg, det_builder, v_builder);
    const dindex vol_idx{vm_builder->vol_index()};

    // Don't build the first gap here, as this will be done in a separate
    // step once all portals of the barrel are constructed
    if (i == 1u) {
      const scalar_t gap_west_z{sign * std::min(std::abs(vol_bounds.upper_z),
                                                std::abs(vol_bounds.lower_z))};

      volume_sizes.push_back({vol_idx, {gap_east_z, gap_west_z}});

      // Update the gap extent for the next gap
      gap_east_z = sign * std::max(std::abs(vol_bounds.upper_z),
                                   std::abs(vol_bounds.lower_z));

      is_gap = !is_gap;
      continue;
    }

    // Every second layer is a gap volume
    is_gap = !is_gap;
    if (is_gap) {
      auto link_east{static_cast<dindex>(static_cast<int>(vol_idx) +
                                         cfg.endcap_config().side() - 2)};
      auto link_west{static_cast<dindex>(static_cast<int>(vol_idx) -
                                         cfg.endcap_config().side() - 2)};

      const scalar_t gap_west_z{sign * std::min(std::abs(vol_bounds.upper_z),
                                                std::abs(vol_bounds.lower_z))};

      const point3_t gap_center{static_cast<scalar_t>(0),
                                static_cast<scalar_t>(0),
                                0.5f * (gap_east_z + gap_west_z)};
      vm_builder->add_volume_placement({gap_center});

      detail::add_cylinder_portals(vm_builder, cfg, gap_west_z, gap_east_z,
                                   inner_radius, outer_radius, link_north,
                                   link_south, link_east, link_west);

      volume_sizes.push_back({vol_idx, {gap_east_z, gap_west_z}});

      // Set the inner gap radius for the next gap volume
      gap_east_z = sign * std::max(std::abs(vol_bounds.upper_z),
                                   std::abs(vol_bounds.lower_z));

      vm_builder->set_name("gap_" + std::to_string(vol_idx));

    } else {
      const dindex j{i / 2u};

      auto link_east{static_cast<dindex>(static_cast<int>(vol_idx) +
                                         cfg.endcap_config().side() + 2)};
      auto link_west{static_cast<dindex>(static_cast<int>(vol_idx) -
                                         cfg.endcap_config().side() + 2)};

      // The first endacp layer needs to link to the connector gap
      // The last endcap layer needs to exit the detector
      if (sign < 0) {
        link_east = (i == 0u) ? connector_link : link_east;
        link_west = (j == cfg.n_edc_layers() - 1u) ? end_of_world : link_west;
      } else {
        link_west = (i == 0u) ? connector_link : link_west;
        link_east = (j == cfg.n_edc_layers() - 1u) ? end_of_world : link_east;
      }

      // Position the volume at the respective endcap layer position
      const scalar_t center_z{sign * cfg.endcap_layer_positions().at(j)};
      const point3_t vol_center{static_cast<scalar_t>(0),
                                static_cast<scalar_t>(0), center_z};
      vm_builder->add_volume_placement({vol_center});

      // Configure the module factory for this layer
      auto &endcap_cfg = cfg.endcap_config();

      endcap_cfg.center(center_z)
          .inner_radius(inner_radius)
          .outer_radius(outer_radius - cfg.envelope());

      // Configure the portal factory
      cylinder_portal_config<scalar_t> portal_cfg{};

      portal_cfg.envelope(cfg.envelope())
          .fixed_inner_radius(inner_radius)
          .fixed_outer_radius(outer_radius)
          // Link the volume portals to their neighbors
          .link_north(link_north)
          .link_south(link_south)
          .link_east(link_east)
          .link_west(link_west);

      // Configure the material
      cfg.material_config().thickness(cfg.module_mat_thickness());
      // Use different material thickness for cylinder portals in endcaps
      const scalar_t t{cfg.cyl_material_map().thickness};
      cfg.cyl_material_map().thickness = 0.15f * unit<scalar_t>::mm;

      // Add a layer of module surfaces (mat have material)
      auto module_mat_factory = decorate_material<detector_t>(
          cfg,
          std::make_unique<endcap_generator<detector_t, trapezoid2D>>(
              endcap_cfg),
          true);

      // Add cylinder and disc portals (may have material, depending on
      // whether material maps are being used or not)
      auto portal_mat_factory = decorate_material<detector_t>(
          cfg,
          std::make_unique<cylinder_portal_generator<detector_t>>(portal_cfg));

      // Reset cylinder material thickness
      cfg.cyl_material_map().thickness = t;

      vm_builder->add_surfaces(module_mat_factory, gctx);
      vm_builder->add_surfaces(portal_mat_factory);

      // Fill the volume extent into 'vol_bounds'
      get_volume_extent(cfg, portal_mat_factory, vol_bounds);

      // Set the new current boundaries to construct the next gap
      volume_sizes.push_back(
          {vol_idx, {vol_bounds.lower_z, vol_bounds.upper_z}});

      vm_builder->set_name("endcap_" + std::to_string(vol_idx));

      if (cfg.use_grids()) {
        // Add a disc grid to every endcap module layer
        add_disc_grid(det_builder, cfg, vol_idx);
      }
    }
  }
  return volume_sizes;
}

/// Helper method for creating the portals of one of endcap sections of the
/// beampipe volume.
///
/// @param beampipe_builder volume builder for the beampipe volume
/// @param cfg config for the toy detector
/// @param neg_edc_lay_sizes indices and z-extent of the endcap volumes of one
///                          detector side (positive or negative)
template <typename detector_builder_t, typename vol_extent_data_t>
inline void add_connector_portals(
    detector_builder_t &det_builder,
    toy_det_config<typename detector_builder_t::detector_type::scalar_type>
        &cfg,
    const dindex beampipe_idx, const vol_extent_data_t &edc_vol_extents,
    const vol_extent_data_t &brl_vol_extents) {
  using detector_t = typename detector_builder_t::detector_type;
  using algebra_t = typename detector_t::algebra_type;
  using transform3_t = dtransform3D<algebra_t>;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  using factory_interface_t = surface_factory_interface<detector_t>;
  using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
  using disc_factory_t = surface_factory<detector_t, ring2D>;

  std::shared_ptr<factory_interface_t> pt_cyl_factory{nullptr};
  std::shared_ptr<factory_interface_t> pt_disc_factory{nullptr};

  if (cfg.use_material_maps()) {
    pt_cyl_factory =
        decorate_material<detector_t>(cfg, std::make_unique<cyl_factory_t>());
    pt_disc_factory =
        decorate_material<detector_t>(cfg, std::make_unique<disc_factory_t>());
  } else {
    pt_cyl_factory = std::make_shared<cyl_factory_t>();
    pt_disc_factory = std::make_shared<disc_factory_t>();
  }

  // Mask volume link for portals that exit the detector
  constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};
  const transform3_t identity{};

  // Set the inner radius of the gap to the radius of the beampipe vol.
  const scalar_t inner_r{cfg.beampipe_vol_radius()};
  const scalar_t outer_r{cfg.outer_radius()};

  // Get the data of the gap volume
  // The connector gap is always the second volume constructed in the endcaps
  const auto &connector_gap_data = edc_vol_extents[1];
  const dindex connector_gap_idx{connector_gap_data.first};

  const scalar_t side{std::copysign(1.f, connector_gap_data.second.lower)};
  const scalar_t gap_east_z{connector_gap_data.second.lower};
  const scalar_t gap_west_z{connector_gap_data.second.upper};
  const scalar_t min_z{math::min(gap_east_z, gap_west_z)};
  const scalar_t max_z{math::max(gap_east_z, gap_west_z)};
  const point3_t gap_center{static_cast<scalar_t>(0), static_cast<scalar_t>(0),
                            0.5f * (gap_east_z + gap_west_z)};

  // Check that the volume under construction is really the connector gap
  assert(std::abs(gap_east_z) == cfg.barrel_config().half_length() ||
         std::abs(gap_west_z) == cfg.barrel_config().half_length());

  volume_builder_interface<detector_t> *connector_builder =
      det_builder[connector_gap_idx];

  connector_builder->set_name("connector_gap_" +
                              std::to_string(connector_gap_idx));
  connector_builder->add_volume_placement({gap_center});

  // Go over all volume extends and build the corresponding disc portal
  std::vector<std::vector<scalar_t>> boundaries{};
  std::vector<nav_link_t> vol_links{};
  for (const auto &e : brl_vol_extents) {
    const scalar_t min_r{math::min(e.second.lower, e.second.upper)};
    const scalar_t max_r{math::max(e.second.lower, e.second.upper)};

    boundaries.push_back(std::vector<scalar_t>{min_r, max_r});
    vol_links.push_back(static_cast<nav_link_t>(e.first));
  }
  assert(boundaries.size() == vol_links.size());

  // Barrel-facing disc portal
  pt_disc_factory->push_back(
      {surface_id::e_portal,
       transform3_t{point3_t{static_cast<scalar_t>(0), static_cast<scalar_t>(0),
                             (side < 0.f) ? max_z : min_z}},
       vol_links, boundaries});

  // Outward-facing disc portal
  pt_disc_factory->push_back(
      {surface_id::e_portal,
       transform3_t{point3_t{static_cast<scalar_t>(0), static_cast<scalar_t>(0),
                             (side < 0.f) ? min_z : max_z}},
       static_cast<nav_link_t>(connector_gap_idx - 1u),
       std::vector<scalar_t>{inner_r, outer_r}});

  // Inner cylinder portal
  pt_cyl_factory->push_back({surface_id::e_portal, identity,
                             static_cast<nav_link_t>(beampipe_idx),
                             std::vector<scalar_t>{inner_r, min_z, max_z}});

  // Outer cylinder portal
  pt_cyl_factory->push_back({surface_id::e_portal, identity, end_of_world,
                             std::vector<scalar_t>{outer_r, min_z, max_z}});

  // Add all portals to the connector gap volume
  connector_builder->add_surfaces(pt_disc_factory);
  connector_builder->add_surfaces(pt_cyl_factory);
}

/// Helper method for creating the portals of the beampipe barrel section.
///
/// @param beampipe_builder volume builder for the beampipe volume
/// @param cfg config for the toy detector
template <typename detector_t>
inline void add_beampipe_portals(
    volume_builder_interface<detector_t> *beampipe_builder,
    toy_det_config<typename detector_t::scalar_type> &cfg) {
  using algebra_t = typename detector_t::algebra_type;
  using transform3_t = dtransform3D<algebra_t>;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  using factory_interface_t = surface_factory_interface<detector_t>;
  using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
  using disc_factory_t = surface_factory<detector_t, ring2D>;

  std::shared_ptr<factory_interface_t> pt_cyl_factory{nullptr};

  if (cfg.use_material_maps()) {
    pt_cyl_factory =
        decorate_material<detector_t>(cfg, std::make_unique<cyl_factory_t>());
  } else {
    pt_cyl_factory = std::make_shared<cyl_factory_t>();
  }

  // Mask volume link for portals that exit the detector
  constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};

  const scalar_t h_z{cfg.barrel_config().half_length()};
  const scalar_t inner_r{0.f};
  const scalar_t outer_r{cfg.beampipe_vol_radius()};

  // If there are no endcaps, add the beampipe disc portals at the boundaries
  // of the barrel section
  if (cfg.n_edc_layers() == 0u) {
    std::shared_ptr<factory_interface_t> pt_disc_factory{nullptr};
    if (cfg.use_material_maps()) {
      pt_disc_factory = decorate_material<detector_t>(
          cfg, std::make_unique<disc_factory_t>());
    } else {
      pt_disc_factory = std::make_shared<disc_factory_t>();
    }

    // Lower dics portal
    pt_disc_factory->push_back(
        {surface_id::e_portal,
         transform3_t{point3_t{static_cast<scalar_t>(0),
                               static_cast<scalar_t>(0), -h_z}},
         end_of_world, std::vector<scalar_t>{inner_r, outer_r}});

    // Outward-facing disc portal
    pt_disc_factory->push_back(
        {surface_id::e_portal,
         transform3_t{
             point3_t{static_cast<scalar_t>(0), static_cast<scalar_t>(0), h_z}},
         end_of_world, std::vector<scalar_t>{inner_r, outer_r}});

    beampipe_builder->add_surfaces(pt_disc_factory);
  }

  // Cylinder portal that leads into the barrel section
  dindex first_barrel_idx{
      cfg.n_brl_layers() == 0u ? end_of_world : 2u * cfg.n_edc_layers() + 2u};
  pt_cyl_factory->push_back(
      {surface_id::e_portal, transform3_t{},
       static_cast<nav_link_t>(first_barrel_idx),
       std::vector<scalar_t>{outer_r,
                             -h_z + std::numeric_limits<scalar_t>::epsilon(),
                             h_z - std::numeric_limits<scalar_t>::epsilon()}});

  beampipe_builder->add_surfaces(pt_cyl_factory);
}

/// Helper method for creating the portals of one of the endcap sections of the
/// beampipe volume.
///
/// @param beampipe_builder volume builder for the beampipe volume
/// @param cfg config for the toy detector
/// @param neg_edc_lay_sizes indices and z-extent of the endcap volumes of one
///                          detector side (positive or negative)
template <typename detector_t, typename layer_size_cont_t>
inline void add_beampipe_portals(
    volume_builder_interface<detector_t> *beampipe_builder,
    toy_det_config<typename detector_t::scalar_type> &cfg,
    const layer_size_cont_t &edc_lay_sizes) {
  using algebra_t = typename detector_t::algebra_type;
  using transform3_t = dtransform3D<algebra_t>;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  using factory_interface_t = surface_factory_interface<detector_t>;
  using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
  using disc_factory_t = surface_factory<detector_t, ring2D>;

  // Mask volume link for portals that exit the detector
  constexpr auto end_of_world{detail::invalid_value<nav_link_t>()};
  // For which endcap side should the portals be constructed?
  const scalar_t side{std::copysign(1.f, edc_lay_sizes.front().second.lower)};

  std::shared_ptr<factory_interface_t> pt_cyl_factory{nullptr};
  std::shared_ptr<factory_interface_t> pt_disc_factory{nullptr};

  if (cfg.use_material_maps()) {
    pt_cyl_factory =
        decorate_material<detector_t>(cfg, std::make_unique<cyl_factory_t>());
    pt_disc_factory =
        decorate_material<detector_t>(cfg, std::make_unique<disc_factory_t>());
  } else {
    pt_cyl_factory = std::make_shared<cyl_factory_t>();
    pt_disc_factory = std::make_shared<disc_factory_t>();
  }

  const scalar_t inner_r{0.f};
  const scalar_t outer_r{cfg.beampipe_vol_radius()};
  scalar_t disc_pos_z{-side * std::numeric_limits<scalar_t>::max()};

  // Go over all volume extends and build the corresponding cylinder portals
  std::vector<std::vector<scalar_t>> boundaries{};
  std::vector<nav_link_t> vol_links{};
  for (const auto &e : edc_lay_sizes) {
    const scalar_t min_z{math::min(e.second.lower, e.second.upper)};
    const scalar_t max_z{math::max(e.second.lower, e.second.upper)};

    if (side < 0.f) {
      disc_pos_z = (disc_pos_z > min_z) ? min_z : disc_pos_z;
    } else {
      disc_pos_z = (disc_pos_z < max_z) ? max_z : disc_pos_z;
    }

    boundaries.push_back(std::vector<scalar_t>{outer_r, min_z, max_z});
    vol_links.push_back(static_cast<nav_link_t>(e.first));
  }
  assert(boundaries.size() == vol_links.size());

  pt_cyl_factory->push_back(
      {surface_id::e_portal, transform3_t{}, vol_links, boundaries});

  // Outward-facing disc portal
  pt_disc_factory->push_back(
      {surface_id::e_portal,
       transform3_t{point3_t{static_cast<scalar_t>(0), static_cast<scalar_t>(0),
                             disc_pos_z}},
       end_of_world, std::vector<scalar_t>{inner_r, outer_r}});

  // Add all portals to the beampipe volume
  beampipe_builder->add_surfaces(pt_disc_factory);
  beampipe_builder->add_surfaces(pt_cyl_factory);
}

}  // namespace detail

/// Builds a detray geometry that contains the innermost tml layers. The number
/// of barrel and endcap layers can be chosen, but all barrel layers should be
/// present when an endcap detector is built to have the barrel region radius
/// match the endcap diameter.
///
/// @param resource vecmem memory resource to use for container allocations
/// @param cfg toy detector configuration
///
/// @returns a complete detector object
template <concepts::algebra algebra_t>
inline auto build_toy_detector(vecmem::memory_resource &resource,
                               toy_det_config<dscalar<algebra_t>> cfg = {}) {
  using scalar_t = dscalar<algebra_t>;
  using transform3_t = dtransform3D<algebra_t>;

  using builder_t = detector_builder<toy_metadata<algebra_t>, volume_builder>;
  using detector_t = typename builder_t::detector_type;
  using nav_link_t = typename detector_t::surface_type::navigation_link;
  using cyl_factory_t = surface_factory<detector_t, concentric_cylinder2D>;
  using vol_extent_container_t =
      std::vector<std::pair<dindex, detail::extent2D<scalar_t>>>;

  // Check config
  if (cfg.n_edc_layers() > cfg.endcap_layer_positions().size()) {
    throw std::invalid_argument(
        "ERROR: Too many endcap layers requested (max " +
        std::to_string(cfg.endcap_layer_positions().size()) + ")!");
  }
  if (cfg.n_brl_layers() > cfg.barrel_layer_radii().size() - 1u) {
    throw std::invalid_argument(
        "ERROR: Too many barrel layers requested (max " +
        std::to_string(cfg.barrel_layer_radii().size() - 1u) + ")!");
  }
  if (cfg.n_edc_layers() > 0 && cfg.n_brl_layers() < 4) {
    throw std::invalid_argument(
        "ERROR: All four barrel layers need to be present in order to add "
        "endcap layers");
  }

  // Toy detector builder
  builder_t det_builder;
  det_builder.set_name("toy_detector");

  // Geometry context object
  typename detector_t::geometry_context gctx{};

  // Add the volume that contains the beampipe
  cfg.material_config().thickness(cfg.beampipe_mat_thickness());
  auto beampipe_builder = detail::decorate_material(
      cfg, det_builder, det_builder.new_volume(volume_id::e_cylinder));

  const dindex beampipe_idx{beampipe_builder->vol_index()};
  beampipe_builder->add_volume_placement(transform3_t{});
  beampipe_builder->set_name("beampipe_" + std::to_string(beampipe_idx));

  // Add the beampipe as a passive material surface
  auto beampipe_factory = detail::decorate_material<detector_t>(
      cfg, std::make_unique<cyl_factory_t>(), true);

  scalar_t max_z{cfg.n_edc_layers() == 0u ? cfg.barrel_config().half_length()
                                          : cfg.endcap_layer_positions().at(
                                                cfg.n_edc_layers() - 1u)};
  scalar_t min_z{-max_z};

  beampipe_factory->push_back(
      {surface_id::e_passive, transform3_t{},
       static_cast<nav_link_t>(beampipe_idx),
       std::vector<scalar_t>{cfg.barrel_layer_radii().at(0), min_z, max_z}});

  beampipe_builder->add_surfaces(beampipe_factory);

  // Build the negative endcap
  vol_extent_container_t neg_edc_vol_extents;
  if (cfg.n_edc_layers() > 0u) {
    cfg.endcap_config().side(-1);

    neg_edc_vol_extents =
        detail::add_endcap_detector(det_builder, gctx, cfg, beampipe_idx);

    // Add the beampipe volume portals for the negative endcap section
    detail::add_beampipe_portals(beampipe_builder, cfg, neg_edc_vol_extents);
  }
  // Build the barrel section
  vol_extent_container_t brl_vol_extents;
  if (cfg.n_brl_layers() > 0u) {
    brl_vol_extents =
        detail::add_barrel_detector(det_builder, gctx, cfg, beampipe_idx);
  }
  // Add the beampipe volume portals for the barrel section
  detail::add_beampipe_portals(beampipe_builder, cfg);

  // Build the positive endcap
  vol_extent_container_t pos_edc_vol_extents;
  if (cfg.n_edc_layers() > 0u) {
    cfg.endcap_config().side(1);

    pos_edc_vol_extents =
        detail::add_endcap_detector(det_builder, gctx, cfg, beampipe_idx);

    // Add the beampipe volume portals for the positive endcap section
    detail::add_beampipe_portals(beampipe_builder, cfg, pos_edc_vol_extents);

    // Add the connection between barrel and both endcaps
    // Negative endcap
    add_connector_portals(det_builder, cfg, beampipe_idx, neg_edc_vol_extents,
                          brl_vol_extents);
    // Positive endcap
    add_connector_portals(det_builder, cfg, beampipe_idx, pos_edc_vol_extents,
                          brl_vol_extents);
  }

  // Build and return the detector and fill the name map
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
