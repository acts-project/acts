// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/surface_factory.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/detectors/nested_volumes_metadata.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/utils/consistency_checker.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/print_detector.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace detray {

/// @brief Configuration of the nested test detector
///
/// The detector is a cuboid world volume with cuboid daughter volumes:
///
/// - volume 0: world, x, y, z in [-200, 200] mm
/// - volume 1: box A, daughter of the world, z in [-140, -60] mm
/// - volume 2: box G, daughter of box A (second hierarchy level), flush with
///             the -z boundary of box A (moved inwards by @c inset )
/// - volume 3: box B, daughter of the world, flush with the +z world boundary
///             (moved inwards by @c inset )
/// - volume 4: box C, daughter of the world, touches box D
/// - volume 5: box D, daughter of the world, touches box C at z = 30 mm
///
/// Every volume contains sensitive rectangle surfaces.
struct nested_detector_config {
  /// Distance by which the boxes B and G are moved away from the boundary of
  /// their mother volume (zero: the boxes are flush with the boundary)
  float m_inset{0.f};
  /// Run the consistency check after building the detector
  bool m_do_check{true};

  /// Setters
  /// @{
  /// @param inset distance between box B and the +z world boundary
  /// @returns this configuration
  nested_detector_config &inset(const float inset) {
    m_inset = inset;
    return *this;
  }
  /// @param check run the consistency check
  /// @returns this configuration
  nested_detector_config &do_check(const bool check) {
    m_do_check = check;
    return *this;
  }
  /// @}

  /// Getters
  /// @{
  /// @returns the distance between box B and the +z world boundary
  float inset() const { return m_inset; }
  /// @returns whether to run the consistency check
  bool do_check() const { return m_do_check; }
  /// @}
};

namespace detail {

/// Box faces in the order of the portal volume links
enum class box_face : unsigned int {
  e_neg_x = 0u,
  e_pos_x = 1u,
  e_neg_y = 2u,
  e_pos_y = 3u,
  e_neg_z = 4u,
  e_pos_z = 5u,
};

/// Add the six rectangle portals of a cuboid volume to a volume builder
///
/// @tparam detector_t the detector type
///
/// @param v_builder the volume builder
/// @param center the center of the box in global coordinates
/// @param half_lengths the half lengths of the box in x, y and z
/// @param links the volume links of the faces (order of @c box_face )
template <typename detector_t>
inline void add_box_portals(
    volume_builder_interface<detector_t> &v_builder,
    const typename detector_t::point3_type &center,
    const std::array<typename detector_t::scalar_type, 3> &half_lengths,
    const std::array<typename detector_t::surface_type::navigation_link, 6>
        &links) {
  using scalar_t = typename detector_t::scalar_type;
  using vector3_t = typename detector_t::vector3_type;
  using transform3_t = typename detector_t::transform3_type;
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;

  const scalar_t h_x{half_lengths[0]};
  const scalar_t h_y{half_lengths[1]};
  const scalar_t h_z{half_lengths[2]};

  auto pt_factory = std::make_shared<rectangle_factory>();

  auto link = [&links](box_face face) {
    return links[static_cast<unsigned int>(face)];
  };

  // x-y plane: no rotation, shift in z
  const vector3_t shift_z{0.f, 0.f, h_z};
  pt_factory->push_back({surface_id::e_portal, transform3_t{center - shift_z},
                         link(box_face::e_neg_z),
                         std::vector<scalar_t>{-h_x, -h_y, h_x, h_y}});
  pt_factory->push_back({surface_id::e_portal, transform3_t{center + shift_z},
                         link(box_face::e_pos_z),
                         std::vector<scalar_t>{-h_x, -h_y, h_x, h_y}});

  // x-z plane: rotate by 90deg around x-axis, shift in y
  const vector3_t shift_y{0.f, h_y, 0.f};
  const vector3_t xz_x{1.f, 0.f, 0.f};
  const vector3_t xz_z{0.f, -1.f, 0.f};
  pt_factory->push_back(
      {surface_id::e_portal, transform3_t{center - shift_y, xz_z, xz_x},
       link(box_face::e_neg_y), std::vector<scalar_t>{-h_x, -h_z, h_x, h_z}});
  pt_factory->push_back(
      {surface_id::e_portal, transform3_t{center + shift_y, xz_z, xz_x},
       link(box_face::e_pos_y), std::vector<scalar_t>{-h_x, -h_z, h_x, h_z}});

  // y-z plane: rotate by 90deg around y-axis, shift in x
  const vector3_t shift_x{h_x, 0.f, 0.f};
  const vector3_t yz_x{0.f, 0.f, -1.f};
  const vector3_t yz_z{1.f, 0.f, 0.f};
  pt_factory->push_back(
      {surface_id::e_portal, transform3_t{center - shift_x, yz_z, yz_x},
       link(box_face::e_neg_x), std::vector<scalar_t>{-h_z, -h_y, h_z, h_y}});
  pt_factory->push_back(
      {surface_id::e_portal, transform3_t{center + shift_x, yz_z, yz_x},
       link(box_face::e_pos_x), std::vector<scalar_t>{-h_z, -h_y, h_z, h_y}});

  v_builder.add_surfaces(pt_factory);
}

/// Add sensitive rectangles in the x-y plane to a volume builder
///
/// @tparam detector_t the detector type
///
/// @param v_builder the volume builder
/// @param z_positions the global z positions of the rectangles
/// @param bounds the rectangle bounds (min x, min y, max x, max y)
/// @param shift global shift of the rectangles in x and y
template <typename detector_t>
inline void add_sensitives(
    volume_builder_interface<detector_t> &v_builder,
    const std::vector<typename detector_t::scalar_type> &z_positions,
    const std::vector<typename detector_t::scalar_type> &bounds,
    const std::array<typename detector_t::scalar_type, 2> &shift = {0.f, 0.f}) {
  using point3_t = typename detector_t::point3_type;
  using transform3_t = typename detector_t::transform3_type;
  using nav_link_t = typename detector_t::surface_type::navigation_link;
  using rectangle_factory = surface_factory<detector_t, rectangle2D>;

  auto sf_factory = std::make_shared<rectangle_factory>();
  const auto vol_link{static_cast<nav_link_t>(v_builder.vol_index())};

  for (const auto z : z_positions) {
    sf_factory->push_back({surface_id::e_sensitive,
                           transform3_t{point3_t{shift[0], shift[1], z}},
                           vol_link, bounds});
  }

  v_builder.add_surfaces(sf_factory);
}

}  // namespace detail

/// Build the nested test detector (@see @c nested_detector_config )
///
/// @tparam algebra_t the algebra type
///
/// @param resource the memory resource for the detector containers
/// @param cfg the detector configuration
///
/// @returns the detector and its volume name map
template <concepts::algebra algebra_t>
inline auto build_nested_detector(vecmem::memory_resource &resource,
                                  const nested_detector_config &cfg = {}) {
  using metadata_t = nested_volumes_metadata<algebra_t>;
  using builder_t = detector_builder<metadata_t>;
  using detector_t = typename builder_t::detector_type;
  using scalar_t = dscalar<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using nav_link_t = typename detector_t::surface_type::navigation_link;

  constexpr nav_link_t end_of_world{detail::invalid_value<nav_link_t>()};

  builder_t det_builder{};
  det_builder.set_name("nested_detector");

  // Volume indices
  constexpr nav_link_t world{0u};
  constexpr nav_link_t box_a{1u};
  constexpr nav_link_t box_g{2u};
  constexpr nav_link_t box_b{3u};
  constexpr nav_link_t box_c{4u};
  constexpr nav_link_t box_d{5u};

  auto add_box = [&det_builder](const std::string &name, const dindex mother,
                                const point3_t &center,
                                const std::array<scalar_t, 3> &half_lengths,
                                const std::array<nav_link_t, 6> &links) {
    auto *v_builder = det_builder.new_volume(volume_id::e_cuboid);
    v_builder->set_name(name);
    v_builder->set_mother(mother);
    v_builder->add_volume_placement(center);
    detail::add_box_portals(*v_builder, center, half_lengths, links);

    return v_builder;
  };

  // World volume
  auto *w_builder =
      add_box("world", dindex_invalid, {0.f, 0.f, 0.f}, {200.f, 200.f, 200.f},
              {end_of_world, end_of_world, end_of_world, end_of_world,
               end_of_world, end_of_world});
  detail::add_sensitives(*w_builder, {-180.f, 100.f},
                         {-190.f, -190.f, 190.f, 190.f});
  // Does not reach into the daughter volumes C and D
  detail::add_sensitives(*w_builder, {0.f}, {-190.f, -190.f, 190.f, 80.f});

  // Box A with a daughter volume (box G)
  auto *a_builder =
      add_box("box_a", world, {0.f, 0.f, -100.f}, {60.f, 60.f, 40.f},
              {world, world, world, world, world, world});
  detail::add_sensitives(*a_builder, {-100.f, -70.f},
                         {-55.f, -55.f, 55.f, 55.f});

  // Box G: flush with the boundary of box A at -z (unless moved inwards)
  const scalar_t z_g{-130.f + cfg.inset()};
  const nav_link_t g_neg_z{cfg.inset() > 0.f ? box_a : world};
  auto *g_builder = add_box("box_g", box_a, {0.f, 0.f, z_g}, {20.f, 20.f, 10.f},
                            {box_a, box_a, box_a, box_a, g_neg_z, box_a});
  detail::add_sensitives(*g_builder, {z_g}, {-15.f, -15.f, 15.f, 15.f});

  // Box B: flush with the world boundary at +z (unless moved inwards)
  const scalar_t z_b{160.f - cfg.inset()};
  const nav_link_t b_pos_z{cfg.inset() > 0.f ? world : end_of_world};
  auto *b_builder = add_box("box_b", world, {0.f, 0.f, z_b}, {60.f, 60.f, 40.f},
                            {world, world, world, world, world, b_pos_z});
  detail::add_sensitives(*b_builder, {z_b}, {-55.f, -55.f, 55.f, 55.f});

  // Boxes C and D: touching siblings
  auto *c_builder =
      add_box("box_c", world, {0.f, 130.f, 0.f}, {40.f, 40.f, 30.f},
              {world, world, world, world, world, box_d});
  detail::add_sensitives(*c_builder, {0.f}, {-35.f, -35.f, 35.f, 35.f},
                         {0.f, 130.f});

  auto *d_builder =
      add_box("box_d", world, {0.f, 130.f, 60.f}, {40.f, 40.f, 30.f},
              {world, world, world, world, box_c, world});
  detail::add_sensitives(*d_builder, {60.f}, {-35.f, -35.f, 35.f, 35.f},
                         {0.f, 130.f});

  static_assert(box_a < box_g && box_g < box_b);
  assert(det_builder.n_volumes() == box_d + 1u);

  // Build and check the detector
  typename detector_t::name_map name_map{};
  auto det = det_builder.build(resource, volume_builder_options{}, name_map);

  if (cfg.do_check()) {
    detray::detail::check_consistency(det, false, name_map);
  }

  DETRAY_DEBUG_HOST("\n" << detray::utils::print_detector(det, name_map));

  return std::make_pair(std::move(det), std::move(name_map));
}

}  // namespace detray
