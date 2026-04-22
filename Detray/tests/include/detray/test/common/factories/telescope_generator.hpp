// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/core/detail/data_context.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/tracks/ray.hpp"
#include "detray/utils/axis_rotation.hpp"
#include "detray/utils/logging.hpp"
#include "detray/utils/unit_vectors.hpp"

// System include(s)
#include <cassert>
#include <limits>

namespace detray {

/// @brief Generates a number of surfaces along a given direction
///
/// @tparam detector_t the type of detector the volume belongs to.
template <typename detector_t, typename mask_shape_t = rectangle2D,
          typename trajectory_t =
              detail::ray<typename detector_t::algebra_type>>
class telescope_generator final : public surface_factory_interface<detector_t> {
  using algebra_t = typename detector_t::algebra_type;
  using scalar_t = dscalar<algebra_t>;
  using transform3_t = dtransform3D<algebra_t>;
  using point3_t = dpoint3D<algebra_t>;
  using vector3_t = dvector3D<algebra_t>;

 public:
  /// Build a surface at with extent given in @param boundaries at every
  /// position in @param positions along the pilot-track @param traj.
  DETRAY_HOST
  telescope_generator(
      std::vector<scalar_t> positions,
      darray<scalar_t, mask_shape_t::boundaries::e_size> boundaries,
      trajectory_t traj)
      : m_traj{traj}, m_positions{positions}, m_boundaries{boundaries} {}

  /// Infer the sensitive surface placement from the telescope @param length
  /// if no concrete positions were given.
  /// @param n_surfaces number of surfaces to be generated
  /// @param boundaries mask boundaries of the surfaces
  /// @param traj pilot track along which to build the telescope
  DETRAY_HOST
  telescope_generator(
      scalar_t length, std::size_t n_surfaces,
      darray<scalar_t, mask_shape_t::boundaries::e_size> boundaries,
      trajectory_t traj)
      : m_traj{traj}, m_positions{}, m_boundaries{boundaries} {
    scalar_t pos{0.f};
    scalar_t dist{n_surfaces > 1u
                      ? length / static_cast<scalar_t>(n_surfaces - 1u)
                      : 0.f};
    for (std::size_t i = 0u; i < n_surfaces; ++i) {
      m_positions.push_back(pos);
      pos += dist;
    }
  }

  /// @returns the number of surfaces this factory will produce
  DETRAY_HOST
  auto size() const -> dindex override {
    return static_cast<dindex>(m_positions.size());
  }

  /// Clear the positions and boundaries of the surfaces.
  DETRAY_HOST
  void clear() override {
    m_positions.clear();
    m_boundaries = {};
  };

  /// This is a surface generator, no external surface data needed
  /// @{
  DETRAY_HOST
  void push_back(surface_data<detector_t> &&) override { /*Do nothing*/ }
  DETRAY_HOST
  auto push_back(std::vector<surface_data<detector_t>> &&)
      -> void override { /*Do nothing*/ }
  /// @}

  /// Create a surface telescope.
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
    DETRAY_VERBOSE_HOST("Generate telescope modules...");
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

    // Create the module centers
    const auto mod_placements = module_positions(m_traj, m_positions);

    // Create geometry data
    for (const auto &mod_placement : mod_placements) {
      auto mask_volume_link{static_cast<nav_link_t>(volume_idx)};

      // Surfaces with the linking into the local containers
      mask_link_t mask_link{mask_id, masks.template size<mask_id>()};
      material_link_t material_link{no_material, dindex_invalid};

      const auto trf_index = transforms.size(ctx);
      surfaces.push_back({trf_index, mask_link, material_link, volume_idx,
                          surface_id::e_sensitive},
                         invalid_src_link);

      // The rectangle bounds for this module
      masks.template emplace_back<mask_id>(empty_context{}, m_boundaries,
                                           mask_volume_link);

      // Build the transform
      // Local z axis is the global normal vector
      vector3_t m_local_z = vector::normalize(mod_placement.dir);

      if constexpr (std::is_same_v<mask_shape_t, detray::line_square> ||
                    std::is_same_v<mask_shape_t, detray::line_circular>) {
        // For a telescope with wires, rotate z axis 90
        // degree around vector on x-y plane
        auto curvi_u =
            unit_vectors<vector3_t>().make_curvilinear_unit_u(m_local_z);
        axis_rotation<algebra_t> axis_rot(curvi_u,
                                          -constant<scalar_t>::pi / 2.f);
        m_local_z = axis_rot(m_local_z);
      }

      // Local x axis is the curvilinear vector with respect to local_z
      auto m_local_x =
          unit_vectors<vector3_t>().make_curvilinear_unit_u(m_local_z);

      // Create the global-to-local transform of the module
      transforms.emplace_back(ctx, mod_placement.pos, m_local_z, m_local_x);
    }

    return {surfaces_offset, static_cast<dindex>(surfaces.size())};
  }

 private:
  /// Where and how to place the telescope modules.
  struct module_placement {
    /// Module position
    point3_t pos;
    /// Module normal
    vector3_t dir;
  };

  /// Helper method for positioning the surfaces.
  ///
  /// @param traj pilot trajectory along which the modules should be placed.
  /// @param steps lengths along the trajectory where surfaces should be
  ///              placed.
  ///
  /// @return a vector of the @c module_placements along the trajectory.
  inline std::vector<module_placement> module_positions(
      const trajectory_t &traj, const std::vector<scalar_t> &steps) const {
    // create and fill the module placements
    std::vector<module_placement> placements;
    placements.reserve(steps.size());

    for (const auto s : steps) {
      placements.push_back({traj.pos(s), traj.dir(s)});
    }

    return placements;
  }

  /// "pilot-track" along which to place the surfaces
  trajectory_t m_traj;
  /// Positions of the surfaces in the telescope along the pilot track
  std::vector<scalar_t> m_positions;
  /// The boundary values for the surface mask
  darray<scalar_t, mask_shape_t::boundaries::e_size> m_boundaries;
};

}  // namespace detray
