// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/material/material_rod.hpp"
#include "detray/material/material_slab.hpp"
#include "detray/navigation/accelerators/brute_force.hpp"

namespace detray {

/// Defines a special telescope test detector type with only rectangle portals
/// and one additional kind of contained module surfaces (@tparam mask_shape_t)
template <concepts::algebra algebra_t, typename mask_shape_t = rectangle2D>
struct telescope_metadata {
  /// Define the algebra type for the geometry and navigation
  using algebra_type = algebra_t;
  using scalar_t = dscalar<algebra_type>;

  /// Index type in the masks to find the adjacent volumes
  using nav_link = std::uint_least16_t;

  /// How to store coordinate transform matrices
  template <template <typename...> class vector_t = dvector>
  using transform_store =
      single_store<dtransform3D<algebra_type>, vector_t, geometry_context>;

  //
  // Surface Primitives
  //

  /// Rectangles are needed for the portals
  using rectangle2D_t = mask<detray::rectangle2D, algebra_type, nav_link>;
  using module_t = mask<mask_shape_t, algebra_type, nav_link>;

  /// Rectangles are always needed as portals (but they have the same type as
  /// module rectangles). Only one additional mask shape is allowed
  enum class mask_id : std::uint_least8_t {
    e_rectangle2D = 0u,
    e_annulus2D = 1u,
    e_cylinder2D = 1u,
    e_ring2D = 1u,
    e_trapezoid2D = 1u,
    e_single1D = 1u,
    e_single2D = 1u,
    e_single3D = 1u,
    e_straw_tube = 1u,
    e_drift_cell = 1u,
    e_unbounded_annulus2 = 1u,
    e_unbounded_cell2 = 1u,
    e_unbounded_cylinder2 = 1u,
    e_unbounded_disc2 = 1u,
    e_unbounded_rectangle2 = 1u,
    e_unbounded_trapezoid2 = 1u,
    e_unbounded_line_circular2 = 1u,
    e_unmasked2 = 1u,
  };

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     mask_id id) {
    switch (id) {
      case mask_id::e_rectangle2D:
        os << "e_rectangle2D";
        break;
      case mask_id::e_annulus2D:
        // All other values are 1u, showing first alphabetically
        os << "telescope module shape";
        break;
      default:
        os << "invalid";
    }
    return os;
  };

  /// How to store masks
  template <template <typename...> class vector_t = dvector>
  using mask_store = std::conditional_t<
      std::same_as<module_t, rectangle2D_t>,
      regular_multi_store<mask_id, empty_context, dtuple, vector_t,
                          rectangle2D_t>,
      regular_multi_store<mask_id, empty_context, dtuple, vector_t,
                          rectangle2D_t, module_t>>;
  //
  // Material Description
  //

  /// Material types (slabs for portals, rods are optional)
  using material_slab_t = material_slab<scalar_t>;
  using material_rod_t = material_rod<scalar_t>;

  enum class material_id : std::uint_least8_t {
    e_material_slab = 0u,
    e_raw_material = 1u,  //< used for homogeneous volume material
    e_material_rod = 2u,
    e_none = 3u,
  };

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     material_id id) {
    switch (id) {
      case material_id::e_material_slab:
        os << "e_material_slab";
        break;
      case material_id::e_raw_material:
        os << "e_raw_material";
        break;
      case material_id::e_material_rod:
        os << "e_material_rod";
        break;
      case material_id::e_none:
        os << "e_none";
        break;
      default:
        os << "invalid";
    }
    return os;
  };

  /// How to store materials
  template <typename container_t = host_container_types>
  using material_store = std::conditional_t<
      concepts::line_shape<mask_shape_t, algebra_type>,
      regular_multi_store<material_id, empty_context, dtuple,
                          container_t::template vector_type, material_slab_t,
                          material<scalar_t>, material_rod_t>,
      regular_multi_store<material_id, empty_context, dtuple,
                          container_t::template vector_type, material_slab_t,
                          material<scalar_t>>>;

  /// How to link to the entries in the data stores
  using transform_link = typename transform_store<>::single_link;
  using mask_link = typename mask_store<>::single_link;
  using material_link = typename material_store<>::single_link;
  /// Surface type used for sensitives, passives and portals
  using surface_type =
      surface_descriptor<mask_link, material_link, transform_link, nav_link>;

  //
  // Acceleration structures
  //

  /// Acceleration data structures
  enum class accel_id : std::uint_least8_t {
    e_surface_brute_force = 0u,
    e_volume_brute_force = 1u,
    e_surface_default = e_surface_brute_force,
    e_volume_default = e_volume_brute_force,
  };

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     accel_id id) {
    switch (id) {
      case accel_id::e_surface_brute_force:
        os << "e_surface_brute_force";
        break;
      case accel_id::e_volume_brute_force:
        os << "e_volume_brute_force";
        break;
      default:
        os << "invalid";
    }
    return os;
  };

  /// How to store the brute force search data structure
  template <typename container_t = host_container_types>
  using accelerator_store =
      multi_store<accel_id, empty_context, dtuple,
                  brute_force_collection<surface_type, container_t>,
                  brute_force_collection<dindex, container_t>>;

  //
  // Volume descriptors
  //

  /// No grids/other acceleration data structure, everything is brute forced
  enum geo_objects : std::uint_least8_t {
    e_sensitive = 0u,
    e_portal = 0u,
    e_passive = 0u,
    e_volume = 1u,
    e_size = 2u,
    e_all = e_size,
  };

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     geo_objects id) {
    switch (id) {
      case geo_objects::e_passive:
        os << "e_sensitive/e_portal/e_passive";
        break;
      case geo_objects::e_volume:
        os << "e_volume";
        break;
      case geo_objects::e_size:
        os << "e_size";
        break;
      default:
        os << "invalid";
    }
    return os;
  };

  /// One link for all surfaces (in the brute force method)
  using object_link_type =
      dmulti_index<dtyped_index<accel_id, dindex>, geo_objects::e_size>;
};

}  // namespace detray
