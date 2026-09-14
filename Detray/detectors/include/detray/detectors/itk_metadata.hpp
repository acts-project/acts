// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/detail/multi_store.hpp"
#include "detray/core/detail/single_store.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/annulus2D.hpp"
#include "detray/geometry/shapes/concentric_cylinder2D.hpp"
#include "detray/geometry/shapes/cylinder3D.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/geometry/shapes/ring2D.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/material/material_map.hpp"
#include "detray/navigation/accelerators/brute_force.hpp"
#include "detray/navigation/accelerators/spatial_grid.hpp"

namespace detray {

/// Assembles the detector type. This metatdata contains all available types
template <concepts::algebra algebra_t>
struct itk_metadata {

  /// Define the algebra type for the geometry and navigation
  using algebra_type = algebra_t;
  using scalar_t = dscalar<algebra_type>;

  /// Mask-to-(next)-volume link (potentially switchable for SoA)
  using nav_link = std::uint_least16_t;

  /// How to store coordinate transform matrices
  template <template <typename...> class vector_t = dvector>
  using transform_store =
      single_store<dtransform3D<algebra_type>, vector_t, geometry_context>;

  //
  // Surface Primitives
  //

  /// Mask types
  using rectangle = mask<rectangle2D, algebra_type, nav_link>;
  using annulus = mask<annulus2D, algebra_type, nav_link>;
  using cylinder_portal = mask<concentric_cylinder2D, algebra_type, nav_link>;
  using disc = mask<ring2D, algebra_type, nav_link>;

  /// Give your mask types a name (needs to be consecutive and has to match
  /// the types position in the mask store!)
  enum class mask_ids : std::uint_least8_t {
    e_rectangle2 = 0u,
    e_annulus2 = 1u,
    e_cylinder2 = 2u,
    e_portal_cylinder2 = 2u,
    e_ring2 = 3u,
    e_portal_ring2 = 3u,
  };

  // Compatibility aliases for newer metadata concept checks.
  using mask_id = mask_ids;

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     mask_ids mid) {

    switch (mid) {
      case mask_ids::e_rectangle2:
        // e_portal_rectangle2 has same value (0u)
        os << "e_rectangle2/e_portal_rectangle2";
        break;
      case mask_ids::e_annulus2:
        os << "e_annulus2";
        break;
      case mask_ids::e_cylinder2:
        os << "e_cylinder2/e_portal_cylinder2";
        break;
      case mask_ids::e_ring2:
        // e_portal_ring2 has same value (3u)
        os << "e_ring2/e_portal_ring2";
        break;
      default:
        os << "invalid";
    }
    return os;
  }

  /// How to store masks
  template <template <typename...> class vector_t = dvector>
  using mask_store =
      regular_multi_store<mask_ids, empty_context, dtuple, vector_t,
                          rectangle, annulus, cylinder_portal, disc>;

  //
  // Material Description
  //

  /// Material grid types (default: closed bounds, regular binning)
  /// @{

  // Disc material grid
  template <typename container_t>
  using disc_map_t = material_map<algebra_type, ring2D, container_t>;

  // Concentric cylindrical material grid
  template <typename container_t>
  using concentric_cylinder2_map_t =
      material_map<algebra_type, concentric_cylinder2D, container_t>;

  /// @}

  /// Give your material types a name (needs to be consecutive and has to
  /// match the types position in the mask store!)
  enum class material_ids : std::uint_least8_t {
    // Material texture (grid) shapes
    e_concentric_cylinder2_map = 0u,
    e_disc2_map = 1u,
    e_none = 2u,
  };

  // Compatibility aliases for newer metadata concept checks.
  using material_id = material_ids;

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     material_ids mid) {

    switch (mid) {
      case material_ids::e_concentric_cylinder2_map:
        os << "e_concentric_cylinder2_map";
        break;
      case material_ids::e_disc2_map:
        // e_annulus2_map has same value (1u)
        os << "e_disc2_map/e_annulus2_map";
        break;
      case material_ids::e_none:
        os << "e_none";
        break;
      default:
        os << "invalid";
    }
    return os;
  }

  /// How to store materials
  template <typename container_t = host_container_types>
  using material_store =
      multi_store<material_ids, empty_context, dtuple,
                  grid_collection<concentric_cylinder2_map_t<container_t>>,
                  grid_collection<disc_map_t<container_t>>>;

  //
  // Acceleration structures
  //

  /// surface grid types (default boundaries: closed binning)
  /// @TODO: Will we need the types for all grid configurations (binning,
  /// bin boundaries, serializers)?
  /// @{

  // surface grid definition: bin-content: darray<surface_type, 9>
  template <typename axes_t, typename bin_entry_t, typename container_t>
  using surface_grid_t =
      spatial_grid<algebra_type, axes_t, bins::dynamic_array<bin_entry_t>,
                   simple_serializer, container_t, false>;

  // 2D cylindrical grid for the barrel layers
  template <typename bin_entry_t, typename container_t>
  using cylinder2D_sf_grid =
      surface_grid_t<axes<concentric_cylinder2D>, bin_entry_t, container_t>;

  // disc grid for the endcap layers
  template <typename bin_entry_t, typename container_t>
  using disc_sf_grid = surface_grid_t<axes<ring2D>, bin_entry_t, container_t>;

  /// @}

  /// How to link to the entries in the data stores
  using transform_link = typename transform_store<>::single_link;
  using mask_link = typename mask_store<>::range_link;
  using material_link = typename material_store<>::single_link;
  /// Surface type used for sensitives, passives and portals
  using surface_type =
      surface_descriptor<mask_link, material_link, transform_link, nav_link>;

  //
  // Volume descriptors
  //

  /// How to index the constituent objects (surfaces) in a volume
  /// If they share the same index value here, they will be added into the
  /// same acceleration data structure (brute force is always at 0)
  enum geo_objects : std::uint_least8_t {
    e_portal = 0u,    // Brute force search
    e_passive = 0u,   // Brute force search
    e_sensitive = 1u, // Grid accelerated search (can be different types)
    e_volume = 2u,    // Daughter volumes
    e_size = 3u,      // Every volume holds two acceleration data structures
    e_all = e_size,   // i.e. the brute force method and one grid type
  };

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     geo_objects gobj) {

    switch (gobj) {
      case geo_objects::e_portal:
        // e_passive has same value (0u)
        os << "e_portal/e_passive";
        break;
      case geo_objects::e_sensitive:
        os << "e_sensitive";
        break;
      case geo_objects::e_volume:
        os << "e_volume";
        break;
      case geo_objects::e_size:
        // e_all has same value (3u)
        os << "e_size/e_all";
        break;
      default:
        os << "invalid";
    }
    return os;
  }

  /// Acceleration data structures
  enum class accel_id : std::uint_least8_t {
    e_surface_brute_force = 0u,
    e_surface_cylinder2_grid = 1u,
    e_surface_disc_grid = 2u,
    e_volume_cylinder3D_grid = 3u,
    e_surface_default = e_surface_brute_force,
    e_volume_default = e_volume_cylinder3D_grid,
  };

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     accel_id aid) {

    switch (aid) {
      case accel_id::e_surface_brute_force:
        os << "e_surface_brute_force";
        break;
      case accel_id::e_surface_cylinder2_grid:
        os << "e_surface_cylinder2_grid";
        break;
      case accel_id::e_surface_disc_grid:
        os << "e_surface_disc_grid";
        break;
      case accel_id::e_volume_cylinder3D_grid:
        os << "e_volume_cylinder3D_grid";
        break;
      default:
        os << "invalid";
    }
    return os;
  }

  /// How a volume links to the acceleration data structures
  /// In this case: One link for portals/passives and one sensitive surfaces
  using object_link_type =
      dmulti_index<dtyped_index<accel_id, dindex>, geo_objects::e_size>;

  //
  // Volume acceleration structure
  //

  /// Volume search grid
  template <typename container_t = host_container_types>
  using volume_accelerator =
      spatial_grid<algebra_type,
                   axes<cylinder3D, axis::bounds::e_open, axis::irregular,
                        axis::regular, axis::irregular>,
                   bins::single<dindex>, simple_serializer, container_t,
                   false>;

  /// The tuple store that hold the acceleration data structures for all
  /// volumes. Every collection of acceleration data structures defines its
  /// own container and view type. Does not make use of conditions data
  /// ( @c empty_context )
  template <typename container_t = host_container_types>
  using accelerator_store = multi_store<
      accel_id, empty_context, dtuple,
      brute_force_collection<surface_type, container_t>,
      grid_collection<cylinder2D_sf_grid<surface_type, container_t>>,
      grid_collection<disc_sf_grid<surface_type, container_t>>,
      grid_collection<volume_accelerator<container_t>>>;
};

} // namespace detray
