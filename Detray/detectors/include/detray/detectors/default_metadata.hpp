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
#include "detray/geometry/shapes/cuboid3D.hpp"
#include "detray/geometry/shapes/cylinder2D.hpp"
#include "detray/geometry/shapes/cylinder3D.hpp"
#include "detray/geometry/shapes/line.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/geometry/shapes/ring2D.hpp"
#include "detray/geometry/shapes/trapezoid2D.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/material/material.hpp"
#include "detray/material/material_map.hpp"
#include "detray/material/material_rod.hpp"
#include "detray/material/material_slab.hpp"
#include "detray/navigation/accelerators/brute_force.hpp"
#include "detray/navigation/accelerators/spatial_grid.hpp"

namespace detray {

template <concepts::algebra algebra_t>
struct default_metadata {
  using algebra_type = algebra_t;
  using scalar_t = dscalar<algebra_type>;

  using nav_link = std::uint_least16_t;

  template <template <typename...> class vector_t = dvector>
  using transform_store =
      single_store<dtransform3D<algebra_type>, vector_t, geometry_context>;

  using rectangle2D_t = mask<detray::rectangle2D, algebra_type, nav_link>;
  using drift_cell_t = mask<detray::line_square, algebra_type, nav_link>;
  using trapezoid2D_t = mask<detray::trapezoid2D, algebra_type, nav_link>;
  using straw_tube_t = mask<detray::line_circular, algebra_type, nav_link>;
  using concentric_cylinder2D_t =
      mask<detray::concentric_cylinder2D, algebra_type, nav_link>;
  using ring2D_t = mask<detray::ring2D, algebra_type, nav_link>;
  using annulus2D_t = mask<detray::annulus2D, algebra_type, nav_link>;
  using cylinder2D_t = mask<detray::cylinder2D, algebra_type, nav_link>;

  enum class mask_id : std::uint_least8_t {
    e_rectangle2D = 0u,
    e_drift_cell = 1u,
    e_trapezoid2D = 2u,
    e_straw_tube = 3u,
    e_concentric_cylinder2D = 4u,
    e_ring2D = 5u,
    e_annulus2D = 6u,
    e_cylinder2D = 7u,
  };

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     mask_id id) {
    switch (id) {
      case mask_id::e_rectangle2D:
        os << "e_rectangle2D";
        break;
      case mask_id::e_drift_cell:
        os << "e_drift_cell";
        break;
      case mask_id::e_trapezoid2D:
        os << "e_trapezoid2D";
        break;
      case mask_id::e_straw_tube:
        os << "e_straw_tube";
        break;
      case mask_id::e_concentric_cylinder2D:
        os << "e_concentric_cylinder2D";
        break;
      case mask_id::e_ring2D:
        os << "e_ring2D";
        break;
      case mask_id::e_annulus2D:
        os << "e_annulus2D";
        break;
      case mask_id::e_cylinder2D:
        os << "e_cylinder2D";
        break;
      default:
        os << "invalid";
    }
    return os;
  };

  template <template <typename...> class vector_t = dvector>
  using mask_store =
      regular_multi_store<mask_id, empty_context, dtuple, vector_t,
                          rectangle2D_t, drift_cell_t, trapezoid2D_t,
                          straw_tube_t, concentric_cylinder2D_t, ring2D_t,
                          annulus2D_t, cylinder2D_t>;

  template <typename container_t>
  using rectangle2D_map_t =
      material_map<algebra_type, detray::rectangle2D, container_t>;
  using material_slab_t = material_slab<scalar_t>;
  template <typename container_t>
  using concentric_cylinder2D_map_t =
      material_map<algebra_type, detray::concentric_cylinder2D, container_t>;
  template <typename container_t>
  using ring2D_map_t = material_map<algebra_type, detray::ring2D, container_t>;
  using raw_material_t = material<scalar_t>;
  using material_rod_t = material_rod<scalar_t>;
  template <typename container_t>
  using cylinder3D_map_t =
      material_map<algebra_type, detray::cylinder3D, container_t>;
  template <typename container_t>
  using cylinder2D_map_t =
      material_map<algebra_type, detray::cylinder2D, container_t>;
  template <typename container_t>
  using cuboid3D_map_t =
      material_map<algebra_type, detray::cuboid3D, container_t>;

  enum class material_id : std::uint_least8_t {
    e_rectangle2D_map = 0u,
    e_trapezoid2D_map = 0u,
    e_material_slab = 1u,
    e_concentric_cylinder2D_map = 2u,
    e_ring2D_map = 3u,
    e_annulus2D_map = 3u,
    e_raw_material = 4u,
    e_material_rod = 5u,
    e_cylinder3D_map = 6u,
    e_cylinder2D_map = 7u,
    e_cuboid3D_map = 8u,
    e_none = 9u,
  };

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     material_id id) {
    switch (id) {
      case material_id::e_trapezoid2D_map:
        os << "e_rectangle2D_map/e_trapezoid2D_map";
        break;
      case material_id::e_material_slab:
        os << "e_material_slab";
        break;
      case material_id::e_concentric_cylinder2D_map:
        os << "e_concentric_cylinder2D_map";
        break;
      case material_id::e_annulus2D_map:
        os << "e_ring2D_map/e_annulus2D_map";
        break;
      case material_id::e_raw_material:
        os << "e_raw_material";
        break;
      case material_id::e_material_rod:
        os << "e_material_rod";
        break;
      case material_id::e_cylinder3D_map:
        os << "e_cylinder3D_map";
        break;
      case material_id::e_cylinder2D_map:
        os << "e_cylinder2D_map";
        break;
      case material_id::e_cuboid3D_map:
        os << "e_cuboid3D_map";
        break;
      case material_id::e_none:
        os << "e_none";
        break;
      default:
        os << "invalid";
    }
    return os;
  };

  template <typename container_t = host_container_types>
  using material_store =
      multi_store<material_id, empty_context, dtuple,
                  grid_collection<rectangle2D_map_t<container_t>>,
                  typename container_t::template vector_type<material_slab_t>,
                  grid_collection<concentric_cylinder2D_map_t<container_t>>,
                  grid_collection<ring2D_map_t<container_t>>,
                  typename container_t::template vector_type<raw_material_t>,
                  typename container_t::template vector_type<material_rod_t>,
                  grid_collection<cylinder3D_map_t<container_t>>,
                  grid_collection<cylinder2D_map_t<container_t>>,
                  grid_collection<cuboid3D_map_t<container_t>>>;

  using transform_link = typename transform_store<>::single_link;
  using mask_link = typename mask_store<>::range_link;
  using material_link = typename material_store<>::single_link;
  using surface_type =
      surface_descriptor<mask_link, material_link, transform_link, nav_link>;

  template <typename axes_t, typename bin_entry_t, typename container_t>
  using dynamic_simple_grid_t =
      spatial_grid<algebra_type, axes_t,
                   detray::bins::dynamic_array<bin_entry_t>,
                   detray::simple_serializer, container_t, false>;

  template <typename container_t>
  using surface_concentric_cylinder2D_grid_t =
      dynamic_simple_grid_t<axes<detray::concentric_cylinder2D>, surface_type,
                            container_t>;
  template <typename container_t>
  using surface_ring2D_grid_t =
      dynamic_simple_grid_t<axes<detray::ring2D>, surface_type, container_t>;
  template <typename container_t>
  using surface_cylinder2D_grid_t =
      dynamic_simple_grid_t<axes<detray::cylinder2D>, surface_type,
                            container_t>;

  template <typename axes_t, typename bin_entry_t, typename container_t>
  using single_simple_grid_t =
      spatial_grid<algebra_type, axes_t, detray::bins::single<bin_entry_t>,
                   detray::simple_serializer, container_t, false>;

  template <typename container_t>
  using volume_cylinder3D_grid_t =
      single_simple_grid_t<axes<detray::cylinder3D>, dindex, container_t>;

  enum class accel_id : std::uint_least8_t {
    e_surface_brute_force = 0u,
    e_surface_concentric_cylinder2D_grid = 1u,
    e_surface_ring2D_grid = 2u,
    e_surface_cylinder2D_grid = 3u,
    e_volume_brute_force = 4u,
    e_volume_cylinder3D_grid = 5u,
    e_surface_default = e_surface_brute_force,
    e_volume_default = e_volume_cylinder3D_grid,
  };

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     accel_id id) {
    switch (id) {
      case accel_id::e_surface_brute_force:
        os << "e_surface_brute_force";
        break;
      case accel_id::e_surface_concentric_cylinder2D_grid:
        os << "e_surface_concentric_cylinder2D_grid";
        break;
      case accel_id::e_surface_ring2D_grid:
        os << "e_surface_ring2D_grid";
        break;
      case accel_id::e_surface_cylinder2D_grid:
        os << "e_surface_cylinder2D_grid";
        break;
      case accel_id::e_volume_brute_force:
        os << "e_volume_brute_force";
        break;
      case accel_id::e_volume_cylinder3D_grid:
        os << "e_volume_cylinder3D_grid";
        break;
      default:
        os << "invalid";
    }
    return os;
  };

  template <typename container_t = host_container_types>
  using accelerator_store = multi_store<
      accel_id, empty_context, dtuple,
      brute_force_collection<surface_type, container_t>,
      grid_collection<surface_concentric_cylinder2D_grid_t<container_t>>,
      grid_collection<surface_ring2D_grid_t<container_t>>,
      grid_collection<surface_cylinder2D_grid_t<container_t>>,
      grid_collection<volume_cylinder3D_grid_t<container_t>>,
      brute_force_collection<dindex, container_t>>;

  enum geo_objects : std::uint_least8_t {
    e_passive = 0u,
    e_portal = 0u,
    e_sensitive = 1u,
    e_volume = 2u,
    e_size = 3u,
    e_all = e_size,
  };

  DETRAY_HOST inline friend std::ostream& operator<<(std::ostream& os,
                                                     geo_objects id) {
    switch (id) {
      case geo_objects::e_sensitive:
        os << "e_sensitive";
        break;
      case geo_objects::e_volume:
        os << "e_volume";
        break;
      case geo_objects::e_portal:
        os << "e_passive/e_portal";
        break;
      case geo_objects::e_size:
        os << "e_size";
        break;
      default:
        os << "invalid";
    }
    return os;
  };

  using object_link_type =
      dmulti_index<dtyped_index<accel_id, dindex>, geo_objects::e_size>;
};

}  // namespace detray
