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
#include "detray/definitions/containers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/mask.hpp"
#include "detray/geometry/shapes/rectangle2D.hpp"
#include "detray/geometry/shapes/trapezoid2D.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/io/backend/detail/type_info.hpp"  // mask_info
#include "detray/material/material_slab.hpp"
#include "detray/navigation/accelerators/brute_force.hpp"

// Linear algebra types
#include "detray/tutorial/types.hpp"

// New geometric shape type
#include "my_square2D.hpp"

/// This example defines a detray geometry type for a detector with cuboid
/// volumes that contain a new surface shape (squares), trapezoids and
/// rectangles for the volume boundary surfaces (portals).
/// For now, all surfaces in the volume(s) are stored in a 'brute force'
/// acceleration data structure which tests all surfaces it contains during
/// navigation.
/// In this example detector design, volumes do not contain other volumes, so
/// the volume lookup is done using a uniform grid.
/// Furthermore, the detector will contain homogeneous material on its surfaces.
///
/// If the new square shape should participate in the file IO, then the type and
/// a corresponding enum entry have to be added to the
/// detray/io/frontend/definitions.hpp header. Preferably appended to the end
/// of the existing structures, so that the written ids in existing files stay
/// valid
namespace detray::tutorial {

/// Defines a detector that contains squares, trapezoids and a bounding portal
/// box.
struct my_metadata {
  /// Define the algebra type for the geometry and navigation
  using algebra_type = detray::tutorial::algebra_t;

  /// Portal link type between volumes
  using nav_link = std::uint_least16_t;

  /// How to store and link transforms. The geometry context allows to resolve
  /// the conditions data for e.g. module alignment
  template <template <typename...> class vector_t = dvector>
  using transform_store = single_store<transform3, vector_t, geometry_context>;

  //
  // Surface Primitives
  //

  /// The mask types for the detector sensitive/passive surfaces
  using square = mask<square2D, algebra_type, nav_link>;
  using trapezoid = mask<trapezoid2D, algebra_type, nav_link>;
  // Types for portals
  using rectangle = mask<rectangle2D, algebra_type, nav_link>;

  /// Assign the mask types to the mask tuple container entries. It may be a
  /// good idea to have the most common types in the first tuple entries, in
  /// order to minimize the depth of the 'unrolling' before a mask is found
  /// in the tuple
  enum class mask_id : std::uint_least8_t {
    e_square2D = 0u,
    e_trapezoid2D = 1u,
    e_rectangle2D = 2u
  };

  friend std::ostream& operator<<(std::ostream& os, const mask_id& id) {
    switch (id) {
      case mask_id::e_square2D:
        os << "e_square2D";
        break;
      case mask_id::e_trapezoid2D:
        os << "e_trapezoid2D";
        break;
      case mask_id::e_rectangle2D:
        os << "e_rectangle2D";
        break;
      default:
        os << "Unknown mask_id";
        break;
    }
    return os;
  }

  /// This is the mask collections tuple (in the detector called 'mask store')
  /// the @c regular_multi_store is a vecemem-ready tuple of vectors of
  /// the detector masks.
  template <template <typename...> class vector_t = dvector>
  using mask_store =
      regular_multi_store<mask_id, empty_context, dtuple, vector_t, square,
                          trapezoid, rectangle>;

  //
  // Material Description
  //

  /// The material types to be mapped onto the surfaces: Here homogeneous
  /// material
  using slab = material_slab<scalar>;

  /// Similar to the mask store, there is a material store, which
  enum class material_id : std::uint_least8_t {
    e_material_slab = 0u,
    e_none = 1u,
  };

  friend std::ostream& operator<<(std::ostream& os, const material_id& id) {
    switch (id) {
      case material_id::e_material_slab:
        os << "e_material_slab";
        break;
      case material_id::e_none:
        os << "e_none";
        break;
      default:
        os << "Unknown material_id";
        break;
    }
    return os;
  }

  /// How to store and link materials. The material does not make use of
  /// conditions data ( @c empty_context )
  template <typename container_t = host_container_types>
  using material_store =
      multi_store<material_id, empty_context, dtuple,
                  typename container_t::template vector_type<slab>>;

  //
  // Acceleration structures
  //

  /// Portals and passives in the brute force search, sensitives in the grids
  enum geo_objects : std::uint_least8_t {
    e_portal = 0u,
    e_passive = 0u,
    e_sensitive = 1u,
    e_volume = 2u,
    e_size = 3u,
    e_all = e_size,
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
        // e_all has same value (2u)
        os << "e_size/e_all";
        break;
    }
    return os;
  }

  /// Surface descriptor type used for sensitives, passives and portals
  /// It holds the indices to the surface data in the detector data stores
  /// that were defined above
  using transform_link = typename transform_store<>::single_link;
  using mask_link = typename mask_store<>::single_link;
  using material_link = typename material_store<>::single_link;
  using surface_type =
      surface_descriptor<mask_link, material_link, transform_link, nav_link>;

  /// The acceleration data structures live in another tuple that needs to
  /// indexed correctly
  enum class accel_id : std::uint_least8_t {
    e_surface_brute_force = 0u,  //< test all surfaces in a volume (brute force)
    e_volume_cylinder3D_grid = 1u,
    e_surface_default = e_surface_brute_force,
    e_volume_default = e_volume_cylinder3D_grid,
  };

  friend std::ostream& operator<<(std::ostream& os, const accel_id& id) {
    switch (id) {
      case accel_id::e_surface_brute_force:
        os << "e_surface_brute_force/e_surface_default";
        break;
      default:
        os << "Unknown accel_id";
        break;
    }
    return os;
  }

  /// One link for portals/passives and one sensitive surfaces
  using object_link_type =
      dmulti_index<dtyped_index<accel_id, dindex>, geo_objects::e_size>;

  /// Data structure that allows to find the current detector volume from a
  /// given position. Here: Uniform grid with a 3D cylindrical shape
  template <typename container_t = host_container_types>
  using volume_accelerator =
      spatial_grid<algebra_type,
                   axes<cylinder3D, axis::bounds::e_open, axis::irregular,
                        axis::regular, axis::irregular>,
                   bins::single<dindex>, simple_serializer, container_t, false>;

  /// The tuple store that hold the acceleration data structures for all
  /// volumes. Every collection of accelerationdata structures defines its
  /// own container and view type. Does not make use of conditions data
  /// ( @c empty_context )
  /// How to store the acceleration data structures
  template <typename container_t = host_container_types>
  using accelerator_store =
      multi_store<accel_id, empty_context, dtuple,
                  brute_force_collection<surface_type, container_t>,
                  grid_collection<volume_accelerator<container_t>>>;
};

}  // namespace detray::tutorial
