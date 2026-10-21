// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/concepts.hpp"
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/core/detail/type_traits.hpp"

// System include(s)
#include <concepts>
#include <type_traits>

namespace detray::concepts {

/// Check if a type can be used as a metadata
template <typename M>
concept metadata = requires {
  requires algebra<typename M::algebra_type>;
  requires std::is_enum_v<typename M::mask_id>;
  requires std::is_enum_v<typename M::accel_id>;
  requires std::is_enum_v<typename M::material_id>;
  requires std::is_enum_v<typename M::geo_objects>;

  typename M::nav_link;
  typename M::transform_link;
  typename M::mask_link;
  typename M::material_link;
  typename M::surface_type;
  typename M::object_link_type;

  typename M::template transform_store<dvector>;
  typename M::template mask_store<dvector>;
  typename M::template material_store<host_container_types>;
  typename M::template accelerator_store<host_container_types>;
};

/// Check if a type can be used as a detector
template <typename D>
concept detector = requires(const D& d) {
  requires viewable<D>;
  requires bufferable<D>;

  requires metadata<typename D::metadata>;
  requires algebra<typename D::algebra_type>;
  requires scalar<typename D::scalar_type>;
  requires point2D<typename D::point2_type>;
  requires point3D<typename D::point3_type>;
  requires vector3D<typename D::vector3_type>;
  requires transform3D<typename D::transform3_type>;

  requires std::is_enum_v<typename D::geo_obj_ids>;

  typename D::name_map;

  typename D::volume_type;
  typename D::volume_container;

  typename D::surface_type;
  typename D::surface_container;
  typename D::surface_lookup_container;

  typename D::transform_container;
  typename D::geometry_context;

  typename D::mask_container;
  typename D::masks;

  typename D::material_container;
  typename D::material_link;
  typename D::material;

  typename D::accelerator_container;
  typename D::accel_link;
  typename D::accel;

  { d.name(typename D::name_map()) } -> std::same_as<std::string>;

  { d.volumes() } -> std::same_as<const typename D::volume_container&>;

  { d.volume(dindex()) } -> std::same_as<const typename D::volume_type&>;

  {
    d.volume(std::string_view(), typename D::name_map())
  } -> std::same_as<const typename D::volume_type&>;

  /// @TODO: Not yet implemented correctly
  /*{
    d.volume(typename D::point3_type())
  } -> std::same_as<const typename D::volume_type&>;*/

  { d.portals() };

  { d.surfaces() } -> std::same_as<const typename D::surface_lookup_container&>;

  {
    d.transform_store(typename D::geometry_context())
  } -> std::same_as<const typename D::transform_container&>;

  { d.mask_store() } -> std::same_as<const typename D::mask_container&>;

  { d.material_store() } -> std::same_as<const typename D::material_container&>;

  {
    d.accelerator_store()
  } -> std::same_as<const typename D::accelerator_container&>;
};

/// Check if a type is a host (mutable) detector type
template <typename D>
concept host_detector =
    detector<D> &&
    std::same_as<typename D::container_types, host_container_types>;

/// Check if a type is a device (immutable) detector type
template <typename D>
concept device_detector =
    detector<D> &&
    std::same_as<typename D::container_types, const_device_container_types>;

/// Check for the the presence of any type of grids in a detector definition
template <class D>
concept has_grids =
    detector<D> && (detail::contains_grids_v<typename D::accel> ||
                    detail::contains_grids_v<typename D::material>);

/// Check for the the presence of surface grids in a detector definition
template <class D>
concept has_surface_grids =
    detector<D> && detail::contains_surface_grids_v<typename D::accel>;

/// Check for the the presence of material slabs in a detector definition
template <class D>
concept has_material_slabs =
    detector<D> && detail::contains_material_slabs_v<typename D::material>;

/// Check for the the presence of material rods in a detector definition
template <class D>
concept has_material_rods =
    detector<D> && detail::contains_material_rods_v<typename D::material>;

/// Check for the the presence of homogeneous material types in a detector
/// definition
template <class D>
concept has_homogeneous_material =
    detector<D> &&
    detail::contains_homogeneous_material_v<typename D::material>;

/// Check for the the presence of material maps in a detector definition
template <class D>
concept has_material_maps =
    detector<D> && detail::contains_material_maps_v<typename D::material>;

/// Check that a type is a Draits instance
template <typename T>
concept detector_traits = requires {
  requires metadata<typename T::metadata_type>;
  requires device_view<typename T::view>;
  requires device_buffer<typename T::buffer>;

  requires detector<typename T::host>;
  requires detector<typename T::device>;
};

}  // namespace detray::concepts
