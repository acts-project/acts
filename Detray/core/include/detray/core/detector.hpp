// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/grid_builder.hpp"
#include "detray/builders/homogeneous_material_builder.hpp"
#include "detray/builders/homogeneous_volume_material_builder.hpp"
#include "detray/builders/material_map_builder.hpp"
#include "detray/builders/volume_builder.hpp"
#include "detray/core/concepts.hpp"
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/core/detail/surface_lookup.hpp"
#include "detray/core/name_map.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/containers.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/geometry/volume_descriptor.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/print_detector.hpp"

// Vecmem include(s)
#include <vecmem/memory/memory_resource.hpp>

// System include(s)
#include <sstream>
#include <string>
#include <string_view>

namespace detray {

namespace detail {
/// Temporary way to manipulate transforms in the transform store
/// @todo Remove as soon as contices can be registered!
template <typename detector_t, concepts::transform3D transform3_t>
void set_transform(detector_t &det, const transform3_t &trf, unsigned int i) {
  DETRAY_WARN_HOST(
      "Modifying transforms in the detector will be deprecated! "
      "Please, use a separate geometry context in this case");
  det._transforms.at(i) = trf;
}

/// @brief The detector definition.
///
/// This class is a heavily templated container aggregation, that owns all data.
/// Its view type is used to move the data between host and device.
/// It is constructed by the @c detector_builder and is afterwards immutable.
///
/// @tparam metadata helper that defines collection and link types centrally
/// @tparam container_t type collection of the underlying containers
template <concepts::metadata metadata_t, typename container_t>
class detector {
  // Allow the building of the detector containers
  friend class detray::volume_builder<detector<metadata_t, container_t>>;
  template <concepts::detector, concepts::grid, typename, typename>
  friend class detray::grid_builder;
  friend class detray::homogeneous_material_builder<
      detector<metadata_t, container_t>>;
  friend class detray::homogeneous_volume_material_builder<
      detector<metadata_t, container_t>>;
  template <concepts::detector, std::size_t, typename>
  friend class detray::material_map_builder;
  /// @todo Remove
  friend void set_transform<detector<metadata_t, container_t>,
                            dtransform3D<typename metadata_t::algebra_type>>(
      detector<metadata_t, container_t> &,
      const dtransform3D<typename metadata_t::algebra_type> &, unsigned int);

  /// Raw container types
  template <typename T>
  using vector_type = typename container_t::template vector_type<T>;

 public:
  /// Main definition of geometry types
  /// Metadata policy defining detector component and link types.
  using metadata = metadata_t;

  /// Algebra types
  /// Algebra backend selected by the metadata policy.
  using algebra_type = typename metadata::algebra_type;
  /// Scalar type associated with the algebra backend.
  using scalar_type = dscalar<algebra_type>;
  /// Two-dimensional point type for local surface coordinates.
  using point2_type = dpoint2D<algebra_type>;
  /// Three-dimensional point type for global positions.
  using point3_type = dpoint3D<algebra_type>;
  /// Three-dimensional vector type for directions.
  using vector3_type = dvector3D<algebra_type>;
  /// Three-dimensional transform type for placements.
  using transform3_type = dtransform3D<algebra_type>;

  /// In case the detector needs to be printed
  /// Lookup table mapping detector names to component indices.
  using name_map = detray::name_map;

  /// The surface takes a mask (defines the local coordinates and the surface
  /// extent), its material, a link to an element in the transform container
  /// to define its placement and a source link to the object it represents.
  /// Surface descriptor type used by this detector.
  using surface_type = typename metadata::surface_type;
  /// Container storing the detector surfaces.
  using surface_container = vector_type<surface_type>;
  /// Lookup structure for resolving surface identifiers and source links.
  using surface_lookup_container = surface_lookup<surface_type, vector_type>;

  /// Forward the alignable transform container (surface placements) and
  /// the geo context (e.g. for alignment)
  /// Transform store used for surface placements and geometry contexts.
  using transform_container =
      typename metadata::template transform_store<vector_type>;
  /// Context type accepted by the transform store.
  using geometry_context = typename transform_container::context_type;

  /// Forward mask types that are present in this detector
  /// Collection of the surface mask types supported by this detector.
  using mask_container = typename metadata::template mask_store<vector_type>;
  /// Variant of mask types supported by this detector.
  using masks = typename mask_container::value_types;

  /// Forward material types that are present in this detector
  /// Collection of the material types supported by this detector.
  using material_container =
      typename metadata::template material_store<container_t>;
  /// Variant of material types supported by this detector.
  using material = typename material_container::value_types;

  /// Surface Finders: structures that enable neighborhood searches in the
  /// detector geometry during navigation. Can be different in each volume
  /// Collection of navigation accelerators supported by this detector.
  using accelerator_container =
      typename metadata::template accelerator_store<container_t>;
  /// Variant of accelerator types supported by this detector.
  using accel = typename accelerator_container::value_types;

  /// Volume type
  /// Identifiers used to classify detector geometry objects.
  using geo_obj_ids = typename metadata::geo_objects;
  /// Link type for entries in the material store.
  using material_link = typename material_container::single_link;
  /// Link type for entries in the accelerator store.
  using accel_link = typename accelerator_container::single_link;
  /// Descriptor type for a tracking volume.
  using volume_type = volume_descriptor<geo_obj_ids, accel_link, material_link>;
  /// Container storing tracking-volume descriptors.
  using volume_container = vector_type<volume_type>;

  /// Detector view types
  /// Non-owning const view over all detector component stores.
  using const_view_type =
      dmulti_view<dvector_view<const volume_type>,
                  typename surface_lookup_container::const_view_type,
                  typename transform_container::const_view_type,
                  typename mask_container::const_view_type,
                  typename material_container::const_view_type,
                  typename accelerator_container::const_view_type>;

  /// View type exposed by the detector.

  using view_type = const_view_type;

  static_assert(concepts::device_view<view_type>,
                "Detector view type ill-formed");
  static_assert(concepts::device_view<const_view_type>,
                "Detector const view type ill-formed");

  /// Detector buffer types
  /// Buffer type aggregating device-transfer buffers for detector stores.
  using buffer_type =
      dmulti_buffer<dvector_buffer<volume_type>,
                    typename surface_lookup_container::buffer_type,
                    typename transform_container::buffer_type,
                    typename mask_container::buffer_type,
                    typename material_container::buffer_type,
                    typename accelerator_container::buffer_type>;

  static_assert(concepts::device_buffer<buffer_type>,
                "Detector buffer type ill-formed");

  detector() = delete;
  // The detector holds a lot of data and should never be copied
  detector(const detector &) = delete;
  detector &operator=(const detector &) = delete;

  /// Allowed constructors
  /// @{
  /// Move constructor
  detector(detector &&) noexcept = default;

  /// Move-assign ownership of the detector stores.
  /// @param other detector whose stores are transferred
  /// @return reference to this detector
  detector &operator=(detector &&other) noexcept = default;

  /// Default construction
  /// @param resource memory resource for the allocation of members
  DETRAY_HOST
  explicit detector(vecmem::memory_resource &resource)
      : _volumes(&resource),
        _surfaces(resource),
        _transforms(resource),
        _masks(resource),
        _materials(resource),
        _accelerators(resource) {}

  /// Constructor from detector data view
  /// @tparam detector_view_t device-compatible detector view type
  /// @param det_data view containing the detector component stores
  template <concepts::device_view detector_view_t>
  DETRAY_HOST_DEVICE explicit detector(detector_view_t &det_data)
      : _volumes(detray::detail::get<0>(det_data.m_view)),
        _surfaces(detray::detail::get<1>(det_data.m_view)),
        _transforms(detray::detail::get<2>(det_data.m_view)),
        _masks(detray::detail::get<3>(det_data.m_view)),
        _materials(detray::detail::get<4>(det_data.m_view)),
        _accelerators(detray::detail::get<5>(det_data.m_view)) {}
  /// @}

  /// Return the detector name from the supplied name map.
  /// @param names map containing the detector name
  /// @returns a string that contains the detector name
  std::string name(const name_map &names) const {
    return names.get_detector_name();
  }

  /// @returns the sub-volumes of the detector - const access
  DETRAY_HOST_DEVICE
  inline auto volumes() const -> const vector_type<volume_type> & {
    return _volumes;
  }

  /// Return a volume by its stored index.
  /// @param volume_index index in the volume container
  /// @returns the volume by @p volume_index - const access
  DETRAY_HOST_DEVICE
  inline const auto &volume(dindex volume_index) const {
    return _volumes[volume_index];
  }

  /// Return a volume by its registered name.
  /// @param volume_name registered volume name
  /// @param names map used to resolve volume names
  /// @returns the volume by @p volume_name - const access
  DETRAY_HOST
  inline const auto &volume(const std::string_view volume_name,
                            const name_map &names) const {
    return _volumes.at(names.at(volume_name));
  }

  /// Find the volume containing a global point.
  /// @param p global Cartesian position to locate
  /// @return the volume by global cartesian @p p - const access
  DETRAY_HOST_DEVICE
  inline const auto &volume(const point3_type &p) const {
    // Allow to call the volume search data structure
    // TODO: Add volume accelerator builder
    volume_type v_desc{};
    v_desc.template set_accel_link<geo_obj_ids::e_volume>(
        accel::id::e_volume_default, 0u);
    tracking_volume world{*this, v_desc};

    dindex volume_index{0u};
    world.template visit_accelerator<geo_obj_ids::e_volume, volume_search>(
        p, &volume_index);
    return _volumes[volume_index];
  }

  /// @returns all portals - const
  /// @note Depending on the detector type, this can also contain other
  /// surfaces
  /// @todo add range filter to skip non-portal surfaces
  DETRAY_HOST_DEVICE
  inline const auto &portals() const {
    // All portals are registered with the brute force search
    return _accelerators.template get<accel::id::e_surface_brute_force>().all();
  }

  /// @returns the sub-volumes of the detector - const access
  DETRAY_HOST_DEVICE
  inline auto surfaces() const -> const surface_lookup_container & {
    return _surfaces;
  }

  /// Find a surface using an index, identifier, or source-link query.
  /// @param q query forwarded to the surface lookup
  /// @returns a surface using a query object @p q. This can be an index,
  /// a identifier or a source link searcher (see @c surface_lookup class)
  template <typename query_t>
  DETRAY_HOST_DEVICE constexpr decltype(auto) surface(query_t &&q) const {
    return _surfaces.search(std::forward<query_t>(q));
  }

  /// @returns detector transform store
  DETRAY_HOST_DEVICE
  inline auto transform_store(const geometry_context & /*ctx*/ = {}) const
      -> const transform_container & {
    return _transforms;
  }

  /// @returns all surface/portal masks in the geometry - const access
  DETRAY_HOST_DEVICE
  inline auto mask_store() const -> const mask_container & { return _masks; }

  /// @returns all materials in the geometry - const access
  DETRAY_HOST_DEVICE
  inline auto material_store() const -> const material_container & {
    return _materials;
  }

  /// @returns access to the surface finder container
  DETRAY_HOST_DEVICE
  inline auto accelerator_store() const -> const accelerator_container & {
    return _accelerators;
  }

  /// @returns view of a detector - vecmem auto-conversion to const view types
  DETRAY_HOST auto get_data() -> const_view_type {
    return const_view_type{detray::get_data(std::as_const(_volumes)),
                           detray::get_data(std::as_const(_surfaces)),
                           detray::get_data(std::as_const(_transforms)),
                           detray::get_data(std::as_const(_masks)),
                           detray::get_data(std::as_const(_materials)),
                           detray::get_data(std::as_const(_accelerators))};
  }

  /// @returns const view of a detector
  DETRAY_HOST auto get_data() const -> const_view_type {
    return const_view_type{
        detray::get_data(_volumes),    detray::get_data(_surfaces),
        detray::get_data(_transforms), detray::get_data(_masks),
        detray::get_data(_materials),  detray::get_data(_accelerators)};
  }

 private:
  /// @returns a string stream that prints the detector details (but without
  /// volumes names)
  DETRAY_HOST
  friend std::ostream &operator<<(std::ostream &os, const detector &d) {
    os << utils::print_detector(d) << std::endl;
    return os;
  }

  /// Volume lookup in the volume acceleration data structures
  struct volume_search {
    /// @todo Move this visitor into a volume-search grid type.
    template <concepts::accelerator_collection accel_coll_t,
              typename accel_index_t>
    DETRAY_HOST_DEVICE inline void operator()(const accel_coll_t &coll,
                                              const accel_index_t index,
                                              const point3_type &p,
                                              dindex *const result) const {
      using accel_type = typename accel_coll_t::value_type;

      if constexpr (concepts::volume_accelerator<accel_type>) {
        const auto volume_accelerator = coll[index];

        // Only one entry per bin
        *result = volume_accelerator.search(p).value();
      }
    }
  };

  /// Contains the detector sub-volumes.
  volume_container _volumes;

  /// Lookup for surfaces from identifiers
  surface_lookup_container _surfaces;

  /// Keeps all of the transform data in contiguous memory
  transform_container _transforms;

  /// Masks of all surfaces in the geometry in contiguous memory
  mask_container _masks;

  /// Materials of all surfaces in the geometry in contiguous memory
  material_container _materials;

  /// All surface finder data structures that are used in the detector volumes
  accelerator_container _accelerators;
};

}  // namespace detail

namespace host {

template <concepts::metadata metadata_t>
using detector = detray::detail::detector<metadata_t, host_container_types>;

}  // namespace host

namespace device {

template <concepts::metadata metadata_t>
using detector =
    detray::detail::detector<metadata_t, const_device_container_types>;

}  // namespace device

template <typename T>
struct detector_traits {};

/// Unified interface to detector types - from a metadata type
template <concepts::metadata metadata_t>
struct detector_traits<metadata_t> {
  /// Metadata type of the detector.
  /// Metadata type used to configure a detector.
  using metadata_type = metadata_t;

  /// Host type of the detector.
  /// Host detector type.
  using host = detray::host::detector<metadata_type>;

  /// Device type of the detector.
  /// Device detector type.
  using device = detray::device::detector<metadata_type>;

  /// View of the detector.
  /// Non-owning detector view type.
  using view = typename host::const_view_type;

  /// Buffer for a detector's data.
  /// Device-transfer buffer type for detector data.
  using buffer = typename host::buffer_type;
};

/// Unified interface to detector types - from any detector type
template <concepts::metadata metadata_t, typename container_t>
struct detector_traits<detail::detector<metadata_t, container_t>> {
 private:
  using T = detector_traits<metadata_t>;

 public:
  /// Metadata type used to configure a detector.
  using metadata_type = typename T::metadata_type;
  /// Host detector type.
  using host = typename T::host;
  /// Device detector type.
  using device = typename T::device;
  /// Non-owning detector view type.
  using view = typename T::view;
  /// Device-transfer buffer type for detector data.
  using buffer = typename T::buffer;
};

/// Some predefined alias types
/// @{

/// Metadata type associated with detector type D.
template <typename D>
using detector_metadata_t = typename detector_traits<D>::metadata_type;

/// Host detector type for metadata type D.
template <typename D>
using detector_host_t = typename detector_traits<D>::host;

/// Device detector type for metadata type D.
template <typename D>
using detector_device_t = typename detector_traits<D>::device;

/// Const view type for detector type D.
template <typename D>
using detector_view_t = typename detector_traits<D>::view;

/// Buffer type for transferring detector type D.
template <typename D>
using detector_buffer_t = typename detector_traits<D>::buffer;

/// @}

}  // namespace detray
