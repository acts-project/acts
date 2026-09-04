// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"

// System include(s)
#include <ostream>
#include <utility>

namespace detray {

/// @brief The detray detector volume descriptor.
///
/// Contains the data and links to describe a detector volume. This is the type
/// that is stored in the detector data stores.
///
/// @tparam ID enum of object types contained in the volume
///         (@see @c default_metadata ).
/// @tparam link_t the type of link to the volumes surfaces finder(s)
///         (accelerator structure, e.g. a grid). The surface finder types
///         cannot be given directly, since the containers differ between host
///         and device. The surface finders reside in an 'unrollable tuple
///         container' and are called per volume in the navigator during local
///         navigation.
template <typename ID, typename acc_link_t = dtyped_index<dindex, dindex>,
          typename mat_link_t = dtyped_index<dindex, dindex>>
class volume_descriptor {
 public:
  /// Ids of objects that can be distinguished by the volume
  using object_id = ID;

  /// How to access the surface ranges in the detector surface lookup
  using sf_link_type =
      dmulti_index<dindex_range,
                   static_cast<std::size_t>(
                       static_cast<std::underlying_type_t<surface_id>>(
                           surface_id::e_all))>;

  /// How to access objects (e.g. sensitives/passives/portals) in this
  /// volume. Keeps one accelerator structure link per object type (by ID):
  ///
  /// acc_link_t : id and index of the accelerator structure in the detector's
  ///          surface store.
  ///
  /// E.g. a 'portal' can be found under @c ID::e_portal in this link,
  /// and will then receive link to the @c brute_force acceleration structure
  /// that holds the portals (the accelerator structure's id and index).
  using accel_link_type = dmulti_index<acc_link_t, ID::e_size>;

  /// How to link to the volume material, if any
  using material_link = mat_link_t;
  using material_id = typename material_link::id_type;

  /// Default constructor builds an ~infinitely long cylinder
  constexpr volume_descriptor() = default;

  /// Constructor from shape id.
  ///
  /// @param id id values that determines how to interpret the bounds.
  explicit constexpr volume_descriptor(const volume_id id) : m_id{id} {}

  /// @returns true if the object ID corresponds to a surface
  static consteval bool is_surface_id(const object_id id) {
    return (id == object_id::e_portal || id == object_id::e_sensitive ||
            id == object_id::e_passive);
  }

  /// @returns true if the object ID corresponds to a [daughter] volume
  static consteval bool is_volume_id(const object_id id) {
    return (id == object_id::e_volume);
  }

  /// @returns the volume shape id, e.g. 'cylinder'
  DETRAY_HOST_DEVICE
  constexpr auto id() const -> volume_id { return m_id; }

  /// @returns the index of the volume in the detector volume container.
  DETRAY_HOST_DEVICE
  constexpr auto index() const -> dindex { return m_index; }

  /// @param index the index of the volume in the detector volume container.
  DETRAY_HOST
  constexpr auto set_index(const dindex index) -> void { m_index = index; }

  /// @returns the index of the volume transform in the transform store.
  DETRAY_HOST_DEVICE
  constexpr auto transform() const -> dindex { return m_transform; }

  /// @param index the index of the volume in the detector volume container.
  DETRAY_HOST
  constexpr auto set_transform(const dindex trf_idx) -> void {
    m_transform = trf_idx;
  }

  /// @returns surface link for all object types - const
  DETRAY_HOST_DEVICE constexpr auto sf_link() const -> const sf_link_type& {
    return m_sf_links;
  }

  /// @returns surface descriptor link for a specific type of object - const
  template <surface_id id>
  DETRAY_HOST_DEVICE constexpr auto sf_link() const -> const
      typename sf_link_type::index_type& {
    return detail::get<static_cast<uint>(id)>(m_sf_links);
  }

  /// @returns surface descriptor link for a specific type of object
  template <surface_id id>
  DETRAY_HOST_DEVICE constexpr auto sf_link() ->
      typename sf_link_type::index_type& {
    return detail::get<static_cast<uint>(id)>(m_sf_links);
  }

  /// @returns surface descriptor link for all surface types
  DETRAY_HOST_DEVICE constexpr auto full_sf_range() const ->
      typename sf_link_type::index_type {
    using idx_range_t = typename sf_link_type::index_type;

    const auto& pt_range = sf_link<surface_id::e_portal>();
    const auto& sen_range = sf_link<surface_id::e_sensitive>();
    const auto& psv_range = sf_link<surface_id::e_passive>();

    // Portal range is never empty
    dindex min{detail::get<0>(pt_range)};
    dindex max{detail::get<1>(pt_range)};

    constexpr idx_range_t empty{};
    if (sen_range != empty) {
      min = detail::get<0>(sen_range) < min ? detail::get<0>(sen_range) : min;
      max = detail::get<1>(sen_range) > max ? detail::get<1>(sen_range) : max;
    }

    if (psv_range != empty) {
      min = detail::get<0>(psv_range) < min ? detail::get<0>(psv_range) : min;
      max = detail::get<1>(psv_range) > max ? detail::get<1>(psv_range) : max;
    }

    return idx_range_t{min, max};
  }

  /// @returns the total number of surfaces contained in the volume
  DETRAY_HOST_DEVICE constexpr dindex n_surfaces() const {
    const auto sf_idx_range = full_sf_range();

    assert(sf_idx_range[1] > sf_idx_range[0]);
    return sf_idx_range[1] - sf_idx_range[0];
  }

  /// @returns a surface index with respect to the volume surface range
  DETRAY_HOST_DEVICE constexpr dindex to_local_sf_index(dindex sf_idx) const {
    auto full_range = full_sf_range();

    assert(full_range[0] <= sf_idx);
    assert(sf_idx < full_range[1]);

    return sf_idx - full_range[0];
  }

  /// @returns a surface index with respect to the global detector containers
  DETRAY_HOST_DEVICE constexpr dindex to_global_sf_index(dindex sf_idx) const {
    auto full_range = full_sf_range();
    dindex glob_index{sf_idx + full_range[0]};

    assert(full_range[0] <= glob_index);
    assert(glob_index < full_range[1]);

    return glob_index;
  }

  /// Set or update the index into a geometry container identified by the
  /// obj_id.
  ///
  /// @note There is no check of overlapping index ranges between the object
  /// types. Use with care!
  ///
  /// @param other Surface index range
  template <surface_id id>
  DETRAY_HOST auto update_sf_link(
      const typename sf_link_type::index_type& other) noexcept -> void {
    auto& rg = sf_link<id>();
    // Range not set yet - initialize
    if (constexpr typename sf_link_type::index_type empty{}; rg == empty) {
      rg = other;
    } else {
      // Update upper border
      assert(detail::get<1>(rg) == detail::get<0>(other));
      rg.set_upper(other.upper());
    }
  }

  /// Set or update the index into a geometry container identified by the
  /// obj_id.
  ///
  /// @note There is no check of overlapping index ranges between the object
  /// types. Use with care!
  ///
  /// @param shift shift of the surface range in a larger container.
  /// @param n_surfaces the number of surfaces in this range.
  template <surface_id id>
  DETRAY_HOST auto update_sf_link(std::size_t shift,
                                  std::size_t n_surfaces = 0) noexcept -> void {
    auto& rg = sf_link<id>();
    // Range not set yet - initialize
    if (constexpr typename sf_link_type::index_type empty{}; rg == empty) {
      rg = {0u, static_cast<dindex>(n_surfaces)};
    }
    // Update
    rg.shift(static_cast<dindex>(shift));
  }

  /// @returns the volume material link
  DETRAY_HOST_DEVICE
  constexpr auto material() const -> const material_link& { return m_mat_link; }

  /// Set the volume material link to @param mat_link
  DETRAY_HOST
  constexpr auto set_material(const material_link& mat_link) -> void {
    m_mat_link = mat_link;
  }

  /// Set the volume material link using the given material type id
  /// @param mat_id and index of the material instance @param mat_idx
  DETRAY_HOST
  constexpr auto set_material(const material_id mat_id, const dindex mat_idx)
      -> void {
    m_mat_link = {mat_id,
                  static_cast<typename material_link::index_type>(mat_idx)};
  }

  /// @returns true if the volume descriptor has a valid material link
  DETRAY_HOST_DEVICE
  constexpr auto has_material() const -> bool {
    return (m_mat_link.id() != material_link::id_type::e_none) &&
           !m_mat_link.is_invalid();
  }

  /// @returns link to all acceleration data structures - const access
  DETRAY_HOST_DEVICE constexpr auto accel_link() const
      -> const accel_link_type& {
    return m_accel_links;
  }

  /// @returns acc data structure link for a specific type of object - const
  /// version with runtime indexing
  DETRAY_HOST_DEVICE constexpr auto accel_link(const ID& id) const
      -> const acc_link_t& {
    static_assert(static_cast<std::size_t>(ID::e_size) >= 1);
    assert(static_cast<std::size_t>(id) < static_cast<std::size_t>(ID::e_size));
    return accel_link_helper<static_cast<ID>(
        static_cast<std::size_t>(ID::e_size) - 1)>(id);
  }

  /// @returns acc data structure link for a specific type of object - const
  template <ID obj_id>
  DETRAY_HOST_DEVICE constexpr auto accel_link() const -> const acc_link_t& {
    return detail::get<obj_id>(m_accel_links);
  }

  /// Set surface finder link from @param link
  template <ID obj_id>
  DETRAY_HOST constexpr auto set_accel_link(const acc_link_t link) -> void {
    detail::get<obj_id>(m_accel_links) = link;
  }

  /// Set link from @param id and @param index of the acceleration data
  /// structure (e.g. type and index of a grid in the accelerator store)
  template <ID obj_id>
  DETRAY_HOST constexpr auto set_accel_link(
      const typename acc_link_t::id_type id,
      const typename acc_link_t::index_type index) -> void {
    detail::get<obj_id>(m_accel_links) = acc_link_t{id, index};
  }

  /// Set link for a type of surfaces ( @param obj_id ) from @param id
  /// and @param index of the acceleration data structure (e.g. type and
  /// index of a grid in the accelerator store)
  DETRAY_HOST constexpr auto set_accel_link(
      const ID obj_id, const typename acc_link_t::id_type accel_id,
      const typename acc_link_t::index_type index) -> void {
    m_accel_links[obj_id] = acc_link_t{accel_id, index};
  }

  /// Equality operator
  ///
  /// @param rhs is the right-hand side to compare against.
  DETRAY_HOST_DEVICE
  constexpr auto operator==(const volume_descriptor& rhs) const -> bool {
    return (m_id == rhs.m_id && m_index == rhs.m_index &&
            m_accel_links == rhs.m_accel_links);
  }

  /// @returns a string stream that prints the volume details
  DETRAY_HOST
  friend std::ostream& operator<<(std::ostream& os,
                                  const volume_descriptor& v_desc) {
    os << "id = " << v_desc.id() << "(" << static_cast<int>(v_desc.id()) << ")";
    os << " | index = " << v_desc.index();
    os << " | trf. = " << v_desc.transform();
    os << " | acc link: " << v_desc.accel_link();
    os << " | sf link: " << v_desc.sf_link();
    os << " | mat link: " << v_desc.material();
    return os;
  }

 private:
  /// @returns acc data structure link for a specific type of object - const
  /// version with compile-time indexing
  template <ID obj_id>
  DETRAY_HOST_DEVICE constexpr auto accel_link_helper(const ID& id) const
      -> const acc_link_t& {
    if (id == obj_id) {
      return accel_link<obj_id>();
    } else if constexpr (static_cast<std::size_t>(obj_id) > 0) {
      return accel_link_helper<static_cast<ID>(
          static_cast<std::size_t>(obj_id) - 1)>(id);
    } else {
      __builtin_unreachable();
    }
  }

  /// How to interpret the boundary values
  volume_id m_id{volume_id::e_unknown};

  /// Volume index in the detector's volume container
  dindex m_index{dindex_invalid};

  /// Volume index in the detector's volume container
  dindex m_transform{dindex_invalid};

  /// Index range for every object type
  sf_link_type m_sf_links{};

  /// Volume material link
  material_link m_mat_link = {};

  /// Links for every object type to an acceleration data structure
  accel_link_type m_accel_links{};
};

}  // namespace detray
